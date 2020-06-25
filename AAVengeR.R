library(yaml)
library(igraph)
library(stringdist)
library(ShortRead)
library(tidyverse)
library(parallel)
library(GenomicRanges)
library(gintools)
library(lubridate)
options(stringsAsFactors = FALSE)


# Read the config file.
# This file contains processing parameters and also points to the sample configuration file
# which contains sample specific parameters.
configFile <- commandArgs(trailingOnly = TRUE)
if(! file.exists(configFile)) stop('Error -- configuration file not found.')
config  <- read_yaml(configFile)
source(file.path(config$softwareDir, 'AAVengeR.lib.R'))


#config <- read_yaml('/home/everett/projects/AAVengeR_runs/configs/config.Sabatino.yml')
#config <- read_yaml('/home/everett/SparkAAV/AAVengeR/configs/config.vector.yml')
#source(file.path(config$softwareDir, 'AAVengeR.lib.R'))



# Prepare run session.
if(dir.exists(config$outputDir)) stop('Error -- output directory already exists.')
dir.create(config$outputDir)
invisible(sapply(c('tmp', 'readIDs', 'readsRemoved', 'seqChunks', 'sampleReads', 'logs', 'fragReads'), 
                 function(x) dir.create(file.path(config$outputDir, x))))
write(capture.output(sessionInfo()), file = file.path(config$outputDir, 'sessionInfo.txt'))
config$startTime <- ymd_hms(format(Sys.time(), "%y-%m-%d %H:%M:%S"))
config$logFile <- file.path(config$outputDir, 'logs', 'log')
write(date(), file = config$logFile)




# Read in the sample data and add default column values if not present.
samples <- read_delim(config$sampleConfigFile, delim  = ',', col_names = TRUE, col_types = cols())
checkConfigFilePaths(samples$refGenomeBLATdb)
if(config$indexReads.rc) samples$index1Seq <- as.character(reverseComplement(DNAStringSet(samples$index1Seq)))
if('vectorSeqFile' %in% names(samples)) checkConfigFilePaths(samples$vectorSeqFile)
if(! 'subject' %in% names(samples))   samples$subject <- 'subject'
if(! 'replicate' %in% names(samples)) samples$replicate <- 1
if(any(grepl('~|\\|', paste(samples$subject, samples$sample, samples$replicate)))) 
  stop('Error -- tildas (~) are reserved characters and can not be used in the subject, sample, or replicate sample configuration columns.')

samples$uniqueSample <- paste0(samples$subject, '~', samples$sample, '~', samples$replicate)
if(any(duplicated(samples$uniqueSample))) stop('Error -- all subject, sample, replicate id combinations are not unique.')


cluster <- makeCluster(config$demultiplexing.CPUs)
clusterExport(cluster, c('config', 'samples'))

# Quality trim virus reads and break reads.
invisible(parLapply(cluster, 
                    list(c(config$breakReadsFile,  config$sequence.chunk.size, 'breakReads',  file.path(config$outputDir, 'seqChunks')),
                         c(config$virusReadsFile,  config$sequence.chunk.size, 'virusReads',  file.path(config$outputDir, 'seqChunks'))), 
                 function(x){
                   library(ShortRead)
                   source(file.path(config$softwareDir, 'AAVengeR.lib.R'))
                   qualTrimReads(x[[1]], x[[2]], x[[3]], x[[4]])
                }))


# Collate trimmed reads.
system(paste('cat', file.path(config$outputDir, 'seqChunks', 'virus*'), ' > ', file.path(config$outputDir, 'trimmedVirusReads.fastq')))
system(paste('cat', file.path(config$outputDir, 'seqChunks', 'break*'), ' > ', file.path(config$outputDir, 'trimmedBreakReads.fastq')))
system(paste('rm',  file.path(config$outputDir, 'seqChunks', '*')))


# Convert reads to DNAStringSets and sync reads.
logMsg(config, 'Converting index reads to DNA strings.', config$logFile)
index1Reads <- readFastq(config$index1ReadsFile)
index1Reads <- Reduce('append', parLapply(cluster, split(index1Reads, ntile(1:length(index1Reads), config$demultiplexing.CPUs)), 
               function(x){source(file.path(config$softwareDir, 'AAVengeR.lib.R')); shortRead2DNAstringSet(x)}))

if(config$correctGolayIndexReads){
  logMsg(config, 'Correcting golay bar code reads.', config$logFile)
  index1Reads <- Reduce('append', parLapply(cluster, split(index1Reads, ntile(1:length(index1Reads), config$demultiplexing.CPUs)), golayCorrection))
}

logMsg(config, 'Converting virus reads to DNA strings.', config$logFile)
virusReads <- readFastq(file.path(config$outputDir, 'trimmedVirusReads.fastq'))
virusReads <- Reduce('append', parLapply(cluster, split(virusReads, ntile(1:length(virusReads), config$demultiplexing.CPUs)), 
                                          function(x){source(file.path(config$softwareDir, 'AAVengeR.lib.R')); shortRead2DNAstringSet(x)}))

logMsg(config, 'Converting break reads to DNA strings.', config$logFile)
breakReads <- readFastq(file.path(config$outputDir, 'trimmedBreakReads.fastq'))
breakReads <- Reduce('append', parLapply(cluster, split(breakReads, ntile(1:length(breakReads), config$demultiplexing.CPUs)), 
                                          function(x){source(file.path(config$softwareDir, 'AAVengeR.lib.R')); shortRead2DNAstringSet(x)}))

invisible(file.remove(file.path(config$outputDir, 'trimmedVirusReads.fastq')))
invisible(file.remove(file.path(config$outputDir, 'trimmedBreakReads.fastq')))


# Slow
logMsg(config, 'Synchronizing trimmed reads.', config$logFile)
reads <- syncReads(index1Reads, virusReads, breakReads)
index1Reads <- reads[[1]];  virusReads  <- reads[[2]];  breakReads  <- reads[[3]]
rm(reads)
gc()


# Split the trimmed reads into chunks for parallel processing.
logMsg(config, 'Distributing read data into chunks ...', config$logFile)
chunkNum <- 1
d <- tibble(i = ntile(1:length(index1Reads), config$demultiplexing.CPUs), n = 1:length(index1Reads))
invisible(lapply(split(d, d$i), function(x){
  index1Reads <- index1Reads[min(x$n):max(x$n)]
  virusReads  <- virusReads[min(x$n):max(x$n)]
  breakReads  <- breakReads[min(x$n):max(x$n)]
  save(index1Reads, virusReads, breakReads, file = file.path(config$outputDir, 'seqChunks', chunkNum))
  chunkNum <<- chunkNum + 1
}))
logMsg(config, paste0(chunkNum-1, ' data chunks created.'), config$logFile)

# Clean up and free up memory. 
rm(d, chunkNum, index1Reads, virusReads, breakReads)
gc()


save(list = ls(all=TRUE), file = file.path(config$outputDir, 'savePoint1.RData'))
logMsg(config, 'Starting sample chunk threads, resetting timer.', config$logFile)
config$startTime <- ymd_hms(format(Sys.time(), "%y-%m-%d %H:%M:%S"))



invisible(parLapply(cluster, list.files(file.path(config$outputDir, 'seqChunks'), full.names = TRUE), function(f){
#invisible(lapply(list.files(file.path(config$outputDir, 'seqChunks'), full.names = TRUE), function(f){
#invisible(lapply("/home/everett/projects/AAVengeR_runs/outputs/Sabatino/seqChunks/1.RData", function(f){
  library(ShortRead)
  library(tidyverse)
  source(file.path(config$softwareDir, 'AAVengeR.lib.R'))
  
  load(f)
  
  # Capture the chunk identifier.
  chunk.n <- unlist(str_match_all(f, '(\\d+)$'))[2]
  logFile <-  file.path(config$outputDir, 'logs', paste0('seqChunk_', chunk.n, '.log'))
  
  # Loop through samples in sample data file.
  invisible(lapply(1:nrow(samples), function(r){
    rowNum <- r
    r <- samples[r,]
    
    # Trim over-read sequences.
    virusReads <- trimOverReadSeq(virusReads, r$virusRead.overReadSeq)
    breakReads <- trimOverReadSeq(breakReads, r$breakRead.overReadSeq)
    
    # Ensure that reads are long enough for upcoming sample specific tests.
    preReadSizeCheck1 <- names(virusReads)
    virusReads <- virusReads[width(virusReads) >= (config$virusReads.minLTRseqLength + config$trimmedRead.minLength)]
    breakReads <- breakReads[width(breakReads) >= (nchar(r$breakReadLinkerSeq) + config$trimmedRead.minLength)]
    
    # Sync becuase reads may of been lost in the previous length filter.
    reads <- syncReads(index1Reads, virusReads, breakReads)
    index1Reads <- reads[[1]];  virusReads  <- reads[[2]];  breakReads  <- reads[[3]]
    
    # Report reads removed from min read size filter1.
    ids <- preReadSizeCheck1[! preReadSizeCheck1 %in% names(virusReads)]
    if(length(ids) > 0) writeLines(ids, file.path(config$outputDir, 'readsRemoved', paste0('readSizeCheck1.', r$uniqueSample, '.', chunk.n)))
    
    
    if(length(virusReads) == 0){
      logMsg(config, paste0('Chunk ', chunk.n, ': (', r$uniqueSample, ') No reads remaining after first read size filter.'), logFile)
      return()
    }
    
    # Create barcode demultiplexing vectors.
    v1 <- rep(TRUE, length(virusReads))
    if('index1Reads.maxMismatch' %in% names(config)){
      v1 <- vcountPattern(r$index1Seq, index1Reads, max.mismatch = config$index1Reads.maxMismatch) > 0
      logMsg(config, paste0('Chunk ', chunk.n, ': (', r$uniqueSample, ') ', sum(v1), ' reads pass barcode filter.'), logFile)
    }
    
    
    # Create break read linker barcode demultiplexing vector.
    v2 <- rep(TRUE, length(breakReads))
    if('breakReads.linkerBarcode.maxMismatch' %in% names(config)){
      testSeq <- substr(r$breakReadLinkerSeq, r$breakReadLinkerBarcode.start, r$breakReadLinkerBarcode.end)
      v2 <- vcountPattern(testSeq, subseq(breakReads, r$breakReadLinkerBarcode.start, r$breakReadLinkerBarcode.end), max.mismatch = config$breakReads.linkerBarcode.maxMismatch) > 0
      logMsg(config, paste0('Chunk ', chunk.n, ': (', r$uniqueSample, ') ', sum(v2), ' reads pass linker code filter.'), logFile)
    }
    
    
    # Test to see if any reads demultiplex to this row of the sample table and then subset reads to this sample.
    i <- base::intersect(which(v1), which(v2))
    if(length(i) == 0){
      logMsg(config, paste0('Chunk ', chunk.n, ': (', r$uniqueSample, ') No reads demultiplexed.'), logFile)
      return()
    } 
    
    reads <- syncReads(index1Reads[i], virusReads[i], breakReads[i])
    index1Reads <- reads[[1]];  virusReads  <- reads[[2]]; breakReads  <- reads[[3]]
    
    if(length(virusReads) == 0){
      logMsg(config, paste0('Chunk ', chunk.n, ': (', r$uniqueSample, ') No reads remaining after post demultiplexing read synchronization.'), logFile)
      return()
    } else {
      logMsg(config, paste0('Chunk ', chunk.n, ': (', r$uniqueSample, ') ', length(virusReads), ' reads pass all demultiplexing filters.'), logFile)
    }
    
    
    
    # Test for break point reads that align to the vector plasmid.
    if(config$filter.removeVectorReadPairs){
      vectorReadIDs <- getVectorReadIDs(breakReads, config, r$vectorSeqFile)
    
      if(length(vectorReadIDs) > 0){
        virusReadsVector <- virusReads[names(virusReads) %in% vectorReadIDs]
        breakReadsVector <- breakReads[names(breakReads) %in% vectorReadIDs]
    
        names(virusReadsVector) <- paste0(names(virusReadsVector), '|', r$uniqueSample)
        names(breakReadsVector) <- paste0(names(breakReadsVector), '|', r$uniqueSample)
      
        writeFasta(virusReadsVector, file = file.path(config$outputDir, 'tmp', paste0(r$uniqueSample, '.virusReadsVector.', chunk.n, '.fasta')))
        writeFasta(breakReadsVector, file = file.path(config$outputDir, 'tmp', paste0(r$uniqueSample, '.breakReadsVector.', chunk.n, '.fasta')))
      }
      
      logMsg(config, paste0('Chunk ', chunk.n, ': (', r$uniqueSample, ') ', length(vectorReadIDs), ' reads aligned to the vector and will be removed.'), logFile)
      
      # Filter break point reads against read ids returned by vectorReadIDs().
      breakReads <- breakReads[! names(breakReads) %in% vectorReadIDs]
    
      
      # Sync reads. 
      reads <- syncReads(index1Reads, virusReads, breakReads)
      index1Reads <- reads[[1]];  virusReads  <- reads[[2]]; breakReads  <- reads[[3]]
      if(length(virusReads) == 0){
        logMsg(config, paste0('Chunk ', chunk.n, ': (', r$uniqueSample, ') No reads remaining after vector read removal and read synchronization.'), logFile)
        return()
      } else {
        logMsg(config, paste0('Chunk ', chunk.n, ': (', r$uniqueSample, ') ', length(virusReads), ' reads remain after removing vector aligned reads and synchronizing.'), logFile)
      }
    }

    
    
    # Test the start of virus reads.
    v1 <- rep(TRUE, length(virusReads))
    if('virusReads.startTest.maxMismatch' %in% names(config)){
      v1 <- Reduce('|', lapply(unlist(strsplit(r$virusLTRseq, ';')), function(x){ 
           testSeq <- substr(unlist(strsplit(x, ','))[2], 1,  config$virusReads.startTest.length)
           vcountPattern(testSeq, subseq(virusReads, 1, config$virusReads.startTest.length), max.mismatch = config$virusReads.startTest.maxMismatch) > 0
        }))
      
      logMsg(config, paste0('Chunk ', chunk.n, ': (', r$uniqueSample, ') ', sum(!v1), ' reads removed by virus read start test filter (', config$virusReads.startTest.length, ' NTs).'), logFile)
    }
    
    
    # Test the entire LTR sequence.
    v2 <- rep(TRUE, length(virusReads))
    if('virusReads.fullTest.maxMismatch' %in% names(config)){
      v2 <- Reduce('|', lapply(unlist(strsplit(r$virusLTRseq, ';')), function(x){ 
        vcountPattern(unlist(strsplit(x, ','))[2], subseq(virusReads, 1, nchar(unlist(strsplit(x, ','))[2])), max.mismatch = config$virusReads.fullTest.maxMismatch) > 0
      }))
      
      logMsg(config, paste0('Chunk ', chunk.n, ': (', r$uniqueSample, ') ', sum(!v2), ' reads removed by virus read full test filter.'), logFile)
    }
    
    
    
    # Test break read common linker.
    v3 <- rep(TRUE, length(virusReads))
    if('breakReads.linkerCommon.maxMismatch' %in% names(config)){
      testSeq <- substr(r$breakReadLinkerSeq, r$breakReadLinkerCommon.start, r$breakReadLinkerCommon.end)
      v3 <- vcountPattern(testSeq, subseq(breakReads, r$breakReadLinkerCommon.start,r$breakReadLinkerCommon.end), max.mismatch = config$breakReads.linkerCommon.maxMismatch) > 0
 
      logMsg(config, paste0('Chunk ', chunk.n, ': (', r$uniqueSample, ') ', sum(!v3), ' reads removed by the break read common linker filter.'), logFile)  
    }
    
    # Test to see which reads pass the last three filters.
    logMsg(config, paste0('Chunk ', chunk.n, ': (', r$uniqueSample, ') ', sum(!(v1 & v2 & v3)), ' reads removed after all filters.'), logFile)
    
    
    i <- base::Reduce(base::intersect, list(which(v1), which(v2), which(v3)))
    if(length(i) == 0){
      logMsg(config, paste0('Chunk ', chunk.n, ': (', r$uniqueSample, ') No reads remaining after user specified read tests and read synchronization.'), logFile)
      return()
    }
    
    # Subset and sync reads.
    reads <- syncReads(index1Reads[i], virusReads[i], breakReads[i])
    index1Reads <- reads[[1]];  virusReads  <- reads[[2]]; breakReads  <- reads[[3]]
  
    
    # Trim leading adapter sequences.
    if(config$virusReads.captureLTRseqs){
      o <- captureLTRseqs(virusReads, r$virusLTRseq, rowNum)
      
      if(length(o) == 0){
        logMsg(config, paste0('Chunk ', chunk.n, ': (', r$uniqueSample, ') No reads returned from LTR/ITR capture.'), logFile)
        return()
      } else {
        logMsg(config, paste0('Chunk ', chunk.n, ': (', r$uniqueSample, ') ', length(o), ' reads remaining after LTR/ITR capture.'), logFile)
      }
      
      o$reads <- o$reads[match(o$LTRs$id, names(o$reads))]                      # Order the reads to match LTR table.
      o$reads <- subseq(o$reads, start = nchar(o$LTRs$LTRseq)+1)                # Remove matched LTR sequences from the reads.
      o$reads <- o$reads[width(o$reads) >= config$trimmedRead.minLength]        # Remove reads which are now too short.
      o$LTRs  <- subset(o$LTRs, id %in% names(o$reads))                         # Now that we have removed some reads, trim the LTR table.
      o$LTRs$readSeq <- as.character(o$reads[match(o$LTRs$id, names(o$reads))]) # Add the trimmed sequences to the LTR table so that we can latter find additional NTs between LTRs and genomic alignments.
      
       saveRDS(o$LTRs, file = file.path(config$outputDir, 'tmp', paste0(r$uniqueSample, '.LTRseqs.', chunk.n, '.rds')))
       
       ids <- names(virusReads)[! names(virusReads) %in% names(o$reads)]
       if(length(ids) > 0) writeLines(ids, file.path(config$outputDir, 'readsRemoved', paste0('LTRcapture.', chunk.n)))
       
       virusReads <- o$reads
       
       reads <- syncReads(index1Reads, virusReads, breakReads)
       index1Reads <- reads[[1]];  virusReads  <- reads[[2]]; breakReads  <- reads[[3]]
       if(length(virusReads) == 0){
         logMsg(config, paste0('Chunk ', chunk.n, ': (', r$uniqueSample, ') No reads remaining after LTR capture and trimming.'), logFile)
         return()
       }
    } else {
      # This route does not support multiple LTR sequences -- add check.
      s <- unlist(strsplit(r$virusLTRseq, ','))[2]
      logMsg(config, paste0('Chunk ', chunk.n, ': (', r$uniqueSample, ') Trimming leading seq: ', s), logFile)
      virusReads <- trimLeadingSeq(virusReads, s)
    }
    
    
    # Capture random ids.
    randomIDs <- data.frame()
    if('breakReadLinkerBarcode.start' %in% names(r) & 'breakReadLinkerBarcode.end' %in% names(r)){
      randomIDs <- as.character(subseq(breakReads, r$breakReadLinkerRandomID.start, r$breakReadLinkerRandomID.end))
      randomIDs <- data.frame(randomSeqID = unname(randomIDs), readID = names(randomIDs))
      logMsg(config, paste0('Chunk ', chunk.n, ': (', r$uniqueSample, ') captured random linker ids.'), logFile)
    }
    
    # Remove leading linker from break point reads.
    breakReads <- trimLeadingSeq(breakReads, r$breakReadLinkerSeq)
    
    # Select reads which have the minimum read lengths post trimming.
    preFilterReads <- names(virusReads)
    virusReads <- virusReads[width(virusReads) >= config$trimmedRead.minLength]
    breakReads <- breakReads[width(breakReads) >= config$trimmedRead.minLength]
    
  
    # Report reads removed from min read size filter2.
    ids <- preFilterReads[! preFilterReads %in% names(virusReads)]
    if(length(ids) > 0) writeLines(ids, file.path(config$outputDir, 'readsRemoved', paste0('readSizeCheck2.', r$uniqueSample, '.', chunk.n)))
    
    # Sync reads after selecting for reads with min. length post-trimming.
    reads <- syncReads(index1Reads, virusReads, breakReads)
    index1Reads <- reads[[1]];  virusReads  <- reads[[2]]; breakReads  <- reads[[3]]
    if(length(virusReads) == 0){
      logMsg(config, paste0('Chunk ', chunk.n, ': (', r$uniqueSample, ') No reads remaining after second read size filter.'), logFile)
      return()
    } else {
      logMsg(config, paste0('Chunk ', chunk.n, ': (', length(virusReads), ') reads remaining after read length filters.'), logFile)
    }
    
    
    # Write out final reads. Add sample names to read IDs.
    names(breakReads) <- paste0(names(breakReads), '|', r$uniqueSample)
    names(virusReads) <- paste0(names(virusReads), '|', r$uniqueSample)
    
    save(randomIDs, file = file.path(config$outputDir, 'tmp', paste0(r$uniqueSample, '.randomIDs.', chunk.n, '.RData')))
    writeFasta(breakReads,  file = file.path(config$outputDir, 'tmp', paste0(r$uniqueSample, '.breakReads.', chunk.n, '.fasta')))
    writeFasta(virusReads,  file = file.path(config$outputDir, 'tmp', paste0(r$uniqueSample, '.virusReads.', chunk.n, '.fasta')))
    
    logMsg(config, paste0('Chunk ', chunk.n, ': (', r$uniqueSample, ') completed with ', length(virusReads), ' reads.'), logFile)
  }))
  
  logMsg(config, paste0('Read data chunk ', chunk.n, ' completed.'), file.path(config$outputDir, 'logs', 'log'))
}))

stopCluster(cluster)


# Collate demultiplexed reads using file name snibets.
collateSampleReads('virusReads')
collateSampleReads('breakReads')


# Collate read pairs where the break point read aligns well to the vector for downstream analysis.
if(config$filter.removeVectorReadPairs){
  logMsg(config, 'Collating read pairs which were removed because the break read aligned to the vector.', config$logFile)
  collateSampleReads('virusReadsVector')
  collateSampleReads('breakReadsVector')
}


# Organize read files by reference genome and read source (virusRead or breakRead).
# This will allow samples to be aligned to different genomes defined in the samples table.
d <- tibble(file = list.files(file.path(config$outputDir, 'sampleReads'), pattern = '\\.virusReads\\.|\\.breakReads\\.'),
            uniqueSample = unlist(lapply(strsplit(file, '\\.'), function(x) paste0(x[1:(length(x) - 2)], collapse = '.'))),
            source = ifelse(grepl('virusRead', file), 'virusReads', 'breakReads')) %>%
     left_join(select(samples, uniqueSample, refGenomeBLATdb), by = 'uniqueSample')



# Align the sample FASTA files defined in the previous data frame to their respective genomes. 
# This is done by grouping reads by type (virusReads / breakReads) and reference geneome 
# and using parLapply() to disrubte read chunks across a number of CPUs. Unique reads 
# are nested and then unnnested after alignment.


logMsg(config, 'Starting refGenome alignments.', config$logFile)

alignments <- bind_rows(lapply(split(d, paste(d$source, d$refGenomeBLATdb)), function(x){
  f <- tmpFile()
  
  # Concatenate all the source / genome fasta files into a single file.
  system(paste('cat ', paste0(file.path(config$outputDir, 'sampleReads', x$file), collapse = ' '), ' > ', 
               paste0(file.path(config$outputDir, 'tmp', paste0(f, '.fasta')))))
  
  fasta <- readFasta(file.path(config$outputDir, 'tmp', paste0(f, '.fasta')))
  invisible(file.remove(file.path(config$outputDir, 'tmp', paste0(f, '.fasta'))))
  
  # Here we create a table of unique sequences to align to the reference geneome.
  # The ids associated with each unique sequnce are stored so that they can be 
  # expanded back to match the original data after the alignemnts. 
  
  d <- tibble(id  = as.character(fasta@id),
              seq = as.character(fasta@sread)) %>%
       group_by(seq) %>%
        summarise(ids = list(id)) %>%
       ungroup() %>%
       mutate(seqID = paste0('s', 1:n()))
  
  
  fasta <- DNAStringSet(d$seq)
  names(fasta) <- d$seqID
  
  db <- x$refGenomeBLATdb[1]
  cluster <- makeCluster(config$genomAlignment.CPUs)
  clusterExport(cluster, c('config', 'db'), envir = environment())
  
  logMsg(config, paste0('Starting ', x$source[1], ' alignments'), config$logFile)
  
  
  r <- bind_rows(parLapply(cluster, split(fasta, ceiling(seq_along(fasta)/config$alignment.chunk.size)), function(y){
         library(ShortRead)
         library(tidyverse)
         source(file.path(config$softwareDir, 'AAVengeR.lib.R'))
         alignReads.BLAT(y, db)
       }))
  
  stopCluster(cluster)
  
  d <- left_join(select(d, -seq), r,  by = c("seqID" = "qName"))
  d <- unnest(d, ids)
  
  d$source <- x$source[1]
  select(d, -seqID)
}))

logMsg(config, 'RefGenome alignments completed.', config$logFile)

save(list = ls(all=TRUE), file = file.path(config$outputDir, 'savePoint2.RData'))


# Apply generic alignment filters.
alignments <- 
  filter(alignments, 
         alignmentPercentID >= config$alignment.genome.minPercentID,
         tNumInsert  <= 1, 
         qNumInsert  <= 1,
         tBaseInsert <= 2,
         qBaseInsert <= 2) 


# Apply direction specific filters.
virusReads.aln <- alignments[which(alignments$source == 'virusReads'),]
breakReads.aln <- alignments[which(alignments$source == 'breakReads'),]


# Break point reads have static 5' adapters which are removed with cutadapt.
# Only consider reads which immediately align to the genome.
breakReads.aln <- filter(breakReads.aln, qStart <= 3, matches >= config$trimmedRead.minLength)


# Virus reads have had over-read sequences removed.
# Only consider virus reads where the end of the read aligns to the genome.
virusReads.aln <- filter(virusReads.aln, (qSize - qEnd) <= 3, matches >= config$trimmedRead.minLength)


# If we are not capturing LTR sequences then a static LTR was removed from viral reads
# and we can apply a filter to the beginning of the trimmed read.
if(! config$virusReads.captureLTRseqs){
  logMsg(config, 'Requiring virus reads to begin alignments near their start because the ITR/LTR should of been fully trimmed off.', config$logFile)
  virusReads.aln <- filter(virusReads.aln, qStart <= 3)
}


alignments <- bind_rows(virusReads.aln, breakReads.aln)
rm(virusReads.aln, breakReads.aln)


# Rename alignment column headers by appending read sources.
alignments <- lapply(split(alignments, alignments$source), function(x){
    names(x) <- paste0(names(x), '.', x$source[1])
    x
})


# Combine read alignments into rows and extract sample names from read ids. 

v <- select(alignments[["virusReads"]][ alignments[["virusReads"]]$ids.virusReads %in% alignments[["breakReads"]]$ids.breakReads,],
            ids.virusReads, qStart.virusReads, strand.virusReads, tName.virusReads, tStart.virusReads, tEnd.virusReads)
v$ids.virusReads <- factor(v$ids.virusReads )

b <- select(alignments[["breakReads"]][ alignments[["breakReads"]]$ids.breakReads %in% alignments[["virusReads"]]$ids.virusReads,],
            ids.breakReads, strand.breakReads, tName.breakReads, tStart.breakReads, tEnd.breakReads)
b$ids.breakReads <- factor(b$ids.breakReads)

rm(alignments)

logMsg(config, 'Starting to create fragments from alignment data.', config$logFile)

### readPerChunk <- 20000
cluster <- makeCluster(3)
clusterExport(cluster, c('config', 'b'))

v$s <- ntile(1:nrow(v), 3)


# Here we build a table of potential fragments by matching up virus and break point reads
# and filtering with basic sanity checks. This is done via a left join of the break reads 
# to the virus reads followed by filtering. This join is split between three groupings 
# because of the potential to exceed R's internal table row limit during the joins.

# capture qStart.virusReads

frags <- bind_rows(parLapply(cluster, split(v, v$s), function(x){
         library(dplyr)
         library(stringr)
  
         frags <- left_join(x, b, by = c('ids.virusReads' = 'ids.breakReads')) %>% tidyr::drop_na()
       
         i <- which(frags$tName.virusReads != frags$tName.breakReads)
         if(length(i) > 0) frags <- frags[-i,]
         if(nrow(frags) == 0) return(data.frame())
         
         i <- which(frags$strand.virusReads == frags$strand.breakReads)
         if(length(i) > 0) frags <- frags[-i,]
         if(nrow(frags) == 0) return(data.frame())
         
         mutate(frags, 
                fragStart = ifelse(strand.virusReads == '+', tStart.virusReads + 1, tStart.breakReads + 1),
                fragEnd   = ifelse(strand.virusReads == '+', tEnd.breakReads + 1,   tEnd.virusReads + 1),
                fragTest  = ifelse(strand.virusReads == '+', tStart.virusReads < tEnd.breakReads, tStart.breakReads < tEnd.virusReads),  
                fragWidth = (fragEnd - fragStart) + 1) %>%
         filter(fragTest == TRUE, 
                fragWidth <= config$fragments.maxLength,
                fragWidth >= config$fragments.minLength) %>%
         mutate(uniqueSample = unlist(lapply(str_split(ids.virusReads, '\\|'), '[[', 2 )),
                readID = sub('\\|.+$', '', ids.virusReads)) 
       }))

stopCluster(cluster)

logMsg(config, 'Initial fragment generation completed.', config$logFile)

frags <- left_join(frags, dplyr::select(samples, refGenomeName, uniqueSample), by = 'uniqueSample')
frags <- mutate(frags, fragID = paste0(uniqueSample, ';', tName.virusReads, ';', strand.virusReads, ';', fragStart, ';', fragEnd))


# Remove read frags which share a random linker id and spand more than one sample.
if(config$removeReadFragsWithSameRandomID) frags <- removeReadFragsWithSameRandomID(frags, config)


# Convert read fragments into grouped fragments.
frags <- group_by(frags, fragID) %>%
  summarise(reads = n(), refGenomeName = refGenomeName[1], virusReadAlignmentStart = qStart.virusReads[1], readIDs = list(readID)) %>%
  ungroup() %>%
  separate(fragID, c('uniqueSample', 'seqnames', 'strand', 'start', 'end'), sep = ';') 

frags <- unpackUniqueSampleID(frags)



# Standardize fragment boundaries.
frags <- standardizedFragments(frags, config)


# Now that the boundaries have been standardized, identical fragments need to be merged.
# This will collapse read id lists to comma delimited strings.
frags <- bind_rows(dplyr::group_by(frags, uniqueSample, seqnames, start, end, strand) %>%
                   dplyr::mutate(reads = sum(reads), readID.list = list(unlist(readIDs)), intSiteRefined = all(intSiteRefined), breakPointRefined = all(breakPointRefined)) %>%
                   dplyr::slice(1) %>%
                   dplyr::ungroup() %>%
                   dplyr::select(-readIDs) %>%
                   dplyr::mutate(posid = paste0(seqnames, strand, ifelse(strand == '+', start, end))))


# Save the fragments for downstream analyses.
saveRDS(frags, file = file.path(config$outputDir, 'readFrags.rds'))

frags$s <- paste0(frags$subject, '_', frags$sample, '_', frags$reads)
createFragUCSCTrack(frags, title = 'AAVengeR_fragments', outputFile = file.path(config$outputDir, 'fragments.ucsc'), label = 's')


posIDclusterTable <- 
  dplyr::group_by(unnest(frags, readID.list), subject, readID.list) %>%
  dplyr::summarise(nPositions = n_distinct(posid), posids = list(posid)) %>%
  dplyr::ungroup() %>%
  dplyr::filter(nPositions > 1) 

# posIDclusterTable -- for each read, provide the multiple posids to which it maps.
# subject readID.list                                   nPositions posids    
# <chr>   <chr>                                              <int> <list>    
#   1 pH19    M03249:365:000000000-C3CH4:1:1101:11087:17666          2 <chr [2]> 
#   2 pH19    M03249:365:000000000-C3CH4:1:1101:11121:16397         13 <chr [13]>
#   3 pH19    M03249:365:000000000-C3CH4:1:1101:11589:17142         13 <chr [13]>

posIDclusters <- do.call(rbind, lapply(split(posIDclusterTable, posIDclusterTable$subject), function(s){
                   nodes <- unique(unlist(s$posids))
                   edges <- bind_rows(lapply(split(s, 1:nrow(s)), function(x){ 
                              data.frame(t(combn(unlist(x$posids), 2)))
                            }))

                   network = graph_from_data_frame(edges, vertices = nodes, directed = F)
                   data.frame(subject = s$subject[1], clusters = I(list(lapply(decompose(network), function(x) V(x)$name))))
                 }))

# posIDclusters
# subject     clusters
# <chr>     <list of posids clusters>
# pH19 c("chr11....
# pHO2 c("chr20....
# pJ60 c("chr2+....

# Clusters for pH19 (first row from posIDclusters)
# > posIDclusters[1,]$clusters[[1]]
# [[1]]
# [1] "chr11-8358680"  "chr11+10996378"
# 
# [[2]]
# [1] "chr1-15117666"  "chr2-56733954"  "chr3-54854044"  "chr4-19729222"  "chr5-6580806"   "chr5+12786336" 
# [7] "chr29-15327160" "chr32+18112660" "chrX+121905654" "chr10-59401573" "chr12+55436282" "chr20-37905558" 


frags <- bind_rows(lapply(split(frags, frags$subject), function(s){
           if(s$subject[1] %in% posIDclusters$subject){
             s$posIDcluster <- sapply(s$posid, findPosIDposIDclusters, posIDclusters[posIDclusters$subject == s$subject[1],]$clusters[[1]])
           } else {
             s$posIDcluster <- NA
           }
           s
         }))

rm(posIDclusters, posIDclusterTable)

# Create a second save point.
save(list = ls(all=TRUE), file = file.path(config$outputDir, 'savePoint3.RData'))
#--------------------------------------------------------------------------------------------------



if(config$virusReads.captureLTRseqs){
  
  # Examine the LTR/ITR sequence remnants (for each read) from each fragment and determine how similiar they are to one another.
  # For each fragment, the most similiar 95% of the reads are considered and representative remnant sequences are determined.
  cluster <- makeCluster(config$demultiplexing.CPUs)
  clusterExport(cluster, c('config'))
  frags$s <- ntile(1:nrow(frags), config$demultiplexing.CPUs)
  
  frags <- bind_rows(parLapply(cluster, split(frags, frags$s), function(a){
              library(dplyr)
              source(file.path(config$softwareDir, 'AAVengeR.lib.R'))
    
              LTR.table <- bind_rows(lapply(list.files(file.path(config$outputDir, 'tmp'), pattern = 'LTRseqs', full.names = TRUE), readRDS))
    
              bind_rows(lapply(1:nrow(a), function(x){
                frag <- a[x,]
                ltrs <- LTR.table[match(unlist(frag$readID.list), LTR.table$id),]
  
                r <- representativeSeq(ltrs$LTRseq)
                frag$maxReadPercentDiff <- r[[1]]
                frag$ltrRepSeq <- r[[2]]
                
                additionalLTRnts <- substr(ltrs$readSeq, 1, frag$virusReadAlignmentStart)
                
                if(all(nchar(additionalLTRnts) == 0)){
                  frag$maxReadPercentDiff2 <- frag$maxReadPercentDiff 
                  frag$ltrRepSeq2 <- frag$ltrRepSeq 
                }else{
                  r <- representativeSeq(paste0(ltrs$LTRseq, frag$additionalLTRnts))
                  frag$maxReadPercentDiff2 <- r[[1]]
                  frag$ltrRepSeq2 <- r[[2]]
                }
                
                frag
              }))
           }))
  
  stopCluster(cluster)
  
  
  # Remove fragments where there is not a clear concensus remnant. 
  config$LTR_remnant_maxDifference <- 0.15
  frags <- subset(frags, frags$maxReadPercentDiff2 <= config$LTR_remnant_maxDifference)
  
  
  # Collapse fragments to individual sites. 
  sites <- bind_rows(lapply(split(frags, paste(frags$subject, frags$sample, frags$posid, frags$posIDcluster)), function(x){

    if(nrow(x) > 1){
      i <- rep(TRUE, nrow(x))
      r2 <- representativeSeq(x$ltrRepSeq2)
      
      if(r2[[1]] > 0.15){
        # There is a conflict, one or more fragments have a markedly different LTR sequence.
        i <- as.vector(stringdist::stringdistmatrix(x$ltrRepSeq2, r2[[2]]) / nchar(x$ltrRepSeq2) < 0.15)
        if(sum(i)/nrow(x) >= 2/3){
          x <- x[i,]
        } else {
          return(tibble())
        }
      }
      
      r <- representativeSeq(x$ltrRepSeq)
      
      return(dplyr::select(x, subject, sample, replicate, uniqueSample, seqnames, start, end, strand, posid, posIDcluster) %>% 
             dplyr::mutate(estAbund = n(), position = ifelse(strand[1] == '+', start[1], end[1]), ltrRepSeq = r[[2]], ltrRepSeq2 = r2[[2]]) %>%
             dplyr::slice(1) %>%
             dplyr::mutate(reads = sum(x$reads), readID.list = list(unlist(x$readID.list)), fragmentsRemoved = sum(!i)) %>%
             dplyr::select(subject, sample, replicate, uniqueSample, seqnames, strand, position, posid, posIDcluster, estAbund, reads, readID.list, fragmentsRemoved, ltrRepSeq, ltrRepSeq2))
    }else{
      return(dplyr::select(x, subject, sample, replicate, uniqueSample, seqnames, start, end, reads, strand, posid, posIDcluster, readID.list, ltrRepSeq, ltrRepSeq2) %>% 
             dplyr::mutate(estAbund = n(), position = ifelse(strand[1] == '+', start[1], end[1]), fragmentsRemoved = 0) %>%
             dplyr::select(subject, sample, replicate, uniqueSample, seqnames, strand, position, posid, posIDcluster, estAbund, reads, readID.list, fragmentsRemoved, ltrRepSeq, ltrRepSeq2))
    }
  }))
} else {
  
  # ...
  
}
  

# Collapse multihits.
if(any(is.numeric(sites$posIDcluster))){
  a <- subset(sites, is.na(posIDcluster))
  b <- group_by(subset(sites, ! is.na(posIDcluster)), posIDcluster) %>%
       dplyr::mutate(estAbund = max(estAbund), reads = sum(reads), readID.list = list(unlist(readID.list)),
                     seqnames = 'chrUn', position = posIDcluster, strand = '*', 
                     fragmentsRemoved = sum(fragmentsRemoved),
                     posid = paste0('chrUn*', posIDcluster, '(', n_distinct(posid), ')')) %>%
       dplyr::slice(1) %>%
       ungroup()
    
  sites <- bind_rows(a, b)
}
  
sites$uniqueSample <- NULL
sites$posIDcluster <- NULL

save(sites, file = file.path(config$outputDir, 'sites.RData'))

d <- dplyr::mutate(sites, subject = sub('^p', '', subject), start = position, end = position, siteLabel = paste0(subject, '_', sample, '_', estAbund)) %>%
     dplyr::filter(seqnames != 'chrUn')
createIntUCSCTrack(d, siteLabel = 'siteLabel', outputFile = file.path(config$outputDir, 'sites.ucsc'))
system(paste0('cat ', file.path(config$outputDir, 'sites.ucsc'), ' ', file.path(config$outputDir, 'fragments.ucsc'), ' > ',file.path(config$outputDir, 'AAVengeR.ucsc')))


logMsg(config, 'done.', config$logFile)
