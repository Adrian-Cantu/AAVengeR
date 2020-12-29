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
#configFile <- commandArgs(trailingOnly = TRUE)
#if(! file.exists(configFile)) stop('Error -- configuration file not found.')
#config  <- read_yaml(configFile)
#source(file.path(config$softwareDir, 'AAVengeR2.lib.R'))

config <- read_yaml('/home/everett/projects/AAVengeR_runs/configs/config.Sabatino.yml2')
#config <- read_yaml('/home/everett/projects/AAVengeR_runs/configs/config.Sabatino.yml')
#config <- read_yaml('/home/everett/SparkAAV/AAVengeR/configs/config.vector.yml')
source(file.path(config$softwareDir, 'AAVengeR2.lib.R'))

config$startTime <- ymd_hms(format(Sys.time(), "%y-%m-%d %H:%M:%S"))


# Prepare run session.
if(dir.exists(config$outputDir)) stop('Error -- output directory already exists.')
dir.create(config$outputDir)
invisible(sapply(c('tmp', 'readIDs', 'readsRemoved', 'seqChunks', 'sampleReads', 'logs', 'fragReads'), 
                 function(x) dir.create(file.path(config$outputDir, x))))
write(capture.output(sessionInfo()), file = file.path(config$outputDir, 'sessionInfo.txt'))
config$logFile <- file.path(config$outputDir, 'logs', 'log')
write(date(), file = config$logFile)

# Max edit distance between ITR/LTR representative and all others for a fragment / length of representative 
if(! 'fragmentProcessing.rep.maxDifference' %in% names(config)) config$fragmentProcessing.rep.maxDifference <- 0.15

# Fragment ITR/LTR assembly conflict resolution.
# When fragments are assembled into sites, their representative ITR/LTR sequences are compared 
# similiar to how fragment reads are compared when building fragments. This is important because 
# two or more closely spaced events may be combined during fragment boundary standardization. 
# When a the ITR/LTR representative sequences of assembled fragments are not similar enough 
# we try to salvage the site by retaining the fragments which are similiar to one another and 
# discarding the disimiliar fragments. This parameter defines the min. threshold percentange of fragments 
# with similiar ITR/LTR sequences to allow a site to progress. Future versions of the software may 
# consider spliting the site.
if(! 'fragmentProcessing.assemblyConflictResolution' %in% names(config)) config$fragmentProcessing.assemblyConflictResolution <- 2/3
if(config$fragmentProcessing.assemblyConflictResolution < 0.5) stop('fragmentProcessing.assemblyConflictResolution parmater must be > 0.5')



# Read in the sample data and add default column values if not present.
samples <- read_delim(config$sampleConfigFile, delim  = ',', col_names = TRUE, col_types = cols())

if(config$indexReads.rc) samples$index1Seq <- as.character(reverseComplement(DNAStringSet(samples$index1Seq)))
if('alignment.removeVectorReadPairs' %in% names(config) &  config$alignment.removeVectorReadPairs == TRUE & !'anchorRead.seqFilter.file' %in% names(samples)){
  stop('A anchorRead.seqFilter.file column must be defined in your sample config file if alignment.removeVectorReadPairs is defined in you config file.')
}

if('anchorRead.seqFilter.file' %in% names(samples)) checkConfigFilePaths(samples$anchorRead.seqFilter.file)
if(! 'subject' %in% names(samples))   samples$subject <- 'subject'
if(! 'replicate' %in% names(samples)) samples$replicate <- 1
if(any(grepl('~|\\|', paste(samples$subject, samples$sample, samples$replicate)))) stop('Error -- tildas (~) are reserved characters and can not be used in the subject, sample, or replicate sample configuration columns.')

samples$uniqueSample <- paste0(samples$subject, '~', samples$sample, '~', samples$replicate)
if(any(duplicated(samples$uniqueSample))) stop('Error -- all subject, sample, replicate id combinations are not unique.')

cluster <- makeCluster(config$demultiplexing.CPUs)
clusterExport(cluster, c('config', 'samples'))

# Quality trim virus reads and break reads.
invisible(parLapply(cluster, 
                    list(c(config$adriftReadsFile,  config$sequence.chunk.size, 'adriftReads',  file.path(config$outputDir, 'seqChunks')),
                         c(config$anchorReadsFile,  config$sequence.chunk.size, 'anchorReads',  file.path(config$outputDir, 'seqChunks'))), 
                 function(x){
                   library(ShortRead)
                   source(file.path(config$softwareDir, 'AAVengeR2.lib.R'))
                   qualTrimReads(x[[1]], x[[2]], x[[3]], x[[4]])
                }))


# Collate trimmed reads.
system(paste('cat', file.path(config$outputDir, 'seqChunks', 'anchor*'), ' > ', file.path(config$outputDir, 'trimmedAnchorReads.fastq')))
system(paste('cat', file.path(config$outputDir, 'seqChunks', 'adrift*'), ' > ', file.path(config$outputDir, 'trimmedAdriftReads.fastq')))
system(paste('rm',  file.path(config$outputDir, 'seqChunks', '*')))


# Convert reads to DNAStringSets and sync reads.
logMsg(config, 'Converting index reads to DNA strings.', config$logFile)
index1Reads <- readFastq(config$index1ReadsFile)
index1Reads <- Reduce('append', parLapply(cluster, split(index1Reads, ntile(1:length(index1Reads), config$demultiplexing.CPUs)), 
               function(x){source(file.path(config$softwareDir, 'AAVengeR2.lib.R')); shortRead2DNAstringSet(x)}))

if(config$correctGolayIndexReads){
  logMsg(config, 'Correcting golay bar code reads.', config$logFile)
  index1Reads <- Reduce('append', parLapply(cluster, split(index1Reads, ntile(1:length(index1Reads), config$demultiplexing.CPUs)), golayCorrection))
}

logMsg(config, 'Converting anchor reads to DNA strings.', config$logFile)
anchorReads <- readFastq(file.path(config$outputDir, 'trimmedAnchorReads.fastq'))
anchorReads <- Reduce('append', parLapply(cluster, split(anchorReads, ntile(1:length(anchorReads), config$demultiplexing.CPUs)), 
                                          function(x){source(file.path(config$softwareDir, 'AAVengeR2.lib.R')); shortRead2DNAstringSet(x)}))

logMsg(config, 'Converting break reads to DNA strings.', config$logFile)
adriftReads <- readFastq(file.path(config$outputDir, 'trimmedAdriftReads.fastq'))
adriftReads <- Reduce('append', parLapply(cluster, split(adriftReads, ntile(1:length(adriftReads), config$demultiplexing.CPUs)), 
                                          function(x){source(file.path(config$softwareDir, 'AAVengeR2.lib.R')); shortRead2DNAstringSet(x)}))

invisible(file.remove(file.path(config$outputDir, 'trimmedAnchorReads.fastq')))
invisible(file.remove(file.path(config$outputDir, 'trimmedAdriftReads.fastq')))


# Slow
logMsg(config, 'Synchronizing trimmed reads.', config$logFile)
reads <- syncReads(index1Reads, anchorReads, adriftReads)
index1Reads <- reads[[1]];  anchorReads  <- reads[[2]];  adriftReads  <- reads[[3]]
rm(reads)
gc()


# Split the trimmed reads into chunks for parallel processing.
logMsg(config, 'Distributing read data into chunks ...', config$logFile)
chunkNum <- 1
d <- tibble(i = ntile(1:length(index1Reads), config$demultiplexing.CPUs), n = 1:length(index1Reads))
invisible(lapply(split(d, d$i), function(x){
  index1Reads <- index1Reads[min(x$n):max(x$n)]
  anchorReads  <- anchorReads[min(x$n):max(x$n)]
  adriftReads  <- adriftReads[min(x$n):max(x$n)]
  save(index1Reads, anchorReads, adriftReads, file = file.path(config$outputDir, 'seqChunks', chunkNum))
  chunkNum <<- chunkNum + 1
}))
logMsg(config, paste0(chunkNum-1, ' data chunks created.'), config$logFile)

# Clean up and free up memory. 
rm(d, chunkNum, index1Reads, anchorReads, adriftReads)
gc()


save(list = ls(all=TRUE), file = file.path(config$outputDir, 'savePoint1.RData'))
logMsg(config, 'Starting sample chunk threads, resetting timer.', config$logFile)
config$startTime <- ymd_hms(format(Sys.time(), "%y-%m-%d %H:%M:%S"))

if(! dir.exists(file.path(config$outputDir, 'logs', 'cutadapt'))) dir.create(file.path(config$outputDir, 'logs', 'cutadapt'))
if(! dir.exists(file.path(config$outputDir, 'tmp', 'cutadapt'))) dir.create(file.path(config$outputDir, 'tmp', 'cutadapt'))


invisible(parLapply(cluster, list.files(file.path(config$outputDir, 'seqChunks'), full.names = TRUE), function(f){
#invisible(lapply(list.files(file.path(config$outputDir, 'seqChunks'), full.names = TRUE), function(f){
  library(ShortRead)
  library(tidyverse)
  source(file.path(config$softwareDir, 'AAVengeR2.lib.R'))
  
  load(f)
  
  # Capture the chunk identifier.
  chunk.n <- unlist(str_match_all(f, '(\\d+)$'))[2]
  logFile <- file.path(config$outputDir, 'logs', paste0('seqChunk_', chunk.n, '.log'))
  
   
  # Loop through samples in sample data file to demultiplex and apply read specific filters.
  invisible(lapply(1:nrow(samples), function(r){
    r <- samples[r,]
    
    # Create barcode demultiplexing vectors.
    v1 <- vcountPattern(r$index1.seq, index1Reads, max.mismatch = config$index1Reads.maxMismatch) > 0
    logMsg(config, paste0('Chunk ', chunk.n, ': (', r$uniqueSample, ') ', sum(v1), ' reads pass barcode filter.'), logFile)
    
    
    # Create break read linker barcode demultiplexing vector.
    v2 <- rep(TRUE, length(adriftReads))
    if('adriftReads.linkerBarcode.maxMismatch' %in% names(config)){
      testSeq <- substr(r$adriftRead.linker.seq, r$adriftRead.linkerBarcode.start, r$adriftRead.linkerBarcode.end)
      v2 <- vcountPattern(testSeq, subseq(adriftReads, r$adriftRead.linkerBarcode.start, r$adriftRead.linkerBarcode.end), max.mismatch = config$adriftReads.linkerBarcode.maxMismatch) > 0
      logMsg(config, paste0('Chunk ', chunk.n, ': (', r$uniqueSample, ') ', sum(v2), ' reads pass linker code filter.'), logFile)
    }
    
    # Test to see if any reads demultiplex to this row of the sample table and then subset reads to this sample.
    i <- base::intersect(which(v1), which(v2))
    if(length(i) == 0){
      logMsg(config, paste0('Chunk ', chunk.n, ': (', r$uniqueSample, ') No reads demultiplexed.'), logFile)
      return()
    } else {
      reads <- syncReads(index1Reads[i], anchorReads[i], adriftReads[i])
      index1Reads <- reads[[1]];  anchorReads  <- reads[[2]]; adriftReads  <- reads[[3]]
    }
    
    # Test the start of anchor reads (static or blast search options only)
    # Not compatible with anchorReads.captureLTRseq.method = 'lentiViralHMM'
    v1 <- rep(TRUE, length(anchorReads))
    if('anchorReads.startTest.maxMismatch' %in% names(config) & config$anchorReads.captureLTRseq.method != 'lentiViralHMM'){
      v1 <- Reduce('|', lapply(unlist(strsplit(r$anchorRead.identification, ';')), function(x){ 
        testSeq <- substr(unlist(strsplit(x, ','))[2], 1,  config$anchorReads.startTest.length)
        vcountPattern(testSeq, subseq(anchorReads, 1, config$anchorReads.startTest.length), max.mismatch = config$anchorReads.startTest.maxMismatch) > 0
      }))
      
      logMsg(config, paste0('Chunk ', chunk.n, ': (', r$uniqueSample, ') ', sum(!v1), ' reads removed by virus read start test filter (', config$anchorReads.startTest.length, ' NTs).'), logFile)
    }
    
    if(! any(v1)){
      logMsg(config, paste0('Chunk ', chunk.n, ': (', r$uniqueSample, ') ', sum(!v1), ' All reads removed by virus read start test filter.'), logFile)
      return()
    }

    
    # Capture the LTR/ITR sequences.
    if(config$anchorReads.captureLTRseq.method == 'blastProvidedTemplates'){
      o <- captureLTRseqs(anchorReads, r$anchorRead.identification)
    } else if(config$anchorReads.captureLTRseq.method == 'lentiViralHMM'){  
      o <- captureLTRseqsLentiHMM(anchorReads, r$anchorRead.identification)
    } else if(config$anchorReads.captureLTRseq.method == 'staticLTRseq'){
      o <- captureStaticLTRseq(anchorReads, r$anchorRead.identification)
    } else {
      stop('Error - No LTR capture method provided.')
    }
    
    # Here we require all ITR / LTR remnants to be at least config$anchorReads.identification.minLength NTs because we will use them as adaptor sequences for cutadapt 
    if(! all(names(o$reads) == o$LTRs$id)) stop('LTRseq capture error.')
    i <- nchar(o$LTRs$LTRseq) >= config$anchorReads.identification.minLength
    
    if(sum(i) == 0){
      logMsg(config, paste0('Chunk ', chunk.n, ': (', r$uniqueSample, ') No anchor reads returned a start sequence >= ', config$anchorReads.identification.minLength, ' NTs'))
      return()
    }
    
    o$reads <- o$reads[i]
    o$LTRs  <- o$LTRs[i,]
    
    
    if(length(o$reads) == 0){
      logMsg(config, paste0('Chunk ', chunk.n, ': (', r$uniqueSample, ') No reads returned from LTR/ITR capture.'), logFile)
      return()
    } else {
      logMsg(config, paste0('Chunk ', chunk.n, ': (', r$uniqueSample, ') ', length(o$reads), ' reads remaining after LTR/ITR capture.'), logFile)
    
      # Here we limit the reads to those with ITRs/LTRs > 12 NT (defined in o)
      anchorReads  <- anchorReads[names(anchorReads) %in% names(o$reads)]
      reads <- syncReads(index1Reads, anchorReads, adriftReads)
      index1Reads <- reads[[1]];  anchorReads  <- reads[[2]]; adriftReads  <- reads[[3]]
    }
    
    # Create the over read sequence from the concensus of the break read common linker and thelast 12 NTs 
    # of the captured virus reads ITR / LTR sequences. 
  
    anchorReadOverReadSeq <- substr(r$adriftRead.linker.seq, nchar(r$adriftRead.linker.seq)-config$anchorReads.identification.minLength, nchar(r$adriftRead.linker.seq))
    anchorReadOverReadSeq <- Biostrings::DNAString(anchorReadOverReadSeq)
    anchorReadOverReadSeq <- as.character(Biostrings::reverseComplement(anchorReadOverReadSeq))
    
    logMsg(config, paste0('Chunk ', chunk.n, ': (', r$uniqueSample, ') ', ' anchor read over read sequence determined to be ', 
                          anchorReadOverReadSeq, ' from parsing provided linker sequence.'), logFile)

    # Trim virus read over-read sequences by looking for the reverse complement of end of the common linker sequence.
    anchorReads <- trimOverReadSeq(anchorReads, anchorReadOverReadSeq, logFile = paste0('seqChunk_', chunk.n, '_', r$uniqueSample, '.anchorReads.trimadapt.log'))
    
    
    # Trim break read over-read sequence by looking for the reverse complement of captured ITR/LTR sequences.
    adriftReads <- trimOverReadSeq2(adriftReads, o, logFile = paste0('seqChunk_', chunk.n, '_', r$uniqueSample, '.adriftReads.trimadapt.log'))
    
    v1 <- width(anchorReads) - nchar(o$LTRs[match(names(anchorReads), o$LTRs$id),]$LTRseq) >= config$trimmedRead.minLength
    v2 <- width(adriftReads) - nchar(r$adriftRead.linker.seq) >= config$trimmedRead.minLength
    
    reads <- syncReads(index1Reads, anchorReads[v1], adriftReads[v2])
    index1Reads <- reads[[1]];  anchorReads  <- reads[[2]]; adriftReads  <- reads[[3]]
    
    if(length(anchorReads) == 0){
      logMsg(config, paste0('Chunk ', chunk.n, ': (', r$uniqueSample, ') No reads remaining after over-read trimming and length check.'), logFile)
      return()
    }
    
    
    # Test for ITR/LTR reads that align to the vector plasmid.
    # This test previously used full adriftReads but here we changed to the ends of virus reads which are 
    # expected to align to genomic sequences.
    
    if(config$alignment.removeVectorReadPairs){
      reads <- subseq(anchorReads, width(anchorReads) - config$trimmedRead.minLength + 1, width(anchorReads))
      vectorReadIDs <- getVectorReadIDs(reads, config, r$anchorRead.seqFilter.file)
    
      if(length(vectorReadIDs) > 0){
        anchorReadsVector <- anchorReads[names(anchorReads) %in% vectorReadIDs]
        adriftReadsVector <- adriftReads[names(adriftReads) %in% vectorReadIDs]
    
        names(anchorReadsVector) <- paste0(names(anchorReadsVector), '|', r$uniqueSample)
        names(adriftReadsVector) <- paste0(names(adriftReadsVector), '|', r$uniqueSample)
        
        writeFasta(anchorReadsVector, file = file.path(config$outputDir, 'tmp', paste0(r$uniqueSample, '.anchorReadsVector.', chunk.n, '.fasta')))
        writeFasta(adriftReadsVector, file = file.path(config$outputDir, 'tmp', paste0(r$uniqueSample, '.adriftReadsVector.', chunk.n, '.fasta')))
      }
      
      logMsg(config, paste0('Chunk ', chunk.n, ': (', r$uniqueSample, ') ', length(vectorReadIDs), ' reads aligned to the vector and will be removed.'), logFile)
      
      # Filter break point reads against read ids returned by vectorReadIDs().
      adriftReads <- adriftReads[! names(adriftReads) %in% vectorReadIDs]
    
      
      # Sync reads. 
      reads <- syncReads(index1Reads, anchorReads, adriftReads)
      index1Reads <- reads[[1]];  anchorReads  <- reads[[2]]; adriftReads  <- reads[[3]]
      
      if(length(anchorReads) == 0){
        logMsg(config, paste0('Chunk ', chunk.n, ': (', r$uniqueSample, ') No reads remaining after vector read removal and read synchronization.'), logFile)
        return()
      } else {
        logMsg(config, paste0('Chunk ', chunk.n, ': (', r$uniqueSample, ') ', length(anchorReads), ' reads remain after removing vector aligned reads and synchronizing.'), logFile)
      }
    }
    
    # Sync the LTR capture objects with the current read list.
    o$reads <- anchorReads
    o$LTRs <- o$LTRs[o$LTRs$id %in% names(o$reads),]
    
    o$reads <- o$reads[match(o$LTRs$id, names(o$reads))]                      # Order the reads to match LTR table.
    o$reads <- subseq(o$reads, start = nchar(o$LTRs$LTRseq)+1)                # Remove matched LTR sequences from the reads.
    o$reads <- o$reads[width(o$reads) >= config$trimmedRead.minLength]        # Remove reads which are now too short.
    o$LTRs  <- subset(o$LTRs, id %in% names(o$reads))                         # Now that we have removed some reads, trim the LTR table.
    o$LTRs$readSeq <- as.character(o$reads[match(o$LTRs$id, names(o$reads))]) # Add the trimmed sequences to the LTR table so that we can latter find additional NTs between LTRs and genomic alignments.
      
    saveRDS(o$LTRs, file = file.path(config$outputDir, 'tmp', paste0(r$uniqueSample, '.LTRseqs.', chunk.n, '.rds')))
       
    # Record reads which were removed.
    ids <- names(anchorReads)[! names(anchorReads) %in% names(o$reads)]
    if(length(ids) > 0) writeLines(ids, file.path(config$outputDir, 'readsRemoved', paste0('LTRcapture.', chunk.n)))
       
    anchorReads <- o$reads
    reads <- syncReads(index1Reads, anchorReads, adriftReads)
    index1Reads <- reads[[1]];  anchorReads  <- reads[[2]]; adriftReads  <- reads[[3]]
      
    if(length(anchorReads) == 0){
      logMsg(config, paste0('Chunk ', chunk.n, ': (', r$uniqueSample, ') No reads remaining after LTR capture and trimming.'), logFile)
      return()
    }

    
    # Capture random ids.
    randomIDs <- data.frame()
    if('breakReadLinkerBarcode.start' %in% names(r) & 'breakReadLinkerBarcode.end' %in% names(r)){
      randomIDs <- as.character(subseq(adriftReads, r$adriftRead.linkerRandomID.start, r$adriftRead.linkerRandomID.end))
      randomIDs <- data.frame(randomSeqID = unname(randomIDs), readID = names(randomIDs))
      logMsg(config, paste0('Chunk ', chunk.n, ': (', r$uniqueSample, ') captured random linker ids.'), logFile)
    }
    
    # Remove leading linker from break point reads.
    adriftReads <- trimLeadingSeq(adriftReads, r$adriftRead.linker.seq)
    adriftReads <- adriftReads[width(adriftReads) >= config$trimmedRead.minLength]
  
    reads <- syncReads(index1Reads, anchorReads, adriftReads)
    index1Reads <- reads[[1]];  anchorReads  <- reads[[2]]; adriftReads  <- reads[[3]]
    
    if(length(anchorReads) == 0){
      logMsg(config, paste0('Chunk ', chunk.n, ': (', r$uniqueSample, ') No reads remaining after second read size filter.'), logFile)
      return()
    } 
    
    # Write out final reads. Add sample names to read IDs.
    names(adriftReads) <- paste0(names(adriftReads), '|', r$uniqueSample)
    names(anchorReads) <- paste0(names(anchorReads), '|', r$uniqueSample)
    
    save(randomIDs, file = file.path(config$outputDir, 'tmp', paste0(r$uniqueSample, '.randomIDs.', chunk.n, '.RData')))
    writeFasta(adriftReads,  file = file.path(config$outputDir, 'tmp', paste0(r$uniqueSample, '.adriftReads.', chunk.n, '.fasta')))
    writeFasta(anchorReads,  file = file.path(config$outputDir, 'tmp', paste0(r$uniqueSample, '.anchorReads.', chunk.n, '.fasta')))
    
    logMsg(config, paste0('Chunk ', chunk.n, ': (', r$uniqueSample, ') completed with ', length(anchorReads), ' reads.'), logFile)
  }))
  
  logMsg(config, paste0('Read data chunk ', chunk.n, ' completed.'), file.path(config$outputDir, 'logs', 'log'))
}))

stopCluster(cluster)



# Collate demultiplexed reads using file name snibets.
collateSampleReads('anchorReads')
collateSampleReads('adriftReads')


# Collate read pairs where the break point read aligns well to the vector for downstream analysis.
if(config$alignment.removeVectorReadPairs){
  logMsg(config, 'Collating read pairs which were removed because the break read aligned to the vector.', config$logFile)
  collateSampleReads('anchorReadsVector')
  collateSampleReads('adriftReadsVector')
}


# Organize read files by reference genome and read source (virusRead or breakRead).
# This will allow samples to be aligned to different genomes defined in the samples table.
df <- tibble(file = list.files(file.path(config$outputDir, 'sampleReads'), pattern = '\\.anchorReads\\.|\\.adriftReads\\.'),
            uniqueSample = unlist(lapply(strsplit(file, '\\.'), function(x) paste0(x[1:(length(x) - 2)], collapse = '.'))),
            source = ifelse(grepl('anchorRead', file), 'anchorReads', 'adriftReads')) %>%
     left_join(select(samples, uniqueSample, refGenome.id), by = 'uniqueSample')

# file                                uniqueSample    source        refGenome.id                                        
# pH19~GTSP2169~1.adriftReads.fasta   pH19~GTSP2169~1  adriftReads  canFam3     
# pH19~GTSP2169~1.anchorReads.fasta   pH19~GTSP2169~1  anchorReads  canFam3     
# pH19~GTSP2170~1.adriftReads.fasta   pH19~GTSP2170~1  adriftReads  canFam3     
# pH19~GTSP2170~1.anchorReads.fasta   pH19~GTSP2170~1  anchorReads  canFam3

# Convert refGenome id to refGenome BLAT database file within AAVenger data directory.
df$refGenome <- file.path(config$softwareDir, 'data', 'blatDBs', paste0(df$refGenome.id, '.2bit'))



# Align the sample FASTA files defined in the previous data frame to their respective genomes. 
# This is done by grouping reads by type (anchorReads / adriftReads) and reference geneome 
# and using parLapply() to disrubte read chunks across a number of CPUs. Unique reads 
# are nested and then unnnested after alignment.


logMsg(config, 'Starting refGenome alignments.', config$logFile)


nSample <- 1

alignments <- bind_rows(lapply(split(df, 1:nrow(df)), function(x){
  
  fasta <- readDNAStringSet(file.path(config$outputDir, 'sampleReads', x$file))
  
  # Here we create a table of unique sequences to align to the reference geneome.
  # The ids associated with each unique sequnce are stored so that they can be 
  # expanded back to match the original data after the alignemnts. 
  
  d <- tibble(id  = names(fasta),
              seq = as.character(fasta)) %>%
       group_by(seq) %>%
        summarise(ids = list(id)) %>%
       ungroup() %>%
       mutate(seqID = paste0('s', 1:n()))
  
  
  fasta <- DNAStringSet(d$seq)
  names(fasta) <- d$seqID
  
  db <- x$refGenome[1]
  cluster <- makeCluster(config$genomAlignment.CPUs)
  clusterExport(cluster, c('config', 'db'), envir = environment())
  
  logMsg(config, paste0('Aligning seq file ', nSample, '/', nrow(df), ' with ', length(fasta), ' reads using ', config$genomAlignment.CPUs, ' cores.'), config$logFile)
  nSample <<- nSample + 1
  
  r <- bind_rows(parLapply(cluster, split(fasta, dplyr::ntile(1:length(fasta), config$genomAlignment.CPUs)), function(y){
  #r <- bind_rows(lapply(split(fasta, dplyr::ntile(1:length(fasta), config$genomAlignment.CPUs)), function(y){
      
         library(ShortRead)
         library(tidyverse)
         source(file.path(config$softwareDir, 'AAVengeR2.lib.R'))
         
         a <- alignReads.BLAT(y, db)
         
         # Apply generic alignment filters.
         if(nrow(a) > 0){
           a <- dplyr::filter(a, alignmentPercentID >= config$alignment.refGenome.minPercentID,
                                 tNumInsert  <= 1, 
                                 qNumInsert  <= 1,
                                 tBaseInsert <= 2,
                                 qBaseInsert <= 2)
         } 
         
         a
       }))
  
  stopCluster(cluster)
  
  if(nrow(r) == 0) return(tibble())
  
  d <- left_join(select(d, -seq), r,  by = c("seqID" = "qName"))
  d <- unnest(d, ids)
  
  d$source <- x$source[1]
  select(d, -seqID)
}))

logMsg(config, 'RefGenome alignments completed.', config$logFile)

save(list = ls(all=TRUE), file = file.path(config$outputDir, 'savePoint2.RData'))



# Create position ids.
alignments$posid <- paste0(alignments$tName, alignments$strand, ifelse(alignments$strand == '+', alignments$tStart, alignments$tEnd))


# Segregate alignments.
anchorReads.aln <- alignments[which(alignments$source == 'anchorReads'),]
adriftReads.aln <- alignments[which(alignments$source == 'adriftReads'),]


# Adrift point reads have static 5' adapters which are removed with cutadapt.
# Only consider reads which immediately align to the genome.
adriftReads.aln <- filter(adriftReads.aln, qStart <= 3, matches >= config$trimmedRead.minLength)


# Virus reads have had over-read sequences removed.
# Only consider virus reads where the end of the read aligns to the genome.
anchorReads.aln <- filter(anchorReads.aln, (qSize - qEnd) <= 3, matches >= config$trimmedRead.minLength)


alignments <- bind_rows(anchorReads.aln, adriftReads.aln)
rm(anchorReads.aln, adriftReads.aln)


# Rename alignment column headers by appending read sources.
alignments <- lapply(split(alignments, alignments$source), function(x){
    names(x) <- paste0(names(x), '.', x$source[1])
    x
})


# Combine read alignments into rows and extract sample names from read ids. 

v <- select(alignments[["anchorReads"]][ alignments[["anchorReads"]]$ids.anchorReads %in% alignments[["adriftReads"]]$ids.adriftReads,],
            ids.anchorReads, qStart.anchorReads, strand.anchorReads, tName.anchorReads, tStart.anchorReads, tEnd.anchorReads)
v$ids.anchorReads <- factor(v$ids.anchorReads )

b <- select(alignments[["adriftReads"]][ alignments[["adriftReads"]]$ids.adriftReads %in% alignments[["anchorReads"]]$ids.anchorReads,],
            ids.adriftReads, strand.adriftReads, tName.adriftReads, tStart.adriftReads, tEnd.adriftReads)
b$ids.adriftReads <- factor(b$ids.adriftReads)

rm(alignments)
logMsg(config, 'Starting to create fragments from alignment data.', config$logFile)


# Here we build a table of potential fragments by matching up virus and break point reads
# and filtering with basic sanity checks. This is done via a left join of the break reads 
# to the virus reads followed by filtering. The virus reads are split into a number 
# of chunks and the correponding break reads are joined. The reads broken into a couple 
# of chunks because a single join may break R's internal able row number limit. We 
# chose three here because too many chunks would require the copying of all the break
# reads to multiple computaion cores.

cluster <- makeCluster(3)
clusterExport(cluster, c('config', 'b'))

v$s <- ntile(1:nrow(v), 3)

frags <- bind_rows(parLapply(cluster, split(v, v$s), function(x){
         library(dplyr)
         library(stringr)
  
         frags <- left_join(x, b, by = c('ids.anchorReads' = 'ids.adriftReads')) %>% tidyr::drop_na()
       
         i <- which(frags$tName.anchorReads != frags$tName.adriftReads)
         if(length(i) > 0) frags <- frags[-i,]
         if(nrow(frags) == 0) return(data.frame())
         
         i <- which(frags$strand.anchorReads == frags$strand.adriftReads)
         if(length(i) > 0) frags <- frags[-i,]
         if(nrow(frags) == 0) return(data.frame())
         
         mutate(frags, 
                fragStart = ifelse(strand.anchorReads == '+', tStart.anchorReads + 1, tStart.adriftReads + 1),
                fragEnd   = ifelse(strand.anchorReads == '+', tEnd.adriftReads + 1,   tEnd.anchorReads + 1),
                fragTest  = ifelse(strand.anchorReads == '+', tStart.anchorReads < tEnd.adriftReads, tStart.adriftReads < tEnd.anchorReads),  
                fragWidth = (fragEnd - fragStart) + 1) %>%
         filter(fragTest == TRUE, 
                fragWidth <= config$fragments.maxLength,
                fragWidth >= config$fragments.minLength) %>%
         mutate(uniqueSample = unlist(lapply(str_split(ids.anchorReads, '\\|'), '[[', 2 )),
                readID = sub('\\|.+$', '', ids.anchorReads)) 
       }))

stopCluster(cluster)

logMsg(config, 'Initial fragment generation completed.', config$logFile)

frags <- left_join(frags, dplyr::select(samples, refGenome.id, uniqueSample), by = 'uniqueSample')
frags <- mutate(frags, fragID = paste0(uniqueSample, ';', tName.anchorReads, ';', strand.anchorReads, ';', fragStart, ';', fragEnd))


# Remove read frags which share a random linker id and spand more than one sample.
if(config$removeReadFragsWithSameRandomID) frags <- removeReadFragsWithSameRandomID(frags, config)


# Convert read fragments into grouped fragments.
frags <- group_by(frags, fragID) %>%
  summarise(reads = n(), refGenome.id = refGenome.id[1], virusReadAlignmentStart = qStart.anchorReads[1], readIDs = list(readID)) %>%
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


posIDclusterTable <- 
  dplyr::group_by(unnest(frags, readID.list), subject, readID.list) %>%
  dplyr::summarise(nPositions = n_distinct(posid), posids = list(posid)) %>%
  dplyr::ungroup() %>%
  dplyr::filter(nPositions > 1) 

# posIDclusterTable -- for each read, provide the multiple posids to which it maps.
#
# subject     readID.list                                       nPositions   posids (list of posids)    
#   1 pH19    M03249:365:000000000-C3CH4:1:1101:11087:17666          2      <chr [2]> 
#   2 pH19    M03249:365:000000000-C3CH4:1:1101:11121:16397         13      <chr [13]>
#   3 pH19    M03249:365:000000000-C3CH4:1:1101:11589:17142         13      <chr [13]>

posIDclusters <- do.call(rbind, lapply(split(posIDclusterTable, posIDclusterTable$subject), function(s){
                   nodes <- unique(unlist(s$posids))
                   edges <- bind_rows(lapply(split(s, 1:nrow(s)), function(x){ 
                              data.frame(t(combn(unlist(x$posids), 2)))
                            }))

                   network = graph_from_data_frame(edges, vertices = nodes, directed = F)
                   data.frame(subject = s$subject[1], clusters = I(list(lapply(decompose(network), function(x) V(x)$name))))
                 }))

# Structure of posIDclusters
# subject     clusters
# <chr>     <list of posids clusters>
# pH19      c("chr11....
# pHO2      c("chr20....
# pJ60      c("chr2+....
#
# Clusters for pH19 (first row from posIDclusters)
#
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



# Examine the LTR/ITR sequence remnants of each read from each fragment and determine how similiar they 
# are to one another. For each fragment, the most similiar 95% of the reads are considered and representative 
# remnant sequences are determined.

# ltrRepSeq is the representative ITR/LTR sequence for a fragment based on examination of its reads.
# ltrRepSeq2 is the  representative ITR/LTR sequence with addition NTs found between the recognizable 
#   ITR/LTR sequence and the start of the alignment to the genome..

cluster <- makeCluster(config$demultiplexing.CPUs)
clusterExport(cluster, c('config'))

#frags$s <- ntile(1:nrow(frags), config$demultiplexing.CPUs)

# Distribute such that frags with several reads are equally spread across groups
frags <- arrange(frags, desc(reads))
frags$s <- rep(1:config$demultiplexing.CPUs, ceiling(nrow(frags)/config$demultiplexing.CPUs))[1:nrow(frags)]

# JKE
#frags <- bind_rows(parLapply(cluster, split(frags, frags$s), function(a){

# Error: `by` can't contain join column `qName` which is missing from RHS


frags <- bind_rows(lapply(split(frags, frags$s), function(a){  
              library(dplyr)
              source(file.path(config$softwareDir, 'AAVengeR2.lib.R'))
    
              logFile <-  file.path(config$outputDir, 'logs', paste0('LTRrep_', a$s[1], '.log'))
             
              # Assemble LTR.table here within the parLapply loop because it is more efficient to
              # assemble it here rather than passing a copy to each computation node.
              LTR.table <- bind_rows(lapply(list.files(file.path(config$outputDir, 'tmp'), pattern = 'LTRseqs', full.names = TRUE), readRDS))
    
              # Examine each fragment individualy and assess the captured ITRs/LTRs defined in LTR.table.
              bind_rows(lapply(1:nrow(a), function(x){
                frag <- a[x,]
                
                logMsg(config, paste0('Processing fragment ', x, '/', nrow(a), ' with ', frag$reads, ' reads.'), logFile)
               
                # Retrieve LTR sequences for the reads that define this fragment.
                ltrs <- LTR.table[match(unlist(frag$readID.list), LTR.table$id),]
                
                # Identify the most common ITR/LTR sequece in the reads supporting this fragment.
                # The returned maxReadPercentDiff metric is the max edit distance between the 
                # returned representative divided by the number of letters in the representative sequence.
         
                r <- representativeSeq(ltrs$LTRseq, logFile = logFile)
                
                frag$maxReadPercentDiff <- r[[1]]
                frag$ltrRepSeq <- r[[2]]
                
                # If the virus read alignments do not start on 1 then there are additional NTs 
                # between the recognized LTR and the following genomic juncture. This may be 
                # the result of alignment error or genuine extensions of the ITR/LTR.
                # substr(string, 1, 0) with return an empty string.
                additionalLTRnts <- substr(ltrs$readSeq, 1, frag$virusReadAlignmentStart)
                
                if(all(nchar(additionalLTRnts) == 0)){
                  frag$maxReadPercentDiff2 <- frag$maxReadPercentDiff 
                  frag$ltrRepSeq2 <- frag$ltrRepSeq 
                }else{
                  logMsg(config, paste0('Reprocessing fragment with additional NTs', x, '/', nrow(a), ' with ', frag$reads, ' reads.'), logFile)
                  r <- representativeSeq(paste0(ltrs$LTRseq, additionalLTRnts), logFile = logFile)
                  frag$maxReadPercentDiff2 <- r[[1]]
                  frag$ltrRepSeq2 <- r[[2]]
                }
                
                frag
              }))
        }))
  
stopCluster(cluster)
  

# JKE1

save(list = ls(all=TRUE), file = file.path(config$outputDir, 'savePoint_dev.RData'))
  
# Remove fragments where there is not a clear concensus remnant. 
frags <- subset(frags, frags$maxReadPercentDiff2 <= config$fragmentProcessing.rep.maxDifference)
  

# Save the fragments for downstream analyses.
saveRDS(frags, file = file.path(config$outputDir, 'readFrags.rds'))

frags$s <- paste0(frags$subject, '_', frags$sample, '_', frags$reads)
createFragUCSCTrack(frags, title = 'AAVengeR_fragments', outputFile = file.path(config$outputDir, 'fragments.ucsc'), label = 's')


# Collapse fragments to individual sites. 
# Fragment ITR/LTR sequences will be evaulated so that fragments with dissimiliar ITR/LTR 
# sequences will not be combined.

sites <- bind_rows(lapply(split(frags, paste(frags$subject, frags$sample, frags$posid, frags$posIDcluster)), function(x){
    if(nrow(x) > 1){
      i <- rep(TRUE, nrow(x))
      r2 <- representativeSeq(x$ltrRepSeq2, logFile = config$logFile)
      
      if(r2[[1]] > config$fragmentProcessing.rep.maxDifference){
        # There is a conflict, one or more fragments have a markedly different LTR sequence representative then the other fragments.
        # Attempt to salvage this site by retaining the majority of fragments with similiar ITR/LTR sequences.
        
        i <- as.vector(stringdist::stringdistmatrix(x$ltrRepSeq2, r2[[2]]) / nchar(x$ltrRepSeq2) < config$fragmentProcessing.rep.maxDifference)
        if(sum(i)/nrow(x) >= config$fragmentProcessing.assemblyConflictResolution){
          x <- x[i,]
        } else {
          return(tibble())
        }
      }
      
      r <- representativeSeq(x$ltrRepSeq, logFile = config$logFile)
      
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

  

# Collapse multihits.
# Multihit sites are represented by being associated with pseudo chromosome chrUn and strand *.
# Multihit posids have the form of chrUn*x(y) where x is the numeric cluster id and y are the
# number of posids in the cluster. Estimated abundance value of a multi hit is the maximum estimated 
# abundance of any site within the cluster and the reads count is the sum of all reads of all sites 
# within a cluster.

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

d <- dplyr::mutate(sites, start = position, end = position, siteLabel = paste0(subject, '_', sample, '_', estAbund)) %>%
     dplyr::filter(seqnames != 'chrUn')
createIntUCSCTrack(d, siteLabel = 'siteLabel', outputFile = file.path(config$outputDir, 'sites.ucsc'))

system(paste0('cat ', file.path(config$outputDir, 'sites.ucsc'), ' ', file.path(config$outputDir, 'fragments.ucsc'), ' > ',file.path(config$outputDir, 'AAVengeR.ucsc')))

logMsg(config, 'done.', config$logFile)
