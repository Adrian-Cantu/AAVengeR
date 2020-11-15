library(RMySQL)
library(dplyr)
library(stringr)

f <- '201004_M03249_0101_000000000-JBCKB.csv'
outputFile <- 'wistarSampleConfig.csv'

refGenomeName <- 'hg38'
refGenomeBLATdb <- '/home/everett/data/sequenceDatabases/BLAT/hg38/hg38.2bit'
vectorSeqFile <- '/home/everett/data/sequenceDatabases/FASTA/HIV-1.ff'
u5_hmm <- '/home/everett/projects/AAVengeR/data/u5_100.hmm'
u3_hmm <- '/home/everett/projects/AAVengeR/data/u3_100_rc.hmm'
u5_breakOverRead <- 'TGCTAGAGATTT'
u3_breakOverRead <- 'TGGATGGAAGCA'

# Read in sample data & subset to the subjects in this report group.
invisible(sapply(dbListConnections(MySQL()), dbDisconnect))
dbConn  <- dbConnect(MySQL(), group='specimen_management')
samples <- dbGetQuery(dbConn, 'select * from gtsp')
dbDisconnect(dbConn)

d <- read.csv(f, header = TRUE)

o <- str_split(d$SampleName, '[\\-_]')

r <- tibble(subject = samples[match(unlist(lapply(o, '[[', 1)), samples$SpecimenAccNum),]$Patient,
                sample = paste0(unlist(lapply(o, '[[', 1)), '_', unlist(lapply(o, '[[', 3))),
                replicate = unlist(lapply(o, '[[', 2)),
                breakReadLinkerBarcode.start = 1,
                breakReadLinkerBarcode.end = nchar(str_match(d$linker, '[^N]+')),
                breakReadLinkerSeq = d$linker,
                index1Seq = d$barcode,
                virusLTRseq = ifelse(grepl('u5', sample), u5_hmm, u3_hmm),
                refGenomeBLATdb = refGenomeBLATdb,
                virusRead.overReadSeq = 'AGTCCCTTAAGCGGAG',
                breakRead.overReadSeq = ifelse(grepl('u5', sample), u5_breakOverRead, u3_breakOverRead),
                breakReadLinkerRandomID.start = breakReadLinkerBarcode.end + 1,
                breakReadLinkerRandomID.end = breakReadLinkerBarcode.end + 13,
                refGenomeName = refGenomeName,
                vectorSeqFile = vectorSeqFile)
                
write.csv(r, file = outputFile, quote = FALSE)
                