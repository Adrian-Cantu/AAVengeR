# AAVengeR

The genomic fragments analyzed by the AAVengeR pipeline begin from fixed positions within ITR / LTR sequences and end with variable positions created when genomic DNA samples are sheared with sonication. Sequencing reads originating from the fixed ITR/LTR positions are labeled as anchored while reads originating from the variable sheared positions are labeled as adrift.   
  
AAVengeR required two configurations files. The first file, the one and only argument provided to the software defines general parameters as well as points to a second sample specific configuration file. Please see the example configuration files in the provided test data set, data/testData. The general configuration file is a YAML formatted file and includes these parameters:
  
## General configuration file.  
  
softwareDir: (Required) Path to the AAVengeR software. eg. /home/everett/projects/AAVengeR  
outputDir: (Required) Path to output directory. This directory is created when the software runs. The software will not continue if the output directory exists when it starts in order not to overwrite previous results. eg.  data/testData/output  
demultiplexing.CPUs: (Required) Number of CPUs to use when demultiplexing samples. eg. 15  
genomAlignment.CPUs: (Required) Number of CPUs to use when aligning reads to reference genomes: eg. 15  
  
adriftReadsFile: (Required) Path to adrift reads. eg. data/testData/R1.test.fastq.gz  
anchorReadsFile: (Required)  Path to anchored reads. eg. data/testData/R2.test.fastq.gz  
index1ReadsFile: (Required)  Path to index reads. eg. data/testData/I1.test.fastq.gz  
sampleConfigFile: (Required) path to sample configuration file (see below). eg. data/testData/sampleConfig.csv  
  
anchorReads.captureLTRseq.method: (Required) Method to capture LTR / ITR sequence. Options include : blastProvidedTemplates, lentiViralHMM, and staticLTRseq.  
  blastProvidedTemplates: list of labeled possible LTRs / ITRs. (see sample configuration file details).   
  lentiViralHMM: name of HMM defined in data/hmms   
  staticLTRseq: sequence of expected LTR / ITR sequence (still under development)  
   
sequence.chunk.size: (Required) Number of reads to pass to each demultiplexing CPU per cycle, eg. 300000  
  
alignment.chunk.size: (Required) Number of reads to pass to each alignment CPU per cycle: eg. 15000  
  
sequence.qualTrim.code: (Required)  FASTQ quality code beneath which sequences will be trimmed using a sliding window approach. eg.'5'  
  
sequence.qualTrim.minLength: (Required) Minimum length of reads allow to proceed after trimming low quality bases. eg. 75  
  
index1Reads.maxMismatch: (Optional) If this parameter is provided, index reads exceeding the number of mismatches defined here will not be used for demultiplexing. eg. 0  

adriftReads.linkerBarcode.maxMismatch: (Optional) If this parameter is provided along with the sample configuration file parameters to extract unique adrift linker sequence , unique linker sequences exceeding the number of mismatches defined here will not be used for demultiplexing. eg. 1  
  
indexReads.rc: (Optional) Boolean controlling if the revere compliment of index reads should be used. eg. FALSE  
  
correctGolayIndexReads: (Optional) Boolean controlling if index reads employ a Golay correction systems and correction should be attempted. eg. TRUE  
  
fragments.maxLength: (Required) Maximum length (NT) of genomic fragments to consider. eg.2500  
fragments.minLength: (Required) Minimum length (NT) of genomic fragments to consider. eg. 25  
trimmedRead.minLength: (Required) Minimum length of genomic NTs required to attempt an alignment to the reference genome. eg. 25    

fragmentProcessing.rep.maxReads (Required):  When combining aligned fragments into groups supporting integration events, large groups of fragments can overload the aligner. Here we define the maximum number of fragment leader sequences to align before we randomly sample. eg. 5000     
fragmentProcessing.rep.maxDifference: (Required) Documentation pending. eg. 0.15  
fragmentProcessing.assemblyConflictResolution: (Required) Documentation pending. eg. 0.66  
  
anchorReads.startTest.length: (Optional) Test anchor reads, from 1 to n (this parameter), to see if they match one of the possible sequences defined in the sample configuration file (anchorRead.identification ). Reads with test mismatched  exceeding the anchorReads.startTest.maxMismatch limit are excluded. eg. 10  
anchorReads.startTest.maxMismatch: (Optional) Defining the maximum number of allowed mismatches for the anchor rads star test here enables the test: 0  
  
The required settings below define how LTR sequences are captured via blast or HMM alignment.  
  Blast approach, anchorReads.captureLTRseq.method = blastProvidedTemplates:  
  anchorReads.identification.maxAlignmentStart: 3  
  anchorReads.identification.minPercentSeqID: 90  
  anchorReads.identification.maxGapOpen: 1  
  anchorReads.identification.repSeqEditDistMaxPercentage: 5  
 
  HMM approach, anchorReads.captureLTRseq.method = lentiViralHMM:  
  anchorReads.identification.HMMmaxStartPos: 3  
  anchorReads.identification.HMMminEval: 0.0001  
   
 anchorReads.identification.minLength: (Required) Minimum length (NT) of captured LTR / ITR required to consider an anchor read. eg. 10  

### Alignment parameters.  
alignment.removeVectorReadPairs: TRUE  
alignment.removeMultiHitReadPairs: TRUE  
alignment.seqFilter.minPercentID: 95  
alignment.seqFilter.minPercentQueryCoverage: 75  
alignment.refGenome.minPercentID: 98  
alignment.refGenom.minPercentQueryCoverage: 90  


### Fragment standardization.  
removeReadFragsWithSameRandomID: TRUE  
standardizeSitesBy: subject  
standardizeBreakPointsBy: replicate  
  
The parameters below point to the primary executables of local third party software installations.  
All software installations are required.  
  
command.blat:        /home/everett/ext/blat  
command.bwa:         /home/everett/ext/bwa  
command.makeblastdb: /home/everett/ext/blast+/bin/makeblastdb  
command.blastn:      /home/everett/ext/blast+/bin/blastn  
command.samtools:    /home/everett/ext/samtools/bin/samtools  
command.cutadapt3:   /usr/bin/cutadapt3  
command.python2:     /usr/bin/python2  
command.muscle:      /home/everett/ext/muscle  
  
 
## Sample configuration file.  
The sample configuration file is a comma delimited file including these data fields.   

subject: Required. Alphanumeric subject identifier.  
sample: Required. Alphanumeric sample identifier.  
replicate: Required. Integer replicate identifier.   
Subject, sample, replicate combinations must be unique.  	
adriftRead.linker.seq: Required. Adrift read linker sequence as encountered by the sequencer. This linker is expected to appear at the beginning of adrift reads and may contain a random ID sequence defined by a series of Ns, eg.  TACGCTATGTCGACCGTGACNNNNNNNNNNNNCTCCGCTTAAGGGACT.   
  
adriftRead.linkerBarcode.start: Optional. If a segment of the adrift liker sequence is sample specific, it can be extracted and used to demultiplex samples along with the barcode sequences. This parameter defines the start of the unique linker sequence.  
  
adriftRead.linkerBarcode.end: Optional. This parameter defines the end of the unique linker sequence.  
  
adriftRead.linkerRandomID.start: Optional. The adrift linker may contain a random sequence which is useful for identifying sample cross-contamination as well as PCR recombination events. This parameter defines the start of the random sequence id.  
adriftRead.linkerRandomID.end: Optional. This parameter defines the end of the random sequence.  
  
index1.seq: Required. Barcode sequence used for demultiplexing. Future versions of the software will support a second index sequence. 	 
anchorRead.seqFilter.file: Optional. FASTA format file against anchor reads are aligned and subsequently removed if alignments metrics exceed those defined in the main configuration file.  
anchorRead.overRead.seq: Required. Short genomic fragments may be completely sequenced resulting in anchor reads which end with its partner adrift reads linker sequence. This parameter is the reverse compliment of the end of the adrift read linker which is used to trim fragment over-reading.   
refGenome.id: Required. Genome identifier such as hg38 or canFam3 for the reference genome against which sequencing reads will be aligned. These identifiers refer to 2bit formatted genomes in the data folder, eg. hg38 -> data/blatDBs/hg38.2bit  

anchorRead.identification:  Required. The format of this parameter is dependent on the anchorReads.captureLTRseq.method parameter defined in the main configuration file.   
  
blastProvidedTemplates option:  This parameter defines the possible sequences one could encounter when sequencing out of LTR / ITR sequences. The parameter is a semicolon delimited list of possibilities where each possibility is defined with a label and expected sequence, eg.   
Label 1,ATGCTAGCTAACGT; Label 2,GCTAGCTAGCTA; Label 3,CGTACAGTCATAC  
  
The beginning of anchor reads are blasted against each possibility and the longest alignment meeting the minimal alignment scores defined in the main configuration file is retained and used to define the recognizable LTR / ITR sequences. Possible sequence may be subsequences of other possible sequence. This parameter needs to be quoted because it may contains commas which would conflict with the commas used to delimit other columns.  
