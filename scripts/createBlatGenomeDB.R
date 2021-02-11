library(BSgenome)

get_reference_genome <- function(reference_genome, type = "full"){
  stopifnot(type %in% c("full", "minimal"))
  pattern <- paste0("\\.", reference_genome, "$")
  
  match_index <- which(grepl(pattern, BSgenome::installed.genomes()))
  
  if (length(match_index) != 1) stop(paste("Cannot find unique genome for", reference_genome))
  
  BS_genome_full_name <- BSgenome::installed.genomes()[match_index]
  
  library(BS_genome_full_name, character.only = TRUE)
  
  gen <- get(BS_genome_full_name)
  
  if (type == "minimal") {
    min_seqs <- grep("_", seqnames(gen), fixed=TRUE, invert=TRUE, value=TRUE)
    gen@user_seqnames <- setNames(min_seqs, min_seqs)
    gen@seqinfo <- gen@seqinfo[min_seqs]
  }
  
  return(gen)
}


export(get_reference_genome('mm9', type = 'minimal'), '../data/blatDBs/mm9.2bit', 'TwoBit')
export(get_reference_genome('hg38', type = 'minimal'), '../data/blatDBs/hg38.2bit', 'TwoBit')
