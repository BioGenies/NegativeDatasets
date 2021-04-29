get_putative_amps <- function(data_path, keywords = c("Antimicrobial", "Antibacterial", "Antiviral", "Fungicide",
                                                      "defense", "Defensin", "Bacteriocin")) {
  seqs <- read_fasta(paste0(data_path, "Data/uniprot-reviewed yes.fasta")) 
  annot <- read.delim(paste0(data_path, "Data/uniprot-reviewed yes.tab"))
  new_seq_names <- sapply(names(seqs), function(x) strsplit(x, "|", fixed = TRUE)[[1]][2])
  names(seqs) <- new_seq_names
  
  putative_amps <- filter_by_annotations(seqs, annot, keywords_vec = keywords, exclude = FALSE) 
  save(putative_amps, file = "data/putative_AMPs.rda")
  putative_amps
}
