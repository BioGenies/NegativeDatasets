read_uniprot_seqs <- function(seq_file) {
  seqs <- read_fasta(seq_file)
  new_seq_names <- sapply(names(seqs), function(x) strsplit(x, "|", fixed = TRUE)[[1]][2])
  names(seqs) <- new_seq_names
  seqs
}