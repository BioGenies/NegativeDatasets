library(dplyr)
library(biogram)

if(Sys.info()[["nodename"]] == "kasia-MACH-WX9") {
  data_path <- "/media/kasia/Data/Dropbox/Projekty/BioNgramProjects/NegativeDatasets/"
}

source("functions/filter_sequences.R")

uniprot_seqs <- read_fasta(paste0(data_path, "Data/uniprot-reviewed yes.fasta"))
uniprot_tab <- read.delim(paste0(data_path, "Data/uniprot-reviewed yes.tab"))
new_seq_names <- sapply(names(uniprot_seqs), function(x) strsplit(x, "|", fixed = TRUE)[[1]][2])
names(uniprot_seqs) <- new_seq_names

# Example - 2016 MLAMP (10.1093/bioinformatics/btw560)
mlamp_dataset <- uniprot_seqs %>% 
  filter_by_lengths(min_len = 5, max_len = 100) %>% 
  filter_by_annotations(uniprot_data = uniprot_tab, 
                        keywords_vec = c("Antimicrobial", "Antibiotic", "Fungicide", "Defensin")) %>% 
  filter_nonstandard_aa() %>% 
  filter_with_cdhit(0.4)
