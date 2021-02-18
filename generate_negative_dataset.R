library(biogram)
library(dplyr)
library(pbapply)

if(Sys.info()[["nodename"]] == "kasia-MACH-WX9") {
  data_path <- "/media/kasia/Data/Dropbox/Projekty/BioNgramProjects/NegativeDatasets/"
}

source("functions/filter_sequences.R")
source("functions/generate_negative_datasets.R")

sequences <- read_fasta(paste0(data_path, "Data/uniprot-reviewed yes.fasta"))
new_seq_names <- sapply(names(sequences), function(x) strsplit(x, "|", fixed = TRUE)[[1]][2])
names(sequences) <- new_seq_names
uniprot_data <- read.delim(paste0(data_path, "Data/uniprot-reviewed yes.tab"))
positive_dataset <- read_fasta("./data/positive_dataset.fa")

exemplary_negative_dataset <- generate_negative_dataset_2(sequences, uniprot_data, positive_dataset) %>% 
  setNames(paste0(names(.), "_AMP=0"))
