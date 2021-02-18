library(biogram)
library(dplyr)

source("functions/filter_sequences.R")

dbaasp <- read.csv("data/dbaasp.csv") %>% 
  filter(Gram. | Gram..1) 
  
dbaasp_seqs <- dbaasp[["sequence"]] %>% 
  strsplit("") %>% 
  setNames(paste0("DBAASP_", dbaasp[["dbaasp_id"]], "_AMP=1")) %>% 
  filter_nonstandard_aa() %>% 
  filter_by_lengths(10, max(lengths(.)))
  

lapply(c(0.5, 0.7, 0.75, 0.80, 0.85, 0.9, 0.95), function(ith_cutoff) {
  data.frame(cutoff = ith_cutoff,
             sequences = length(filter_with_cdhit(dbaasp_seqs, ith_cutoff)))
}) %>% bind_rows()


positive_dataset <- filter_with_cdhit(dbaasp_seqs, 0.95)
write_fasta(positive_dataset, "./data/positive_dataset.fa")
