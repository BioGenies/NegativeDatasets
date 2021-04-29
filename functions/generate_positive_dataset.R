generate_positive_dataset <- function(sequence_file) {
  dbaasp <- read.csv("data/dbaasp.csv") %>% 
    filter(Gram. | Gram..1) 
  dbaasp_seqs <- dbaasp[["sequence"]] %>% 
    strsplit("") %>% 
    setNames(paste0("DBAASP_", dbaasp[["dbaasp_id"]], "_AMP=1")) %>% 
    filter_nonstandard_aa() %>% 
    filter_by_lengths(5, max(lengths(.)))
  positive_dataset <- filter_with_cdhit(dbaasp_seqs, 0.9, l = 4)
  write_fasta(positive_dataset, "./data/positive_dataset.fa")
  positive_dataset
}