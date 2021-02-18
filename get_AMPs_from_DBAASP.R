library(tidysq)
library(dplyr)

dbaasp <- read.csv("./data/dbaasp.csv")

filter(dbaasp, Gram. | Gram..1) %>% 
  mutate(name = paste0("dbaasp", dbaasp_id)) %>% 
  select(name, sequence) %>% 
  filter(!duplicated(sequence)) %>% 
  na.omit %>%
  as_tibble() %>% 
  mutate(sequence = sq(sequence, alphabet = "ami_bsc")) %>% 
  mutate(sequence = remove_na(sequence)) %>% 
  filter(get_sq_lengths(sequence) > 0) %>% {
    write_fasta(x = .[["sequence"]], name = .[["name"]], file = "./data/dbaasp.fasta")
  }
