library(dplyr)
library(biogram)
library(pbapply)

if(Sys.info()[["nodename"]] == "kasia-MACH-WX9") {
  data_path <- "/media/kasia/Data/Dropbox/Projekty/BioNgramProjects/NegativeDatasets/"
}

if(Sys.info()[["nodename"]] == "ryzen") {
  data_path <- "~/Dropbox/Projekty/BioNgramProjects/NegativeDatasets/"
}


source("functions/filter_sequences.R")
source("functions/generate_holdouts.R")
source("functions/generate_negative_datasets.R")
load("./data/putative_AMPs.rda")

uniprot_seqs <- read_fasta(paste0(data_path, "Data/uniprot-reviewed yes.fasta"))
uniprot_tab <- read.delim(paste0(data_path, "Data/uniprot-reviewed yes.tab"))
new_seq_names <- sapply(names(uniprot_seqs), function(x) strsplit(x, "|", fixed = TRUE)[[1]][2])
names(uniprot_seqs) <- new_seq_names

positive_dataset <- read_fasta("./data/positive_dataset.fa")
holdouts <- generate_holdout_groups(positive_dataset)
positive_traintest <- positive_dataset[unlist(sapply(holdouts, function(i) i[["traintest"]]), use.names = FALSE)]
positive_benchmark <- positive_dataset[unlist(sapply(holdouts, function(i) i[["benchmark"]]), use.names = FALSE)]


for (i in 1:5) {
  # Sampling algorithm that need only sequences and uniprot data
  s1 <- generate_negative_dataset_1(sequences = uniprot_seqs,
                                    uniprot_data = uniprot_tab)
  write_fasta(c(positive_traintest, s1), paste0(data_path, "Datasets/Training_method1_rep", i, ".fasta"))
  
  # Sampling algorithms that need only sequences and positive dataset
  benchmark_2_3 <- lapply(paste0("generate_negative_dataset_", 2:3), function(ith_fun) {
    s <- match.fun(ith_fun)(sequences = uniprot_seqs, 
                            positive_dataset = positive_traintest,
                            n_threads = 12)
    write_fasta(c(positive_traintest, s), paste0(data_path, "Datasets/Training_method", last(strsplit(ith_fun, "_")[[1]]), "_rep", i, ".fasta"))
    match.fun(ith_fun)(sequences = uniprot_seqs, 
                       positive_dataset = positive_traintest)
  })
  
  # Sampling algorithms that need sequences, uniprot data and positive dataset
  benchmark_4_11 <- lapply(paste0("generate_negative_dataset_", 4:11), function(ith_fun) {
    s <- match.fun(ith_fun)(sequences = uniprot_seqs, 
                            uniprot_data = uniprot_tab,
                            positive_dataset = positive_traintest)
    write_fasta(c(positive_traintest, s), paste0(data_path, "Datasets/Training_method", last(strsplit(ith_fun, "_")[[1]]), "_rep", i, ".fasta"))
  })
  
  # Sampling algorithm that needs sequences, uniprot data, positive dataset and putative AMPs
  s12 <- generate_negative_dataset_12(sequences = uniprot_seqs,
                                      uniprot_data = uniprot_tab, 
                                      positive_dataset = positive_traintest,
                                      potential_AMPs = putative_amps)
  write_fasta(c(positive_traintest, s12), paste0(data_path, "Datasets/Training_method12_rep", i, ".fasta"))
  
  benchmark_1 <- generate_negative_dataset_1(sequences = uniprot_seqs,
                                             uniprot_data = uniprot_tab) 
  benchmark_12 <- generate_negative_dataset_12(sequences = uniprot_seqs,
                                               uniprot_data = uniprot_tab,
                                               positive_dataset = positive_traintest,
                                               potential_AMPs = putative_amps)
  
  benchmark_all <- c(positive_benchmark,
                     benchmark_1, 
                     unlist(benchmark_2_3, recursive = FALSE),
                     unlist(benchmark_4_11, recursive = FALSE),
                     benchmark_12)
  write_fasta(benchmark_all, paste0(data_path, "Datasets/Benchmark_rep", i, ".fasta"))
}
