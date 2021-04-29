generate_all_traintest_datasets <- function(data_path, seed_vector, n_rep, uniprot_seqs, uniprot_tab, positive_traintest, putative_amps, positive_benchmark) {
  for (i in 1:n_rep) {
    set.seed(seed_vector[[i]])
    # Sampling algorithms that need only sequences and uniprot data - iAMP-2L, AmPEP
    benchmark_seq_dat <- lapply(paste0("generate_negative_dataset_"), c("iAMP2L", "AmPEP"), function(ith_fun) {
      s <- match.fun(ith_fun)(sequences = uniprot_seqs,
                              uniprot_data = uniprot_tab)
      s_train <- sample(s, 0.8*length(s), replace = FALSE)
      write_fasta(c(positive_traintest, s_train), paste0(data_path, "Datasets/Training_method_", last(strsplit(ith_fun, "_")[[1]]), "_rep", i, ".fasta"))
      s[which(!(names(s) %in% names(s_train)))]
    })
    
    # Sampling algorithms that need only sequences and positive dataset - ampir
    benchmark_seq_pos <- lapply(paste0("generate_negative_dataset_", c("ampir-precursor", "ampir-mature")), function(ith_fun) {
      s <- match.fun(ith_fun)(sequences = uniprot_seqs, 
                              positive_dataset = positive_traintest,
                              n_threads = 12)
      write_fasta(c(positive_traintest, s), paste0(data_path, "Datasets/Training_method_", last(strsplit(ith_fun, "_")[[1]]), "_rep", i, ".fasta"))
      match.fun(ith_fun)(sequences = uniprot_seqs, 
                         positive_dataset = positive_benchmark)
    })
    
    # Sampling algorithms that need sequences, uniprot data and positive dataset
    benchmark_seq_dat_pos <- lapply(paste0("generate_negative_dataset_", c("CSAMPPred", "Wang", "AmpGram", "dbAMP", "Witten", "AMPScannerV2", "GabereNoble", "AMAP")), 
                                    function(ith_fun) {
                                      s <- match.fun(ith_fun)(sequences = uniprot_seqs, 
                                                              uniprot_data = uniprot_tab,
                                                              positive_dataset = positive_traintest)
                                      write_fasta(c(positive_traintest, s), paste0(data_path, "Datasets/Training_method_", last(strsplit(ith_fun, "_")[[1]]), "_rep", i, ".fasta"))
                                      match.fun(ith_fun)(sequences = uniprot_seqs,
                                                         uniprot_data = uniprot_tab,
                                                         positive_dataset = positive_benchmark)
                                    })
    
    # Sampling algorithm that needs sequences, uniprot data, positive dataset and putative AMPs - AMPlify
    s_seq_dat_pos_putative <- generate_negative_dataset_AMPlify(sequences = uniprot_seqs,
                                                                uniprot_data = uniprot_tab, 
                                                                positive_dataset = positive_traintest,
                                                                potential_AMPs = putative_amps)
    write_fasta(c(positive_traintest, s_seq_dat_pos_putative), paste0(data_path, "Datasets/Training_method_AMPlify_rep", i, ".fasta"))
    
    benchmark_seq_dat_pos_putative <- generate_negative_dataset_AMPlify(sequences = uniprot_seqs,
                                                                        uniprot_data = uniprot_tab,
                                                                        positive_dataset = positive_benchmark,
                                                                        potential_AMPs = putative_amps)
    
    benchmark_all <- c(positive_benchmark,
                       unlist(benchmark_seq_dat, recursive = FALSE),
                       unlist(benchmark_seq_pos, recursive = FALSE),
                       unlist(benchmark_seq_dat_pos, recursive = FALSE),
                       benchmark_seq_dat_pos_putative)
    write_fasta(benchmark_all, paste0(data_path, "Datasets/Benchmark_rep", i, ".fasta"))
  }
}

