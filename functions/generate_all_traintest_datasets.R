generate_all_traintest_datasets <- function(data_path, seed_vector, n_rep, uniprot_seqs, uniprot_tab, positive_traintest, putative_amps, positive_benchmark) {
  for (i in 1:n_rep) {
    set.seed(seed_vector[[i]])
    # Sampling algorithms that need only sequences and uniprot data - iAMP-2L, AmPEP, dbAMP
    benchmark_seq_dat <- lapply(paste0("generate_negative_dataset_", c("iAMP2L", "AmPEP", "dbAMP")), function(ith_fun) {
      s <- match.fun(ith_fun)(sequences = uniprot_seqs,
                              uniprot_data = uniprot_tab)
      s_train <- sample(s, 0.8*length(s), replace = FALSE)
      write_fasta(c(positive_traintest, s_train), paste0(data_path, "Datasets/Training_method_", last(strsplit(ith_fun, "_")[[1]]), "_rep", i, ".fasta"))
      s[which(!(names(s) %in% names(s_train)))]
    })
    
    # Sampling algorithm that need only sequences and positive dataset - ampir precursor
    s_ampir_precursor <- generate_negative_dataset_ampir_precursor(sequences = uniprot_seqs,
                                                                   positive_dataset = positive_traintest,
                                                                   sequences_to_filter_out = positive_traintest,
                                                                   n_threads = 20)
    write_fasta(c(positive_traintest, s_ampir_precursor),  paste0(data_path, "Datasets/Training_method_ampir-precursor_rep", i, ".fasta"))
    benchmark_ampir_precursor <- generate_negative_dataset_ampir_precursor(sequences = uniprot_seqs,
                                                                           positive_dataset = positive_traintest,
                                                                           sequences_to_filter_out = c(positive_traintest, s_ampir_precursor),
                                                                           n_threads = 20)
    
    # Sampling algorithm that need only sequences and positive dataset - ampir mature
    all_ampir_mature <- filter_nonstandard_aa(
      filter_by_lengths(
        filter_with_cdhit(sequences = uniprot_seqs, threshold = 0.9, n_threads = n_threads),
        11, 39))
    s_ampir_mature <- sample(all_ampir_mature, 0.8*length(all_ampir_mature), replace = FALSE) %>% 
      filter_out_positive_sequences(positive_traintest)
    write_fasta(c(positive_traintest, s_ampir_mature), paste0(data_path, "Datasets/Training_method_ampir-mature_rep", i, ".fasta"))
    benchmark_ampir_mature <- all_ampir_mature[which(!(names(all_ampir_mature) %in% names(s_ampir_mature)))] %>% 
      filter_out_positive_sequences(positive_benchmark)
    
    # Sampling algorithms that need sequences, uniprot data and positive dataset - CSAMPpred
    csamppred <- generate_negative_dataset_CSAMPPred(sequences = uniprot_seqs, 
                                                     uniprot_data = uniprot_tab,
                                                     positive_dataset = c(positive_traintest, positive_benchmark))
    s_csamppred <- sample(csamppred, length(positive_traintest), replace = FALSE)
    write_fasta(c(positive_traintest, s_csamppred), paste0(data_path, "Datasets/Training_method_CSAMPPred_rep", i, ".fasta"))
    benchmark_csamppred <- csamppred[which(!(names(csamppred) %in% names(s_csamppred)))]
    
    
    # Sampling algorithms that need sequences, uniprot data and positive dataset
    benchmark_seq_dat_pos <- lapply(paste0("generate_negative_dataset_", c("Wang", "AmpGram", "Witten", "AMPScannerV2", "GabereNoble", "AMAP")), 
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
                       benchmark_ampir_precursor,
                       benchmark_ampir_mature,
                       benchmark_csamppred,
                       unlist(benchmark_seq_dat_pos, recursive = FALSE),
                       benchmark_seq_dat_pos_putative)
    write_fasta(benchmark_all, paste0(data_path, "Datasets/Benchmark_rep", i, ".fasta"))
  }
}

