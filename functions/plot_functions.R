encode_seq <- function(x, property) {
  sapply(x, function(ith_seq) {
    mean(aaprop[property, tolower(ith_seq)])
  })
}

get_prop_df <- function(methods, n_rep, data_path) {
  lapply(methods, function(ith_method) {
    lapply(1:n_rep, function(j) {
      ds <- read_fasta(paste0(data_path, "Training_method_", ith_method, "_rep", j, ".fasta"))
      ds_neg <- ds[which(grepl("AMP=0", names(ds)))]
      calculate_properties(ds_neg, ith_method, j)
    }) %>% bind_rows()
  }) %>% bind_rows()
}

calculate_properties <- function(ds_neg, method, rep) {
  data.frame(prot = names(ds_neg),
             method = method,
             rep = rep,
             len = lengths(ds_neg),
             BIGC670101 = encode_seq(ds_neg, "BIGC670101"),
             ARGP820101 = encode_seq(ds_neg, "ARGP820101"),
             CHAM820101 = encode_seq(ds_neg, "CHAM820101"),
             CHOP780201 = encode_seq(ds_neg, "CHOP780201"),
             CHOP780202 = encode_seq(ds_neg, "CHOP780202"),
             CHOP780203 = encode_seq(ds_neg, "CHOP780203"),
             FASG760101 = encode_seq(ds_neg, "FASG760101"),
             FASG760104 = encode_seq(ds_neg, "FASG760104"),
             FASG760105 = encode_seq(ds_neg, "FASG760105"),
             FAUJ880103 = encode_seq(ds_neg, "FAUJ880103"),
             KLEP840101 = encode_seq(ds_neg, "KLEP840101"),
             KYTJ820101 = encode_seq(ds_neg, "KYTJ820101"),
             ZIMJ680103 = encode_seq(ds_neg, "ZIMJ680103"),
             ENGD860101 = encode_seq(ds_neg, "ENGD860101"),
             FASG890101 = encode_seq(ds_neg, "FASG890101"))
}


calculate_aa_comp_datasets <- function(methods, n_rep, data_path) {
  lapply(methods, function(i) {
    lapply(1:n_rep, function(j) {
      ds <- read_fasta(paste0(data_path, "Training_method_", i, "_rep", j, ".fasta"))
      aa <- unlist(ds[which(grepl("AMP=0", names(ds)))], use.names = FALSE)
      as.data.frame(table(aa)/length(aa)) %>% 
        mutate(method = i,
               rep = j)
    }) %>% bind_rows()
  }) %>% bind_rows()
}


calculate_aa_comp_peptides <- function(methods, n_rep, data_path) {
  lapply(methods, function(i) {
    lapply(1:n_rep, function(j) {
      ds <- read_fasta(paste0(data_path, "Training_method_", i, "_rep", j, ".fasta"))
      neg <- ds[which(grepl("AMP=0", names(ds)))]
      lapply(names(neg), function(ith_prot) {
        data.frame(table(factor(neg[[ith_prot]], levels = c("A", "C", "D", "E", "F", "G", "H",
                                                            "I", "K", "L", "M", "N", "P", "Q",
                                                            "R", "S", "T", "V", "W", "Y")))/length(neg[[ith_prot]])) %>%
          setNames(c("Amino acid", "Frequency"))  %>%
          mutate(method = i,
                 rep = j,
                 prot = names(neg))
      }) %>% bind_rows()
    }) %>% bind_rows()
  }) %>% bind_rows()
} 

get_aa_comp_peptides_positive <- function(positive) {
  lapply(names(positive), function(ith_prot) {
    data.frame(table(factor(positive[[ith_prot]], levels = c("A", "C", "D", "E", "F", "G", "H",
                                                             "I", "K", "L", "M", "N", "P", "Q",
                                                             "R", "S", "T", "V", "W", "Y")))/length(positive[[ith_prot]])) %>%
      setNames(c("Amino acid", "Frequency"))  %>%
      mutate(method = "Positive",
             rep = 1)
  }) %>% bind_rows()
}

get_aa_comp_pos <- function(positive) {
  aa_pos <- unlist(positive, use.names = FALSE)
  as.data.frame(table(aa_pos)/length(aa_pos)) %>% 
    setNames(c("aa", "Freq")) %>% 
    mutate(method = "Positive",
           rep = "1")
}

get_ngram_counts_sum <- function(methods, n_rep, data_path) {
  lapply(methods, function(i) {
    lapply(1:n_rep, function(j) {
      ds <- read_fasta(paste0(data_path, "Training_method_", i, "_rep", j, ".fasta"))
      seqs <- ds[which(grepl("AMP=0", names(ds)))]
      ngrams <- count_multimers(seqs,
                                k_vector = c(rep(2, 4), rep(3, 4)),
                                kmer_gaps_list = list(NULL, 1, 2, 3, c(0, 0), c(1, 0), c(0, 1), c(1, 1)),
                                alphabet = toupper(colnames(biogram::aaprop))) %>%
        as.matrix() %>% 
        colSums()/length(seqs) 
      ngrams %>% 
        t() %>% 
        as.data.frame() %>% 
        mutate(method = i,
               rep = j)
    }) %>% bind_rows()
  }) %>% bind_rows()
}

get_ngram_counts_pos <- function(positive) {
  count_multimers(positive,
                  k_vector = c(rep(2, 4), rep(3, 4)),
                  kmer_gaps_list = list(NULL, 1, 2, 3, c(0, 0), c(1, 0), c(0, 1), c(1, 1)),
                  alphabet = toupper(colnames(biogram::aaprop))) %>%
    as.matrix() %>% 
    colSums()/length(positive) 
  ngram_counts_sum_pos <- ngram_counts_sum_pos %>% 
    t() %>% 
    as.data.frame() %>% 
    mutate(method = "Positive",
           rep = 1)
}