encode_seq <- function(x, property) {
  sapply(x, function(ith_seq) {
    mean(aaprop[property, tolower(ith_seq)])
  })
}

calculate_properties <- function(methods, n_rep, data_path) {
  lapply(methods, function(ith_method) {
    lapply(1:n_rep, function(j) {
      ds <- read_fasta(paste0(data_path, "Training_method_", ith_method, "_rep", j, ".fasta"))
      ds_neg <- ds[which(grepl("AMP=0", names(ds)))]
      data.frame(prot = names(ds_neg),
                 method = ith_method,
                 rep = j,
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
    }) %>% bind_rows()
  }) %>% bind_rows()
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
                 rep = j)
      }) %>% bind_rows()
    }) %>% bind_rows()
  }) %>% bind_rows()
} 