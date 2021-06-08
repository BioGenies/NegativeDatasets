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
                 rep = j)
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

get_ngram_counts_sum_pos <- function(positive) {
  ngram_counts_sum_pos <- count_multimers(positive,
                  k_vector = c(rep(2, 4), rep(3, 4)),
                  kmer_gaps_list = list(NULL, 1, 2, 3, c(0, 0), c(1, 0), c(0, 1), c(1, 1)),
                  alphabet = toupper(colnames(biogram::aaprop))) %>%
    as.matrix() %>% 
    colSums()/length(positive) 
  ngram_counts_sum_pos %>% 
    t() %>% 
    as.data.frame() %>% 
    mutate(method = "Positive",
           rep = 1)
}


get_aa_comp_heatmap <- function(aa_comp_all) {
  aa_comp_heatmap_dat <- aa_comp_all %>% 
    pivot_wider(names_from = "aa", values_from = "Freq", values_fill = 0) %>% 
    mutate(Dataset = ifelse(Dataset == "Positive", "Positive", paste0(method, "_rep", rep))) %>% 
    select(-c(method, rep)) 

  clustering_methods <- as.dendrogram(hclust(dist(as.matrix(aa_comp_heatmap_dat[,2:21]))))
  methods_order <- order.dendrogram(clustering_methods)
  
  dendro <- clustering_methods %>%
    dendro_data %>%
    segment %>%
    ggplot(aes(x = x, y = y, xend = xend, yend = yend)) +
    geom_segment() +
    scale_y_continuous("") +
    scale_x_discrete("",
                     limits = factor(1L:nobs(clustering_methods))) + 
    theme_void() + 
    coord_flip()
  
  aa_comp_clustered_dat <- pivot_longer(aa_comp_heatmap_dat, 2:21, names_to = "Amino acid", values_to = "Frequency")
  aa_comp_clustered_dat[["Dataset"]] <- factor(aa_comp_clustered_dat[["Dataset"]],
                                               levels = aa_comp_heatmap_dat[["Dataset"]][methods_order], ordered = TRUE)
  
  heatmap <- aa_comp_clustered_dat %>% 
    ggplot(aes(x = `Amino acid`, y = Dataset)) +
    geom_tile(aes(fill = Frequency)) +
    scale_fill_gradientn(colors = c("#ffffff", "#ffe96b", "#ff4242", "#630000"), values = rescale(c(0, 0.03, 0.08, 0.14), to = c(0, 1))) +
    theme_bw() +
    theme(legend.position = "bottom",
          legend.key.width = unit(2, "lines"))
  
  max_height <- unit.pmax(ggplotGrob(heatmap)[["heights"]],
                          ggplotGrob(dendro)[["heights"]])
  
  grob_list <- list(heatmap = ggplotGrob(heatmap),
                    dendrogram = ggplotGrob(dendro))
  grob_list[["heatmap"]][["heights"]] <-
    grob_list[["dendrogram"]][["heights"]] <-
    max_height
  
  grid.arrange(grob_list[["heatmap"]], grob_list[["dendrogram"]], widths = c(0.8, 0.2))
}


get_pca_aa_comp_plot <- function(aa_comp_all, dataset_colors) {
  
  pca_data_all <- aa_comp_all %>% 
    pivot_wider(names_from = aa, values_from = Freq, values_fill = 0) %>% 
    mutate(method = factor(method, levels = names(dataset_colors)))
  
  pca_res_all <- prcomp(pca_data_all[, 4:ncol(pca_data_all)], center = TRUE, scale = TRUE)
  
  ggbiplot(pca_res_all, choices = 1:2) +
    theme_bw() +
    geom_point(aes(color = pca_data_all[["method"]], shape = pca_data_all[["method"]]), size = 3) +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_color_manual("Dataset", values = dataset_colors, labels = names(dataset_colors)) +
    scale_shape_manual("Dataset", values = c(17, rep(16, 13)), labels = names(dataset_colors))
}


get_pca_prop_plot <- function(df_all, methods, dataset_colors) {
  props <- lapply(methods, function(i) {
    lapply(1:5, function(j) {
      data.frame(t(colMeans(filter(df_all, method == i, rep == j)[, 5:(ncol(df_all)-1)]))) %>% 
        mutate(method = factor(i),
               repetition = factor(j),
               dataset = "Negative")
    }) %>% bind_rows()
  }) %>% bind_rows()
  props_pos <- data.frame(t(colMeans(filter(df_all, method == "Positive")[, 5:(ncol(df_all)-1)]))) %>% 
    mutate(method = "Positive", 
           repetition = "1",
           dataset = "Positive")
  props_all <- bind_rows(props, props_pos) %>% 
    mutate(method = factor(method, levels = names(dataset_colors)))
  
  pca_prop_res_all <- prcomp(props_all[, 1:(ncol(props_all)-3)], center = TRUE, scale = TRUE)
  
  ggbiplot(pca_prop_res_all, choices = 1:2) +
    #ggtitle("PCA on means of physicochemical properties with positive dataset") +
    geom_point(aes(color = props_all[["method"]], shape = props_all[["method"]]), size = 3) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_color_manual("Dataset", values = dataset_colors, labels = names(dataset_colors)) +
    scale_shape_manual("Dataset", values = c(17, rep(16, 13)), labels = names(dataset_colors)) +
    xlim(c(-2.5, 6))
}


get_aa_comp_barplot <- function(aa_comp_all, dataset_colors) {
  aa_comp_all %>% 
    group_by(Dataset, aa, method) %>% 
    dplyr::summarise(Frequency = mean(Freq), sd = sd(Freq)) %>% 
    mutate(method = factor(method, levels = names(dataset_colors))) %>% 
    ggplot(aes(x = method, y = Frequency, fill = Dataset)) +
    geom_col(color = "black", size = 0.25) +
    facet_wrap(~aa, nrow = 4) +
    scale_fill_manual("Dataset", values = c(Negative = "#76bef2", Positive = "#ff4242")) +
    geom_errorbar(aes(ymin = Frequency - sd, ymax = Frequency + sd), width = 0.2,
                  position = position_dodge(0.9)) + 
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90)) +
    #ggtitle("Differences in amino acid frequency using different sampling methods") +
    theme(legend.position = "none")
}

get_pca_res_ngrams <- function(ngram_counts_sum_all) {
  ngram_counts_sum_all[is.na(ngram_counts_sum_all)] <- 0
  pca_ngram_res <- prcomp(select(ngram_counts_sum_all, -c("method", "rep", "Dataset")), center = TRUE, scale = TRUE)
}

plot_pca_res_ngrams <- function(pca_ngram_res, ngram_counts_sum_all, dataset_colors) {
  ggbiplot(pca_ngram_res, choices = 1:2, var.axes = FALSE) +
    geom_point(aes(color = ngram_counts_sum_all[["method"]], shape = ngram_counts_sum_all[["method"]]), size = 3) +
    theme_bw() +
    scale_color_manual("Dataset", values = dataset_colors, labels = names(dataset_colors)) +
    scale_shape_manual("Dataset", values = c(17, rep(16, 13)), labels = names(dataset_colors)) 
}

plot_pca_res_ngrams_zoom <- function(pca_ngram_res, ngram_counts_sum_all, dataset_colors) {
  ggbiplot(pca_ngram_res, choices = 1:2, var.axes = FALSE) +
    geom_point(aes(color = ngram_counts_sum_all[["method"]], shape = ngram_counts_sum_all[["method"]]), size = 3) +
    theme_bw() +
    scale_color_manual("Dataset", values = dataset_colors, labels = names(dataset_colors)) +
    scale_shape_manual("Dataset", values = c(17, rep(16, 13)), labels = names(dataset_colors)) +
    xlim(c(-0.7, -0.25)) +
    ylim(c(-0.85, -0.35)) + 
    theme(legend.position = "none", 
          axis.title.y = element_blank(),
          axis.title.x = element_blank())
}

get_sequence_length_plot <- function(df_all) {
  
  plist <- lapply(unique(df_all[["method"]]), function(ith_method) {
    df_all %>% 
      filter(method == ith_method) %>% 
      ggplot(aes(x = rep, y = len, group = rep, fill = Dataset)) +
      geom_violin() +
      theme_bw() +
      facet_wrap(~method) +
      theme(legend.position = "none",
            axis.title.x = element_blank(),
            axis.title.y = element_blank()) +
      scale_fill_manual("Dataset", values = c(Negative = "#76bef2", Positive = "#ff4242")) 
  }) %>% setNames(unique(df_all[["method"]]))
  
  blank <- ggplot() + theme_void()
  
  p <- plot_grid(plotlist = list(plist[["Positive"]], plist[["AMAP"]], plist[["GabereNoble"]], plist[["CSAMPPred"]], plist[["AmPEP"]],
                 blank, plist[["AmpGram"]], plist[["Wang"]], plist[["dbAMP"]], plist[["ampir-mature"]],
                 blank, plist[["AMPlify"]], blank, plist[["iAMP2L"]], plist[["ampir-precursor"]],
                 blank, plist[["AMPScannerV2"]], blank, blank, blank, blank, plist[["Witten"]]),
            nrow = 5, ncol = 5, rel_widths = c(1, 2.5, 2.5, 2.5, 2.5)) 
  
  grid.arrange(arrangeGrob(p, left = textGrob("Length", rot = 90), bottom = textGrob("Replication")))

}

get_statistical_analysis_plot_aa_comp_replicates <- function(aa_comp_peptides) {
  combns <- combn(unique(aa_comp_peptides[["rep"]]), 2, simplify = FALSE)
  test_res <- lapply(unique(aa_comp_peptides[["method"]]), function(ith_method) {
    lapply(seq_along(combns), function(ith_combn) {
      test_dat <- filter(aa_comp_peptides, method == ith_method, rep %in% combns[[ith_combn]])
      lapply(unique(aa_comp_peptides[["Amino acid"]]), function(ith_aa) {
        data.frame(method = ith_method,
                   comparison = paste0(combns[[ith_combn]][1], "_", combns[[ith_combn]][2]),
                   aa = ith_aa,
                   pval = wilcox.test(x = filter(test_dat, `Amino acid` == ith_aa, rep == combns[[ith_combn]][1])[["Frequency"]],
                                      y = filter(test_dat, `Amino acid` == ith_aa, rep == combns[[ith_combn]][2])[["Frequency"]],
                                      exact = FALSE)[["p.value"]])
      }) %>% bind_rows() 
    }) %>% bind_rows() %>%
      mutate(pval_adjusted = p.adjust(pval))
  }) %>% bind_rows()
  
  test_res %>%
    select(-pval) %>%
    group_by(aa, method) %>%
    dplyr::summarise(n_signif = as.factor(sum(pval_adjusted < 0.05))) %>%
    ggplot(aes(x = aa, y = method, fill = n_signif)) +
    geom_tile(color = "white") +
    scale_fill_manual("n of significant comparisons", values = "#76bef2", na.value = "grey90") +
    theme_bw() +
    theme(legend.position = "bottom")
}


get_sequence_length_table <- function(df_all, data_path) {
  df_all %>% 
    group_by(method, rep) %>% 
    dplyr::summarise(n = n()) %>% 
    pivot_wider(names_from = "rep", values_from = "n") %>% 
    xtable %>% 
    print(file = paste0(data_path, "Publication_results/sequence_length_table.txt"),
          include.rownames = FALSE)
}