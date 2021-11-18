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
  }) %>% bind_rows() %>% 
    change_method_names()
}


change_method_names <- function(df) {
  mutate(df, method = case_when(method == "iAMP2L" ~ "iAMP-2L",
                                method == "CSAMPPred" ~ "CS-AMPPred",
                                method == "Wang" ~ "Wang et al.",
                                method == "Witten" ~ "Witten&Witten",
                                method == "GabereNoble" ~ "Gabere&Noble",
                                method %in% c("dbAMP", "ampir-precursor", "ampir-mature", "AmpGram", "AMAP", 
                                              "AMPlify", "AMPScannerV2", "Positive") ~ method))
}

modify_labels <- function(df, types, shortcuts) {
  for(i in seq_along(types)) {
    df[[types[i]]] <- sapply(df[[types[i]]], function(j) paste0(shortcuts[i], ":", j))
  }
  df
}

calculate_properties <- function(ds_neg, method, rep) {
  data.frame(prot = names(ds_neg),
             method = method,
             rep = rep,
             len = lengths(ds_neg),
             `Residue volume (Bigelow, 1967)` = encode_seq(ds_neg, "BIGC670101"),
             `Polarizability parameter (Charton-Charton, 1982)` = encode_seq(ds_neg, "CHAM820101"),
             `Normalized frequency of alpha-helix (Chou-Fasman, 1978b)` = encode_seq(ds_neg, "CHOP780201"),
             `Normalized frequency of beta-sheet (Chou-Fasman, 1978b)` = encode_seq(ds_neg, "CHOP780202"),
             `Normalized frequency of beta-turn (Chou-Fasman, 1978b)` = encode_seq(ds_neg, "CHOP780203"),
             `Normalized van der Waals volume (Fauchere et al., 1988)` = encode_seq(ds_neg, "FAUJ880103"),
             `Net charge (Klein et al., 1984)` = encode_seq(ds_neg, "KLEP840101"),
             `Hydropathy index (Kyte-Doolittle, 1982)` = encode_seq(ds_neg, "KYTJ820101"),
             `Polarity (Zimmerman et al., 1968)` = encode_seq(ds_neg, "ZIMJ680103"),
             check.names = FALSE)
}

# This is a modified version of ggbiplot function from ggbiplot package
# https://github.com/vqv/ggbiplot/blob/master/R/ggbiplot.r
# Copyright 2011 Vincent Q. Vu.
my_ggbiplot <- function(pcobj, choices = 1:2, scale = 1, pc.biplot = TRUE, 
                        obs.scale = 1 - scale, var.scale = scale, 
                        groups = NULL, ellipse = FALSE, ellipse.prob = 0.68, 
                        labels = NULL, labels.size = 3, alpha = 1, 
                        var.axes = TRUE, 
                        circle = FALSE, circle.prob = 0.69, 
                        varname.size = 3, varname.adjust = 1.5, 
                        varname.abbrev = FALSE, ...) {
  
  stopifnot(length(choices) == 2)
  
  # Recover the SVD
  if(inherits(pcobj, 'prcomp')){
    nobs.factor <- sqrt(nrow(pcobj$x) - 1)
    d <- pcobj$sdev
    u <- sweep(pcobj$x, 2, 1 / (d * nobs.factor), FUN = '*')
    v <- pcobj$rotation
  } else if(inherits(pcobj, 'princomp')) {
    nobs.factor <- sqrt(pcobj$n.obs)
    d <- pcobj$sdev
    u <- sweep(pcobj$scores, 2, 1 / (d * nobs.factor), FUN = '*')
    v <- pcobj$loadings
  } else if(inherits(pcobj, 'PCA')) {
    nobs.factor <- sqrt(nrow(pcobj$call$X))
    d <- unlist(sqrt(pcobj$eig)[1])
    u <- sweep(pcobj$ind$coord, 2, 1 / (d * nobs.factor), FUN = '*')
    v <- sweep(pcobj$var$coord,2,sqrt(pcobj$eig[1:ncol(pcobj$var$coord),1]),FUN="/")
  } else if(inherits(pcobj, "lda")) {
    nobs.factor <- sqrt(pcobj$N)
    d <- pcobj$svd
    u <- predict(pcobj)$x/nobs.factor
    v <- pcobj$scaling
    d.total <- sum(d^2)
  } else {
    stop('Expected a object of class prcomp, princomp, PCA, or lda')
  }
  
  # Scores
  choices <- pmin(choices, ncol(u))
  df.u <- as.data.frame(sweep(u[,choices], 2, d[choices]^obs.scale, FUN='*'))
  
  # Directions
  v <- sweep(v, 2, d^var.scale, FUN='*')
  df.v <- as.data.frame(v[, choices])
  
  names(df.u) <- c('xvar', 'yvar')
  names(df.v) <- names(df.u)
  
  if(pc.biplot) {
    df.u <- df.u * nobs.factor
  }
  
  # Scale the radius of the correlation circle so that it corresponds to 
  # a data ellipse for the standardized PC scores
  r <- sqrt(qchisq(circle.prob, df = 2)) * prod(colMeans(df.u^2))^(1/4)
  
  # Scale directions
  v.scale <- rowSums(v^2)
  df.v <- r * df.v / sqrt(max(v.scale))
  
  # Change the labels for the axes
  if(obs.scale == 0) {
    u.axis.labs <- paste('standardized PC', choices, sep='')
  } else {
    u.axis.labs <- paste('PC', choices, sep='')
  }
  
  # Append the proportion of explained variance to the axis labels
  u.axis.labs <- paste(u.axis.labs, 
                       sprintf('(%0.1f%% explained var.)', 
                               100 * pcobj$sdev[choices]^2/sum(pcobj$sdev^2)))
  
  # Score Labels
  if(!is.null(labels)) {
    df.u$labels <- labels
  }
  
  # Grouping variable
  if(!is.null(groups)) {
    df.u$groups <- groups
  }
  
  # Variable Names
  if(varname.abbrev) {
    df.v$varname <- abbreviate(rownames(v))
  } else {
    df.v$varname <- rownames(v)
  }
  
  # Variables for text label placement
  df.v$angle <- with(df.v, (180/pi) * atan(yvar / xvar))
  df.v$hjust = with(df.v, (1 - varname.adjust * sign(xvar)) / 2)
  
  # Base plot
  g <- ggplot(data = df.u, aes(x = xvar, y = yvar)) + 
    xlab(u.axis.labs[1]) + ylab(u.axis.labs[2]) + coord_equal()
  
  if(var.axes) {
    # Draw circle
    if(circle) 
    {
      theta <- c(seq(-pi, pi, length = 50), seq(pi, -pi, length = 50))
      circle <- data.frame(xvar = r * cos(theta), yvar = r * sin(theta))
      g <- g + geom_path(data = circle, color = muted('white'), 
                         size = 1/2, alpha = 1/3)
    }
    
    # Draw directions
    g <- g +
      geom_segment(data = df.v,
                   aes(x = 0, y = 0, xend = xvar, yend = yvar),
                   arrow = arrow(length = unit(1/2, 'picas'), type = "closed"), 
                   color = 'grey90', linetype = "dashed")
  }
  
  # Draw either labels or points
  if(!is.null(df.u$labels)) {
    if(!is.null(df.u$groups)) {
      g <- g + geom_text(aes(label = labels, color = groups), 
                         size = labels.size)
    } else {
      g <- g + geom_text(aes(label = labels), size = labels.size)      
    }
  } else {
    if(!is.null(df.u$groups)) {
      # g <- g + geom_point(aes(color = groups), alpha = alpha)
    } else {
      #  g <- g + geom_point(alpha = alpha)      
    }
  }
  
  # Overlay a concentration ellipse if there are groups
  if(!is.null(df.u$groups) && ellipse) {
    theta <- c(seq(-pi, pi, length = 50), seq(pi, -pi, length = 50))
    circle <- cbind(cos(theta), sin(theta))
    
    ell <- ddply(df.u, 'groups', function(x) {
      if(nrow(x) <= 2) {
        return(NULL)
      }
      sigma <- var(cbind(x$xvar, x$yvar))
      mu <- c(mean(x$xvar), mean(x$yvar))
      ed <- sqrt(qchisq(ellipse.prob, df = 2))
      data.frame(sweep(circle %*% chol(sigma) * ed, 2, mu, FUN = '+'), 
                 groups = x$groups[1])
    })
    names(ell)[1:2] <- c('xvar', 'yvar')
    g <- g + geom_path(data = ell, aes(color = groups, group = groups))
  }
  
  # Label the variable axes
  if(var.axes) {
    g <- g + 
      geom_text(data = df.v, 
                aes(label = varname, x = xvar, y = yvar, 
                    angle = angle, hjust = hjust), 
                color = 'gray50', size = varname.size)
  }
  
  return(g)
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
  }) %>% bind_rows() %>% 
    change_method_names()
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
  }) %>% bind_rows() %>% 
    change_method_names()
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
                                kmer_alphabet = toupper(colnames(biogram::aaprop))) %>%
        as.matrix() %>% 
        colSums()/length(seqs) 
      ngrams %>% 
        t() %>% 
        as.data.frame() %>% 
        mutate(method = i,
               rep = j)
    }) %>% bind_rows()
  }) %>% bind_rows() %>% 
    change_method_names()
}

get_ngram_counts_sum_pos <- function(positive) {
  ngram_counts_sum_pos <- count_multimers(positive,
                                          k_vector = c(rep(2, 4), rep(3, 4)),
                                          kmer_gaps_list = list(NULL, 1, 2, 3, c(0, 0), c(1, 0), c(0, 1), c(1, 1)),
                                          kmer_alphabet = toupper(colnames(biogram::aaprop))) %>%
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
    select(-c(method, rep)) %>% 
    modify_labels(types = "Dataset", shortcuts = "TSM")
  
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
          legend.key.width = unit(2, "lines")) +
    ylab("Data set")
  
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
  
  my_ggbiplot(pca_res_all, choices = 1:2) +
    theme_bw() +
    geom_point(aes(color = pca_data_all[["method"]], shape = pca_data_all[["method"]]), size = 3) +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_color_manual("Data set", values = dataset_colors, labels = sapply(names(dataset_colors), function(i) ifelse(i != "Positive", paste0("TSM:", i), i))) +
    scale_shape_manual("Data set", values = c(17, rep(16, 11)), labels = sapply(names(dataset_colors), function(i) ifelse(i != "Positive", paste0("TSM:", i), i)))
}


get_pca_prop_plot <- function(df_all, dataset_colors) {
  props <- lapply(unique(filter(df_all, Dataset == "Negative")[["method"]]), function(i) {
    lapply(1:5, function(j) {
      data.frame(t(colMeans(filter(df_all, method == i, rep == j)[, 5:(ncol(df_all)-1)])),
                 check.names = FALSE) %>% 
        mutate(method = factor(i),
               repetition = factor(j),
               dataset = "Negative")
    }) %>% bind_rows()
  }) %>% bind_rows()
  props_pos <- data.frame(t(colMeans(filter(df_all, method == "Positive")[, 5:(ncol(df_all)-1)])),
                          check.names = FALSE) %>% 
    mutate(method = "Positive", 
           repetition = "1",
           dataset = "Positive")
  props_all <- bind_rows(props, props_pos) %>% 
    mutate(method = factor(method, levels = names(dataset_colors))) 
  
  pca_prop_res_all <- prcomp(props_all[, 1:(ncol(props_all)-3)], center = TRUE, scale = TRUE)
  
  my_ggbiplot(pca_prop_res_all, choices = 1:2, varname.size = 2, varname.adjust = 1.1) +
    geom_point(aes(color = props_all[["method"]], shape = props_all[["method"]]), size = 1.5) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_color_manual("Data set", values = dataset_colors, labels = sapply(names(dataset_colors), function(i) ifelse(i != "Positive", paste0("TSM:", i), i))) +
    scale_shape_manual("Data set", values = c(17, rep(16, 11)), labels = sapply(names(dataset_colors), function(i) ifelse(i != "Positive", paste0("TSM:", i), i))) +
    xlim(c(-5.5, 6.2)) +
    ylim(c(-5.5, 4))
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
    theme(legend.position = "none") +
    xlab("Sampling method used for generation of training negative data set (TSM)")
}

get_pca_res_ngrams <- function(ngram_counts_sum_all) {
  ngram_counts_sum_all[is.na(ngram_counts_sum_all)] <- 0
  pca_ngram_res <- prcomp(select(ngram_counts_sum_all, -c("method", "rep", "Dataset")), center = TRUE, scale = TRUE)
}

plot_pca_res_ngrams <- function(pca_ngram_res, ngram_counts_sum_all, dataset_colors) {
  ggbiplot(pca_ngram_res, choices = 1:2, var.axes = FALSE) +
    geom_point(aes(color = ngram_counts_sum_all[["method"]], shape = ngram_counts_sum_all[["method"]]), size = 3) +
    theme_bw() +
    scale_color_manual("Data set", values = dataset_colors, labels = sapply(names(dataset_colors), function(i) ifelse(i != "Positive", paste0("TSM:", i), i))) +
    scale_shape_manual("Data set", values = c(17, rep(16, 11)), labels = sapply(names(dataset_colors), function(i) ifelse(i != "Positive", paste0("TSM:", i), i)))
}



get_sequence_length_plot <- function(df_all) {
  df_all <- mutate(df_all, method = ifelse(method == "Positive", "Positive", paste0("TSM:", method)))
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
  
  p <- plot_grid(plotlist = list(plist[["Positive"]], plist[["TSM:AMAP"]], plist[["TSM:Gabere&Noble"]], plist[["TSM:CS-AMPPred"]], plist[["TSM:ampir-mature"]], 
                                 blank, plist[["TSM:AmpGram"]], plist[["TSM:Wang et al."]], plist[["TSM:dbAMP"]], blank,
                                 blank, plist[["TSM:AMPlify"]], blank, plist[["TSM:iAMP-2L"]], blank,
                                 blank, plist[["TSM:AMPScannerV2"]], blank, blank, blank, blank, plist[["TSM:Witten&Witten"]]),
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
    scale_fill_manual("Number of significant comparisons", values = "#76bef2", na.value = "grey90") +
    theme_bw() +
    theme(legend.position = "bottom") +
    xlab("Amino acid") +
    ylab("Sampling method used for generation\nof training negative data set (TSM)")
}


get_train_dataset_size_table <- function(df_all, data_path) {
  df_all %>% 
    group_by(method, rep) %>% 
    dplyr::summarise(n = n()) %>% 
    pivot_wider(names_from = "rep", values_from = "n") %>% 
    xtable %>% 
    print(file = paste0(data_path, "Publication_results/train_dataset_size_table.txt"),
          include.rownames = FALSE)
}

get_benchmark_dataset_size_table <- function(data_path, methods) {
  df <- lapply(1:5, function(ith_rep) {
    s <- read_fasta(paste0(data_path, "Datasets/Benchmark_rep", ith_rep, ".fasta"))
    lapply(methods, function(ith_method) {
      data.frame(method = ith_method,
                 rep = ith_rep,
                 n = as.character(sum(grepl(ith_method, names(s)))))
    }) %>% bind_rows()
  }) %>% bind_rows()
  pos <- read_fasta(paste0(data_path, "Datasets/Benchmark_rep1.fasta"))
  pos_df <- data.frame(method = "Positive",
                       rep = 1:5,
                       n = c(as.character(sum(grepl("AMP=1", names(pos)))), rep("", 4)))
  bind_rows(pos_df, df) %>% 
    pivot_wider(names_from = "rep", values_from = "n") %>% 
    xtable %>% 
    print(file = paste0(data_path, "Publication_results/benchmark_dataset_size_table.txt"),
          include.rownames = FALSE)
}

get_statistical_analysis_plot_aa_comp_methods <- function(aa_comp_peptides_all) {
  m <- as.character(unique(filter(aa_comp_peptides_all, method != "Positive")[["method"]]))
  aa_comp_peptides_all <- mutate(aa_comp_peptides_all, method = factor(ifelse(method == "Positive", "Positive", paste0("TSM:", method)),
                                                                       levels = c(paste0("TSM:", m), "Positive")))
  combns <- combn(unique(aa_comp_peptides_all[["method"]]), 2, simplify = FALSE)
  test_res_methods <- lapply(seq_along(combns), function(ith_combn) {
    test_dat <- filter(aa_comp_peptides_all, method  %in% combns[[ith_combn]])
    lapply(unique(aa_comp_peptides_all[["Amino acid"]]), function(ith_aa) {
      data.frame(method1 = combns[[ith_combn]][1], 
                 method2 = combns[[ith_combn]][2],
                 aa = ith_aa,
                 pval = wilcox.test(x = filter(test_dat, `Amino acid` == ith_aa, method == combns[[ith_combn]][1])[["Frequency"]],
                                    y = filter(test_dat, `Amino acid` == ith_aa, method == combns[[ith_combn]][2])[["Frequency"]],
                                    exact = FALSE)[["p.value"]])
    }) %>% bind_rows() 
  }) %>% bind_rows() %>%
    mutate(pval_adjusted = p.adjust(pval)) %>%
    select(-pval) %>%
    group_by(aa, method1, method2) %>%
    dplyr::summarise(is_signif = pval_adjusted < 0.05) 
  
  ggplot(test_res_methods, aes(x = method1, y = method2, fill = is_signif)) +
    geom_tile(color = "white") +
    facet_wrap(~aa) +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          axis.text.x = element_text(angle = 90)) +
    scale_fill_manual("Is significant", values = c(`FALSE` = "#76bef2", `TRUE` = "#ff4242")) +
    xlab("First method") +
    ylab("Second method")
}

get_results_plot_mean_auc_sd <- function(detailed_stats_mean) {
  dat <- detailed_stats_mean %>% 
    mutate(ident = seq_source == method) %>% 
    modify_labels(types = c("architecture"), shortcuts = c("A"))
  ggplot(dat, aes(x = method, y = seq_source, fill = mean_AUC)) +
    geom_tile(size = 1) +
    geom_tile(data = dat[dat[["ident"]] == TRUE, ], aes(color = ident), size = 1) +
    geom_point(data = dat, aes(x = method, y = seq_source, size = sd),
               color = "black") +
    facet_wrap(~architecture, ncol = 4) +
    scale_fill_gradient("Mean AUC", low =  "#ffe96b",  high = "#ff4242",
                        trans = scales::trans_new("square_exp", function(x) exp(x)^2, function(x) log(sqrt(x)))) +
    scale_size_continuous("Standard deviation") +
    scale_color_manual(guide = FALSE, values = c(`TRUE` = "black")) +
    theme_bw(base_size = 16) +
    theme(axis.text.x = element_text(angle = 90), legend.position = "bottom",
          legend.key.width = unit(2, "cm")) +
    xlab("Sampling method used for generation of training negative data set (TSM)") +
    ylab("Sampling method used for generation of benchmark negative data set (BSM)")
}

get_reference_auc_df <- function(detailed_stats_mean) {
  detailed_stats_mean %>%
    dplyr::filter(method == seq_source) %>% 
    select(architecture, seq_source, reference_AUC = mean_AUC)
}

get_reference_nonreference_AUC <- function(detailed_stats_mean) {
  reference_auc_df <- get_reference_auc_df(detailed_stats_mean)
  
  inner_join(reference_auc_df %>% 
               dplyr::group_by(architecture) %>% 
               dplyr::summarise(reference_mean_AUC = mean(reference_AUC)),
             detailed_stats_mean %>%
               filter(method != seq_source) %>% 
               select(architecture, method, seq_source, nonreference_AUC = mean_AUC) %>% 
               dplyr::group_by(architecture) %>% 
               dplyr::summarise(nonreference_mean_AUC = mean(nonreference_AUC))) 
}

get_reference_nonreference_AUC_table <- function(detailed_stats_mean) {
  df <- get_reference_nonreference_AUC(detailed_stats_mean)
  colnames(df) <- c("Architecture", "Mean reference AUC", "Mean nonreference AUC")
  print(xtable(df),
        file = paste0(data_path, "Publication_results/reference_vs_nonreference_auc_table.txt"),
        include.rownames = FALSE)
}

plot_reference_vs_nonreference <- function(detailed_stats_mean, architecture_colors) {
  names(architecture_colors) <- paste0("A:", names(architecture_colors))
  detailed_stats_mean %>% 
    modify_labels(types = c("architecture"), shortcuts = c("A")) %>% 
    get_reference_nonreference_AUC() %>% 
    ggplot(aes(x = reference_mean_AUC, y = nonreference_mean_AUC,
               color = architecture, label = architecture)) +
    geom_point(size = 3) +
    geom_abline(slope = 1, intercept = 0) +
    geom_label_repel() +
    scale_x_continuous("Mean AUC if trained and tested using\ndata sampled with the same method", 
                       limits = c(0.5, 1)) +
    scale_y_continuous("Mean AUC if trained and tested using\ndata sampled with different methods", 
                       limits = c(0.5, 1)) +
    scale_color_manual("Architecture", values = architecture_colors) +
    coord_equal() +
    theme_bw() +
    labs(tag = "A") +
    theme(plot.tag = element_text(size = 24),
          legend.position = "none")
}

plot_reference_vs_nonreference_by_train_method <- function(detailed_stats_mean, architecture_colors) {
  names(architecture_colors) <- paste0("A:", names(architecture_colors))
  reference_auc_df <- detailed_stats_mean %>% 
    get_reference_auc_df() 
  
  inner_join(reference_auc_df %>% 
               dplyr::group_by(architecture) %>% 
               dplyr::summarise(reference_mean_AUC = mean(reference_AUC)),
             detailed_stats_mean %>%
               filter(method != seq_source) %>% 
               select(architecture, method, seq_source, nonreference_AUC = mean_AUC) %>% 
               dplyr::group_by(architecture, method) %>% 
               dplyr::summarise(nonreference_mean_AUC = mean(nonreference_AUC)))  %>% 
    modify_labels(types = c("architecture", "method"), shortcuts = c("A", "TSM")) %>% 
    ggplot(aes(x = reference_mean_AUC, y = nonreference_mean_AUC,
               color = architecture, label = architecture)) +
    geom_point(size = 3) +
    facet_wrap(~method) +
    geom_abline(slope = 1, intercept = 0) +
    geom_label_repel() +
    scale_x_continuous("Mean AUC if trained and tested using data sampled with the same method", 
                       limits = c(0.5, 1)) +
    scale_y_continuous("Mean AUC if trained and tested using data sampled with different methods", 
                       limits = c(0.5, 1)) +
    scale_color_manual("Architecture", values = architecture_colors) +
    coord_equal() +
    theme_bw() +
    theme(legend.position = "none")
}


plot_auc_boxplot <- function(detailed_stats, type) {
  ggplot(detailed_stats, aes(x = get(type), y = AUC, fill = get(type))) +
    geom_boxplot(alpha = 0.25) +
    theme_bw() + 
    coord_flip() +
    theme(plot.title = element_text(hjust = 0.5),
          plot.tag = element_text(size = 24),
          aspect.ratio = 1,
          legend.position = "none")
}

plot_effect_boxplots <- function(detailed_stats, architecture_colors, dataset_colors) {
  a <- plot_auc_boxplot(detailed_stats, "architecture") +
    labs(tag = "B") +
    xlab("Architecture (A)") +
    scale_fill_manual(values = architecture_colors)
  m <- plot_auc_boxplot(detailed_stats, "method") +
    labs(tag = "C") +
    xlab("Training data set sampling method (TSM)") +
    scale_fill_manual(values = dataset_colors)
  s <- plot_auc_boxplot(detailed_stats, "seq_source") +
    labs(tag = "D") +
    xlab("Benchmark data set sampling method (BSM)")  +
    scale_fill_manual(values = dataset_colors)
  list(a, m, s)
}

plot_ref_vs_nonref_and_effects <- function(detailed_stats, detailed_stats_mean, architecture_colors, dataset_colors) {
  boxplots <- plot_effect_boxplots(detailed_stats, architecture_colors, dataset_colors)
  ref_vs_nonref <- plot_reference_vs_nonreference(detailed_stats_mean, architecture_colors) +
    theme(legend.position = "none") 
  (ref_vs_nonref + boxplots[[1]]) / (boxplots[[2]] + boxplots[[3]])
}

get_ref_vs_nonref_test_table <- function(detailed_stats) {
  lapply(unique(detailed_stats[["architecture"]]), function(i) {
    dat <- filter(detailed_stats, architecture == i) %>% 
      mutate(type = ifelse(method == seq_source, "same", "diff"))
    data.frame(Architecture = i,
               pval = kruskal.test(dat[["AUC"]], dat[["type"]])[["p.value"]])
  }) %>% bind_rows() %>% 
    mutate(`Bonferroni corrected p-value` = p.adjust(pval, method = "bonferroni")) %>% 
    arrange(Architecture) %>% 
    select(-pval) %>% 
    xtable(digits = -2) %>% 
    print(file = paste0(data_path, "Publication_results/reference_vs_nonreference_test_table.txt"),
          include.rownames = FALSE)
}

get_pairwise_paired_wilcox_test_table <- function(detailed_stats_mean, type, outfile) {
  res <- pairwise.wilcox.test(detailed_stats_mean[["mean_AUC"]], 
                       as.factor(detailed_stats_mean[[type]]), 
                       p.adjust.method="bonferroni", 
                       paired = TRUE) 
  res[["p.value"]] %>% 
    xtable(digits = -2)  %>% 
    print(file = outfile)
}

get_mean_sd_table <- function(detailed_stats_mean, outfile) {
  detailed_stats_mean %>% 
    group_by(architecture) %>% 
    summarise(`Mean SD` = mean(sd),
              `Min SD` = min(sd),
              `Max SD` = max(sd)) %>% 
    xtable(digits = 5) %>% 
    print(file = outfile)
}
