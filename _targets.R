library(dplyr)
library(biogram)
library(pbapply)
library(targets)
library(ggplot2)
library(tidyr)
library(ggdendro)
library(grid)
library(gridExtra)
library(ggbiplot)
library(cowplot)
library(xtable)
library(seqR)


if(Sys.info()[["nodename"]] == "kasia-MACH-WX9") {
  data_path <- "/media/kasia/Data/Dropbox/Projekty/BioNgramProjects/NegativeDatasets/"
}

if(Sys.info()[["nodename"]] == "ryzen") {
  data_path <- "~/Dropbox/Projekty/BioNgramProjects/NegativeDatasets/"
}

source("functions/generate_positive_dataset.R")
source("functions/read_uniprot_seqs.R")
source("functions/filter_sequences.R")
source("functions/generate_holdouts.R")
source("functions/generate_negative_datasets.R")
source("functions/generate_all_traintest_datasets.R")
source("functions/get_putative_amps.R")
source("functions/plot_functions.R")


list(
  tar_target(
    uniprot_seqs_file,
    paste0(data_path, "Data/uniprot-reviewed yes.fasta"),
    format = "file"),
  tar_target(
    uniprot_tab_file,
    paste0(data_path, "Data/uniprot-reviewed yes.tab"),
    format = "file"),
  tar_target(
    positive_dataset,
    generate_positive_dataset("data/dbaasp.csv")),
  tar_target(
    uniprot_seqs,
    read_uniprot_seqs(uniprot_seqs_file)),
  tar_target(
    uniprot_tab,
    read.delim(uniprot_tab_file)),
  tar_target(
    putative_amps,
    get_putative_amps(data_path)),
  tar_target(
    holdouts,
    generate_holdout_groups(positive_dataset)),
  tar_target(
    positive_traintest,
    positive_dataset[unlist(sapply(holdouts, function(i) i[["traintest"]]), use.names = FALSE)]),
  tar_target(
    positive_benchmark,
    positive_dataset[unlist(sapply(holdouts, function(i) i[["benchmark"]]), use.names = FALSE)]),
  tar_target(
    datasets,
    generate_all_traintest_datasets(data_path = data_path,
                                    seed_vector = c(74983, 10846, 26542, 90183, 76351),
                                    n_rep = 5,
                                    uniprot_seqs = uniprot_seqs,
                                    uniprot_tab = uniprot_tab,
                                    positive_traintest = positive_traintest,
                                    putative_amps = putative_amps, 
                                    positive_benchmark = positive_benchmark)
  ),
  tar_target(
    methods,
    c("AMAP", "AmpGram", "ampir-mature", "ampir-precursor", "AMPlify",
      "AMPScannerV2", "CSAMPPred", "dbAMP", "GabereNoble", "iAMP2L", "Wang", "Witten")
  ),
  tar_target(
    dataset_colors,
    c(Positive = "#ff4242", AMAP = "#f2b176", AmpGram = "#f27676", `ampir-mature` = "#76bef2",
      `ampir-precursor` = "#76eff2", AMPlify = "#f2d676", AMPScannerV2 = "#b976f2", `CS-AMPPred` = "#76f2be", 
      dbAMP = "#f276e8", `Gabere&Noble` = "#7688f2", `iAMP-2L` = "#eff275", `Wang et. al` = "#80f276", `Witten&Witten` = "#ccf276")
  ),
  tar_target(
    aa_comp_peptides,
    calculate_aa_comp_peptides(methods, 5, paste0(data_path, "Datasets/"))
  ),
  tar_target(
    aa_comp_peptides_pos,
    get_aa_comp_peptides_positive(positive_traintest)
  ),
  tar_target(
    aa_comp_peptides_all,
    bind_rows(aa_comp_peptides,
              aa_comp_peptides_pos) %>% 
      mutate(method = factor(method, levels = c(methods, "Positive")))
  ),
  tar_target(
    aa_comp,
    calculate_aa_comp_datasets(methods, 5, paste0(data_path, "Datasets/"))
  ),
  tar_target(
    aa_comp_pos,
    get_aa_comp_pos(positive_traintest)
  ),
  tar_target(
    aa_comp_all,
    bind_rows(mutate(aa_comp, 
                     rep = as.character(rep),
                     Dataset = "Negative"), 
              mutate(aa_comp_pos,
                     Dataset = "Positive"))
  ),
  tar_target(
    df_neg,
    get_prop_df(methods, 5, paste0(data_path, "Datasets/"))
  ),
  tar_target(
    df_pos,
    calculate_properties(positive_traintest, "Positive", "1")
  ),
  tar_target(
    df_all,
    bind_rows(mutate(df_pos, Dataset = "Positive"),
              mutate(df_neg, Dataset = "Negative",
                     rep = as.character(rep))) %>% 
      mutate(method = factor(method, levels = names(dataset_colors)))
  ),
  tar_target(
    ngram_counts_sum,
    get_ngram_counts_sum(methods, 5, paste0(data_path, "Datasets/"))
  ),
  tar_target(
    ngram_counts_sum_pos,
    get_ngram_counts_sum_pos(positive_traintest)
  ),
  tar_target(
    ngram_counts_sum_all,
    bind_rows(mutate(ngram_counts_sum, Dataset = "Negative"),
              mutate(ngram_counts_sum_pos, Dataset = "Positive")) %>% 
      mutate(method = factor(method, levels = names(dataset_colors))) 
  ),
  tar_target(
    aa_comp_heatmap,
    ggsave(filename = "aa_comp_heatmap.eps",
           plot = get_aa_comp_heatmap(aa_comp_all),
           path = paste0(data_path, "Publication_results/"),
           width = 10, height = 10)
  ),
  tar_target(
    pca_aa_comp,
    ggsave(filename = "pca_aa_comp.eps",
           get_pca_aa_comp_plot(aa_comp_all, dataset_colors),
           path = paste0(data_path, "Publication_results/"),
           width = 10, height = 10)
  ),
  tar_target(
    pca_prop,
    ggsave(filename = "pca_prop.eps",
           get_pca_prop_plot(df_all, dataset_colors),
           path = paste0(data_path, "Publication_results/"),
           width = 10, height = 10)
  ),
  tar_target(
    ngram_pca,
    get_pca_res_ngrams(ngram_counts_sum_all)
  ),
  tar_target(
    ngram_pca_plot,
    plot_pca_res_ngrams(ngram_pca, ngram_counts_sum_all, dataset_colors)
  ),
  tar_target(
    ngram_pca_plot_zoom,
    plot_pca_res_ngrams_zoom(ngram_pca, ngram_counts_sum_all, dataset_colors)
  ),
  tar_target(
    ngram_pca_combined,
    ggsave(filename = "pca_ngrams.eps",
           ngram_pca_plot + annotation_custom(ggplotGrob(ngram_pca_plot_zoom), 
                                              xmin = 1, xmax = 3, ymin = 0, ymax = 1.75),
           path = paste0(data_path, "Publication_results/"),
           width = 10, height = 10)
  ),
  tar_target(
    sequence_length_plot,
    ggsave(filename = "sequence_length.eps",
           get_sequence_length_plot(df_all),
           path = paste0(data_path, "Publication_results/"),
           width = 12, height = 8)
  ),
  tar_target(
    aa_comp_peptides_test_plot,
    ggsave(filename = "aa_comp_test_replicates.eps",
           get_statistical_analysis_plot_aa_comp_replicates(aa_comp_peptides),
           path = paste0(data_path, "Publication_results/"),
           width = 6, height = 4)
  ),
  tar_target(
    aa_comp_barplot,
    ggsave(filename = "aa_comp_barplot.eps",
           get_aa_comp_barplot(aa_comp_all, dataset_colors),
           path = paste0(data_path, "Publication_results/"),
           width = 10, height = 8)
  ),
  tar_target(
    sequence_length_table,
    get_sequence_length_table(df_all, data_path) 
  ),
  tar_target(
    aa_comp_methods_test_plot,
    ggsave(filename = "aa_comp_test_methods.eps",
          get_statistical_analysis_plot_aa_comp_methods(aa_comp_peptides_all),
          path = paste0(data_path, "Publication_results/"),
          width = 10, height = 9))
)

