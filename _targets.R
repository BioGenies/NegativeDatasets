library(biogram)
library(pbapply)
library(targets)
library(ggplot2)
library(tidyr)
library(ggdendro)
library(grid)
library(gridExtra)
library(ggbiplot)
library(xtable)
library(seqR)
library(scales)
library(mlr3measures)
library(hmeasure)
library(ggrepel)
library(dplyr)
library(patchwork)


if(Sys.info()[["nodename"]] == "kasia-MACH-WX9") {
  data_path <- "/media/kasia/Data/Dropbox/Projekty/BioNgramProjects/NegativeDatasets/"
}

if(Sys.info()[["nodename"]] == "ryzen") {
  data_path <- "~/Dropbox/Projekty/BioNgramProjects/NegativeDatasets/"
}

if(Sys.info()[["nodename"]] == "huawei") {
  data_path <- "/home/michal/Dropbox/BioNgramProjects/NegativeDatasets/"
}

source("functions/generate_positive_dataset.R")
source("functions/read_uniprot_seqs.R")
source("functions/filter_sequences.R")
source("functions/generate_holdouts.R")
source("functions/generate_negative_datasets.R")
source("functions/generate_all_traintest_datasets.R")
source("functions/get_putative_amps.R")
source("functions/plot_functions.R")
source("functions/analyse_results.R")


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
    c("AMAP", "AmpGram", "ampir-mature", "AMPlify", "AMPScannerV2", "CSAMPPred", 
      "dbAMP", "GabereNoble", "iAMP2L", "Wang", "Witten")
  ),
  tar_target(
    dataset_colors,
    c(Positive = "#ff4242", AMAP = "#f2b176", AmpGram = "#f27676", `ampir-mature` = "#76bef2",
      AMPlify = "#f2d676", AMPScannerV2 = "#b976f2", `CS-AMPPred` = "#76f2be", dbAMP = "#f276e8", 
      `Gabere&Noble` = "#7688f2", `iAMP-2L` = "#eff275", `Wang et al.` = "#80f276", `Witten&Witten` = "#ccf276")
  ),
  tar_target(
    architecture_colors,
    c(AMAP = "#b67f49", AmPEP = "#6cb649", AmPEPpy = "#33803f", AmpGram = "#bc5658", Ampir = "#497db6", AMPScannerV2 = "#7f49b6", 
      `CS-AMPPred` = "#49b5b6", `Deep-AmPEP30` = "#b6498b", `iAMP-2L` = "#e1df81", MACREL = "#81b6e1", MLAMP = "#8e8e8e", `SVM-LZ` = "#d0ad2f")
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
      mutate(method = factor(method, levels = c(names(dataset_colors[2:length(dataset_colors)]), "Positive")))
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
    ggsave(filename = "pca_ngrams.eps",
           plot_pca_res_ngrams(ngram_pca, ngram_counts_sum_all, dataset_colors),
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
    train_dataset_size_table,
    get_train_dataset_size_table(df_all, data_path) 
  ),
  tar_target(
    benchmark_dataset_size_table,
    get_benchmark_dataset_size_table(data_path, 
                                     c("AMAP", "AmpGram", "ampir-mature", "AMPlify", "AMPScannerV2", 
                                       "CS-AMPPred", "dbAMP", "Gabere&Noble", "iAMP-2L", "Wang", "Witten&Witten"))
  ),
  tar_target(
    aa_comp_methods_test_plot,
    ggsave(filename = "aa_comp_test_methods.eps",
           get_statistical_analysis_plot_aa_comp_methods(aa_comp_peptides_all),
           path = paste0(data_path, "Publication_results/"),
           width = 10, height = 9)),
  tar_target(
    architectures,
    c("AMAP", "AmPEP", "AmPEPpy", "Ampir", "AMPScannerV2", "CS-AMPPred", 
      "Deep-AmPEP30", "iAMP-2L", "AmpGram", "MACREL", "MLAMP", "SVM-LZ")
  ),
  tar_target(
    all_results,
    aggregate_all_results(results_path = paste0(data_path, "Results/"),
                          methods = methods,
                          architectures = architectures)
  ),
  tar_target(
    seqtype_all_results,
    get_seq_source(all_results)
  ),
  tar_target(
    detailed_stats,
    get_detailed_stats(seqtype_all_results, architectures)
  ),
  tar_target(
    detailed_stats_mean,
    get_detailed_stats_mean(detailed_stats)
  ),
  tar_target(
    mean_auc_sd_plot,
    ggsave(filename = "results_mean_AUC+SD.eps",
           get_results_plot_mean_auc_sd(detailed_stats_mean),
           path = paste0(data_path, "Publication_results/"),
           width = 18, height = 16)
  ),
  tar_target(
    ref_vs_nonref_and_effects_plot,
    ggsave(filename = "reference_vs_nonreference+effects.pdf",
           plot_ref_vs_nonref_and_effects(detailed_stats, detailed_stats_mean, architecture_colors, dataset_colors),
           path = paste0(data_path, "Publication_results/"),
           width = 14, height = 14)
  ),
  tar_target(
    ref_vs_nonref_by_train_method_plot,
    ggsave(filename = "reference_vs_nonreference_by_train_method.eps",
           plot_reference_vs_nonreference_by_train_method(detailed_stats_mean, architecture_colors),
           path = paste0(data_path, "Publication_results/"),
           width = 24, height = 22)
  ),
  tar_target(
    ref_vs_nonref_table,
    get_reference_nonreference_AUC_table(detailed_stats_mean)
  ),
  tar_target(
    ref_vs_nonref_stat_test_table,
    get_ref_vs_nonref_test_table(detailed_stats)
  ),
  tar_target(
    wilcox_test_architectures,
    get_pairwise_paired_wilcox_test_table(detailed_stats_mean, 
                                          "architecture",
                                          paste0(data_path, "Publication_results/pairwise_wilcoxon_table_architectures.txt"))
  ),
  tar_target(
    wilcox_test_tsm,
    get_pairwise_paired_wilcox_test_table(detailed_stats_mean, 
                                          "method",
                                          paste0(data_path, "Publication_results/pairwise_wilcoxon_table_TSM.txt"))
  ),
  tar_target(
    wilcox_test_bsm,
    get_pairwise_paired_wilcox_test_table(detailed_stats_mean, 
                                          "seq_source",
                                          paste0(data_path, "Publication_results/pairwise_wilcoxon_table_BSM.txt"))
  ),
  tar_target(
    mean_sd_table,
    get_mean_sd_table(detailed_stats_mean, 
                      paste0(data_path, "Publication_results/mean_sd_table.txt"))
  )
)

