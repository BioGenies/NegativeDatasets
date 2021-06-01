library(dplyr)
library(biogram)
library(pbapply)
library(targets)

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
    c("AMAP", "AmPEP", "AmpGram", "ampir-mature", "ampir-precursor", "AMPlify",
      "AMPScannerV2", "CSAMPPred", "dbAMP", "GabereNoble", "iAMP2L", "Wang", "Witten")
  ),
  tar_target(
    dataset_colors,
    c(Positive = "#ff4242", AMAP = "#f2b176", AmPEP = "#b4b4b4", AmpGram = "#f27676", `ampir-mature` = "#76bef2",
      `ampir-precursor` = "#76eff2", AMPlify = "#f2d676", AMPScannerV2 = "#b976f2", CSAMPPred = "#76f2be", 
      dbAMP = "#f276e8", GabereNoble = "#7688f2", iAMP2L = "#eff275", Wang = "#80f276", Witten = "#ccf276")
  ),
  tar_target(
    aa_comp_peptides_with_names,
    calculate_aa_comp_peptides(methods, 5, paste0(data_path, "Datasets/"))
  ),
  tar_target(
    aa_comp_peptides,
    select(aa_comp_peptides_with_names, -prot)
  ),
  tar_target(
    aa_comp_peptides_positive,
    get_aa_comp_positive(positive_traintest)
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
  )
)

