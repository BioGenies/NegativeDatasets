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
    get_positive_dataset("data/dbaasp.csv")),
  tar_target(
    uniprot_seqs,
    read_uniprot_seqs(uniprot_seqs_file)),
  tar_target(
    uniprot_tab,
    read.delim(uniprot_tab_file)),
  tar_target(
    putative_amps,
    get_putative_amps()),
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
  )
)

