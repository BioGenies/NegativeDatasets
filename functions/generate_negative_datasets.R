# CS-AMPPred
generate_negative_dataset_1 <- function(sequences, uniprot_data, positive_dataset) {
  filter_random_sequences(
    filter_with_cdhit(
      filter_by_lengths(
        filter_by_annotations(sequences = sequences, 
                              uniprot_data = uniprot_data, 
                              keywords_vec = "Antimicrobial", 
                              exclude = TRUE),
        16, 90),
      0.4),
    length(positive_dataset))
}


# SVM-LZ 2nd dataset
generate_negative_dataset_2 <- function(sequences, positive_dataset) {
  filter_with_cdhit(
    filter_nonstandard_aa(
      generate_cutted_sequences(
        filter_random_sequences(
          filter_by_annotations(sequences = sequences, 
                                uniprot_data = uniprot_data, 
                                keywords_vec = c("Antimicrobial", "Secreted"),
                                exclude = TRUE),
          10000),
        sample(lengths(positive_dataset), 10000, replace = TRUE))), 
    0.7)
}


# MLAMP
generate_negative_dataset_3 <- function(sequences, uniprot_data) {
  filter_with_cdhit(
    filter_nonstandard_aa(
      filter_by_annotations(
        filter_by_lengths(sequences, 5, 100),
        uniprot_data, 
        c("Antimicrobial", "Antibiotic", "Fungicide", "Defensin"))),
    0.4)
}


# AmpGram
generate_negative_dataset_4 <- function(sequences, uniprot_data, positive_dataset) {
  generate_cutted_sequences(
    filter_nonstandard_aa(
      filter_by_annotations(sequences = sequences, 
                            uniprot_data = uniprot_data, 
                            keywords_vec = c("Antimicrobial", "Transit peptide", "Antibacterial", "Antiviral", "Fungicide", "Secreted"), 
                            exclude = TRUE)),
    lengths(positive_dataset))
}


# dbAMP
generate_negative_dataset_5 <- function(sequences, uniprot_data, positive_dataset) {
  filter_with_cdhit(
    filter_by_lengths(
      filter_by_annotations(sequences = sequences, 
                            uniprot_data = uniprot_data, 
                            keywords_vec = c("Antimicrobial", "Antibiotic", "Antiviral", "Fungicide", "Secreted", "Defensin", "Toxin", "Transmembrane"), 
                            exclude = TRUE),
      10, 100),
    0.4)
}


# Deep learning regression model for antimicrobial peptide design 
generate_negative_dataset_6 <- function(sequences, uniprot_data, positive_dataset) {
  generate_cutted_sequences(
    filter_with_cdhit(
      filter_by_location(
        filter_by_annotations(sequences = sequences, 
                              uniprot_data = uniprot_data, 
                              keywords_vec = c("Antimicrobial", "Antibiotic", "Antiviral", "Fungicide", "Secreted"), 
                              exclude = TRUE),
        uniprot_data, 
        location = "Cytoplasm", 
        evidence = "ECO:0000269"),
      0.4),
    lengths(positive_dataset), 
    excluded_aa = "C")
}


# ampir (precursor)
generate_negative_dataset_7 <- function(sequences, positive_dataset) {
  filter_random_sequences(
    filter_by_lengths(
      filter_nonstandard_aa(
        filter_out_positive_sequences(
          filter_with_cdhit(sequences, 0.9),
          positive_dataset)),
      51, 500),
    10*length(positive_dataset))
}


# ampir (mature)
generate_negative_dataset_8 <- function(sequences, positive_dataset) {
  filter_nonstandard_aa(
    filter_by_lengths(
      filter_out_positive_sequences(
        filter_with_cdhit(sequences, 0.9),
        positive_dataset),
      11, 39))
}


# AMPlify
generate_negative_dataset_9 <- function(sequences, uniprot_data, positive_dataset, potential_AMPs) {
  generate_cutted_sequences(
    filter_out_positive_sequences(
      filter_nonstandard_aa(
        filter_by_annotations(
          filter_by_lengths(sequences, 1, 200),
          uniprot_data,
          c("Antimicrobial", "Antibiotic", "defense", "Defensin", "Bacteriocin", "Fungicide"))), 
      c(positive_dataset, potential_AMPs)),
    lengths(positive_dataset))
}


# ACEP/AMPScanner V2
generate_negative_dataset_10 <- function(sequences, uniprot_data, positive_dataset) {
  generate_cutted_sequences(
    filter_with_blat(
      filter_by_lengths(
        filter_by_annotations(
          filter_by_location(sequences, uniprot_data, "Cytoplasm"),
          uniprot_dat,
          c("Antimicrobial", "Antiviral", "Antibiotic", "Fungicide", "Secreted")),
        10, max(lengths(sequences))),
      positive_dataset),
    lengths(positive_dataset))
}
