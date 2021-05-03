# iAMP-2L
generate_negative_dataset_iAMP2L <- function(sequences, uniprot_data) {
  seqs <- filter_with_cdhit(
    filter_nonstandard_aa(
      filter_by_annotations(
        filter_by_lengths(sequences, 5, 100),
        uniprot_data, 
        c("Antimicrobial", "Antibiotic", "Fungicide", "Defensin"))),
    0.4)
  names(seqs) <- paste0("Seq", 1:length(seqs), "_sampling_method=iAMP-2L_AMP=0")
  seqs
}

# AmPEP
generate_negative_dataset_AmPEP <- function(sequences, uniprot_data) {
  seqs <- filter_nonstandard_aa(
    filter_by_location(
      filter_by_annotations(
        filter_by_lengths(sequences, 5, 255),
        uniprot_data,
        c("Antimicrobial","Antibiotic", "Antiviral", "Fungicide", "Secreted", "Toxin", "Defensin", "defense")),
      uniprot_data, "Membrane|membrane", exclude = TRUE))
  names(seqs) <- paste0("Seq", 1:length(seqs), "_sampling_method=AmPEP_AMP=0")
  seqs
}

# ampir (precursor)
generate_negative_dataset_ampir_precursor <- function(sequences, positive_dataset, sequences_to_filter_out, n_threads = 1) {
  seqs <- filter_random_sequences(
    filter_by_lengths(
      filter_nonstandard_aa(
        filter_out_positive_sequences(
          filter_with_cdhit(sequences, 0.9, n_threads = n_threads),
          sequences_to_filter_out)),
      51, 500),
    10*length(positive_dataset))
  names(seqs) <- paste0("Seq", 1:length(seqs), "_sampling_method=ampir-precursor_AMP=0")
  seqs
}


# ampir (mature)
generate_negative_dataset_ampir_mature <- function(sequences, positive_dataset, n_threads = 1) {
  seqs <- filter_nonstandard_aa(
    filter_by_lengths(
      filter_out_positive_sequences(
        filter_with_cdhit(sequences, 0.9, n_threads = n_threads),
        positive_dataset),
      11, 39))
  names(seqs) <- paste0("Seq", 1:length(seqs), "_sampling_method=ampir-mature_AMP=0")
  seqs
}


# CS-AMPPred
generate_negative_dataset_CSAMPPred <- function(sequences, uniprot_data, positive_dataset) {
  seqs <- filter_random_sequences(
    filter_with_cdhit(
      filter_nonstandard_aa(
        filter_by_lengths(
          filter_by_annotations(sequences = sequences, 
                                uniprot_data = uniprot_data, 
                                keywords_vec = "Antimicrobial", 
                                exclude = TRUE),
          16, 90)),
      0.4),
    length(positive_dataset))
  names(seqs) <- paste0("Seq", 1:length(seqs), "_sampling_method=CS-AMPPred_AMP=0")
  seqs
}


# Wang et al.
generate_negative_dataset_Wang <- function(sequences, uniprot_data, positive_dataset) {
  seqs <- filter_with_cdhit(
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
  names(seqs) <- paste0("Seq", 1:length(seqs), "_sampling_method=Wang-et-al_AMP=0")
  seqs
}


# AmpGram
generate_negative_dataset_AmpGram <- function(sequences, uniprot_data, positive_dataset) {
  seqs <- generate_cutted_sequences(
    filter_nonstandard_aa(
      filter_by_annotations(sequences = sequences, 
                            uniprot_data = uniprot_data, 
                            keywords_vec = c("Antimicrobial", "Transit peptide", "Antibacterial", "Antiviral", "Fungicide", "Secreted"), 
                            exclude = TRUE)),
    lengths(positive_dataset))
  names(seqs) <- paste0("Seq", 1:length(seqs), "_sampling_method=AmpGram_AMP=0")
  seqs
}


# dbAMP
generate_negative_dataset_dbAMP <- function(sequences, uniprot_data) {
  seqs <- filter_with_cdhit(
    filter_nonstandard_aa(
      filter_by_lengths(
        filter_by_annotations(sequences = sequences, 
                              uniprot_data = uniprot_data, 
                              keywords_vec = c("Antimicrobial", "Antibiotic", "Antiviral", "Fungicide", "Secreted", "Defensin", "Toxin", "Transmembrane"), 
                              exclude = TRUE),
        10, 100)),
    0.4)
  names(seqs) <- paste0("Seq", 1:length(seqs), "_sampling_method=dbAMP_AMP=0")
  seqs
}


# Witten & Witten 
generate_negative_dataset_Witten <- function(sequences, uniprot_data, positive_dataset) {
  seqs <- generate_cutted_sequences(
    filter_with_cdhit(
      filter_nonstandard_aa(
        filter_by_location(
          filter_by_annotations(sequences = sequences, 
                                uniprot_data = uniprot_data, 
                                keywords_vec = c("Antimicrobial", "Antibiotic", "Antiviral", "Fungicide", "Secreted"), 
                                exclude = TRUE),
          uniprot_data, 
          location = "Cytoplasm", 
          evidence = "ECO:0000269")),
      0.4),
    lengths(positive_dataset), 
    excluded_aa = "C")
  names(seqs) <- paste0("Seq", 1:length(seqs), "_sampling_method=Witten&Witten_AMP=0")
  seqs
}


# AMPScanner V2
generate_negative_dataset_AMPScannerV2 <- function(sequences, uniprot_data, positive_dataset) {
  seqs <- generate_cutted_sequences(
    filter_with_blat(
      filter_nonstandard_aa(
        filter_by_lengths(
          filter_by_annotations(
            filter_by_location(sequences, uniprot_data, "Cytoplasm"),
            uniprot_data,
            c("Antimicrobial", "Antiviral", "Antibiotic", "Fungicide", "Secreted")),
          10, max(lengths(sequences)))),
      positive_dataset),
    lengths(positive_dataset))
  names(seqs) <- paste0("Seq", 1:length(seqs), "_sampling_method=AMPScannerV2_AMP=0")
  seqs
}


# Gabere&Noble
generate_negative_dataset_GabereNoble <- function(sequences, uniprot_data, positive_dataset) {
  seqs <- generate_cutted_sequences(
    filter_nonstandard_aa(
      filter_by_location(
        filter_by_annotations(sequences = sequences,
                              uniprot_data = uniprot_data,
                              keywords_vec = "Antimicrobial",
                              exclude = TRUE),
        uniprot_data = uniprot_data,
        location = "Cytoplasm|Endoplasmic reticulum|Mitochondrion|Golgi")),
    6*lengths(positive_dataset))
  names(seqs) <- paste0("Seq", 1:length(seqs), "_sampling_method=Gabere&Noble_AMP=0")
  seqs
}


# AMAP
generate_negative_dataset_AMAP <- function(sequences, uniprot_data, positive_dataset) {
  seqs <- filter_with_cdhit(
    generate_negative_dataset_10(sequences = sequences,
                                 uniprot_data = uniprot_data,
                                 positive_dataset = positive_dataset),
    0.4)
  names(seqs) <- paste0("Seq", 1:length(seqs), "_sampling_method=AMAP_AMP=0")
  seqs
}


# AMPlify
generate_negative_dataset_AMPlify <- function(sequences, uniprot_data, positive_dataset, potential_AMPs) {
  seqs <- generate_cutted_sequences(
    filter_out_positive_sequences(
      filter_nonstandard_aa(
        filter_by_annotations(
          filter_by_lengths(sequences, 1, 200),
          uniprot_data,
          c("Antimicrobial", "Antibiotic", "defense", "Defensin", "Bacteriocin", "Fungicide"))), 
      c(positive_dataset, potential_AMPs)),
    lengths(positive_dataset))
  names(seqs) <- paste0("Seq", 1:length(seqs), "_sampling_method=AMPlify_AMP=0")
  seqs
}

