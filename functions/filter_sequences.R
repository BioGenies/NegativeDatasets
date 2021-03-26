filter_nonstandard_aa <- function(sequences) {
  standard <- toupper(biogram:::return_elements(seq_type = "prot"))
  is_standard <- vapply(sequences, function(seq) all(seq %in% standard), c(FALSE))
  sequences[is_standard]
}

filter_by_annotations <- function(sequences, uniprot_data, keywords_vec, exclude = TRUE) {
  has_annot <- sapply(keywords_vec, function(x) uniprot_data[["Entry"]][which(grepl(x, uniprot_data[["Keywords"]]))]) %>% 
    unlist() %>% 
    unique()
  if(exclude == TRUE) {
    sequences[which(!(names(sequences) %in% has_annot))]
  } else {
    sequences[which(names(sequences) %in% has_annot)]
  }
}

filter_by_lengths <- function(sequences, min_len, max_len) {
  sequences[which(lengths(sequences) >= min_len & lengths(sequences) <= max_len)]
}

# From CD-HIT manual
# Choose of word size:
# -n 5 for thresholds 0.7 ~ 1.0
# -n 4 for thresholds 0.6 ~ 0.7
# -n 3 for thresholds 0.5 ~ 0.6
# -n 2 for thresholds 0.4 ~ 0.5
filter_with_cdhit <- function(sequences, threshold, word_length = 2, cdhit_path = "./third-party", l = 10, n_threads = 1) {
  input <- tempfile(tmpdir = getwd())
  output <- tempfile(tmpdir = getwd())
  cdhit <- paste0(cdhit_path, "/cdhit -i ", input,  " -o ", output, " -c ", threshold, " -n ", word_length, " -l ", l, " -T ", n_threads)
  write_fasta(sequences, input)
  system(cdhit)
  res <- read_fasta(output)
  file.remove(input, output, paste0(output, ".clstr"))
  #names(res)
  res
}

filter_with_blat <- function(sequences, database, blat_path = "./third-party") {
  input <- tempfile(tmpdir = getwd())
  output <- tempfile(tmpdir = getwd())
  database_file <- tempfile(tmpdir = getwd())
  blat <- paste0(blat_path, "/blat ", database_file, " ", input, " -prot -out=blast8 ", output)
  write_fasta(sequences, input)
  write_fasta(database, database_file)
  system(blat)
  res <- unique(read.delim(output, header = FALSE)[["V1"]])
  file.remove(input, output, database_file)
  sequences[which(!(names(sequences) %in% res))]
}


generate_cutted_sequences <- function(sequences, lens, excluded_aa = NULL) {
  if(is.null(excluded_aa)) {
    seq_vec <- unlist(sequences, use.names = FALSE)
  } else {
    seq_vec <- unlist(sequences, use.names = FALSE)
    seq_vec <- seq_vec[which(seq_vec != excluded_aa)]
  }

  end_vec <- unname(cumsum(sample(lens)))
  start_vec <- c(1, end_vec[-length(end_vec)] + 1)
  
  seq_parts <- floor(1L:length(seq_vec)/max(end_vec)) + 1
  max_part <- length(seq_vec)%/%max(end_vec) 
  
  splits <- split(seq_vec, seq_parts)
  pos_df <- data.frame(start = start_vec, end = end_vec, part = sample(1L:max_part, size = length(lens), replace = TRUE))
  cutted_sequences <- pblapply(1L:nrow(pos_df), function(ith_row) {
    unname(splits[[pos_df[ith_row, "part"]]][pos_df[ith_row, "start"]:pos_df[ith_row, "end"]])
  })
  
  names(cutted_sequences) <- paste0("CUTTED", 1L:length(cutted_sequences))
  cutted_sequences
}


filter_out_positive_sequences <- function(sequences, positive_dataset) {
  sequences[which(!(sequences %in% positive_dataset))]
}


filter_random_sequences <- function(sequences, number_of_sequences) {
  sample(sequences, number_of_sequences)
}


filter_by_location <- function(sequences, uniprot_data, location, evidence = NULL, exclude = FALSE) {
  has_annot <- uniprot_data[["Entry"]][which(grepl(location, uniprot_data[["Subcellular.location..CC."]]))] 
  if(!is.null(evidence)) {
    x <- filter(uniprot_data, Entry %in% has_annot)
    has_annot <- x[["Entry"]][which(grepl(evidence, x[["Subcellular.location..CC."]]))]
  }
  has_annot <- unique(unlist(has_annot))
  if(exclude == TRUE) {
    sequences[which(!(names(sequences) %in% has_annot))]
  } else {
    sequences[which(names(sequences) %in% has_annot)]
  }
}
