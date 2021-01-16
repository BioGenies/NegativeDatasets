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
filter_with_cdhit <- function(sequences, threshold, word_length = 2, cdhit_path = "./third-party") {
  input <- tempfile(tmpdir = getwd())
  output <- tempfile(tmpdir = getwd())
  cdhit <- paste0(cdhit_path, "/cdhit -i ", input,  " -o ", output, " -c ", threshold, " -n ", word_length)
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

