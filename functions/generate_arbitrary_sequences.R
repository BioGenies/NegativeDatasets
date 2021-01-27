generate_arbitrary_sequences <- function(lens) {
 lapply(lens, function(ith_len) {
   sample(toupper(biogram:::return_elements(seq_type = "prot")),
          ith_len, replace = TRUE)
 }) %>% setNames(paste0("ARBITRARY", 1:length(lens)))
}
