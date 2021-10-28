generate_holdout_groups <- function(sequences, seed) {
  set.seed(seed = seed)
  seq_length_groups <- cut(lengths(sequences), 
                           breaks = as.numeric(quantile(lengths(sequences), probs = seq(0, 1, 0.2))),
                           include.lowest = TRUE)
  
  names(seq_length_groups) <- names(sequences)
  
  holdout_list <- lapply(levels(seq_length_groups), function(ith_group) {
    peptides_in_group <- names(seq_length_groups)[seq_length_groups == ith_group]
    group_benchmark <- sample(peptides_in_group, round(length(peptides_in_group)*0.20, 0))
    list(benchmark = group_benchmark,
         traintest = setdiff(peptides_in_group, group_benchmark))
  }) 
  
  names(holdout_list) <- levels(seq_length_groups)
  holdout_list
}