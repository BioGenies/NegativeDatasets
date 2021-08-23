library(dplyr)
library(pbapply)
library(broom)

detailed_stats <- readRDS("./drafts/detailed_stats.rds")

models <- list("AUC ~ architecture",
               "AUC ~ method",
               "AUC ~ seq_source",
               "AUC ~ architecture + method + seq_source",
               "AUC ~ method*seq_source",
               "AUC ~ architecture + method + seq_source + method*seq_source",
               "AUC ~ architecture + method + seq_source + architecture*method + architecture*seq_source + method*seq_source",
               "AUC ~ architecture*method*seq_source")

lm_res <- lapply(models, function(i) {
  lm_model <- lm(eval(parse(text = i)), data = detailed_stats)
  summ <- summary(lm_model)
  n_signif <- data.frame(summ[["coefficients"]]) %>% 
    filter(grepl(":", rownames(.)) & `Pr...t..` < 0.05) %>% 
    nrow()
  n_int <- data.frame(summ[["coefficients"]]) %>% 
    filter(grepl(":", rownames(.))) %>% 
    nrow()
  data.frame(
    model = i,
    adjusted_R_squared = summ[["adj.r.squared"]],
    n_interactions = n_int,
    n_significant_interactions = n_signif
  )
}) %>% bind_rows()

big_model <- lm(AUC ~ architecture*method*seq_source, data = detailed_stats)

smaller_model <- update(big_model, . ~ . -`architectureMLAMP:methodGabere&Noble:seq_sourceGabere&Noble`)


all_terms <- names(coef(big_model))

only_interactions <- all_terms[grepl(pattern = ":", all_terms, fixed = TRUE)]

all_interactions_list <- pblapply(only_interactions, function(ith_interaction) {
  interaction_term <- paste0(". ~ . -`", ith_interaction, "`")
  smaller_model <- update(big_model, interaction_term)
  test_res <- anova(big_model, smaller_model)
  tidy(test_res) %>% 
    mutate(term = ith_interaction)
})

