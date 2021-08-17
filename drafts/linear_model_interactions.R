library(dplyr)

detailed_stats <- readRDS("./detailed_stats.rds")

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
  data.frame(
    model = i,
    adjusted_R_squared = summ[["adj.r.squared"]],
    n_significant_interactions = n_signif
  )
}) %>% bind_rows()

