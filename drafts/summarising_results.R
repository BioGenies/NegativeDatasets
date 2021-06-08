library(dplyr)
library(mlr3measures)
library(tidyr)
library(ggplot2)
library(hmeasure)

results_path <- "/media/kasia/Data/Dropbox/Projekty/BioNgramProjects/NegativeDatasets/Results/"

architectures <- list.files(results_path)

methods <- c("iAMP2L", "AmPEP", "dbAMP", "ampir-precursor", "ampir-mature", "CSAMPPred", 
             "Wang", "AmpGram", "Witten", "AMPScannerV2", "GabereNoble", "AMAP", "AMPlify")
methods_seqammes <- c("iAMP-2L", "AmPEP", "dbAMP", "ampir-precursor", "ampir-mature", "CS-AMPPred",
                      "Wang", "AmpGram", "Witten&Witten", "AMPScannerV2", "Gabere&Noble", "AMAP", "AMPlify")

# Results on whole benchmark datasets
all_results <- lapply(architectures, function(ith_architecture) {
  lapply(methods, function(ith_method) {
    lapply(1:5, function(ith_rep) {
      res <- read.csv(paste0(results_path, ith_architecture, "/training_method_", ith_method, "_rep", ith_rep, ".csv")) %>% 
        mutate(rep = ith_rep)
    }) %>% bind_rows() %>% 
      mutate(method = ith_method)
  }) %>% bind_rows %>% 
    mutate(architecture = ith_architecture)
}) %>% bind_rows() 

# stats <- all_results %>% 
#   group_by(architecture, method, rep) %>% 
#   summarise(AUC = auc(truth = as.factor(target), prob = probability, positive = "1"),
#             Accuracy = acc(truth = as.factor(target), response = as.factor(prediction)),
#             MCC = mcc(truth = as.factor(target), response = as.factor(prediction), positive = "1"))

# stats %>% 
#   select(-c(Accuracy, MCC)) %>% 
#   ggplot(aes(x = factor(method), y = AUC)) +
#   geom_jitter(width = 0.2) +
#   facet_wrap(~architecture) +
#   theme_bw() +
#   theme(axis.text.x = element_text(angle = 90))


# Results on each method separately 
seqtype_all_results <- all_results %>% 
  mutate(seq_source = sapply(all_results[["ID"]], function(i) gsub("method=", "", strsplit(i, "_")[[1]][3])))
detailed_stats <- lapply(architectures, function(ith_architecture) {
  lapply(methods, function(ith_method) {
    lapply(1:5, function(ith_rep) {
      lapply(unique(seqtype_all_results[["seq_source"]])[2:14], function(ith_seq_source) {
        dat <- filter(seqtype_all_results, architecture == ith_architecture, method == ith_method, 
                      seq_source %in% c(ith_seq_source, "AMP=1"), rep == ith_rep)
        data.frame(architecture = ith_architecture,
                   method = ith_method,
                   rep = ith_rep,
                   seq_source = ith_seq_source,
                   AUC = ifelse(all(is.na(dat[["probability"]])),
                                HMeasure(dat[["target"]], dat[["prediction"]])[["metrics"]][["AUC"]],
                                HMeasure(dat[["target"]], dat[["probability"]])[["metrics"]][["AUC"]]),
                   probs = ifelse(all(is.na(dat[["probability"]])), TRUE, FALSE))
      }) %>% bind_rows()
    }) %>% bind_rows()
  }) %>% bind_rows()
}) %>% bind_rows()


detailed_stats <- readRDS("detailed_stats.rds")

detailed_stats_mean <- detailed_stats %>% 
  group_by(architecture, method, seq_source) %>% 
  summarise(mean_AUC = mean(AUC),
            sd = sd(AUC))
ggplot(detailed_stats_mean, aes(x = method, y = seq_source, fill = mean_AUC)) +
  geom_tile() +
  geom_point(data = detailed_stats_mean, aes(x = method, y = seq_source, size = sd)) +
  facet_wrap(~architecture, ncol = 3) +
  scale_fill_gradient2("Mean AUC", low = "#ffffff", mid = "#ffe96b",  high = "#ff4242", midpoint = 0.5) +
  scale_size_continuous("Standard deviation") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90)) +
  xlab("Sampling method used for generation of training negative dataset") +
  ylab("Sampling method used for generation of test negative dataset")

