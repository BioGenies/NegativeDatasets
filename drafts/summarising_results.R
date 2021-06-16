library(dplyr)
library(mlr3measures)
library(tidyr)
library(ggplot2)
library(hmeasure)
library(pbapply)

results_path <- "/media/kasia/Data/Dropbox/Projekty/BioNgramProjects/NegativeDatasets/Results/"


if(Sys.info()[["nodename"]] == "huawei")
  results_path <- "/home/michal/Dropbox/BioNgramProjects/NegativeDatasets/Results/"

architectures <- list.files(results_path)

methods <- c("iAMP2L", "AmPEP", "dbAMP", "ampir-precursor", "ampir-mature", "CSAMPPred", 
             "Wang", "AmpGram", "Witten", "AMPScannerV2", "GabereNoble", "AMAP", "AMPlify")
methods_seqammes <- c("iAMP-2L", "AmPEP", "dbAMP", "ampir-precursor", "ampir-mature", "CS-AMPPred",
                      "Wang", "AmpGram", "Witten&Witten", "AMPScannerV2", "Gabere&Noble", "AMAP", "AMPlify")



# Results on whole benchmark datasets
all_results <- pblapply(architectures, function(ith_architecture) {
  lapply(methods, function(ith_method) {
    lapply(1:5, function(ith_rep) {
      res <- try(data.table::fread(file = paste0(results_path, ith_architecture, "/training_method_", ith_method, "_rep", ith_rep, ".csv"), data.table = FALSE) %>% 
                   mutate(rep = ith_rep,
                          method = ith_method,
                          architecture = ith_architecture), silent = TRUE)
      
      if(inherits(res, "try-error"))
        res <- data.frame(ID = c(), target = c(), prediction = c(), probability = c(), rep = c(), method = c(), 
                          architecture = c())
      
      res
    }) %>% 
      bind_rows()
  }) %>% 
    bind_rows()
})%>% 
  bind_rows()


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


detailed_stats <- readRDS("./drafts/detailed_stats.rds")

detailed_stats_mean <- detailed_stats %>% 
  group_by(architecture, method, seq_source) %>% 
  summarise(mean_AUC = mean(AUC),
            sd = sd(AUC))

pdf("./drafts/architecture-benchmark.pdf", height = 10, width = 9)
ggplot(detailed_stats_mean, aes(x = method, y = seq_source, fill = mean_AUC)) +
  geom_tile() +
  geom_point(data = detailed_stats_mean, aes(x = method, y = seq_source, size = sd)) +
  facet_wrap(~architecture, ncol = 3) +
  scale_fill_gradient("Mean AUC", low =  "#ffe96b",  high = "#ff4242",
                       trans = scales::trans_new("square_exp", function(x) exp(x)^2, function(x) log(sqrt(x)))) +
  scale_size_continuous("Standard deviation") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90), legend.position = "bottom") +
  xlab("Sampling method used for generation of training negative dataset") +
  ylab("Sampling method used for generation of test negative dataset")
dev.off()

filter(detailed_stats, architecture == "Deep-AmPEP30",
       method %in% c("iAMP2L", "Wang")) %>% 
  ggplot(aes(x = seq_source, y = AUC)) +
  geom_point() +
  facet_wrap(~ method) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90), legend.position = "bottom")


filter(detailed_stats, architecture == "Deep-AmPEP30") %>% 
  ggplot(aes(x = seq_source, y = AUC)) +
  geom_point() +
  facet_wrap(~ method) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90), legend.position = "bottom")


library(ggbeeswarm)
filter(detailed_stats_mean, architecture == "Deep-AmPEP30") %>% 
  ggplot(aes(x = seq_source, y = mean_AUC, color = method)) +
  geom_quasirandom(method = "smiley") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90), legend.position = "bottom")

