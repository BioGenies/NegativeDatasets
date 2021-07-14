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


part_dat <- filter(all_results, rep == 1) 

proper_pred_df <- group_by(part_dat, ID, target, method) %>% 
  summarise(mean_pred = mean(prediction))

group_by(proper_pred_df, target, mean_pred, method) %>% 
  summarise(n = length(mean_pred)) %>% 
  ungroup() %>% 
  group_by(target, mean_pred) %>% 
  mutate(n_frac = n/sum(n)) %>% 
  ggplot(aes(x = mean_pred, y = n_frac, fill = as.factor(target))) +
  geom_col(width = 0.05, position = "dodge") +
  facet_wrap( ~ method, labeller = label_both, nrow = 2)

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


detailed_stats <- readRDS("./drafts/detailed_stats.rds") %>%
  mutate(method = factor(method, levels = sort(unique(method))),
         seq_source = factor(seq_source, levels = sort(unique(seq_source)), labels = levels(method)))


detailed_stats_mean <- detailed_stats %>% 
  group_by(architecture, method, seq_source) %>% 
  summarise(mean_AUC = mean(AUC),
            sd = sd(AUC)) %>% 
  ungroup()

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
filter(detailed_stats_mean, architecture == "AMPScannerV2") %>% 
  ggplot(aes(x = seq_source, y = mean_AUC, color = method)) +
  geom_quasirandom(method = "smiley") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90), legend.position = "bottom")

filter(detailed_stats, architecture == "AMPScannerV2") %>% 
  ggplot(aes(x = seq_source, y = AUC)) + 
  geom_point() +
  facet_wrap(~method)

# method = train
# seq_source = test

reference_auc_df <- detailed_stats_mean %>%
  filter(method == seq_source) %>% 
  select(architecture, seq_source, reference_AUC = mean_AUC)

library(ggrepel)

pdf("./drafts/mad_mean_plot.pdf")
inner_join(detailed_stats_mean, reference_auc_df) %>% 
  mutate(delta_AUC = abs(reference_AUC - mean_AUC)) %>% 
  group_by(architecture) %>% 
  summarise(median_mean_AUC = median(mean_AUC),
            median_mean_AUC_min = median(mean_AUC) - mad(mean_AUC),
            median_mean_AUC_max = median(mean_AUC) + mad(mean_AUC),
            mad_delta_AUC = mad(delta_AUC)) %>% 
  ggplot(aes(x = mad_delta_AUC, y = median_mean_AUC, 
             ymin = median_mean_AUC_min, ymax = median_mean_AUC_max,
             color = architecture, label = architecture)) +
  geom_point(size = 3) +
  geom_errorbar() +
  geom_label_repel()
dev.off()

head(detailed_stats)

pdf("./drafts/sd_mean_plot.pdf")
group_by(detailed_stats, architecture, method, seq_source) %>% 
  summarise(mean_AUC = mean(AUC)) %>% 
  ungroup() %>% 
  group_by(architecture) %>% 
  summarise(sd_AUC = sd(mean_AUC), mean_AUC = mean(mean_AUC)) %>% 
  ggplot(aes(x = sd_AUC, y = mean_AUC,
             color = architecture, label = architecture)) +
  geom_point(size = 3) +
  geom_label_repel()
dev.off()

p <- inner_join(reference_auc_df %>% 
             group_by(architecture) %>% 
             summarise(reference_mean_AUC = mean(reference_AUC)),
           detailed_stats_mean %>%
             filter(method != seq_source) %>% 
             select(architecture, seq_source, nonreference_AUC = mean_AUC) %>% 
             group_by(architecture) %>% 
             summarise(nonreference_mean_AUC = mean(nonreference_AUC))) %>% 
  ggplot(aes(x = reference_mean_AUC, y = nonreference_mean_AUC,
             color = architecture, label = architecture)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0) +
  geom_label_repel() +
  coord_equal()

png("./drafts/reference_vs_nonreference.png", width = 800, height = 800)
p
dev.off()


group_by(detailed_stats, architecture, seq_source) %>% 
  summarise(mean_AUC = mean(AUC)) %>% 
  ggplot(aes(x = seq_source, y = mean_AUC)) +
  geom_col() +
  facet_wrap(~architecture) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90), legend.position = "bottom")
