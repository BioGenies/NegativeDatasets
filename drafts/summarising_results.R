library(dplyr)
library(mlr3measures)
library(tidyr)
library(ggplot2)
library(hmeasure)
library(pbapply)
library(biogram)
library(scales)

results_path <- "/media/kasia/Data/Dropbox/Projekty/BioNgramProjects/NegativeDatasets/Results/"


if(Sys.info()[["nodename"]] == "huawei")
  results_path <- "/home/michal/Dropbox/BioNgramProjects/NegativeDatasets/Results/"

architectures <- list.files(results_path)

methods <- c("iAMP2L", "dbAMP", "ampir-mature", "CSAMPPred", 
             "Wang", "AmpGram", "Witten", "AMPScannerV2", "GabereNoble", "AMAP", "AMPlify")
methods_seqnames <- c("iAMP-2L", "dbAMP", "ampir-mature", "CS-AMPPred",
                      "Wang", "AmpGram", "Witten&Witten", "AMPScannerV2", "Gabere&Noble", "AMAP", "AMPlify")

change_method_names <- function(df) {
  mutate(df, method = case_when(method == "iAMP2L" ~ "iAMP-2L",
                                method == "CSAMPPred" ~ "CS-AMPPred",
                                method == "Wang" ~ "Wang et al.",
                                method == "Witten" ~ "Witten&Witten",
                                method == "GabereNoble" ~ "Gabere&Noble",
                                method %in% c("dbAMP", "ampir-precursor", "ampir-mature", "AmpGram", "AMAP", "AMPlify", "AMPScannerV2", "Positive") ~ method))
}

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
      
      res %>% 
        mutate(target = ifelse(grepl("AMP=1", ID), 1, 0))
    }) %>% 
      bind_rows()
  }) %>% 
    bind_rows()
})%>% 
  bind_rows() %>% 
  change_method_names()


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
  mutate(seq_source = sapply(all_results[["ID"]], function(i) gsub("method=", "", strsplit(i, "_")[[1]][3]))) %>% 
  filter(!(seq_source %in% c("AmPEP", "ampir-precursor"))) %>% 
  mutate(seq_source = ifelse(seq_source == "Wang-et-al", "Wang et. al", seq_source))

detailed_stats <- lapply(architectures, function(ith_architecture) {
  lapply(unique(seqtype_all_results[["method"]]), function(ith_method) {
    lapply(1:5, function(ith_rep) {
      lapply(unique(seqtype_all_results[["seq_source"]])[2:12], function(ith_seq_source) {
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
  coord_equal() +
  theme_bw() +
  xlim(c(0.5, 1)) +
  ylim(c(0.5, 1))

png("./drafts/reference_vs_nonreference.png", width = 800, height = 800)
p
dev.off()

# reference - nonreference for different train datasets
inner_join(reference_auc_df %>% 
             group_by(architecture) %>% 
             summarise(reference_mean_AUC = mean(reference_AUC)),
           detailed_stats_mean %>%
             filter(method != seq_source) %>% 
             select(architecture, method, seq_source, nonreference_AUC = mean_AUC) %>% 
             group_by(architecture, method) %>% 
             summarise(nonreference_mean_AUC = mean(nonreference_AUC))) %>% 
  ggplot(aes(x = reference_mean_AUC, y = nonreference_mean_AUC,
             color = architecture, label = architecture)) +
  geom_point() +
  facet_wrap(~method) +
  geom_abline(slope = 1, intercept = 0) +
  geom_label_repel() +
  coord_equal() +
  theme_bw() +
  xlim(c(0.5, 1)) +
  ylim(c(0.5, 1))
  


group_by(detailed_stats, architecture, seq_source) %>% 
  summarise(mean_AUC = mean(AUC)) %>% 
  ggplot(aes(x = seq_source, y = mean_AUC)) +
  geom_col() +
  facet_wrap(~architecture) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90), legend.position = "bottom")



# Pairwise comparisons of training and testing datasets vs. AUC
# Sequence length
seq_path <- "/media/kasia/Data/Dropbox/Projekty/BioNgramProjects/NegativeDatasets/Datasets/"
test_res_length <- lapply(unique(detailed_stats_mean[["architecture"]]), function(ith_architecture) {
  lapply(1:5, function(ith_rep) {
    lapply(methods, function(ith_method) {
      train <- read_fasta(paste0(seq_path, "Training_method_", ith_method, "_rep", ith_rep, ".fasta"))
      test <- read_fasta(paste0(seq_path, "Benchmark_rep", ith_rep, ".fasta"))
      train_filtered <- train[which(grepl("AMP=0", names(train)))]
      lapply(methods_seqnames, function(ith_seqtype) {
        test_filtered <- test[which(grepl(ith_seqtype, names(test)))]
        data.frame(architecture = ith_architecture,
                   method = ith_method,
                   seq_source = ith_seqtype,
                   pval = ks.test(x = lengths(train_filtered), y = lengths(test_filtered))[["p.value"]])
      }) %>% bind_rows()
    }) %>% bind_rows()
  }) %>% bind_rows()
}) %>% bind_rows()

ihs <- function(x) log(x + sqrt(x^2+1))

inner_join(test_res_length, detailed_stats) %>% 
  ggplot(aes(x = ihs(pval), y = AUC, color = architecture)) +
  geom_jitter(size = 3, alpha = 0.1) 

inner_join(test_res_length, detailed_stats) %>% 
  ggplot(aes(x = pval, y = AUC, color = architecture)) +
  geom_jitter(size = 3, alpha = 0.1) +
  facet_wrap(~architecture) 


# Distribution of replication mean AUC 
detailed_stats_mean %>% 
  ggplot(aes(y = mean_AUC, x = architecture, color = architecture)) +
  geom_violin(aes(fill = architecture), alpha = 0.25) +
  geom_point() +
  facet_wrap(~method)

### Fractions of problematic sequences

seqtype_all_results %>% 
  filter(target == 0) %>% 
  group_by(rep, ID, seq_source) %>% 
  summarise(mean_pred = mean(prediction)) %>% 
  group_by(mean_pred, seq_source) %>% 
  summarise(n = length(mean_pred)) %>% 
  ungroup() %>% 
  ggplot(aes(x = mean_pred, y = n)) +
  geom_col() +
  facet_wrap(~seq_source, ncol = 2, scales = "free_y")

seqtype_all_results %>% 
  filter(target == 0) %>% 
  group_by(rep, ID, seq_source, method) %>% 
  summarise(mean_pred = mean(prediction)) %>% 
  group_by(mean_pred, seq_source, method) %>% 
  summarise(n = length(mean_pred)) %>% 
  ungroup() %>% 
  ggplot(aes(x = mean_pred, y = n)) +
  geom_col() +
  facet_grid(seq_source~method, scales = "free_y")

### Problematic sequences vs. AUC
# on test dataset
mean_architecture_method_independent_auc <- detailed_stats_mean %>% 
  group_by(seq_source) %>% 
  summarise(mean_AUC = mean(mean_AUC))

seqtype_all_results %>% 
  filter(target == 0) %>% 
  group_by(rep, ID, seq_source) %>% 
  summarise(mean_pred = mean(prediction)) %>% 
  group_by(mean_pred, seq_source) %>% 
  summarise(n = length(mean_pred)) %>% 
  ungroup() %>% 
  left_join(mean_architecture_method_independent_auc) %>% 
  ggplot(aes(x = mean_pred, y = n, fill = mean_AUC)) +
  geom_col() +
  scale_fill_gradientn(colors = c("#ffe96b", "#ff4242", "#630000"), values = rescale(c(0, 0.08, 0.14), to = c(0, 1))) +
  theme_bw() +
  facet_wrap(~seq_source, ncol = 3, scales = "free_y")

mean_architecture_independent_auc <- detailed_stats_mean %>% 
  group_by(seq_source, method) %>% 
  summarise(mean_AUC = mean(mean_AUC))

seqtype_all_results %>% 
  filter(target == 0) %>% 
  group_by(rep, ID, seq_source, method) %>% 
  summarise(mean_pred = mean(prediction)) %>% 
  group_by(mean_pred, seq_source, method) %>% 
  summarise(n = length(mean_pred)) %>% 
  ungroup() %>% 
  left_join(mean_architecture_independent_auc, by = c("seq_source", "method")) %>% 
  ggplot(aes(x = mean_pred, y = n, fill = mean_AUC)) +
  geom_col() +
  scale_fill_gradientn(colors = c("#ffe96b", "#ff4242", "#630000"), values = rescale(c(0, 0.08, 0.14), to = c(0, 1))) +
  theme_bw() +
  facet_grid(seq_source~method, scales = "free_y") +
  ggtitle("Training dataset sampling method (top), testing dataset sampling method (right)")

# on methods
mean_seqsource_independent_auc <- detailed_stats_mean %>% 
  group_by(method) %>% 
  summarise(mean_AUC = mean(mean_AUC))

seqtype_all_results %>% 
  filter(target == 0) %>% 
  group_by(rep, ID, method) %>% 
  summarise(mean_pred = mean(prediction)) %>% 
  group_by(mean_pred, method) %>% 
  summarise(n = length(mean_pred)) %>% 
  ungroup() %>% 
  left_join(mean_seqsource_independent_auc) %>% 
  ggplot(aes(x = mean_pred, y = n, fill = mean_AUC)) +
  geom_col() +
  scale_fill_gradientn(colors = c("#ffe96b", "#ff4242", "#630000"), values = rescale(c(0, 0.08, 0.14), to = c(0, 1))) +
  theme_bw() +
  facet_wrap(~method, ncol = 3, scales = "free_y")

# problematic sequences vs. AUC
seqtype_all_results %>% 
  filter(target == 0) %>% 
  group_by(rep, ID, method) %>% 
  summarise(mean_pred = mean(prediction)) %>% 
  group_by(mean_pred, method) %>% 
  summarise(n = length(mean_pred)) %>% 
  ungroup() %>% 
 # filter(mean_pred <= 8/11) %>% 
  group_by(method) %>% 
  summarise(frac_problematic_seqs = sum(mean_pred <= 8/10)/n()) %>% 
  left_join(mean_seqsource_independent_auc) %>% 
  ggplot(aes(x = mean_AUC, y = frac_problematic_seqs, color = method, label = method)) +
  geom_point() +
  geom_label_repel() +
  theme_bw() +
  ylab("Fraction of problematic sequences") +
  ggtitle("Sequence is problematic if predicted correctly by max 8 architectures")


### AUC vs. differences in amino acid composition

targets::tar_load(aa_comp_all)
data_path <- "/media/kasia/Data/Dropbox/Projekty/BioNgramProjects/NegativeDatasets/Datasets/"
aa_comp_traintest <- lapply(1:5, function(j) {
  ds <- read_fasta(paste0(data_path, "Benchmark_rep", j, ".fasta"))
  pos <- unlist(ds[which(grepl("AMP=1", names(ds)))], use.names = FALSE)
  neg_df <- lapply(1:length(methods_seqnames), function(ith_method) {
    aa <- unlist(ds[which(grepl(methods_seqnames[ith_method], names(ds), fixed = TRUE))], use.names = FALSE)
    as.data.frame(table(aa)/length(aa)) %>% 
      mutate(method = ifelse(methods_seqnames[ith_method] == "Wang", "Wang et al.", methods_seqnames[ith_method]))
  }) %>% bind_rows() %>% 
    mutate(rep = as.character(j),
           Dataset = "Negative")
  pos_df <- mutate(setNames(as.data.frame(table(pos)/length(pos)), c("aa", "Freq")),
                        method = "Positive", rep = as.character(j), Dataset = "Positive")
  bind_rows(neg_df, pos_df)
}) %>% bind_rows() %>% 
  mutate(Type = "Testing") %>% 
  bind_rows(change_method_names(mutate(filter(aa_comp_all, !(method %in% c("AmPEP", "ampir-precursor"))), Type = "Training")))

# Negative train vs positive train
dat <- filter(aa_comp_traintest, Type == "Training" & Dataset == "Negative")
neg_pos_diffs <- lapply(1:5, function(ith_rep) {
    lapply(unique(dat[["method"]]), function(ith_method) {
      lapply(unique(dat[["aa"]]), function(ith_aa) {
        neg <- filter(dat, method == ith_method, rep == ith_rep, aa == ith_aa)[["Freq"]]
        pos <- filter(aa_comp_traintest, Type == "Training", Dataset == "Positive", aa == ith_aa)[["Freq"]]
        n <- ifelse(length(neg) == 0, 0, neg)
        p <- ifelse(length(pos) == 0, 0, pos)
        data.frame(method = ith_method,
                   aa = ith_aa,
                   rep = ith_rep,
                   diff = abs(n-p),
                   norm_diff = abs(n-p)/p,
                   square_diff = (n-p)^2,
                   norm_square_diff = ((n-p)^2)/p)
      }) %>% bind_rows()
    }) %>% bind_rows()
  }) %>% bind_rows()

# adding mean AUC for a given method, replication and architecture (mean AUC from results for different test datasets)
neg_pos_diffs_plot_dat <- detailed_stats %>% 
  group_by(method, rep, architecture) %>% 
  summarise(mean_AUC = mean(AUC)) %>% 
  left_join(neg_pos_diffs, by = c("method", "rep"))

# Differences for each amino acid vs. mean AUC 
ggplot(neg_pos_diffs_plot_dat, aes(x = diff, y = mean_AUC, color = method)) +
         geom_point() +
         facet_wrap(~aa) +
  scale_color_brewer(palette = "Paired")

# Aggregated differences vs. AUC
neg_pos_diffs %>% 
  group_by(method, rep) %>% 
  summarise(sum_diff = sum(diff)) %>% 
  left_join(., summarise(
    group_by(detailed_stats, method, rep, architecture),
    mean_AUC = mean(AUC))) %>% 
  ggplot(aes(x = sum_diff, y = mean_AUC, color = method)) +
  geom_point() +
  scale_color_brewer(palette = "Paired")

# Distance 
neg_pos_diffs %>% 
  group_by(method, rep) %>% 
  summarise(distance = sqrt(sum(square_diff))) %>% 
  left_join(., summarise(
    group_by(detailed_stats, method, rep, architecture),
    mean_AUC = mean(AUC))) %>% 
  ggplot(aes(x = distance, y = mean_AUC, color = method)) +
  geom_point() +
  scale_color_brewer(palette = "Paired")

# Negative test vs negative train
test_train_diffs <- lapply(1:5, function(ith_rep) {
  lapply(unique(aa_comp_traintest[["method"]])[1:11], function(ith_method) {
    lapply(unique(dat[["aa"]]), function(ith_aa) {
      train <- filter(aa_comp_traintest, method == ith_method, rep == ith_rep, aa == ith_aa, Type == "Training")[["Freq"]]
      test <- filter(aa_comp_traintest, method == ith_method, rep == ith_rep, aa == ith_aa, Type == "Testing")[["Freq"]]
      tr <- ifelse(length(train) == 0, 0, train)
      te <- ifelse(length(test) == 0, 0, test)
      data.frame(method = ith_method,
                 aa = ith_aa,
                 rep = ith_rep,
                 diff = abs(te-tr),
                 norm_diff = abs(te-tr)/tr,
                 square_diff = (te-tr)^2,
                 norm_square_diff = ((te-tr)^2)/tr)
    }) %>% bind_rows()
  }) %>% bind_rows()
}) %>% bind_rows()

test_train_diffs_plot_dat <- detailed_stats %>% 
  group_by(method, rep, architecture) %>% 
  summarise(mean_AUC = mean(AUC)) %>% 
  left_join(test_train_diffs, by = c("method", "rep"))

# Differences for each amino acid
ggplot(test_train_diffs_plot_dat, aes(x = diff, y = mean_AUC, color = method)) +
  geom_point() +
  facet_wrap(~aa) +
  scale_color_brewer(palette = "Paired")

# Aggregated differences
test_train_diffs %>% 
  mutate(norm_diff = ifelse(is.na(norm_diff), 0, norm_diff)) %>% 
  group_by(method, rep) %>% 
  summarise(sum_norm_diff = sum(norm_diff)) %>% 
  mutate(sum_diff = ifelse(is.na(sum_norm_diff), 0, sum_norm_diff)) %>% 
  left_join(., summarise(
    group_by(detailed_stats, method, rep, architecture),
    mean_AUC = mean(AUC))) %>% 
  ggplot(aes(x = sum_norm_diff, y = mean_AUC, color = method)) +
  geom_point() +
  scale_color_brewer(palette = "Paired")

# Distance
test_train_diffs %>% 
  mutate(norm_square_diff = ifelse(is.na(norm_square_diff), 0, norm_square_diff)) %>% 
  group_by(method, rep) %>% 
  summarise(norm_distance = sqrt(sum(norm_square_diff))) %>% 
  left_join(., summarise(
    group_by(detailed_stats, method, rep, architecture),
    mean_AUC = mean(AUC))) %>% 
  ggplot(aes(x = norm_distance, y = mean_AUC, color = method)) +
  geom_point() +
  scale_color_brewer(palette = "Paired")
