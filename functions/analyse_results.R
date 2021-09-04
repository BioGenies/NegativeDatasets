aggregate_all_results <- function(results_path, architectures, methods) {
  pblapply(architectures, function(ith_architecture) {
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
}

get_seq_source <- function(all_results) {
  all_results %>% 
    mutate(seq_source = sapply(all_results[["ID"]], function(i) gsub("method=", "", strsplit(i, "_")[[1]][3]))) %>% 
    filter(!(seq_source %in% c("AmPEP", "ampir-precursor"))) %>% 
    mutate(seq_source = ifelse(seq_source == "Wang-et-al", "Wang et al.", seq_source))
}

get_detailed_stats <- function(seqtype_all_results, architectures) {
  lapply(architectures, function(ith_architecture) {
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
}

get_detailed_stats_mean <- function(detailed_stats) {
  detailed_stats %>% 
    dplyr::group_by(architecture, method, seq_source) %>% 
    dplyr::summarise(mean_AUC = mean(AUC),
                     sd = sd(AUC)) %>% 
    ungroup()
}