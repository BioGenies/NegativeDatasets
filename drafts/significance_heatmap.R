library(biogram)
library(dplyr)
library(reshape2)
library(ggplot2)

set.seed(15390)


before_test <- lapply(toupper(colnames(aaprop)), function(ith_aa) {
  data.frame(aa = ith_aa, occ = runif(5), repl = paste0("repl", 1L:5))
}) %>% 
  bind_rows()


dat <- lapply(paste0("method", 1L:12), function(ith_method) {
  lapply(combn(paste0("repl", 1L:5), 2, simplify = FALSE), function(ith_repl) {
    lapply(toupper(colnames(aaprop)), function(ith_aa) {
      part_dat <- filter(before_test, aa == ith_aa, repl %in% ith_repl)
      data.frame(aa = ith_aa, 
                 occ1 = part_dat[1, "occ"], occ2 = part_dat[2, "occ"],
                 repl_id1 = part_dat[1, "repl"], repl_id2 = part_dat[2, "repl"],
                 pval = runif(1, max = 0.2))
    }) %>% bind_rows()
  }) %>% bind_rows() %>% 
    mutate(signif = pval < 0.05,
           method = ith_method)
}) %>% bind_rows()
  
  
# 
# select(dat, aa, repl_id1, repl_id2, pval) %>% 
#   ggplot(aes(x = repl_id1, y = repl_id2, fill = pval)) +
#   geom_tile() +
#   facet_wrap(~ aa)

select(dat, aa, repl_id1, repl_id2, pval, method) %>% 
  ggplot(aes(x = aa, y = paste0(repl_id1, repl_id2), fill = pval < 0.05)) +
  geom_tile() +
  facet_wrap(~ method)

select(dat, aa, repl_id1, repl_id2, pval, method) %>% 
  group_by(aa, method) %>% 
  summarise(n_signif = sum(pval < 0.05)) %>% 
  ggplot(aes(x = aa, y = method, fill = n_signif)) +
  geom_tile()
