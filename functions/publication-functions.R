library(dplyr)
library(ggplot2)
library(gargle)
library(googlesheets4)
library(xtable)
dat <- read_sheet("https://docs.google.com/spreadsheets/d/1kKV11R3LOF7-vpeaCc4LgmGBDo3UaqLin_dc5Wl6GI8/edit?usp=sharing")

dat[, c("Software", "DOI", "Implementation of the negative dataset sampling algorithm?", 
        "Negative dataset sampling algorithm comment", "Implementation of the model architecture?", 
        "Model architecture comment")] %>% 
  xtable() %>% 
  print.xtable() %>% 
  cat(sep = "\n", file = "./publication-results/sampling_architectures.tex", append = FALSE)

