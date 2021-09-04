library(dplyr)
library(tidyr)
library(ggplot2)

square_bracketize <- function(x)
  paste0("[", x, "]")

all_lines <- readLines("./data/amp-software.txt")

year <- strsplit(all_lines, "):", fixed = TRUE) %>% 
  sapply(first) %>% 
  strsplit(" (", fixed = TRUE) %>% 
  sapply(last)
  

field_names <- c("repository", "web server", "standalone software")

field_df <- lapply(field_names, function(ith_field) {
  grepl(square_bracketize(ith_field), all_lines, fixed = TRUE)
}) %>% 
  setNames(gsub(" ", "", field_names)) %>% 
  do.call(cbind, .)

pub_plot <- data.frame(year = year) %>% 
  cbind(field_df) %>% 
  group_by(year) %>% 
  arrange(year) %>% 
  summarise(n_repository = sum(repository),
            n_ws = sum(webserver),
            n_ss = sum(standalonesoftware),
            n_published = length(year)) %>% 
  mutate(n_model = n_ws + n_ss) %>% 
  select(-n_ws, -n_ss) %>% 
  pivot_longer(cols = c(n_repository, n_model, n_published)) %>% 
  mutate(name = factor(name, 
                       levels = c("n_published", "n_model", "n_repository"),
                       labels = c("Published models", "Available models", 
                                  "Code repositories"))) %>% 
  ggplot(aes(x = year, y = value, label = value)) +
  geom_col() +
  geom_text(vjust = -0.5, size = 6) +
  coord_cartesian(ylim = c(0, 11)) +
  scale_x_discrete("Year") +
  scale_y_continuous("Value") +
  facet_wrap(~ name, nrow = 3) +
  theme_bw(base_size = 16)

png(filename = "./reports/pub-plot.png", width = 500, height = 700)
pub_plot
dev.off()
