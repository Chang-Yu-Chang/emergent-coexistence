library(tidyverse)
library(cowplot)
source(here::here("analysis/00-metadata.R"))

accuracy <- read_csv(paste0(folder_data, "temp/91-accuracy.csv"), show_col_types = F)
accuracy_count <- accuracy %>%
    group_by(PairType) %>%
    count(name = "Count")

p1a <- accuracy %>%
    ggplot() +
    geom_histogram(aes(x = Accuracy), color = 1, binwidth = 0.01, breaks = seq(0.6,1,.01)) +
    geom_text(data = accuracy_count, x = -Inf, y = Inf, aes(label = paste0("N=",Count)), vjust = 2, hjust = -1) +
    geom_vline(xintercept = 0.9, color = "red", linetype = 2) +
    facet_grid(PairType~.) +
    theme_classic() +
    theme(panel.border = element_rect(color = 1, fill = NA)) +
    labs(x = "Accuracy", y = "Count")

p1b <- accuracy %>%
    ggplot() +
    geom_histogram(aes(x = Accuracy), color = 1, binwidth = 0.01, breaks = seq(0.6,1,.01)) +
    geom_text(data = accuracy_count, x = -Inf, y = Inf, aes(label = paste0("N=",Count)), vjust = 2, hjust = -1) +
    geom_vline(xintercept = 0.9, color = "red", linetype = 2) +
    scale_y_log10() +
    facet_grid(PairType~.) +
    theme_classic() +
    theme(panel.border = element_rect(color = 1, fill = NA)) +
    labs(x = "Accuracy", y = "Count")

p1 <- plot_grid(p1a, p1b, nrow = 1, axis = "tb", align = "h", scale = 0.9, labels = c("count", "log-scale")) + theme(plot.background = element_rect(fill = "white", color = NA))
ggsave(paste0(folder_data, "temp/24a-01-random_forest-accuracy.png"), p1, width = 6, height = 6)
