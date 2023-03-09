#' Script for supplement figures

library(tidyverse)
library(cowplot)
source(here::here("analysis/00-metadata.R"))

pairs <- read_csv(paste0(folder_data, "output/pairs.csv"), show_col_types = F)
pairs_freq <- read_csv(paste0(folder_data, "temp/93-pairs_freq.csv"), show_col_types = F)
load(paste0(folder_data, "temp/95-communities_network.Rdata"))
communities_hierarchy <- read_csv(paste0(folder_data, "temp/95-communities_hierarchy.csv"), show_col_types = F)

# Clean up the pairs data ----
pairs <- pairs %>%
    # Remove no-colony pairs
    unite(col = "Pair", Community, Isolate1, Isolate2, sep = "_", remove = F) %>%
    filter(!(Pair %in% pairs_no_colony)) %>%
    # Remove low-accuracy model pairs
    filter(AccuracyMean > 0.9)


# Figure S5 machine vs. human ----
pairs_to_include <- pairs %>%
    select(Community, Isolate1, Isolate2) %>%
    mutate(Include = T)
pairs_T8_combined <- read_csv(paste0(folder_data, "temp/92-pairs_T8_combined.csv"), show_col_types = F) %>%
    left_join(pairs_to_include) %>%
    filter(!is.na(Include))

p1 <- pairs_T8_combined %>%
    ggplot() +
    geom_abline(slope = 1, intercept = 0, color = "red", linetype = 2) +
    geom_point(aes(x = TotalCount_human, y = TotalCount_machine), shape = 21, size = 2) +
    geom_text(x = -Inf, y = Inf, label = paste0("N=", nrow(pairs_T8_combined)), vjust = 2, hjust = -1) +
    scale_x_log10() +
    scale_y_log10() +
    theme_classic() +
    labs(x = "Segmentation CFU count", y = "Manual CFU count")


p2 <- pairs_T8_combined %>%
    filter(!is.na(Isolate1Count_human)) %>%
    ggplot() +
    geom_abline(slope = 1, intercept = 0, color = "red", linetype = 2) +
    geom_hline(yintercept = c(0,1), color = gray(.8), linetype = 2) +
    geom_vline(xintercept = c(0,1), color = gray(.8), linetype = 2) +
    geom_point(aes(x = Isolate1CFUFreq_human, y = Isolate1CFUFreq_machine), shape = 21, size = 2, stroke = .4) +
    geom_text(x = 0.5, y = 0.9, label = paste0("N=", pairs_T8_combined %>% filter(!is.na(Isolate1Count_human)) %>% nrow())) +
    scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 21)) +
    theme_classic() +
    labs(x = "Random Forest CFU frequency", y = "Manual CFU frequency")

p <- plot_grid(p1, p2, nrow = 1, axis = "tblr", align = "h", scale = .9, labels = c("A", "B")) +
    paint_white_background()

ggsave(here::here("plots/FigS5-human_machine_comparison.png"), p, width = 8, height = 4)

## cocultures without human results
pairs_T8_combined %>%
    filter(is.na(Isolate1Count_human)) %>%
    nrow()
## R-squared
pairs_T8_combined %>%
    filter(!is.na(Isolate1Count_human)) %>%
    lm(TotalCount_human ~ TotalCount_machine, data = .) %>%
    summary()
pairs_T8_combined %>%
    filter(!is.na(Isolate1Count_human)) %>%
    lm(Isolate1CFUFreq_human ~ Isolate1CFUFreq_machine, data = .) %>%
    summary()
