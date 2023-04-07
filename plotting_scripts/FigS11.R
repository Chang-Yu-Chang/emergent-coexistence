library(tidyverse)
library(cowplot)
source(here::here("processing_scripts/00-metadata.R"))

isolates <- read_csv(paste0(folder_data, "output/isolates.csv"), show_col_types = F) # 68 isolates
accuracy <- read_csv(paste0(folder_data, "temp/24-accuracy.csv"), show_col_types = F)

# Random forest model accuracy
isolates_removal <- isolates$ExpID[which(is.na(isolates$BasePairMismatch))] # Isolates that do not match ESV
accuracy_to_plot <- accuracy %>%
    # Remove pairs that have cocultures with no colony
    unite(col = "Pair", Community, Isolate1, Isolate2, sep = "_", remove = F) %>%
    filter(!(Pair %in% pairs_no_colony)) %>% # 186-6=180 pairs. 180*3 = 540 cocultures
    # Remove the pairs containing the isolates that do not match ESV
    left_join(select(isolates, Community, Isolate1 = Isolate, ExpID1 = ExpID)) %>%
    left_join(select(isolates, Community, Isolate2 = Isolate, ExpID2 = ExpID)) %>%
    filter(!(ExpID1 %in% isolates_removal) & !(ExpID2 %in% isolates_removal))
nrow(accuracy_to_plot) # 180-20 = 160 pairs. (160-6)*3 = 154*3=462 cocultures


p1 <- accuracy_to_plot %>%
    ggplot() +
    geom_histogram(aes(x = Accuracy), color = 1, binwidth = 0.01, breaks = seq(0.6,1,.01), fill = NA) +
    geom_text(x = -Inf, y = Inf, label = paste0("N=", nrow(accuracy_to_plot)), vjust = 2, hjust = -1) +
    geom_vline(xintercept = 0.9, color = "red", linetype = 2) +
    theme_classic() +
    labs(x = "Accuracy", y = "Count")

#
accuracy_to_plot_count <- accuracy_to_plot %>%
    group_by(AccuracyPassThreshold) %>%
    count(name = "Count")

p2 <- accuracy_to_plot_count %>%
    ggplot() +
    geom_col(aes(x = AccuracyPassThreshold, y = Count), color = 1, fill = NA, width = .8) +
    geom_text(aes(x = AccuracyPassThreshold, y = Count, label = paste0("n=", Count)), vjust = -1) +
    geom_text(x = -Inf, y = Inf, label = paste0("N=", nrow(accuracy_to_plot)), vjust = 2, hjust = -1) +
    scale_x_discrete(labels = c("FALSE" = "Accuracy<0.9", "TRUE" = "Accuracy>0.9")) +
    scale_y_continuous(limits = c(0, 560)) +
    theme_classic() +
    labs(x = "")
p <- plot_grid(p1, p2, nrow = 1, labels = LETTERS[1:2], scale = 0.9) + paint_white_background()
ggsave(here::here("plots/FigS11-random_forest_accuracy.png"), p, width = 8, height = 4)

