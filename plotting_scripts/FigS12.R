library(tidyverse)
library(cowplot)
source(here::here("processing_scripts/00-metadata.R"))

isolates <- read_csv(paste0(folder_data, "output/isolates.csv"), show_col_types = F) # 65 isolates
pairs_freq_machine_human <- read_csv(paste0(folder_data, "temp/27-pairs_freq_machine_human.csv"), show_col_types = F) # 459
accuracy <- read_csv(paste0(folder_data, "temp/24-accuracy.csv"), show_col_types = F)

# Check numbers
pairs_machine <- pairs_freq_machine_human %>%
    select(Community, Isolate1, Isolate2, Isolate1InitialODFreq, Isolate1CFUFreq_machine) %>%
    pivot_wider(names_from = Isolate1InitialODFreq, names_prefix = "F", values_from = Isolate1CFUFreq_machine) %>%
    filter(!is.na(F5), !is.na(F50), !is.na(F95)) %>%
    select(Community, Isolate1, Isolate2) %>%
    mutate(ContainMachine = T) %>%
    left_join(filter(accuracy, Isolate1InitialODFreq == 5)) %>% filter(Accuracy > 0.9)
nrow(pairs_machine) # 144 pairs have machine result. 159 total - 6 pairs that have no colony - 9 with low accurary = 144

pairs_human <- pairs_freq_machine_human %>%
    select(Community, Isolate1, Isolate2, Isolate1InitialODFreq, Isolate1CFUFreq_human) %>%
    pivot_wider(names_from = Isolate1InitialODFreq, names_prefix = "F", values_from = Isolate1CFUFreq_human) %>%
    filter(!is.na(F5), !is.na(F50), !is.na(F95)) %>%
    select(Community, Isolate1, Isolate2) %>%
    mutate(ContainHuman = T)
nrow(pairs_human) # 128 pairs have human result. 159 total - 6 pairs that have no colony - 25 hard to distinguish by eyes = 130

pairs_machine_human <- pairs_machine %>%
    left_join(pairs_human) %>%
    left_join(filter(accuracy, Isolate1InitialODFreq == 5)) %>%
    filter(ContainMachine, ContainHuman, Accuracy > 0.9)
nrow(pairs_machine_human) # 127 pairs that have both human and machine data, and the machine accurarcy is > 0.9

pairs_freq_machine_human_cleaned <- pairs_freq_machine_human %>% # 153*3 = 459 cocultures
    filter(paste0(Community, Isolate1, Isolate2) %in% paste0(pairs_machine_human$Community, pairs_machine_human$Isolate1, pairs_machine_human$Isolate2))
nrow(pairs_freq_machine_human_cleaned) # 127 pairs that have human result. 127*3 = 381 cocultures



#
p1 <- pairs_freq_machine_human_cleaned %>%
    ggplot() +
    geom_abline(slope = 1, intercept = 0, color = "red", linetype = 2) +
    geom_point(aes(x = TotalCount_human, y = TotalCount_machine), shape = 21, size = 2) +
    geom_text(x = -Inf, y = Inf, label = paste0("N=", nrow(pairs_freq_machine_human_cleaned)), vjust = 2, hjust = -1) +
    scale_x_log10() +
    scale_y_log10() +
    theme_classic() +
    labs(x = "Segmentation CFU count", y = "Manual CFU count")


p2 <- pairs_freq_machine_human_cleaned %>%
    ggplot() +
    geom_abline(slope = 1, intercept = 0, color = "red", linetype = 2) +
    geom_hline(yintercept = c(0,1), color = gray(.8), linetype = 2) +
    geom_vline(xintercept = c(0,1), color = gray(.8), linetype = 2) +
    geom_point(aes(x = Isolate1CFUFreq_human, y = Isolate1CFUFreq_machine), shape = 21, size = 2, stroke = .4) +
    geom_text(x = 0.5, y = 0.9, label = paste0("N=", nrow(pairs_freq_machine_human_cleaned))) +
    scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 21)) +
    theme_classic() +
    labs(x = "Random Forest CFU frequency", y = "Manual CFU frequency")

p <- plot_grid(p1, p2, nrow = 1, axis = "tblr", align = "h", scale = .9, labels = c("A", "B")) +
    paint_white_background()

ggsave(here::here("plots/FigS12.png"), p, width = 8, height = 4)

# statistics
## colony count
model <- pairs_freq_machine_human_cleaned %>%
    filter(!is.na(Isolate1Count_human)) %>%
    lm(TotalCount_human ~ TotalCount_machine, data = .)
summary(model)$r.squared # R2 = 0.85
sqrt(mean(model$residuals^2)) # RMSE=17.67

## frequency
model <- pairs_freq_machine_human_cleaned %>%
    filter(!is.na(Isolate1Count_human)) %>%
    lm(Isolate1CFUFreq_human ~ Isolate1CFUFreq_machine, data = .)
summary(model)$r.squared # R2 = 0.87
sqrt(mean(model$residuals^2)) # RMSE=0.17


