library(tidyverse)
library(cowplot)
source(here::here("processing_scripts/00-metadata.R"))

isolates <- read_csv(paste0(folder_data, "output/isolates.csv"), show_col_types = F) # 68 isolates
pairs_freq_machine_human <- read_csv(paste0(folder_data, "temp/27-pairs_freq_machine_human.csv"), show_col_types = F) # 558 pairs
accuracy <- read_csv(paste0(folder_data, "temp/24-accuracy.csv"), show_col_types = F)

# Step 1: remove the pairs containing the four isolates
isolates_removal <- isolates$ExpID[which(is.na(isolates$BasePairMismatch))] # Isolates that do not match ESV
pairs_freq_machine_human <- pairs_freq_machine_human %>%
    left_join(accuracy) %>%
    # Remove the pairs containing the isolates that do not match ESV
    left_join(select(isolates, Community, Isolate1 = Isolate, ExpID1 = ExpID)) %>%
    left_join(select(isolates, Community, Isolate2 = Isolate, ExpID2 = ExpID)) %>%
    filter(!(ExpID1 %in% isolates_removal) & !(ExpID2 %in% isolates_removal))

nrow(pairs_freq_machine_human) # 186-26 = 160 pairs. 160*3=480 cocultures

# Check numbers
pairs_machine <- pairs_freq_machine_human %>%
    select(Community, Isolate1, Isolate2, Isolate1InitialODFreq, Isolate1CFUFreq_machine) %>%
    pivot_wider(names_from = Isolate1InitialODFreq, names_prefix = "F", values_from = Isolate1CFUFreq_machine) %>%
    filter(!is.na(F5), !is.na(F50), !is.na(F95)) %>%
    select(Community, Isolate1, Isolate2) %>%
    mutate(ContainMachine = T)
nrow(pairs_machine) # 154 pairs have machine result. 160 total - 6 pairs that have no colony = 154

pairs_human <- pairs_freq_machine_human %>%
    select(Community, Isolate1, Isolate2, Isolate1InitialODFreq, Isolate1CFUFreq_human) %>%
    pivot_wider(names_from = Isolate1InitialODFreq, names_prefix = "F", values_from = Isolate1CFUFreq_human) %>%
    filter(!is.na(F5), !is.na(F50), !is.na(F95)) %>%
    select(Community, Isolate1, Isolate2) %>%
    mutate(ContainHuman = T)
nrow(pairs_human) # 130 pairs have human result. 160 total - 6 pairs that have no colony - 24 hard to distinguish by eyes = 130

pairs_human_machine <- pairs_freq_machine_human %>%
    select(Community, Isolate1, Isolate2, Isolate1InitialODFreq, machine = Isolate1CFUFreq_machine, human = Isolate1CFUFreq_human) %>%
    pivot_wider(names_from = Isolate1InitialODFreq, names_prefix = "F", values_from = c(human, machine)) %>%
    filter(!is.na(human_F5), !is.na(human_F50), !is.na(human_F95), !is.na(machine_F5), !is.na(machine_F50), !is.na(machine_F95)) %>%
    select(Community, Isolate1, Isolate2) %>%
    mutate(ContainHumanMachine = T)
nrow(pairs_human_machine) # 130 pairs have both human and machine data

pairs_machine %>%
    left_join(pairs_human) %>%
    left_join(filter(accuracy, Isolate1InitialODFreq == 5)) %>%
    filter(ContainMachine, ContainHuman, Accuracy > 0.9) %>%
    nrow() # 129 pairs that have both human and machine data, and the machine accurarcy is > 0.9


#
pairs_freq_machine_human_cleaned1 <- pairs_freq_machine_human %>% # 160*3 = 480 cocultures
    left_join(pairs_machine) %>%
    filter(ContainMachine) %>% # 154 pairs that have machine result. 154*3 = 462 cocultures
    # Remove those with low accuracy
    left_join(accuracy) %>%
    filter(Accuracy > 0.9) # 145 pairs that have human result. 145*3=435 cocultures
nrow(pairs_freq_machine_human_cleaned1) # 145*3 = 135 cocultures that have high accurarcy machine result

pairs_freq_machine_human_cleaned2 <- pairs_freq_machine_human_cleaned1 %>% # 145*3=435 cocultures
    # Remove pairs that have no human results
    left_join(pairs_human_machine) %>%
    select(Community, Isolate1, Isolate2, Isolate1InitialODFreq,
           Isolate1Count_machine, Isolate1Count_human,
           TotalCount_machine, TotalCount_human,
           Isolate1CFUFreq_machine, Isolate1CFUFreq_human,
           ContainHumanMachine) %>%
    filter(ContainHumanMachine)
nrow(pairs_freq_machine_human_cleaned2) # 129*3 = 387 cocultures

pairs_freq_machine_human_cleaned3 <- pairs_freq_machine_human %>% # 160*3 = 480 cocultures
    left_join(pairs_human) %>%
    filter(ContainHuman)
nrow(pairs_freq_machine_human_cleaned3) # 130 pairs that have human result. 130*3 = 390 cocultures


#
p1 <- pairs_freq_machine_human_cleaned2 %>%
    ggplot() +
    geom_abline(slope = 1, intercept = 0, color = "red", linetype = 2) +
    geom_point(aes(x = TotalCount_human, y = TotalCount_machine), shape = 21, size = 2) +
    geom_text(x = -Inf, y = Inf, label = paste0("N=", nrow(pairs_freq_machine_human_cleaned2)), vjust = 2, hjust = -1) +
    scale_x_log10() +
    scale_y_log10() +
    theme_classic() +
    labs(x = "Segmentation CFU count", y = "Manual CFU count")


p2 <- pairs_freq_machine_human_cleaned2 %>%
    ggplot() +
    geom_abline(slope = 1, intercept = 0, color = "red", linetype = 2) +
    geom_hline(yintercept = c(0,1), color = gray(.8), linetype = 2) +
    geom_vline(xintercept = c(0,1), color = gray(.8), linetype = 2) +
    geom_point(aes(x = Isolate1CFUFreq_human, y = Isolate1CFUFreq_machine), shape = 21, size = 2, stroke = .4) +
    geom_text(x = 0.5, y = 0.9, label = paste0("N=", nrow(pairs_freq_machine_human_cleaned2))) +
    scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 21)) +
    theme_classic() +
    labs(x = "Random Forest CFU frequency", y = "Manual CFU frequency")

p <- plot_grid(p1, p2, nrow = 1, axis = "tblr", align = "h", scale = .9, labels = c("A", "B")) +
    paint_white_background()

ggsave(here::here("plots/FigS12-human_machine_comparison.png"), p, width = 8, height = 4)

# R-squared
pairs_freq_machine_human_cleaned2 %>%
    filter(!is.na(Isolate1Count_human)) %>%
    lm(TotalCount_human ~ TotalCount_machine, data = .) %>%
    summary() # adjusted R-squared:  0.8478
pairs_freq_machine_human_cleaned2 %>%
    filter(!is.na(Isolate1Count_human)) %>%
    lm(Isolate1CFUFreq_human ~ Isolate1CFUFreq_machine, data = .) %>%
    summary() # adjusted R-squared:  0.8665
