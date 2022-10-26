#' This script compare the random forest prediction to the human eye counts

library(tidyverse)
library(cowplot)
source(here::here("analysis/00-metadata.R"))

pairs_freq_ID <- read_csv(paste0(folder_data, "temp/00c-pairs_freq_ID.csv"), show_col_types = F)
pairs_T0_boots <- read_csv(paste0(folder_data, "temp/06-pairs_T0_boots.csv"), show_col_types = F) # bootstraps using T0 mean and sd
pairs_T8_boots <- read_csv(paste0(folder_data, "temp/06-pairs_T8_boots.csv"), show_col_types = F) # bootstraps using random forest object probabilities
pairs_T8 <- read_csv(paste0(folder_data, "temp/06-pairs_T8.csv"), show_col_types = F) # random forest classification
# Human manual results
pairs_freq <- read_csv(paste0(folder_pipeline, "result_pairwise_competition_arranged.csv"), show_col_types = F) # human-eye results
isolates_duplicate <- tibble(ID = c(462, 355, 356, 461, 452, 446, 305, 435, 444, 348, 460, 454), Duplicated = T)
isolates_ID_match <- read_csv(paste0(folder_data, "raw/pairwise_competition/isolates1.csv"), col_types = cols()) %>%
    mutate(Assembly = "self_assembly") %>%
    select(ID, Community, Isolate) %>%
    left_join(isolates_duplicate, by = "ID") %>%
    replace_na(list(Duplicated = F))


# 1. Clean up column names ----
pairs_freq_renamed <- pairs_freq %>%
    mutate(Experiment = str_replace(Experiment, "Transitivity_", "")) %>%
    select(Batch = Experiment, Community, Isolate1, Isolate2, Isolate1InitialODFreq = Isolate1Freq, Isolate1Count = ColonyCount1, TotalCount = ColonyCount) %>%
    group_by(Batch, Community, Isolate1, Isolate2, Isolate1InitialODFreq) %>%
    mutate(Isolate1CFUFreq = Isolate1Count / TotalCount) %>%
    mutate(Type = "human") %>%
    ungroup()

pairs_T8_renamed <- pairs_T8 %>%
    select(Batch, Community, Isolate1, Isolate2, Isolate1InitialODFreq, Isolate1Count, TotalCount, Isolate1CFUFreq) %>%
    mutate(Type = "machine")

pairs_T8_boots_renamed <- pairs_T8_boots %>%
    #select(Batch, Community, Isolate1, Isolate2, Isolate1InitialODFreq, Isolate1Count, TotalCount, Isolate1CFUFreq) %>%
    mutate(Type = "machine-bootstrap") %>%
    group_by(Batch, Community, Isolate1, Isolate2, Isolate1InitialODFreq) %>%
    # Mean and sd over 1000 bootstraps
    summarize(Isolate1CFUFreqMean = mean(Isolate1CFUFreq),
              Isolate1CFUFreqSd = sd(Isolate1CFUFreq)) %>%
    #mutate(Type = "machine-bootstrap") %>%
    ungroup() %>%
    rename(Isolate1CFUFreqMean_machine = Isolate1CFUFreqMean, Isolate1CFUFreqSd_machine = Isolate1CFUFreqSd)

pairs_T8_boots_combined <- pairs_freq_renamed %>%
    select(-Type, -Isolate1Count, -TotalCount) %>%
    rename(Isolate1CFUFreq_human = Isolate1CFUFreq) %>%
    left_join(pairs_T8_boots_renamed) %>%
    left_join(pairs_freq_ID) %>%
    select(-image_name_isolate1, -image_name_isolate2) %>%
    # Label pairs containing duplicated
    left_join(rename(isolates_ID_match, ID1 = ID, Isolate1 = Isolate, Duplicated1 = Duplicated), by = c("Community", "Isolate1")) %>%
    left_join(rename(isolates_ID_match, ID2 = ID, Isolate2 = Isolate, Duplicated2 = Duplicated), by = c("Community", "Isolate2")) %>%
    mutate(PairType = case_when(
        Duplicated1 == F & Duplicated2 == F ~ "clean",
        Duplicated1 == T & Duplicated2 == F ~ "one duplicate",
        Duplicated1 == F & Duplicated2 == T ~ "one duplicate",
        Duplicated1 == T & Duplicated2 == T ~ "both duplicate"
    )) %>%
    mutate(PairType = factor(PairType, c("clean", "one duplicate", "both duplicate"))) %>%
    select(-ID1, -ID2) %>%
    select(image_name_pair, everything())


pairs_T8_combined <- pairs_T8_renamed %>%
    bind_rows(pairs_freq_renamed) %>%
    pivot_wider(names_from = Type, values_from = c(Isolate1Count, TotalCount, Isolate1CFUFreq)) %>%
    left_join(pairs_freq_ID) %>%
    select(-image_name_isolate1, -image_name_isolate2) %>%
    # Label pairs containing duplicated
    left_join(rename(isolates_ID_match, ID1 = ID, Isolate1 = Isolate, Duplicated1 = Duplicated), by = c("Community", "Isolate1")) %>%
    left_join(rename(isolates_ID_match, ID2 = ID, Isolate2 = Isolate, Duplicated2 = Duplicated), by = c("Community", "Isolate2")) %>%
    mutate(PairType = case_when(
        Duplicated1 == F & Duplicated2 == F ~ "clean",
        Duplicated1 == T & Duplicated2 == F ~ "one duplicate",
        Duplicated1 == F & Duplicated2 == T ~ "one duplicate",
        Duplicated1 == T & Duplicated2 == T ~ "both duplicate"
    )) %>%
    mutate(PairType = factor(PairType, c("clean", "one duplicate", "both duplicate"))) %>%
    select(-ID1, -ID2) %>%
    select(image_name_pair, everything())

write_csv(pairs_T8_combined, paste0(folder_data, "temp/92-pairs_T8_combined.csv"))

# 2. Plots ----
# 2.1 Total count

p1a <- pairs_T8_combined %>%
    ggplot() +
    geom_abline(slope = 1, intercept = 0, color = "red", linetype = 2) +
    geom_point(aes(x = TotalCount_human, y = TotalCount_machine), shape = 21, size = 2) +
    theme_classic()

p1b <- pairs_T8_combined %>%
    ggplot() +
    geom_abline(slope = 1, intercept = 0, color = "red", linetype = 2) +
    geom_point(aes(x = TotalCount_human, y = TotalCount_machine), shape = 21, size = 2) +
    scale_x_log10() +
    scale_y_log10() +
    theme_classic() +
    labs(x = "log(TotalCount_human)", y = "log(TotalCount_machine)")
p1 <- plot_grid(p1a, p1b, nrow = 1, axis = "tblr", align = "h", scale = .9) +
    theme(plot.background = element_rect(fill = "white", color = NA))
ggsave(paste0(folder_data, "temp/92-comparison-total_count.png"), p1, width = 8, height = 4)


# 2.2 Frequency
## bootstraps combined
p2 <- pairs_T8_combined %>%
    ggplot() +
    geom_abline(slope = 1, intercept = 0, color = "red", linetype = 2) +
    geom_hline(yintercept = c(0,1), color = gray(.8), linetype = 2) +
    geom_vline(xintercept = c(0,1), color = gray(.8), linetype = 2) +
    geom_point(aes(x = Isolate1CFUFreq_human, y = Isolate1CFUFreq_machine),
               shape = 21, size = 2, stroke = .4) +
    scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 21)) +
    theme_classic() +
    ggtitle("")

## bootstraps facets by duplicate pairs
p3 <- pairs_T8_boots_combined %>%
    ggplot() +
    geom_abline(slope = 1, intercept = 0, color = "red", linetype = 2) +
    geom_hline(yintercept = c(0,1), color = gray(.8), linetype = 2) +
    geom_vline(xintercept = c(0,1), color = gray(.8), linetype = 2) +
    #geom_smooth(aes(x = Isolate1CFUFreq_human, y = Isolate1CFUFreq_machine), method = "lm") +
    geom_point(aes(x = Isolate1CFUFreq_human, y = Isolate1CFUFreqMean_machine),
               shape = 21, size = 2, stroke = .4) +
    geom_segment(aes(x = Isolate1CFUFreq_human, xend = Isolate1CFUFreq_human,
                     y = Isolate1CFUFreqMean_machine + 1* Isolate1CFUFreqSd_machine,
                     yend = Isolate1CFUFreqMean_machine - 1* Isolate1CFUFreqSd_machine),
                 size = .1) +
    scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 21)) +
    facet_grid(.~PairType) +
    theme_classic() +
    ggtitle("")

lm(Isolate1CFUFreqMean_machine ~ Isolate1CFUFreq_human, data = pairs_T8_boots_combined) %>%
    summary()


ggsave(paste0(folder_data, "temp/92-comparison-coculture_frequency.png"), p2, width = 4, height = 4)
ggsave(paste0(folder_data, "temp/92-comparison-coculture_frequency_facet.png"), p3, width = 10, height = 4)


