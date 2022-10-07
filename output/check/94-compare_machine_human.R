#' This script compare the random forest prediction to the human eye counts
library(tidyverse)
library(cowplot)
folder_script <- "~/Desktop/lab/emergent-coexistence/output/check/"
folder_main <- "~/Dropbox/lab/emergent-coexistence/data/raw/plate_scan/emergent_coexistence_plate_scan_check/"

# object_prediction_batch <- write_csv(paste0(folder_main, "meta/object_prediction_batch.csv")) # Object probability for all pairs
# accuracy_validation_batch <- write_csv(paste0(folder_main, "meta/accuracy_validation_batch.csv")) # Accuracy of random forest for all pairs
# object <- read_csv(paste0(folder_main, "meta/object.csv"), show_col_types = F) # Object probability after cleaning up pairs
# accuracy <- read_csv(paste0(folder_main, "meta/accuracy.csv"), show_col_types = F) # Accuracy after cleaning up pairs
# boots <- read_csv(paste0(folder_main, "meta/bootstraps.csv"), show_col_types = F) # Bootstrapped results using the object probability
pairs_freq_ID <- read_csv(paste0(folder_main, "meta/00-pairs_freq_ID.csv"), show_col_types = F)
pairs_T0_boots <- read_csv(paste0(folder_main, "meta/93-pairs_T0_boots.csv"), show_col_types = F) # bootstraps using T0 mean and sd
pairs_T8_boots <- read_csv(paste0(folder_main, "meta/93-pairs_T8_boots.csv"), show_col_types = F) # bootstraps using random forest object probabilities
pairs_T8 <- read_csv(paste0(folder_main, "meta/93-pairs_T8.csv"), show_col_types = F) # random forest classification
pairs_freq <- read_csv("~/Dropbox/lab/emergent-coexistence/data/raw/plate_scan/emergent_coexistence_plate_scan_check/result_pairwise_competition_arranged.csv", show_col_types = F) # human-eye results
isolates_duplicate <- tibble(ID = c(462, 355, 356, 461, 452, 446, 305, 435, 444, 348, 460, 454), Duplicated = T)
isolates_ID_match <- read_csv("~/Dropbox/lab/emergent-coexistence/data/raw/pairwise_competition/isolates1.csv", col_types = cols()) %>%
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
ggsave(paste0(folder_main, "meta/94-comparison-total_count.png"), p1, width = 8, height = 4)


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


# p2 <- plot_grid(p2a, p2b, nrow = 2, axis = "tblr", align = "v", scale = .9) +
#     theme(plot.background = element_rect(fill = "white", color = NA))
ggsave(paste0(folder_main, "meta/94-comparison-coculture_frequency.png"), p2, width = 4, height = 4)
ggsave(paste0(folder_main, "meta/94-comparison-coculture_frequency_facet.png"), p3, width = 10, height = 4)













if (FALSE) {

isolates_duplicate <- tibble(ID = c(462, 355, 356, 461, 452, 446, 305, 435, 444, 348, 460, 454), Duplicated = T)
isolates_ID_match <- read_csv("~/Dropbox/lab/emergent-coexistence/data/raw/pairwise_competition/isolates1.csv", col_types = cols()) %>%
    mutate(Assembly = "self_assembly") %>%
    select(ID, Community, Isolate) %>%
    left_join(isolates_duplicate, by = "ID") %>%
    replace_na(list(Duplicated = F))
pairs_similar <- tibble(
    Community = c(rep("C1R4", 3),
                  "C2R6", "C2R8", "C4R1",
                  rep("C7R1", 3),
                  "C8R4",
                  rep("C11R1", 5),
                  rep("C11R2", 7),
                  rep("C11R5", 4)),
    Isolate1 = c(2, 2, 4,
                 1, 3, 1,
                 1, 1, 3,
                 1,
                 1,1,2,2,8,
                 2,3,3,3,5,7,8,
                 1,1,2,3),
    Isolate2 = c(4, 5, 5,
                 3, 4, 3,
                 3, 4, 4,
                 3,
                 3,4,8,9,9,
                 8,5,7,12,6,12,9,
                 3,5,4,5),
    SimilarByEye = T
)


# 2. Clean up column names and result ----
pairs_freq_renamed <- pairs_freq %>%
    mutate(Experiment = str_replace(Experiment, "Transitivity_", "")) %>%
    select(Batch = Experiment, Community, Isolate1, Isolate2, Isolate1Freq, Isolate1Count = ColonyCount1, TotalCount = ColonyCount) %>%
    group_by(Batch, Community, Isolate1, Isolate2, Isolate1Freq) %>%
    mutate(Isolate1FinalFraction_human = Isolate1Count / TotalCount) %>%
    select(-Isolate1Count) %>%
    rename(TotalCount_human = TotalCount)
boots_summary <- boots %>%
    mutate(Isolate1FinalFraction = Isolate1Count / TotalCount) %>%
    group_by(image_name_pair) %>%
    # Mean and sd over 1000 bootstraps
    summarize(Isolate1FinalFractionMean_machine = mean(Isolate1FinalFraction),
              Isolate1FinalFractionSd_machine = sd(Isolate1FinalFraction),
              TotalCount_machine = unique(TotalCount)) %>%
    separate(image_name_pair, into = c("Batch", "temp", "Community", "Isolate2Freq", "Isolate1Freq", "Isolate1", "Isolate2"), convert = T, remove = F) %>%
    #filter(Isolate1 < Isolate2) %>%
    select(image_name_pair, Batch, Community, Isolate1, Isolate2, Isolate1Freq, Isolate2Freq, Isolate1FinalFractionMean_machine, Isolate1FinalFractionSd_machine, TotalCount_machine) %>%
    rowwise() %>%
    mutate(Isolate1Freq = ifelse(Isolate1 > Isolate2, 5, Isolate1Freq),
           Isolate2Freq = ifelse(Isolate1 > Isolate2, 95, Isolate2Freq),
           Isolate1FinalFractionMean_machine = ifelse(Isolate1 > Isolate2, 1 - Isolate1FinalFractionMean_machine, Isolate1FinalFractionMean_machine)
           ) %>%
    mutate(temp = min(Isolate1,Isolate2), Isolate2 = max(Isolate1, Isolate2), Isolate1 = temp) %>%
    select(-temp) %>%
    mutate(Community = factor(Community, paste0("C", rep(1:12, each = 8), "R", rep(1:8, 12)))) %>%
    arrange(Community, Isolate1, Isolate2, Isolate1Freq)

pairs_freq_combined <- pairs_freq_renamed %>%
    left_join(boots_summary) %>%
    left_join(select(accuracy, image_name_pair, Accuracy, AccuracyPassThreshold)) %>%
    # # Remove C11R2 isolate 13, which is a Staph. It's also not included in isolates_ID_match
    filter(!((Community == "C11R2" & Isolate1 == 13) | ((Community == "C11R2" & Isolate2 == 13)))) %>%
    # Remove batch B2 C11R1 pairs that contains isolate 1
    filter(!(Batch == "B2" & Community == "C11R1" & (Isolate1 == 1 | Isolate2 == 1))) %>%
    # Keep batch C C11R1 pairs that contains isolate 1
    filter(!(Batch == "C" & Community == "C11R1" & (Isolate1 != 1 & Isolate2 != 1))) %>%
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
    select(image_name_pair, Batch, Community, Isolate1, Isolate2, Isolate1Freq,
           Isolate1FinalFraction_human, Isolate1FinalFractionMean_machine, Isolate1FinalFractionSd_machine,
           TotalCount_human, TotalCount_machine,
           Accuracy, AccuracyPassThreshold, PairType, Duplicated1, Duplicated2)



#
pairs_freq_combined %>%
    filter(Isolate1FinalFraction_human == 0) %>%
    mutate(Isolate1FinalFractionMean_machine = factor(Isolate1FinalFractionMean_machine, Isolate1FinalFractionMean_machine)) %>%
    arrange(Isolate1FinalFractionMean_machine) %>%

    ggplot() +
    geom_point(aes(x = image_name_pair, y = Isolate1FinalFractionMean_machine),
               shape = 21, size = 2, stroke = .2) +
    # geom_segment(aes(x = Isolate1FinalFraction_human, xend = Isolate1FinalFraction_human,
    #                  y = Isolate1FinalFractionMean_machine + 1* Isolate1FinalFractionSd_machine,
    #                  yend = Isolate1FinalFractionMean_machine - 1* Isolate1FinalFractionSd_machine),
    #              size = .1) +
    #scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 21)) +
    facet_grid(.~PairType, scales = "free_x") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, hjust = 0.5)) +
    labs()



# See those that are off
pairs_freq_combined %>%
    filter(Isolate1FinalFraction_human > 0.9,
           Isolate1FinalFractionMean_machine < 0.1) %>%
    view



pairs_freq_combined %>%
    filter(Isolate1FinalFraction_human < 0.1,
           Isolate1FinalFractionMean_machine > 0.8) %>%
    view




}
