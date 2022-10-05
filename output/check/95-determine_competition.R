#' This scripts reads the T0 vs. T8 frequencies and determine the competition outcomes
#' for each unique species pair by computing the fitness functions
#' 1. Combine T0 and T8 results
#' 2. Map the frequency changes to the fitness function
#' 3. Determine the competition outcomes

library(tidyverse)
library(cowplot)
folder_script <- "~/Desktop/lab/emergent-coexistence/output/check/"
folder_main <- "~/Dropbox/lab/emergent-coexistence/data/raw/plate_scan/emergent_coexistence_plate_scan_check/"

pairs_freq_ID <- read_csv(paste0(folder_main, "meta/00-pairs_freq_ID.csv"), show_col_types = F)
pairs_T0_boots <- read_csv(paste0(folder_main, "meta/93-pairs_T0_boots.csv"), show_col_types = F) # bootstraps using T0 mean and sd
pairs_T8_boots <- read_csv(paste0(folder_main, "meta/93-pairs_T8_boots.csv"), show_col_types = F) # bootstraps using random forest object probabilities
pairs_T8 <- read_csv(paste0(folder_main, "meta/93-pairs_T8.csv"), show_col_types = F) # random forest classification
pairs_ID <- distinct(pairs_freq_ID, Batch, Community, Isolate1, Isolate2)

pairs_no_colony <- c(
    "C11R1_2_8",
    "C11R1_2_9",
    "C11R1_8_9",
    "C11R2_2_10"
)
plates_no_colony <- c(
    "B2_T8_C11R1_5-95_2_8",
    "B2_T8_C11R1_5-95_2_9",
    "B2_T8_C11R1_5-95_8_2",
    "B2_T8_C11R1_5-95_9_8",
    "B2_T8_C11R1_50-50_2_8",
    "B2_T8_C11R1_50-50_2_9",
    "C2_T8_C11R2_50-50_2_10",
    "C2_T8_C11R2_50-50_9_13",
    "C_T8_C11R1_50-50_1_2", # no plate
    "C_T8_C11R1_50-50_1_3" # no plate
)


# 1. Calculate the significance of frequency change
read_pairs_boots_table <- function (T8_freq_type = "bootstrapped") {
    temp <- rep(list(NA), nrow(pairs_ID))
    for (i in 1:nrow(pairs_ID)) {
        pair_name <- paste0(pairs_ID$Community[i], "_", pairs_ID$Isolate1[i], "_", pairs_ID$Isolate2[i])
        if (pair_name %in% pairs_no_colony) {cat("\nT8 has no colony, skip pair\t", pair_name); next}
        temp[[i]] <- read_csv(paste0(folder_main, "check/meta-93-T0_T8_frequencies/", T8_freq_type, "/", pair_name, ".csv"), show_col_types = F)
        cat("\n", i, "/", nrow(pairs_ID))
    }
    pairs_boots_table <- bind_rows(temp[which(!is.na(temp))])
}
compute_pairs_fitness <- function (pairs_boots_table) {
    pairs_boots_table %>%
        group_by(Community) %>%
        mutate(Community = factor(Community, paste0("C", rep(1:12, each = 8), "R", rep(1:8, 12)))) %>%
        arrange(Community, Isolate1, Isolate2, Isolate1InitialODFreq) %>%
        group_by(Community, Isolate1, Isolate2, Isolate1InitialODFreq, .drop = F) %>%
        mutate(Fraction = Count / sum(Count)) %>%
        # C11R1 isolate 6
        filter(!is.na(FreqChangeSign)) %>%
        pivot_wider(id_cols = c(Community, Isolate1, Isolate2, Isolate1InitialODFreq), names_from = FreqChangeSign, values_from = Fraction) %>%
        mutate(FitnessChange = case_when(
            increase > 0.95 ~ 1,
            decrease > 0.95 ~ -1,
            increase < 0.95 & decrease < 0.95 ~ 0
        ))
}
make_interaction_type <- function () {
    interaction_type_three <- tibble(
        FromRare = rep(c(1, -1, 0), each = 9),
        FromMedium = rep(rep(c(1, -1, 0), each = 3), 3),
        FromAbundant = rep(c(1, -1, 0), 9),
        InteractionType = NA,
        InteractionTypeFiner = NA
    )
    interaction_type_two <- tibble(
        FromRare = rep(c(1, -1, 0), each = 3),
        FromMedium = rep(NA, 9),
        FromAbundant = rep(c(1, -1, 0), 3),
        InteractionType = NA,
        InteractionTypeFiner = NA
    )

    interaction_type <- bind_rows(interaction_type_three, interaction_type_two)

    ## Assign interaction types to combinations of frequency changes signs
    interaction_type$InteractionType[c(1, 14, 28, 32, 10, 13, 31)] <- "exclusion"
    interaction_type$InteractionType[c(2, 3, 5, 8, 9, 23, 26, 29, 30, 33, 4, 11, 12, 15, 17, 20, 34, 35)] <- "coexistence"
    interaction_type$InteractionType[c(27, 36)] <- "coexistence"

    ## Assign finer interaction types to combinations of frequency changes signs
    interaction_type$InteractionTypeFiner[c(1, 14, 28, 32)] <- "competitive exclusion"
    interaction_type$InteractionTypeFiner[c(10, 13, 31)] <- "mutual exclusion"
    interaction_type$InteractionTypeFiner[c(2, 3, 5, 8, 9, 23, 26, 29, 30, 33)] <- "stable coexistence"
    interaction_type$InteractionTypeFiner[c(4, 11, 12, 15, 17, 20, 34, 35)] <- "frequency-dependent coexistence"
    interaction_type$InteractionTypeFiner[c(27,36)] <- "neutrality"
    interaction_type <- interaction_type %>%  mutate(FitnessFunction = paste(FromRare, FromMedium, FromAbundant, sep = "_"))
}
append_pairs_outocme <- function (pairs_fitness) {
    pairs_fitness %>%
    pivot_wider(id_cols = c(Community, Isolate1, Isolate2), names_from = Isolate1InitialODFreq, values_from = FitnessChange) %>%
    rename(FromRare = `5`, FromMedium = `50`, FromAbundant = `95`) %>%
    left_join(interaction_type) %>%
    select(Community, Isolate1, Isolate2, InteractionType, InteractionTypeFiner, FitnessFunction)

}
interaction_type <- make_interaction_type()

pairs_outcome_classifiedT8 <- read_pairs_boots_table("classified") %>%
    compute_pairs_fitness() %>%
    append_pairs_outocme()

pairs_outcome_bootstrappedT8 <- read_pairs_boots_table("bootstrapped") %>%
    compute_pairs_fitness() %>%
    append_pairs_outocme()

write_csv(pairs_outcome_classifiedT8, paste0(folder_main, "meta/95-pairs_outcome_classifiedT8.csv"))
write_csv(pairs_outcome_bootstrappedT8, paste0(folder_main, "meta/95-pairs_outcome_bootstrappedT8.csv"))

# 2. Combine T0 and T8 frequencies ----
pairs_freq_T0_boots <- pairs_T0_boots %>%
    group_by(Community, Isolate1, Isolate2, Isolate1InitialODFreq, Time) %>%
    summarize(Isolate1CFUFreqMean = mean(Isolate1CFUFreq),
              Isolate1CFUFreqSd = sd(Isolate1CFUFreq))
pairs_freq_T8_boots <- pairs_T8_boots %>%
    group_by(Community, Isolate1, Isolate2, Isolate1InitialODFreq, Time) %>%
    summarize(Isolate1CFUFreqMean = mean(Isolate1CFUFreq),
              Isolate1CFUFreqSd = sd(Isolate1CFUFreq))
pairs_freq_T8 <- pairs_T8 %>%
    select(Community, Isolate1, Isolate2, Isolate1InitialODFreq, Time,
           Isolate1CFUFreqMean = Isolate1CFUFreq) %>%
    mutate(Isolate1CFUFreqSd = NA)


# T8 frequency bootstrapping so there is SD at T8
pairs_freq_boots <- bind_rows(pairs_freq_T0_boots, pairs_freq_T8_boots)

# No bootstrapping, the random forest predicted classification, so there is no SD at T8
pairs_freq <- bind_rows(pairs_freq_T0_boots, pairs_freq_T8)

write_csv(pairs_freq_boots, paste0(folder_main, "meta/95-pairs_freq_boots.csv"))
write_csv(pairs_freq, paste0(folder_main, "meta/95-pairs_freq.csv"))









