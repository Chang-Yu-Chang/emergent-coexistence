#' This scripts uses the updated criterion of defining pairwise competition outcomes
#'
#' 1. Combine T0 and T8 results into a pairs_freq table
#' 2. Determine the competition outcomes: match the significance of frequency changes to the fitness function table
#' 3. Calculate isolate tournament
#'
#' This scrip is meant to replace 93-determine_competition.R
#'
#' There are five categories
#' 1) Loser is extinct: EXTINCT
#' 2) Loser declines in all three competitions, not extinct in all three.
#' 3) ALL THREE have both species, both invade when rare.
#' 4) ALL THREE have both species, they do not BOTH invade when rare
#' 5) 2/3 replicates OR 1/3 have the loser going extinct. Not conclusive

library(tidyverse)
library(cowplot)
source(here::here("analysis/00-metadata.R"))

communities <- read_csv(paste0(folder_data, "temp/00c-communities.csv"), show_col_types = F)
pairs_ID <- read_csv(paste0(folder_data, "temp/00c-pairs_ID.csv"), show_col_types = F)
pairs_T0_boots <- read_csv(paste0(folder_data, "temp/06-pairs_T0_boots.csv"), show_col_types = F) # bootstraps using T0 mean and sd
pairs_T8_boots <- read_csv(paste0(folder_data, "temp/06-pairs_T8_boots.csv"), show_col_types = F) # bootstraps using random forest object probabilities=

# 1. Bind the pair frequency data ----
pairs_freq_T0_boots <- pairs_T0_boots %>%
    group_by(Batch, Community, Isolate1, Isolate2, Isolate1InitialODFreq, Time) %>%
    summarize(Isolate1CFUFreqMean = mean(Isolate1CFUFreq),
              Isolate1CFUFreqSd = sd(Isolate1CFUFreq)) %>%
    mutate(Community = factor(Community, paste0("C", rep(1:12, each = 8), "R", rep(1:8, 12)))) %>%
    arrange(Community, Isolate1, Isolate2, Isolate1InitialODFreq) %>%
    ungroup() %>%
    left_join(pairs_ID) %>%
    # PairFreqID
    mutate(PairFreqID = 1:n()) %>%
    select(PairID, PairFreqID, everything())
pairs_freq_T8_boots <- pairs_T8_boots %>%
    group_by(Batch, Community, Isolate1, Isolate2, Isolate1InitialODFreq, Time) %>%
    summarize(Isolate1CFUFreqMean = mean(Isolate1CFUFreq),
              Isolate1CFUFreqSd = sd(Isolate1CFUFreq)) %>%
    mutate(Community = factor(Community, paste0("C", rep(1:12, each = 8), "R", rep(1:8, 12)))) %>%
    arrange(Community, Isolate1, Isolate2, Isolate1InitialODFreq) %>%
    left_join(pairs_ID) %>%
    left_join(select(pairs_freq_T0_boots, PairFreqID, Community, Isolate1, Isolate2, Isolate1InitialODFreq))

pairs_freq <- bind_rows(pairs_freq_T0_boots, pairs_freq_T8_boots) %>%
    select(PairID, PairFreqID, everything())

write_csv(pairs_freq, paste0(folder_data, "temp/93a-pairs_freq.csv"))


# 2. Determine pairwise competition ----
if (FALSE) {

read_pairs_boots_table <- function () {
    #' This function batchly reads the table of boostrapped frequency change per pair
    temp <- rep(list(NA), nrow(pairs_ID))
    for (i in 1:nrow(pairs_ID)) {
        pair_name <- paste0(pairs_ID$Community[i], "_", pairs_ID$Isolate1[i], "_", pairs_ID$Isolate2[i])
        if (pair_name %in% pairs_no_colony) {cat("\nT8 has no colony, skip pair\t", pair_name); next}
        temp[[i]] <- read_csv(paste0(folder_pipeline, "images/", pairs_ID$Batch[i],"-09-bootstrap/", pair_name, ".csv"), show_col_types = F)
        cat("\n", i, "/", nrow(pairs_ID))
    }
    pairs_boots_table <- bind_rows(temp[which(!is.na(temp))])
}
}
read_pairs_boots_table <- function () {
    #' This function binds the table of bootstrapped frequency change per pair
    pairs_boots_table_list <- rep(list(NA), nrow(pairs_ID))
    for (i in 1:nrow(pairs_ID)) {
        community <- pairs_ID$Community[i]
        isolate1 <- pairs_ID$Isolate1[i]
        isolate2 <- pairs_ID$Isolate2[i]
        pair_name <- paste0(community, "_", isolate1, "_", isolate2)
        if (pair_name %in% pairs_no_colony) {cat("\nT8 has no colony, skip pair\t", pair_name); next}

        # Increase or decrease sign
        pair_boots <- bind_rows(
            filter(pairs_T0_boots, Community == community, Isolate1 == isolate1, Isolate2 == isolate2),
            filter(pairs_T8_boots, Community == community, Isolate1 == isolate1, Isolate2 == isolate2),
        ) %>%
            select(Community, Isolate1, Isolate2, Isolate1InitialODFreq, Time, BootstrapID, Isolate1CFUFreq) %>%
            pivot_wider(names_from = Time, values_from = Isolate1CFUFreq) %>%
            mutate(FreqChange = T8-T0) %>%
            mutate(FreqChangeSign = case_when(
                FreqChange > 0 ~ "increase",
                FreqChange == 0 ~ "same",
                FreqChange < 0 ~ "decrease",
            )) %>%
            pivot_longer(cols = c(T0, T8), names_to = "Time", values_to = "Isolate1CFUFreq")

        # Table of increase and decrease
        pairs_boots_table_list[[i]] <- pair_boots %>%
            filter(Time == "T0") %>%
            mutate(FreqChangeSign = factor(FreqChangeSign, c("increase", "same", "decrease"))) %>%
            group_by(Community, Isolate1, Isolate2, Isolate1InitialODFreq, FreqChangeSign, .drop = F) %>%
            count(name = "Count") %>%
            group_by(Community, Isolate1, Isolate2, Isolate1InitialODFreq) %>%
            mutate(Fraction = Count / sum(Count))
    }
    pairs_boots_table <- bind_rows(pairs_boots_table_list[!is.na(pairs_boots_table_list)])
    return(pairs_boots_table)
}
compute_pairs_fitness <- function (pairs_boots_table) {
    #' This script takes the bootstrap samples and compute the significance of frequency changes
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
    #' This function generates the fitness function table.
    #' There are a total of 27 possibilities
    interaction_type <- tibble(
        FromRare = rep(c(1, -1, 0), each = 9),
        FromMedium = rep(rep(c(1, -1, 0), each = 3), 3),
        FromAbundant = rep(c(1, -1, 0), 9),
        InteractionType = NA,
        InteractionTypeFiner = NA
    )
    ## Assign interaction types to combinations of frequency changes signs
    interaction_type$InteractionType[c(1,10,13,14)] <- "exclusion"
    interaction_type$InteractionType[c(2:6,8,11,20,23, 9,18,21,24,25,26,27)] <- "coexistence"
    interaction_type$InteractionType[c(7,12,15:17,19,22)] <- "unknown"

    ## Assign finer interaction types to combinations of frequency changes signs
    interaction_type$InteractionTypeFiner[c(1,14)] <- "competitive exclusion"
    interaction_type$InteractionTypeFiner[c(10,13)] <- "mutual exclusion"
    interaction_type$InteractionTypeFiner[c(2,5,8)] <- "stable coexistence"
    interaction_type$InteractionTypeFiner[c(4,6,11,20)] <- "frequency-dependent coexistence"
    interaction_type$InteractionTypeFiner[c(3)] <- "coexistence at 95%"
    interaction_type$InteractionTypeFiner[c(23)] <- "coexistence at 5%"
    interaction_type$InteractionTypeFiner[c(9,18,21,24:26)] <- "neutrality"
    interaction_type$InteractionTypeFiner[c(9,18,21,24:26)] <- "2-freq neutrality"
    interaction_type$InteractionTypeFiner[c(27)] <- "3-freq neutrality"
    interaction_type$InteractionTypeFiner[c(7,12,15:17,19,22)] <- "unknown"

    interaction_type <- interaction_type %>%  mutate(FitnessFunction = paste(FromRare, FromMedium, FromAbundant, sep = "_"))
}
append_pairs_outcome <- function (pairs_fitness) {
    # This function appends the bootstrap result to the fitness function table
    pairs_fitness %>%
        pivot_wider(id_cols = c(Community, Isolate1, Isolate2), names_from = Isolate1InitialODFreq, values_from = FitnessChange) %>%
        rename(FromRare = `5`, FromMedium = `50`, FromAbundant = `95`) %>%
        left_join(interaction_type) %>%
        select(Community, Isolate1, Isolate2, InteractionType, InteractionTypeFiner, FitnessFunction)

}
interaction_type <- make_interaction_type()

pairs_boots_table <- read_pairs_boots_table()
pairs_fitness <- pairs_boots_table %>%
    compute_pairs_fitness() %>%
    append_pairs_outcome() %>%
    ungroup() %>%
    select(Community, Isolate1, Isolate2, FitnessFunction)

pairs_values <- pairs_freq %>%
    group_by(PairID) %>%
    select(-PairFreqID) %>%
    filter(Time == "T8") %>%
    left_join(select(pairs_fitness, Community, Isolate1, Isolate2, FitnessFunction)) %>%
    pivot_wider(id_cols = -Isolate1CFUFreqSd, names_from = Isolate1InitialODFreq, names_prefix = "F", values_from = Isolate1CFUFreqMean) %>%
    ungroup() %>%
    select(-Time)


pairs_interaction <- pairs_values %>%
    mutate(InteractionType = case_when(
        (F5 == 0 & F50 == 0 & F95 == 0) ~ "exclusion",
        (F5 == 1 & F50 == 1 & F95 == 1) ~ "exclusion",
        FitnessFunction == "-1_-1_-1" ~ "exclusion",
        FitnessFunction == "1_1_1" ~ "exclusion",
        (F5 != 0 & F50 != 0 & F95 != 0) & (F5 != 1 & F50 != 1 & F95 != 1) & FitnessFunction %in% c("1_1_-1", "1_0_-1", "1_-1_-1") ~ "coexistence",
        (F5 != 0 & F50 != 0 & F95 != 0) & (F5 != 1 & F50 != 1 & F95 != 1) ~ "coexistence",
        T ~ "inconclusive"
    )) %>%
    mutate(InteractionTypeFiner = case_when(
        (F5 == 0 & F50 == 0 & F95 == 0) ~ "1-exclusion",
        (F5 == 1 & F50 == 1 & F95 == 1) ~ "1-exclusion",
        FitnessFunction == "-1_-1_-1" ~ "2-exclusion",
        FitnessFunction == "1_1_1" ~ "2-exclusion",
        (F5 != 0 & F50 != 0 & F95 != 0) & (F5 != 1 & F50 != 1 & F95 != 1) & FitnessFunction %in% c("1_1_-1", "1_0_-1", "1_-1_-1") ~ "3-coexistence",
        (F5 != 0 & F50 != 0 & F95 != 0) & (F5 != 1 & F50 != 1 & F95 != 1) ~ "4-coexistence",
        T ~ "5-inconclusive"
    )) %>%
    mutate(FlipLoser = case_when(
        (F5 == 0 & F50 == 0 & F95 == 0) | FitnessFunction == "-1_-1_-1"  ~ F,
        (F5 == 1 & F50 == 1 & F95 == 1) | FitnessFunction == "1_1_1" ~ T,
        FitnessFunction %in% c("1_1_-1", "1_0_-1", "1_-1_-1") ~ F,
        InteractionType %in% c("3-coexistence", "4-coexistence", "5-inconclusive") ~ F
    )) %>%
    # From and To
    mutate(From = case_when(
        FitnessFunction %in% c("1_1_1") ~ Isolate1, # Isolate1 wins
        FitnessFunction %in% c("-1_-1_-1") ~ Isolate2, # Isolate2 wins
        TRUE ~ Isolate1
    )) %>%
    mutate(To = case_when(
        FitnessFunction %in% c("1_1_1") ~ Isolate2,
        FitnessFunction %in% c("-1_-1_-1") ~ Isolate1,
        TRUE ~ Isolate2
    ))

write_csv(pairs_interaction, paste0(folder_data, "temp/93a-pairs_interaction.csv"))





















