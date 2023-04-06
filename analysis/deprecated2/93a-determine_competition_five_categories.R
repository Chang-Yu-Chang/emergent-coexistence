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

pairs_boots <- bind_rows(pairs_T0_boots, pairs_T8_boots) %>%
    select(Community, Isolate1, Isolate2, Isolate1InitialODFreq, Time, BootstrapID, Isolate1CFUFreq)

write_csv(pairs_boots, paste0(folder_data, "temp/93a-pairs_boots.csv"))

# Find 5% and 95%
pairs_boots_percentile <- pairs_boots %>%
    group_by(Community, Isolate1, Isolate2, Isolate1InitialODFreq, Time) %>%
    arrange(Isolate1CFUFreq) %>%
    mutate(Percentile = paste0("Percentile", 0.1 * (1:1000))) %>%
    filter(Percentile %in% c("Percentile5", "Percentile95")) %>%
    ungroup() %>%
    select(Community, Isolate1, Isolate2, Isolate1InitialODFreq, Time, Isolate1CFUFreq, Measure = Percentile) %>%
    pivot_wider(names_from = Measure, names_prefix = "Isolate1CFUFreq", values_from = Isolate1CFUFreq)

# Find mean and median
pairs_boots_mean <- pairs_boots %>%
    group_by(Community, Isolate1, Isolate2, Isolate1InitialODFreq, Time) %>%
    summarize(Isolate1CFUFreqMean = mean(Isolate1CFUFreq),
              Isolate1CFUFreqMedian = median(Isolate1CFUFreq)) %>%
    ungroup()


pairs_freq <- left_join(pairs_boots_mean, pairs_boots_percentile) %>%
    mutate(Community = factor(Community, paste0("C", rep(1:12, each = 8), "R", rep(1:8, 12)))) %>%
    arrange(Community, Isolate1, Isolate2, Isolate1InitialODFreq) %>%
    ungroup() %>%
    left_join(pairs_ID) %>%
    arrange(Time, PairID, Isolate1InitialODFreq) %>%
    #mutate(PairFreqID = 1:n()) %>%
    select(Batch, PairID, everything())

write_csv(pairs_freq, paste0(folder_data, "temp/93a-pairs_freq.csv"))

# Pair outcomes from Djordje's scripts ----
pairs_outcome <- read_csv(paste0(folder_data, "raw/pairs_outcome.csv"), show_col_types = F) %>%
    left_join(pairs_ID)

# Determine direction of exclusion
pairs_exclusion <- pairs_outcome %>%
    filter(outcome %in% c("1-exclusion", "2-exclusion"))
pairs_freq_exclusion <- pairs_freq %>%
    select(PairID, Community, Isolate1, Isolate2, Isolate1InitialODFreq, Time, Isolate1CFUFreqMean) %>%
    pivot_wider(names_from = Time, values_from = Isolate1CFUFreqMean) %>%
    drop_na(T8) %>%
    filter(PairID %in% pairs_exclusion$PairID) %>%
    mutate(Sign = sign(T8-T0)) %>%
    group_by(PairID) %>%
    summarize(Sign = case_when(all(Sign==1) ~ 1, all(Sign == -1) ~-1))

pairs_ID_winner <- pairs_freq_exclusion$PairID[pairs_freq_exclusion$Sign == 1]
pairs_ID_loser <- pairs_freq_exclusion$PairID[pairs_freq_exclusion$Sign == -1]
pairs_ID_others <- pairs_outcome$PairID[pairs_outcome$outcome %in% c("3-coexistence", "4-coexistence", "5-inconclusive")]

pairs_outcome <- pairs_outcome %>%
    mutate(From = case_when( # For network
        PairID %in% pairs_ID_winner ~ Isolate1,
        PairID %in% pairs_ID_loser ~ Isolate2,
        PairID %in% pairs_ID_others ~ Isolate1,
    )) %>%
    mutate(To = case_when( # For network
        PairID %in% pairs_ID_winner ~ Isolate2,
        PairID %in% pairs_ID_loser ~ Isolate1,
        PairID %in% pairs_ID_others ~ Isolate2,
    )) %>% # For plotting the frequency
    mutate(Isolate1IsLoser = case_when(
        PairID %in% pairs_ID_winner ~ F,
        PairID %in% pairs_ID_loser ~ T,
        PairID %in% pairs_ID_others ~ NA,
    ))


write_csv(pairs_outcome, paste0(folder_data, "temp/93a-pairs_outcome.csv"))






