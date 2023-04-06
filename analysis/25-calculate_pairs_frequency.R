#' This scripts read the boostrapped pairs T0 and T8 results into a pairs_freq table
#' and calculate the mean, 5th and 95th percentile of the pairwise data

library(tidyverse)
library(cowplot)
source(here::here("analysis/00-metadata.R"))

pairs_ID <- read_csv(paste0(folder_data, "temp/00c-pairs_ID.csv"), show_col_types = F)
pairs_freq_ID <- read_csv(paste0(folder_data, "temp/00c-pairs_freq_ID.csv"), show_col_types = F) %>%
    mutate(PairFreqID = 1:n()) %>%
    select(PairFreqID, Community, Isolate1, Isolate2, Isolate1InitialODFreq)
pairs_boots <- read_csv(paste0(folder_data, "temp/07-pairs_boots.csv"), show_col_types = F)


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
    left_join(pairs_freq_ID) %>%
    arrange(PairFreqID, PairID, Isolate1InitialODFreq) %>%
    #mutate(PairFreqID = 1:n()) %>%
    select(PairFreqID, Batch, PairID, everything())

nrow(pairs_freq) # 1107 coclutures. 558 at T0 and 549 at T8
write_csv(pairs_freq, paste0(folder_data, "temp/25-pairs_freq.csv"))

