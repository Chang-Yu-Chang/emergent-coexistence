#' This script reads the measured T0 and T8 pairwise data and bootstrap it
#' 1. Bootstrap T0 CFU frequency. Output temp/07-pairs_T0_boots.csv
#' 2. Bootstrap T8 CFU frequencies. Output temp/07-pairs_T8_boots.csv

library(tidyverse)
library(cowplot)
source(here::here("processing_scripts/00-metadata.R"))

pairs_T0 <- read_csv(paste0(folder_data, "temp/06-pairs_T0.csv"), show_col_types = F)
pairs_T8 <- read_csv(paste0(folder_data, "temp/06-pairs_T8.csv"), show_col_types = F)

# 1. Bootstrap T0 freq_A from Poisson ----
n_bootstraps = 1000
set.seed(9)
pairs_T0_boots <- pairs_T0 %>%
    select(Batch, Community, Isolate1, Isolate2, Isolate1InitialODFreq, Isolate2InitialODFreq, cfu_A, cfu_B) %>%
    mutate(Time = "T0", RawDataType = "ODtoCFU") %>%
    rowwise() %>%
    mutate(bootstrap = list(
        tibble(BootstrapID = 1:n_bootstraps,
               n_A = rpois(n_bootstraps, cfu_A),
               n_B = rpois(n_bootstraps, cfu_B),
               Isolate1CFUFreq = n_A / (n_A + n_B))
    )) %>%
    unnest(cols = c(bootstrap)) %>%
    select(Batch, Community, Isolate1, Isolate2, Isolate1InitialODFreq,
           Time, RawDataType, BootstrapID, Isolate1CFUFreq)


# 2. Bootstrap T8 freq_A from Poisson -----
pairs_T8_boots <- pairs_T8 %>%
    mutate(Isolate2Count = TotalCount - Isolate1Count) %>%
    select(Batch, Community, Isolate1, Isolate2, Isolate1InitialODFreq, Isolate1Count, Isolate2Count) %>%
    mutate(Time = "T8", RawDataType = "CFU") %>%
    rowwise() %>%
    mutate(bootstrap = list(
        tibble(BootstrapID = 1:n_bootstraps,
               n_A = rpois(n_bootstraps, Isolate1Count*2), # To avoid drawing 0 for small colony count
               n_B = rpois(n_bootstraps, Isolate2Count*2),
               Isolate1CFUFreq = n_A / (n_A + n_B))
    )) %>%
    unnest(cols = c(bootstrap)) %>%
    select(Batch, Community, Isolate1, Isolate2, Isolate1InitialODFreq,
           Time, RawDataType, BootstrapID, Isolate1CFUFreq)


pairs_boots <- bind_rows(pairs_T0_boots, pairs_T8_boots) %>%
    select(Community, Isolate1, Isolate2, Isolate1InitialODFreq, Time, BootstrapID, Isolate1CFUFreq)

nrow(pairs_T0_boots)/1000 # 159*3=477
nrow(pairs_T8_boots)/1000 # 477-18=459
nrow(pairs_boots)/1000 # 477+459=936

write_csv(pairs_boots, paste0(folder_data, "temp/07-pairs_boots.csv"))




