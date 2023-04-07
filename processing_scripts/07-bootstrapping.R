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
set.seed(1)
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
temp <- rep(list(NA), nrow(pairs_T8))
for (i in 1:nrow(pairs_T8)) {
    bootstrap_freq <- function (total_count, cfu_freq) sum(runif(total_count, 0, 1) < cfu_freq)/total_count
    cfu_freqs <- NULL
    for (k in 1:n_bootstraps)  cfu_freqs[k] <- bootstrap_freq(pairs_T8$TotalCount[i], pairs_T8$Isolate1CFUFreq[i])

    temp[[i]] <- tibble(
        Batch = pairs_T8$Batch[i],
        Community = pairs_T8$Community[i],
        Isolate1 = pairs_T8$Isolate1[i],
        Isolate2 = pairs_T8$Isolate2[i],
        Isolate1InitialODFreq = pairs_T8$Isolate1InitialODFreq[i],
        Time = "T8",
        RawDataType = "CFU",
        BootstrapID = 1:n_bootstraps,
        Isolate1CFUFreq = cfu_freqs
    )
    cat("\t", i)
}

pairs_T8_boots <- bind_rows(temp)

pairs_boots <- bind_rows(pairs_T0_boots, pairs_T8_boots) %>%
    select(Community, Isolate1, Isolate2, Isolate1InitialODFreq, Time, BootstrapID, Isolate1CFUFreq)

nrow(pairs_T0_boots)/1000 # 186*3=558
nrow(pairs_T8_boots)/1000 # 558-9=549 because 9 images do not have colony
nrow(pairs_boots)/1000 # 558+549=1107

write_csv(pairs_boots, paste0(folder_data, "temp/07-pairs_boots.csv"))




