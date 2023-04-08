#' This script compares the random forest prediction to the human eye counts

library(tidyverse)
library(cowplot)
source(here::here("processing_scripts/00-metadata.R"))

pairs_freq_ID <- read_csv(paste0(folder_data, "temp/00c-pairs_freq_ID.csv"), show_col_types = F)
isolates_duplicate <- tibble(ID = as.character(c(462, 355, 356, 461, 452, 446, 305, 435, 444, 348, 460, 454)), Duplicated = T)
isolates_ID <- read_csv(paste0(folder_data, "temp/00c-isolates_ID.csv"), show_col_types = F) %>%
    left_join(isolates_duplicate, by = "ID") %>%
    select(ID, Community, Isolate, Duplicated) %>%
    replace_na(list(Duplicated = F))

# Human result
pairs_freq_human <- read_csv(paste0(folder_data, "raw/pairs_freq_human.csv"), show_col_types = F) # human-eye results
pairs_freq_human <- pairs_freq_human %>%
    mutate(Experiment = str_replace(Experiment, "Transitivity_", "")) %>%
    select(Batch = Experiment, Community, Isolate1, Isolate2, Isolate1InitialODFreq = Isolate1Freq, Isolate1Count = ColonyCount1, TotalCount = ColonyCount) %>%
    group_by(Batch, Community, Isolate1, Isolate2, Isolate1InitialODFreq) %>%
    mutate(Isolate1CFUFreq = Isolate1Count / TotalCount) %>%
    mutate(Type = "human") %>%
    ungroup()

# Machine result
pairs_freq_machine <- read_csv(paste0(folder_data, "temp/06-pairs_T8.csv"), show_col_types = F) # random forest classification
pairs_freq_machine <- pairs_freq_machine %>%
    select(Batch, Community, Isolate1, Isolate2, Isolate1InitialODFreq, Isolate1Count, TotalCount, Isolate1CFUFreq) %>%
    mutate(Type = "machine")

pairs_freq_machine_human <- pairs_freq_machine %>%
    bind_rows(pairs_freq_human) %>%
    pivot_wider(names_from = Type, values_from = c(Isolate1Count, TotalCount, Isolate1CFUFreq)) %>%
    left_join(pairs_freq_ID) %>%
    select(-image_name_isolate1, -image_name_isolate2) %>%
    # Label pairs containing duplicated
    left_join(rename(isolates_ID, ID1 = ID, Isolate1 = Isolate, Duplicated1 = Duplicated)) %>%
    left_join(rename(isolates_ID, ID2 = ID, Isolate2 = Isolate, Duplicated2 = Duplicated)) %>%
    mutate(PairType = case_when(
        Duplicated1 == F & Duplicated2 == F ~ "clean",
        Duplicated1 == T & Duplicated2 == F ~ "one duplicate",
        Duplicated1 == F & Duplicated2 == T ~ "one duplicate",
        Duplicated1 == T & Duplicated2 == T ~ "both duplicate"
    )) %>%
    mutate(PairType = factor(PairType, c("clean", "one duplicate", "both duplicate"))) %>%
    select(-ID1, -ID2, -starts_with("folder")) %>%
    select(image_name_pair, everything())

write_csv(pairs_freq_machine_human, paste0(folder_data, "temp/27-pairs_freq_machine_human.csv"))
