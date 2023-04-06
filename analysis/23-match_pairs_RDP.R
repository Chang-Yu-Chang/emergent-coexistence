#' This script match the isolates' 16S and RDP to pairs

library(tidyverse)
source(here::here("analysis/00-metadata.R"))

isolates_ID <- read_csv(paste0(folder_data, "temp/00c-isolates_ID.csv"), show_col_types = F)
pairs_ID <- read_csv(paste0(folder_data, "temp/00c-pairs_ID.csv"), show_col_types = F)
isolates_RDP <- read_csv(paste0(folder_data, "temp/12-isolates_RDP.csv"), show_col_types = F)


# 1. Match the isolates taxonomic information to pairs ----
pairs_RDP <- pairs_ID %>%
    left_join(rename_with(isolates_ID, ~ paste0(., "1"), !contains("Community"))) %>%
    left_join(rename_with(isolates_ID, ~ paste0(., "2"), !contains("Community"))) %>%
    left_join(rename_with(isolates_RDP, ~ paste0(., "1"), !contains("Community"))) %>%
    left_join(rename_with(isolates_RDP, ~ paste0(., "2"), !contains("Community"))) %>%
    # Fermenter
    mutate(PairFermenter = case_when(
        Fermenter1 == T & Fermenter2 == T ~ "FF",
        Fermenter1 == F & Fermenter2 == F ~ "RR",
        (Fermenter1 == T & Fermenter2 == F) | (Fermenter1 == F & Fermenter2 == T) ~"FR"
    )) %>%
    # Family
    mutate(PairFamily = case_when(
        Family1 == "Enterobacteriaceae" & Family2 == "Pseudomonadaceae"~ "EE",
        Family1 == "Pseudomonadaceae" & Family2 == "Pseudomonadaceae" ~ "PP",
        (Family1 == "Enterobacteriaceae" & Family2 == "Pseudomonadaceae") | (Family1 == "Pseudomonadaceae" & Family2 == "Enterobacteriaceae") ~"FR",
        T ~ "others"
    )) %>%
    select(-PairID)

write_csv(pairs_RDP, paste0(folder_data, "temp/23-pairs_RDP.csv"))
















