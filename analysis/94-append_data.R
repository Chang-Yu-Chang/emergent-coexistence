#' This script append the competition outcome with other isolate features to
#' generate the meta table
#'
#' 0. communities
#' 1. isolates
#' 2. pairs
#' 3. example pairs for plotting the frequencies

library(tidyverse)
source(here::here("analysis/00-metadata.R"))

# 0. Communities ----
communities <- read_csv(paste0(folder_data, "temp/00c-communities.csv"), show_col_types = F)

# 1. Isolate metadata ----
isolates_ID <- read_csv(paste0(folder_data, "temp/00c-isolates_ID.csv"), show_col_types = F)
isolates_RDP <- read_csv(paste0(folder_data, "temp/12-isolates_RDP.csv"), show_col_types = F)
isolates_epsilon <- read_csv(paste0(folder_data, "temp/06-isolates_epsilon.csv"), show_col_types = F)
isolates_abundance <- read_csv(paste0(folder_data, "temp/32-isolates_abundance.csv"), show_col_types = F)
isolates_tournament <- read_csv(paste0(folder_data, "temp/93-isolates_tournament.csv"), show_col_types = F)

isolates <- isolates_ID %>%
    left_join(select(isolates_RDP, -ID), by = c("ExpID", "Community", "Isolate")) %>%
    left_join(isolates_epsilon, by = c("Community", "Isolate")) %>%
    left_join(isolates_tournament, by = c("Community", "Isolate")) %>%
    left_join(isolates_abundance, by = join_by(ExpID, ID, Community, Isolate, Sequence, Family, Genus)) %>%
    mutate(Community = ordered(Community, levels = communities$Community))

write_csv(isolates, paste0(folder_data, "output/isolates.csv"))

# 2. pairs metadata ----
pairs_ID <- read_csv(paste0(folder_data, "temp/00c-pairs_ID.csv"), show_col_types = F)
pairs_mismatch <- read_csv(paste0(folder_data, "temp/15-pairs_mismatch.csv"), show_col_types = F)
pairs_RDP <- read_csv(paste0(folder_data, "temp/16-pairs_RDP.csv"), show_col_types = F)
pairs_accuracy <- read_csv(paste0(folder_data, "temp/91-pairs_accuracy.csv"), show_col_types = F)
#pairs_interaction <- read_csv(paste0(folder_data, "temp/93a-pairs_interaction.csv"), show_col_types = F)
pairs_outcome <- read_csv(paste0(folder_data, "temp/93a-pairs_outcome.csv"), show_col_types = F)
pairs_freq <- read_csv(paste0(folder_data, "temp/93a-pairs_freq.csv"), show_col_types = F)

pairs <- pairs_ID %>%
    left_join(pairs_mismatch, by = c("Community", "Isolate1", "Isolate2")) %>%
    left_join(pairs_RDP, by = c("Batch", "Community", "Isolate1", "Isolate2", "ID1", "ID2")) %>%
    left_join(pairs_accuracy, by = c("Community", "Isolate1", "Isolate2")) %>%
    left_join(pairs_outcome, by = c("PairID", "Community", "Isolate1", "Isolate2")) %>%
    select(PairID, Community, Isolate1, Isolate2, From, To,
           ExpID1, ID1, Fermenter1, GramPositive1, Family1, Genus1, GenusScore1, Sequence1,
           ExpID2, ID2, Fermenter2, GramPositive2, Family2, Genus2, GenusScore2, Sequence2,
           PairFermenter, PairFamily, Mismatch, AccuracyMean, AccuracySd,
           outcome, Isolate1IsLoser)

write_csv(pairs, paste0(folder_data, "output/pairs.csv"))
