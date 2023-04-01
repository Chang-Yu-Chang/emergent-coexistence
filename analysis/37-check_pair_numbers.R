library(tidyverse)
library(cowplot)
source(here::here("analysis/00-metadata.R"))

isolates <- read_csv(paste0(folder_data, "output/isolates.csv"), show_col_types = F)
pairs <- read_csv(paste0(folder_data, "output/pairs.csv"), show_col_types = F)
communities <- read_csv(paste0(folder_data, "temp/00c-communities.csv"), show_col_types = F)
#communities_remained <- read_csv(paste0(folder_data, "output/communities_remained.csv"), show_col_types = F)

# Check number of pairs after removing the four isolates that do not match to ESV
nrow(pairs) # 186 total number of pairs
isolates_removal <- isolates$ExpID[which(is.na(isolates$BasePairMismatch))] # Isolates that do not match to ESVs
pairs_remained <- pairs %>%
    filter(!(ExpID1 %in% isolates_removal) & !(ExpID2 %in% isolates_removal))
nrow(pairs_remained) # 160 pairs


# Expected number of pairs
count_pairs <- function(x) choose(x,2)
isolates %>%
    filter(!is.na(BasePairMismatch)) %>%
    group_by(Community) %>%
    count(name = "Count") %>%
    mutate(PairSize = count_pairs(Count)) %>%
    ungroup() %>%
    summarize(sum(PairSize)) # Expected 160 pairs given the number of isolates per communities


# Check the number of pairs afte removing the siz no-colony pairs and 9 pairs of low arracuary
pairs_remained <- pairs_remained %>% # 160 pairs
    arrange(outcome, PairID) %>%
    mutate(PairID = factor(PairID, unique(PairID))) %>%
    # Remove no-colony pairs. six pairs
    drop_na(outcome) %>% # Remaining 154
    # Remove low-accuracy model pairs. nine pairs
    filter(AccuracyMean > 0.9) # remaining 145 pairs

# Check the number of pairs per communities
pairs_remained %>%
    left_join(communities) %>%
    group_by(Community, CommunityLabel) %>%
    count(name = "Count") %>%
    arrange(CommunityLabel)
