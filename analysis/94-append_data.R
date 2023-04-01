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
#isolates_tournament <- read_csv(paste0(folder_data, "temp/93-isolates_tournament.csv"), show_col_types = F)

isolates <- isolates_ID %>%
    left_join(select(isolates_RDP, -ID), by = c("ExpID", "Community", "Isolate")) %>%
    left_join(isolates_epsilon, by = c("Community", "Isolate")) %>%
    #left_join(isolates_tournament, by = c("Community", "Isolate")) %>%
    left_join(mutate(isolates_abundance, ID = as.character(ID))) %>%
    mutate(Community = ordered(Community, levels = communities$Community))

write_csv(isolates, paste0(folder_data, "output/isolates.csv"))

# 2. pairs metadata ----
pairs_ID <- read_csv(paste0(folder_data, "temp/00c-pairs_ID.csv"), show_col_types = F)
pairs_mismatch <- read_csv(paste0(folder_data, "temp/15-pairs_mismatch.csv"), show_col_types = F) %>% select(-ID1, -ID2)#mutate(ID1 = as.character(ID1), ID2 = as.character(ID2))
pairs_RDP <- read_csv(paste0(folder_data, "temp/16-pairs_RDP.csv"), show_col_types = F)
pairs_accuracy <- read_csv(paste0(folder_data, "temp/91-pairs_accuracy.csv"), show_col_types = F)
#pairs_interaction <- read_csv(paste0(folder_data, "temp/93a-pairs_interaction.csv"), show_col_types = F)
pairs_outcome <- read_csv(paste0(folder_data, "temp/93a-pairs_outcome.csv"), show_col_types = F) # Djorje's result with link direction
pairs_freq <- read_csv(paste0(folder_data, "temp/93a-pairs_freq.csv"), show_col_types = F)

pairs <- pairs_ID %>%
    right_join(pairs_RDP, by = c("Batch", "Community", "Isolate1", "Isolate2")) %>%
    left_join(pairs_mismatch, by = c("Community", "Isolate1", "Isolate2")) %>%
    left_join(pairs_accuracy, by = c("Community", "Isolate1", "Isolate2")) %>%
    left_join(pairs_outcome, by = c("PairID", "Community", "Isolate1", "Isolate2")) %>%
    select(PairID, Community, Isolate1, Isolate2, From, To,
           ExpID1, ID1, Fermenter1, GramPositive1, Family1, Genus1, GenusScore1, Sequence1,
           ExpID2, ID2, Fermenter2, GramPositive2, Family2, Genus2, GenusScore2, Sequence2,
           PairFermenter, PairFamily, Mismatch, AccuracyMean, AccuracySd,
           outcome, Isolate1IsLoser)

write_csv(pairs, paste0(folder_data, "output/pairs.csv"))

# 3. Remove pairs containing the 4 isolates with bad ESV-Sanger alignments ----
# WE need this section to get updated pairs and to compute the isolate rank with correct number of isolates per community
nrow(pairs) # 186 pairs in pairwise competition
isolates_removal <- isolates$ExpID[which(is.na(isolates$BasePairMismatch))] # Isolates that do not match to ESVs
pairs_remained <- pairs %>%
    filter(!(ExpID1 %in% isolates_removal) & !(ExpID2 %in% isolates_removal))
nrow(pairs_remained) # 160 pairs
pairs_remained <- pairs_remained %>%
    # Remove no-colony pairs. six pairs
    drop_na(outcome) %>% # 154 pairs
    # Remove low-accuracy model pairs. nine pairs
    filter(AccuracyMean > 0.9) # 145 pairs
nrow(pairs_remained) # 145 pairs

write_csv(pairs_remained, paste0(folder_data, "output/pairs_remained.csv"))

# Remove the four isolates with bad sabger fgie  ----
tournament_rank <- function(pairs_comm) {
    isolate_name <- pairs_comm %>% select(Isolate1, Isolate2) %>% unlist %>% unique %>% sort()
    # Isolates' ranks in the tournament
    count_isolates <- function (pairs, outcomes = c("1-exclusion", "2-exclusion"), column = "From") {
        pull(filter(pairs_comm, outcome %in% outcomes), {{column}}) %>%
            factor(isolate_name) %>%
            table() %>%
            as.vector
    }
    tour_rank <- tibble(
        Isolate = isolate_name,
        # Win
        Win = count_isolates(pairs_comm, outcomes = c("1-exclusion", "2-exclusion"), column = "From"),
        # Draw; Note that I consider neturality and bistability as draw in the tournament
        Draw = count_isolates(pairs_comm, outcomes = c("3-coexistence", "4-coexistence"), column = "From") +
            count_isolates(pairs_comm, outcomes = c("3-coexistence", "4-coexistence"), column = "To"),
        # Lose
        Lose = count_isolates(pairs_comm, outcomes = c("1-exclusion", "2-exclusion"), column = "To"),
        # Inconclusive
        Inconclusive = count_isolates(pairs_comm, outcomes = c("5-inconclusive"), column = "From") +
            count_isolates(pairs_comm, outcomes = c("5-inconclusive"), column = "To"),
    )

    # Competition score
    tour_rank <- tour_rank %>%
        mutate(Game = Win + Lose + Draw + Inconclusive,
               Score = (Win - Lose + 0 * Draw)/Game) %>%
        arrange(desc(Score))

    # Calculate rank by score; same scores means the same ranks
    temp_score <- ordered(tour_rank$Score, levels = sort(unique(tour_rank$Score), decreasing = T))
    temp_score_table <- table(temp_score)
    temp <- NULL; temp_counter = 1
    for (i in 1:length(temp_score_table)) {
        temp <- c(temp, rep(temp_counter, temp_score_table[i]))
        temp_counter <- temp_counter + temp_score_table[i]
    }

    tour_rank$Rank <- temp
    tour_rank$PlotRank <- 1:nrow(tour_rank)
    return(tour_rank)
}
# pairs_comm <- pairs_remained %>%
#     filter(Community == "C2R6")
# tournament_rank(pairs_comm)

isolates_tournament <- communities %>%
    select(comm = Community, everything()) %>%
    rowwise() %>%
    mutate(pairs_comm = pairs_remained %>% filter(Community == comm) %>% list()) %>%
    mutate(tournaments_comm = pairs_comm %>% tournament_rank() %>% list()) %>%
    select(Community = comm, tournaments_comm) %>%
    unnest(cols = tournaments_comm)

# pairs_outcome %>%
#     filter(Community == "C11R2")  %>%
#     tournament_rank()

write_csv(isolates_tournament, paste0(folder_data, "temp/94-isolates_tournament.csv"))

# Remove the 4 isolates with bad ESV-Sanger aligments
nrow(isolates) # 68 pairs in pairwise competition
isolates_remained <- isolates %>%
    filter(!is.na(BasePairMismatch)) %>%
    left_join(isolates_tournament)
write_csv(isolates_remained, paste0(folder_data, "output/isolates_remained.csv"))



# Update the community size
count_pairs <- function(x) choose(x,2)
communities_remained <- isolates_remained %>%
    left_join(select(communities, Community, CommunityLabel)) %>%
    group_by(Community, CommunityLabel) %>%
    summarize(CommunitySize = n()) %>%
    mutate(CommunityPairSize = count_pairs(CommunitySize)) %>%
    arrange(CommunityLabel)
write_csv(communities_remained, paste0(folder_data, "output/communities_remained.csv"))







