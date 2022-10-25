#' This script append the competition outcome with other isolate features to
#' generate the meta table
#' 1. isolates
#' 2. pairs
#' 3. example pairs for plotting the frequencies

library(tidyverse)
communities <- read_csv("~/Dropbox/lab/emergent-coexistence/data/output/communities.csv", col_types = cols())

folder_script <- "~/Desktop/lab/emergent-coexistence/analysis/"
folder_main <- "~/Dropbox/lab/emergent-coexistence/plate_scan_pipeline/"

pairs_freq_ID <- read_csv(paste0(folder_main, "meta/00-pairs_freq_ID.csv"), show_col_types = F)
pairs_ID <- distinct(pairs_freq_ID, Batch, Community, Isolate1, Isolate2)


# 1. Isolate metadata ----
## From 01A
isolates_ID_match <- read_csv("~/Dropbox/lab/emergent-coexistence/data/temp/isolates_ID_match.csv", col_types = cols())

## From 01B
isolates_RDP <- read_csv("~/Dropbox/lab/emergent-coexistence/data/temp/isolates_RDP.csv", col_types = cols())
#load("~/Dropbox/lab/emergent-coexistence/data/temp/isolates_sanger_seq.Rdata")) # it has four objects: aln, aln2, tree, and isolates_seq

## From 93-
isolates_epsilon <- read_csv(paste0(folder_main, "meta/93-isolates_epsilon.csv"), col_types = cols())

## From 01D
isolates_growth_traits <- read_csv("~/Dropbox/lab/emergent-coexistence/data/temp/isolates_growth_traits.csv", col_types = cols()) %>%
    select(-Assembly, -ExpID)

## From 01E
isolates_abundance <- read_csv("~/Dropbox/lab/emergent-coexistence/data/temp/isolates_abundance.csv", col_types = cols())

## From 95-
#' Note that `isolates_tournament` is from 02D
isolates_tournament <- read_csv(paste0(folder_main, "meta/95-isolates_tournament.csv"), col_types = cols())

# Match isolates' information of 68 isolates
isolates <- isolates_ID_match %>%
    select(Assembly, ExpID, ID, Community, Isolate) %>%
    left_join(isolates_RDP, by = c("ExpID")) %>%
    left_join(isolates_epsilon, by = c("Community", "Isolate")) %>%
    left_join(isolates_tournament, by = c("Community", "Isolate")) %>%
    left_join(isolates_abundance, by = c("Assembly", "ExpID", "Community", "Isolate")) %>%
    #select(-OD620) %>%
    #left_join(isolates_growth_traits, by = c("ID", "Community", "Isolate")) %>%
    mutate(Community = ordered(Community, levels = communities$Community)) %>%
    filter(!(!grepl("JVN", ExpID) & is.na(Community))) %>%
    select(Assembly, everything()) %>%
    as_tibble()

write_csv(isolates, paste0(folder_main, "meta/97-isolates.csv"))

# 2. pairs metadata ----

#' Read and combine data for pairwise competition
#' 1. `pairs.csv` pairwise competition only
#' 2. `pairs_meta.csv` append the isolate' metabolic trait
communities <- read_csv("~/Dropbox/lab/emergent-coexistence/data/output/communities.csv", col_types = cols())

## From 95
pairs_freq <- read_csv(paste0(folder_main, "meta/95-pairs_freq.csv"), show_col_types = F)
pairs_interaction <- read_csv(paste0(folder_main, "meta/95-pairs_interaction.csv"), show_col_types = F)

## From 02E
pairs_taxonomy <- read_csv("~/Dropbox/lab/emergent-coexistence/data/temp/pairs_taxonomy.csv", col_types = cols())
pairs_16S <- read_csv("~/Dropbox/lab/emergent-coexistence/data/temp/pairs_16S.csv", col_types = cols())

## Combine data by assigning pairs information.186 pairs
pairs <- pairs_taxonomy %>%
    #left_join(pairs_ID, by = c("Community", "Isolate1", "Isolate2")) %>%
    left_join(pairs_interaction, by = c("Community", "Isolate1", "Isolate2")) %>%
    left_join(pairs_16S, by = c("PairID", "Community", "Isolate1", "Isolate2")) %>%
    select(PairID, Community, Isolate1, Isolate2, From, To,
           ExpID1, ID1, Fermenter1, GramPositive1, Family1, Genus1, GenusScore1, Sequence1,
           ExpID2, ID2, Fermenter2, GramPositive2, Family2, Genus2, GenusScore2, Sequence2,
           PairFermenter, PairFamily, InteractionType, InteractionTypeFiner, FitnessFunction,
           SequenceDifference, SequenceLength)

# Flip the sign so that Isolate1 is the dominant and Isolate2 is subdominant
#' 50-50 frequency changes. To determine the dominance in coexistence pairs.
#' The isolate that increases frequency in the 50-50 competition is the dominant strain
pairs_coexist_dominant <- pairs %>%
    filter(InteractionType == "coexistence") %>%
    select(Community, Isolate1, Isolate2) %>%
    left_join(pairs_freq) %>%
    filter(Isolate1InitialODFreq == 50) %>%
    select(Community, Isolate1, Isolate2, Time, Isolate1CFUFreqMean) %>%
    pivot_wider(names_from = "Time", values_from = "Isolate1CFUFreqMean") %>%
    mutate(Isolate1Dominant = T8 > T0) %>%
    select(Community, Isolate1, Isolate2, Isolate1Dominant)

pairs <- left_join(pairs, pairs_coexist_dominant)

#
write_csv(pairs, paste0(folder_main, "meta/97-pairs.csv"))


# 3. example pairs ----
#' Subset the example pairs
communities <- read_csv("~/Dropbox/lab/emergent-coexistence/data/output/communities.csv", col_types = cols())
pairs_freq <- read_csv(paste0(folder_main, "meta/95-pairs_freq.csv"), show_col_types = F)
pairs <- read_csv(paste0(folder_main, "meta/97-pairs.csv"), show_col_types = F)

# Example for plotting interaction types ----
## Stable coexistence and coexistence ain examples
pairs_example_outcomes <- pairs %>%
    filter((Community == "C1R2" & Isolate1 == 1 & Isolate2 == 2) |
               (Community == "C1R7" & Isolate1 == 3 & Isolate2 == 6) |
               (Community == "C2R8" & Isolate1 == 1 & Isolate2 == 2))


## Examples of outcomes in a finer scale
pairs_example_outcomes_finer <- pairs %>%
    filter((Community == "C1R2" & Isolate1 == 1 & Isolate2 == 2) |
               (Community == "C1R7" & Isolate1 == 3 & Isolate2 == 6) |
               (Community == "C2R8" & Isolate1 == 1 & Isolate2 == 2) |
               (Community == "C11R1" & Isolate1 == 1 & Isolate2 == 5) |
               (Community == "C2R6" & Isolate1 == 1 & Isolate2 == 4))

#
write_csv(pairs_example_outcomes, paste0(folder_main, "meta/97-pairs_example_outcomes.csv"))
write_csv(pairs_example_outcomes_finer, paste0(folder_main, "meta/97-pairs_example_outcomes_finer.csv"))

