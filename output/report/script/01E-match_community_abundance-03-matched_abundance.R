#' Pick the best-matched ESV-isolate pairs and plot the relative abundances
library(tidyverse)

isolates_ID_match <- read_csv("~/Dropbox/lab/emergent-coexistence/data/temp/isolates_ID_match.csv", col_types = cols())
communities <- read_csv("~/Dropbox/lab/emergent-coexistence/data/output/communities.csv", col_types = cols())

# Sanger to ESV alignment
sequences_alignment <- read_csv("~/Dropbox/lab/emergent-coexistence/data/temp/sequences_alignment.csv", col_types = cols()) %>%
    mutate(Community = factor(Community, levels = communities$Community))

# Find the match with highest alignment score; 68 rows
sequences_abundance_list <- rep(list(NA), 3)
allow_mismatch <- c(0:2, Inf)

for (i in 1:4) {
    sequences_abundance_list[[i]] <- sequences_alignment %>%
        filter(AlignmentType == "local") %>%
        # Filter for BasePairMatch
        filter(BasePairMismatch <= allow_mismatch[i]) %>%
        # For each Sanger, find the Sanger-ESV match with highest alignment score
        group_by(AlignmentType, ExpID) %>%
        arrange(desc(AlignmentScore)) %>%
        slice(1) %>%
        ungroup() %>%
        # For the scenario that two (or more) Sangers match to the same ESV, choose the match with the higher alignment score
        group_by(AlignmentType, Community, CommunityESVID) %>%
        arrange(desc(AlignmentScore)) %>%
        slice(1) %>%
        ungroup() %>%
        arrange(Community) %>%
        # Specify mismatch allowed
        mutate(AllowMismatch = allow_mismatch[i])
}

isolates_abundance <- bind_rows(sequences_abundance_list) %>%
    mutate(AlignmentType = factor(AlignmentType, levels = c("global", "local", "overlap", "global-local", "local-global"))) %>%
    arrange(AlignmentType, AllowMismatch, Community) %>%
    filter(AlignmentType == "local", AllowMismatch == "Inf") %>%
    # Select only necessary variables
    select(Community, ExpID, RelativeAbundance, CommunityESVID) %>%
    # Match it to isolate
    right_join(isolates_ID_match, by = c("Community", "ExpID")) %>%
    select(Assembly, Community, Isolate, ExpID, RelativeAbundance) %>%
    group_by(Community) %>%
    mutate(RankRelativeAbundance = rank(-RelativeAbundance))


write_csv(isolates_abundance, "~/Dropbox/lab/emergent-coexistence/data/temp/isolates_abundance.csv")







