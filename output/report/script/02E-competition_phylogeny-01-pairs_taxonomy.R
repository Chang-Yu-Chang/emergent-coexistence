#' Read and match isolates' taxonomy (family, genus, and fermenter) to pairs
library(tidyverse)

isolates_RDP <- read_csv("~/Dropbox/lab/emergent-coexistence/data/temp/isolates_RDP.csv", col_types = cols())
isolates_ID_match <- read_csv("~/Dropbox/lab/emergent-coexistence/data/temp/isolates_ID_match.csv", col_types = cols())
pairs_ID <- read_csv("~/Dropbox/lab/emergent-coexistence/data/temp/pairs_ID.csv", col_types = cols())

# Match the isolates taxonomic information to pairs
pairs_taxonomy <- pairs_ID %>%
    left_join(rename_with(isolates_ID_match, ~ paste0(., "1"), !contains("Community"))) %>%
    left_join(rename_with(isolates_ID_match, ~ paste0(., "2"), !contains("Community"))) %>%
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
        (Family1 == "Enterobacteriaceae" & Family2 == "Pseudomonadaceae") | (Family1 == "Pseudomonadaceae" & Family2 == "Enterobacteriaceae") ~"FR"
    ))

#
write_csv(pairs_taxonomy, "~/Dropbox/lab/emergent-coexistence/data/temp/pairs_taxonomy.csv")

