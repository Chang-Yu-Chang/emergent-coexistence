#' Read and match isolates' information

library(tidyverse)
communities <- read_csv("~/Dropbox/lab/emergent-coexistence/data/output/communities.csv", col_types = cols())

## From 01A
isolates_ID_match <- read_csv("~/Dropbox/lab/emergent-coexistence/data/temp/isolates_ID_match.csv", col_types = cols())

## From 01B
isolates_RDP <- read_csv("~/Dropbox/lab/emergent-coexistence/data/temp/isolates_RDP.csv", col_types = cols())
#load("~/Dropbox/lab/emergent-coexistence/data/temp/isolates_sanger_seq.Rdata")) # it has four objects: aln, aln2, tree, and isolates_seq

## From 01C
isolates_epsilon <- read_csv("~/Dropbox/lab/emergent-coexistence/data/temp/isolates_epsilon.csv", col_types = cols())

## From 01D
isolates_growth_traits <- read_csv("~/Dropbox/lab/emergent-coexistence/data/temp/isolates_growth_traits.csv", col_types = cols()) %>%
    select(-Assembly, -ExpID, -Family, -Genus)

## From 02D
#' Note that `isolates_tournament` is from 02D
isolates_tournament <- read_csv("~/Dropbox/lab/emergent-coexistence/data/temp/isolates_tournament.csv", col_types = cols())

# Match isolates' information of 68 isolates
isolates <- isolates_ID_match %>%
    select(Assembly, ExpID, ID, Community, Isolate) %>%
    left_join(isolates_RDP, by = c("ExpID")) %>%
    left_join(isolates_epsilon, by = c("Community", "Isolate")) %>%
    left_join(isolates_tournament, by = c("Community", "Isolate")) %>%
    select(-OD620) %>%
    left_join(isolates_growth_traits, by = c("ID", "Community", "Isolate")) %>%
    mutate(Community = ordered(Community, levels = communities$Community)) %>%
    filter(!(!grepl("JVN", ExpID) & is.na(Community))) %>%
    select(Assembly, everything()) %>%
    as_tibble()

write_csv(isolates, "~/Dropbox/lab/emergent-coexistence/data/output/isolates.csv")

