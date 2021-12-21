#' Read and match isolates' information
#' Note that `isolates_tournament` is from 02D
library(tidyverse)
communities <- read_csv(here::here("data/output/communities.csv"))

# Read data ----
## From 01A
isolates_ID_match <- read_csv(here::here("data/temp/isolates_ID_match.csv"))

## From 01B
isolates_RDP <- read_csv(here::here("data/temp/isolates_RDP.csv"))
#load(here::here("data/temp/isolates_sanger_seq.Rdata")) # it has four objects: aln, aln2, tree, and isolates_seq

## From 01C
isolates_epsilon <- read_csv(here::here("data/temp/isolates_epsilon.csv"))

## From 01D
isolates_growth_traits <- read_csv(here::here("data/temp/isolates_growth_traits.csv")) %>%
    select(-Assembly, -ExpID, -Family, -Genus)

## From 02D
isolates_tournament <- read_csv(here::here("data/temp/isolates_tournament.csv"))

# Match isolates' information ----
## 68 isolates
isolates <- isolates_ID_match %>%
    select(Assembly, ExpID, ID, Community, Isolate) %>%
    left_join(isolates_RDP, by = c("ID")) %>%
    left_join(isolates_epsilon, by = c("Community", "Isolate")) %>%
    left_join(isolates_tournament, by = c("Community", "Isolate")) %>%
    select(-OD620) %>%
    left_join(isolates_growth_traits, by = c("ID", "Community", "Isolate")) %>%
    mutate(Community = ordered(Community, levels = communities$Community)) %>%
    filter(!(!grepl("JVN", ExpID) & is.na(Community))) %>%
    select(Assembly, everything()) %>%
    as_tibble()

write_csv(isolates, here::here("data/output/isolates.csv"))

