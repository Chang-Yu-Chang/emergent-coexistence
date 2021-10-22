#' Read and match isolates' information
library(tidyverse)
library(data.table)
communities <- read_csv(here::here("data/output/communities.csv"))
communities_name <- communities$Community

# Read data ----
## From 01A
isolates_ID_match <- read_csv(here::here("data/temp/isolates_ID_match.csv"))

## From 01B
isolates_RDP <- read_csv(here::here("data/temp/isolates_RDP.csv"))
#load(here::here("data/temp/isolates_sanger_seq.Rdata")) # it has four objects: aln, aln2, tree, and isolates_seq

## From 01C
isolates_epsilon <- read_csv(here::here("data/temp/isolates_epsilon.csv"))

if (FALSE) {
## From 01D
isolates_growth_rate_ID <- read_csv(here::here("data/temp/isolates_growth_rate_ID.csv"))
isolates_growth_rate_spread_ID <- read_csv(here::here("data/temp/isolates_growth_rate_spread_ID.csv"))

## From 02D
isolates_tournament <- read_csv(here::here("data/temp/isolates_tournament.csv"))

}


# Match isolates' information ----
## 100 isolates. 68 from djordje, 16 from jean
isolates <- isolates_ID_match %>%
  left_join(isolates_RDP, by = c("ID")) %>%
  left_join(isolates_epsilon, by = c("Community", "Isolate")) %>%
  # left_join(isolates_tournament, by = c("Community", "Isolate")) %>%
  # left_join(isolates_growth_rate_spread_ID, by = c("ExpID", "ID", "Community", "Isolate")) %>%
  mutate(Community = ordered(Community, levels = communities_name)) %>%
  filter(!(!grepl("JVN", ExpID) & is.na(Community))) %>%
  as_tibble()

# Melted isolates df for growth curve plot
isolates_melted <- isolates_ID_match %>%
  left_join(isolates_RDP, by = c("ID")) %>%
  left_join(isolates_epsilon, by = c("Community", "Isolate")) %>%
  # left_join(isolates_growth_rate_ID, by = c("ExpID", "ID", "Community", "Isolate")) %>%
  mutate(Community = ordered(Community, levels = communities_name)) %>%
  as_tibble()


fwrite(isolates, here::here("data/output/isolates.csv"))
fwrite(isolates_melted, here::here("data/output/isolates_melted.csv"))


