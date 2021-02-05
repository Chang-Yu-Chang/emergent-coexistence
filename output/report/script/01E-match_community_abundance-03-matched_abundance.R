#' Filter for matched ESV-isolate and plot the relative abundances
library(tidyverse)
library(data.table)

communities <- fread(here::here("data/temp/communities.csv"))
communities_name <- communities$Community
communities_size <- communities$CommunitySize
communities_name_pool <- c(paste0("C", 1:12, "Rpool"), paste0("C", rep(1:12, each = 8), "R", rep(1:8, 12)))
families_name <- c("Aeromonadaceae", "Alcaligenaceae", "Bradyrhizobiaceae", "Brucellaceae", "Burkholderiaceae", "Caulobacteraceae", "Cellulomonadaceae", "Chitinophagaceae", "Chthoniobacteraceae", "Comamonadaceae", "Cryomorphaceae", "Enterobacteriaceae", "Enterococcaceae", "Flavobacteriaceae", "Hyphomicrobiaceae", "Listeriaceae", "Microbacteriaceae", "Moraxellaceae", "Nocardiaceae", "Obscuribacterales.17", "Oxalobacteraceae", "Paenibacillaceae", "Phyllobacteriaceae", "Porphyromonadaceae", "Pseudomonadaceae", "Rhizobiaceae", "Sanguibacteraceae", "Sphingobacteriaceae", "Sphingomonadaceae", "Xanthomonadaceae")

# Read data
## Sanger to ESV alignment
sequences_alignment_syl <- fread(here::here("data/temp/sequences_alignment_syl.csv")) %>%
  mutate(Community = ordered(Community, levels = communities_name_pool))
  #mutate(Family = factor(Family, families_name)) %>% as_tibble()

# ESV abundance in community
communities_abundance_syl <- fread(here::here("data/temp/communities_abundance_syl.csv")) %>%
  mutate(Community = ordered(Community, levels = communities_name_pool)) %>%
  mutate(Family = factor(Family, families_name)) %>% as_tibble()

# Find the match with highest alignment score; 68 rows
sequences_abundance_list <- rep(list(NA), 3)
allow_mismatch <- c(0:2, Inf)

for (i in 1:4) {
  sequences_abundance_list[[i]] <-
    sequences_alignment_syl %>%
    filter(AlignmentType == "local") %>%
    # Filter for BasePairMatch
    filter(BasePairMismatch <= allow_mismatch[i]) %>%
    # For each Sanger, find the Sanger-ESV match with highest alignment score
    group_by(AlignmentType, ExpID) %>%
    arrange(desc(AlignmentScore)) %>%
    dplyr::slice(1) %>%
    ungroup() %>%
    # Remove duplicates that matches two Sangers to one ESV
    group_by(AlignmentType, Community, CommunityESVID) %>%
    distinct(RelativeAbundance, .keep_all = T) %>%
    arrange(Community) %>%
    # Specify mismatch allowed
    mutate(AllowMismatch = allow_mismatch[i])
}

sequences_abundance <- rbindlist(sequences_abundance_list) %>% as_tibble


# Remove unessential variables
sequences_abundance <- sequences_abundance %>%
  mutate(AlignmentType = ordered(AlignmentType, levels = c("global", "local", "overlap", "global-local", "local-global"))) %>%
  select(AlignmentType, AllowMismatch, SampleID, Community, ExpID, IsolateGenus,
    RelativeAbundance, CommunityESVID, Family, ConsensusLength, BasePairGap, BasePairMismatch, AlignmentScore) %>%
  arrange(AlignmentType, AllowMismatch, Community)


# Isolate abundances
isolates_abundance <-
  sequences_abundance %>%
  filter(AlignmentType == "local", AllowMismatch == 2) %>%
  select(Community, IsolateGenus, RelativeAbundance) %>%
  right_join(isolates_ID_match)

fwrite(sequences_abundance, file = here::here("data/temp/sequences_abundance.csv"))
fwrite(isolates_abundance, file = here::here("data/temp/isolates_abundance.csv"))







