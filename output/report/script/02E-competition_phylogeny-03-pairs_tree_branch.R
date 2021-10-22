#' Compute the pairwise distances between two branch tips
library(tidyverse)
library(data.table)

isolates_ID_match <- read_csv(here::here("data/temp/isolates_ID_match.csv"))
pairs_ID <- read_csv(here::here("data/temp/pairs_ID.csv"))
load(here::here("data/temp/isolates_sanger_seq.Rdata"))
communities <- read_csv(here::here("data/output/communities.csv"))


# Read the pairwise distances between the pairs of tips from a phylogenetic tree using its branch lengths
pairs_tree_distance_matrix <- ape::cophenetic.phylo(tree)

# Extract the pairwise distances within community
pairs_tree_distance <- pairs_tree_distance_matrix %>%
  # Melt the pairwise distance matrix
  as.table() %>% as.data.frame() %>%
  setNames(c("ExpID1", "ExpID2", "PairTreeDistance")) %>%
  filter(ExpID1 %in% isolates_ID_match$ExpID, ExpID2 %in% isolates_ID_match$ExpID) %>%
  # Match isolates's identity
  left_join(select(mutate(isolates_ID_match, ExpID1 = ExpID, Community1 = Community, Isolate1 = Isolate), ExpID1, Community1, Isolate1)) %>%
  left_join(select(mutate(isolates_ID_match, ExpID2 = ExpID, Community2 = Community, Isolate2 = Isolate), ExpID2, Community2, Isolate2)) %>%
  as_tibble()


pairs_tree_distance <- pairs_tree_distance %>%
  filter(Community1 == Community2, Isolate1 != Isolate2, Isolate1 < Isolate2) %>%
  # Configurate variables
  mutate(Community = ordered(Community1, levels = communities$Community),
    Isolate1 = ordered(Isolate1, 1:12),
    Isolate2 = ordered(Isolate2, 1:12)) %>%
  select(Community, Isolate1, Isolate2, PairTreeDistance) %>%
  arrange(Community, Isolate1, Isolate2)


## Match all pairs, including the 6 pairs that are missing because of 2 unsequenced isolates
pairs_tree_distance <- pairs_ID %>%
  mutate(Community = ordered(Community, levels = communities$Community),
    Isolate1 = ordered(Isolate1, 1:12),
    Isolate2 = ordered(Isolate2, 1:12)) %>%
  left_join(pairs_tree_distance) %>%
  as_tibble()

#
fwrite(pairs_tree_distance, file = here::here("data/temp/pairs_tree_distance.csv"))
