#' Read and combine data for pairwise competition
library(tidyverse)
library(data.table)
communities <- read_csv(here::here("data/output/communities.csv"))

# Read pairwise data ----
pairs_ID <- read_csv(here::here("data/temp/pairs_ID.csv"))

## From 02D
pairs_freq <- read_csv(here::here("data/output/pairs_freq.csv"))
pairs_interaction <- read_csv(here::here("data/temp/pairs_interaction.csv"))

## From 02E
pairs_taxonomy <- read_csv(here::here("data/temp/pairs_taxonomy.csv"))
pairs_16S <- read_csv(here::here("data/temp/pairs_16S.csv"))
pairs_tree_distance <- read_csv(here::here("data/temp/pairs_tree_distance.csv"))

# Combine data by assigning pairs information ----
## 186 pairs
pairs <- pairs_ID %>%
  left_join(pairs_interaction, by = c("Community", "Isolate1", "Isolate2")) %>%
  left_join(pairs_taxonomy, by = c("Community", "Isolate1", "Isolate2")) %>%
  left_join(pairs_16S, by = c("Community", "Isolate1", "Isolate2")) %>%
  #left_join(pairs_tree_distance, by = c("Community", "Isolate1", "Isolate2")) %>%
  as_tibble() %>%
  mutate(
    Community = ordered(Community, levels = communities$Community),
    Isolate1 = ordered(Isolate1, levels = 1:12),
    Isolate2 = ordered(Isolate2, levels = 1:12)
  )


## Pairs melt by T0/T8 and three frequencies
pairs_melted <- pairs_freq %>%
  left_join(pairs_taxonomy, by = c("Community", "Isolate1", "Isolate2")) %>%
  left_join(pairs_16S, by = c("Community", "Isolate1", "Isolate2")) %>%
  #left_join(pairs_tree_distance, by = c("Community", "Isolate1", "Isolate2")) %>%
  as_tibble() %>%
  mutate(
    Community = ordered(Community, levels = communities$Community),
    Isolate1 = ordered(Isolate1, levels = 1:12),
    Isolate2 = ordered(Isolate2, levels = 1:12)
  )

#
fwrite(pairs, here::here("data/output/pairs.csv"))
fwrite(pairs_melted, here::here("data/output/pairs_melted.csv"))


