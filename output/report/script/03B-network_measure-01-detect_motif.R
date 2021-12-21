#' Detect the motifs in an invasion network
#' 1. Detect motif counts in observed network
#' 2. Detect motif counts in randomized networks
#' 3. Compute 5% and 95% percentiles of motif counts in randomized networks
library(tidyverse)
source(here::here("plotting_scripts/network_functions.R"))

# Data
isolates <- read_csv(here::here("data/output/isolates.csv"))
pairs <- read_csv(here::here("data/output/pairs.csv"))
communities <- read_csv(here::here("data/output/communities.csv"))
load(here::here("data/output/network_community.Rdata")) # net_list
load(here::here("data/output/network_randomized.Rdata")) # net_randomized_list

# Count motifs
## Count the motifs in observed networks
networks_motif <- tibble(Community = names(net_list), Graph = net_list) %>%
    rowwise() %>%
    mutate(Motif = list(count_motif(Graph))) %>%
    unnest_longer(col = Motif, indices_to = "Motif", values = "Count") %>%
    group_by(Community) %>%
    mutate(Fraction = Count / sum(Count)) %>%
    select(Community, Motif, Count, Fraction)


# Count the motifs in randomized networks
networks_motif_randomized <- tibble(Community = names(net_list), Graph = net_randomized_list) %>%
    unnest_longer(col = Graph, indices_to = "Replicate") %>%
    rowwise() %>%
    mutate(Motif = list(count_motif(Graph)))

networks_motif_randomized <- networks_motif_randomized %>%
    unnest_longer(col = Motif, indices_to = "Motif", values = "Count") %>%
    group_by(Community, Replicate) %>%
    mutate(Fraction = Count / sum(Count)) %>%
    select(Community, Replicate, Motif, Count, Fraction)

# Find 5th and 95th percentiles for each motif within a community
b = 1000
networks_motif_randomized_percentile <-
  networks_motif_randomized %>%
  group_by(Community, Motif) %>%
  arrange(Count) %>%
  select(-Replicate) %>%
  slice(c(b*0.05, b*0.95)) %>%
  mutate(Percentile = c("p5", "p95")) %>%
  {.}

write_csv(networks_motif, file = here::here("data/output/networks_motif.csv"))
write_csv(networks_motif_randomized, file = here::here("data/output/networks_motif_randomized.csv"))
write_csv(networks_motif_randomized_percentile, here::here("data/output/networks_motif_randomized_percentile.csv"))

