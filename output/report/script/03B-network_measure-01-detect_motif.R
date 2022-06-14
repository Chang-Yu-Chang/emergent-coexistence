#' Detect the motifs in an invasion network
#' 1. Detect motif counts in observed network
#' 2. Detect motif counts in randomized networks

library(tidyverse)
source(here::here("plotting_scripts/misc.R"))

communities <- read_csv("~/Dropbox/lab/emergent-coexistence/data/output/communities.csv", col_types = cols())
load("~/Dropbox/lab/emergent-coexistence/data/output/communities_network.Rdata") # net_list
load("~/Dropbox/lab/emergent-coexistence/data/output/communities_network_randomized.Rdata") # net_randomized_list

# Count the motifs in observed networks: 13 communities x 7 motifs = 91 rows
communities_motif <- communities_network %>%
    rowwise() %>%
    mutate(Motif = list(count_motif(Network))) %>%
    unnest_longer(col = Motif, indices_to = "Motif", values = "Count") %>%
    group_by(Community) %>%
    mutate(Fraction = Count / sum(Count)) %>%
    select(Community, Motif, Count, Fraction)


# Count the motifs in randomized networks: 13 communities x 7 motifs x 1000 replicate = 91000 rows
communities_motif_randomized <- communities_network_randomized %>%
    select(Community, NetworkRandomized) %>%
    unnest_longer(col = NetworkRandomized, indices_to = "Replicate", values_to = "Network") %>%
    rowwise() %>%
    mutate(Motif = list(count_motif(Network))) %>%
    unnest_longer(col = Motif, indices_to = "Motif", values = "Count") %>%
    group_by(Community, Replicate) %>%
    mutate(Fraction = Count / sum(Count)) %>%
    select(Community, Replicate, Motif, Count, Fraction)


write_csv(communities_motif, "~/Dropbox/lab/emergent-coexistence/data/output/communities_motif.csv")
write_csv(communities_motif_randomized, "~/Dropbox/lab/emergent-coexistence/data/output/communities_motif_randomized.csv")

