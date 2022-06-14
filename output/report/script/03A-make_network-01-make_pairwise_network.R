#' Make pairwise network from outcome of pairwise competitions
library(tidyverse)
library(tidygraph)
library(ggraph)
source(here::here("plotting_scripts/misc.R"))

communities <- read_csv("~/Dropbox/lab/emergent-coexistence/data/output/communities.csv", col_types = cols())
isolates_ID_match <- read_csv("~/Dropbox/lab/emergent-coexistence/data/temp/isolates_ID_match.csv", col_types = cols())
isolates_tournament <- read_csv("~/Dropbox/lab/emergent-coexistence/data/temp/isolates_tournament.csv", col_types = cols())
isolates_temp <- left_join(isolates_ID_match, isolates_tournament)
pairs_interaction <- read_csv("~/Dropbox/lab/emergent-coexistence/data/output/pairs_interaction.csv", col_types = cols()) %>%
    mutate(InteractionType = ifelse(InteractionType == "neutrality", "coexistence", InteractionType)) %>%
    filter(Set == "CFUandCASEU")

communities_network <- communities %>%
    rename(comm = Community) %>%
    rowwise() %>%
    mutate(Isolates = filter(isolates_temp, Community == comm) %>% list) %>%
    mutate(Pairs = filter(pairs_interaction, Community == comm) %>% list) %>%
    mutate(Network = make_network(Isolates, Pairs) %>% list) %>%
    rename(Community = comm) %>%
    select(Community, CommunityLabel, CommunitySize, CommunityPairSize, Network)

save(communities_network, file = "~/Dropbox/lab/emergent-coexistence/data/output/communities_network.Rdata")
