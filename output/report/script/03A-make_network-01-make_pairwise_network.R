#' Make pairwise network from outcome of pairwise competitions
library(tidyverse)
library(tidygraph)
library(ggraph)
source(here::here("plotting_scripts/misc.R"))

communities <- read_csv("~/Dropbox/lab/emergent-coexistence/data/output/communities.csv", col_types = cols())
isolates <- read_csv("~/Dropbox/lab/emergent-coexistence/data/output/isolates.csv", col_types = cols())
pairs <- read_csv("~/Dropbox/lab/emergent-coexistence/data/output/pairs.csv", col_types = cols()) %>%
    mutate(InteractionType = ifelse(InteractionType == "neutrality", "coexistence", InteractionType))

communities_network <- communities %>%
    rename(comm = Community) %>%
    rowwise() %>%
    mutate(Isolates = filter(isolates, Community == comm) %>% list) %>%
    mutate(Pairs = filter(pairs, Community == comm) %>% list) %>%
    mutate(Network = make_network(Isolates, Pairs) %>% list) %>%
    mutate(NetworkPlot = plot_competitive_network(Network) %>% list) %>%
    rename(Community = comm)

save(communities_network, file = "~/Dropbox/lab/emergent-coexistence/data/output/communities_network.Rdata")
