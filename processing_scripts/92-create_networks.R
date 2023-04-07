#' This script read the isolate and pair data and generate network objects
#' 1. generate network objects
library(tidyverse)
library(tidygraph)
source(here::here("processing_scripts/00-metadata.R"))

communities <- read_csv(paste0(folder_data, "output/communities_remained.csv"), show_col_types = F)
isolates <- read_csv(paste0(folder_data, "output/isolates_remained.csv"), show_col_types = F)
pairs <- read_csv(paste0(folder_data, "output/pairs_remained.csv"), show_col_types = F)

make_network <- function(isolates, pairs) {
    # Nodes
    nodes <- isolates %>%
        # To avoid the tbl_graph error, add dummy isolates
        bind_rows(tibble(Isolate = c(1:12)[-isolates$Isolate])) %>%
        arrange(Isolate)

    # Edges
    edges <- pairs %>%
        filter(outcome %in% c("1-exclusion", "2-exclusion", "3-coexistence", "4-coexistence", "5-inconclusive")) %>%
        mutate(from = From, to = To) %>% select(from, to, outcome)
    edges_coext <- edges[edges$outcome %in% c("3-coexistence","4-coexistence"),]
    edges_coext[,c("from", "to")] <- edges_coext[,c("to", "from")] # Add the mutual edges for coexistence links
    edges <- rbind(edges, edges_coext)

    # Network
    graph <- tbl_graph(nodes = nodes, edges = edges, directed = T) %>%
        # Drop the dummy isolates
        activate(nodes) %>%
        filter(!is.na(Community))

    return(graph)
}

# 1. Make network object from isolate and pair ----
communities_network <- communities %>%
    rename(comm = Community) %>%
    rowwise() %>%
    mutate(Isolates = isolates %>% filter(Community == comm) %>% list) %>%
    mutate(Pairs = pairs %>% filter(Community == comm) %>% list) %>%
    mutate(Network = make_network(Isolates, Pairs) %>% list) %>%
    rename(Community = comm) %>%
    select(Community, CommunityLabel, CommunitySize, CommunityPairSize, Network)

save(communities_network, file = paste0(folder_data, "output/communities_network.Rdata"))
