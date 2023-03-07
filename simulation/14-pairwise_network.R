#' This script analyze calculate the network hierarchy for both
#' 1. monocutlure sets + pool pairs
#' 2. communities + within community pairs
#'
#' by making the network objects

library(tidyverse)
library(tidygraph)
source(here::here("analysis/00-metadata.R"))

# 0. parameters ----
input_parameters <- read_csv(here::here("simulation/01-input_parameters.csv"), col_types = cols())
input_poolPairs <- read_csv(here::here("simulation/03a-input_poolPairs.csv"), col_types = cols())
input_withinCommunityPairs <- read_csv(here::here("simulation/03b-input_withinCommunityPairs.csv"), col_types = cols())

# Generate family-species and class-resource table for matching
sa <- input_parameters$sa[1]
ma <- input_parameters$ma[1]
sal <- tibble(Family = paste0("F", c(rep(0, sa), rep(1, sa))), Species = paste0("S", 0:(sa * 2 - 1)))
mal <- tibble(Class = paste0("T", c(rep(0, ma), rep(1, ma))), Resource = paste0("R", 0:(ma * 2 - 1)))

# 0.1 Functions ----
make_network <- function(isolates, pairs) {
    # Nodes
    if (nrow(isolates) == 1) {
        nodes <- isolates %>% select(Species)
    } else {
        nodes <- isolates %>% select(Species, Rank, PlotRank)
    }

    # Edges
    ## Remove no-growth
    edges <- pairs %>%
        filter(InteractionType %in% c("coexistence", "exclusion", "unknown", "no colony or low accuracy")) %>%
        mutate(from=From, to=To) %>% select(from, to, InteractionType)
    edges_coext <- edges[edges$InteractionType == "coexistence",]
    edges_coext[,c("from", "to")] <- edges_coext[,c("to", "from")] # Add the mutual edges for coexistence links
    edges <- rbind(edges, edges_coext)

    # Network
    graph <- tbl_graph(nodes = nodes, edges = edges, directed = T)

    return(graph)
}
#' Violation of ranks ----
## Compute fraction of pairs that follow the ranks (computed by the number of wins)
tournament_rank <- function(pairs) {

    # Compute the rank of an isolate based on its number of wins and lose in pairwise competition
    if (igraph::is.igraph(pairs)) {
        net <- pairs
        pairs <- as_tibble(get.edgelist(net)) %>%
            setNames(c("From", "To")) %>%
            mutate(InteractionType = get.edge.attribute(net)$InteractionType) %>%
            rowwise() %>%
            mutate(Species1 = min(From, To), Species2 = max(From, To))
    }
    isolate_name <- pairs %>% select(Species1, Species2) %>% unlist %>% unique %>% sort()
    # Isolates' ranks in the tournament
    tour_rank <- data.frame(
        Species = isolate_name,
        # Win
        Win = filter(pairs, InteractionType == "exclusion") %>%
            select(From) %>% unlist() %>% factor(isolate_name) %>% table() %>% as.vector(),
        # Lose
        Lose = filter(pairs, InteractionType == "exclusion") %>%
            select(To) %>% unlist() %>% factor(isolate_name) %>% table() %>% as.vector(),
        # Draw; Note that I consider neturality and bistability as draw in the tournament
        Draw = filter(pairs, InteractionType %in% c("coexistence", "neutrality", "bistability")) %>%
            select(From, To) %>% unlist() %>% factor(isolate_name) %>% table() %>% as.vector())

    # Arrange the df by score
    tour_rank <- tour_rank %>%
        mutate(Score = Win - Lose + 0 * Draw, Game = Win + Lose + Draw) %>%
        arrange(desc(Score))

    # Calculate rank by score; same scores means the same ranks
    temp_score <- ordered(tour_rank$Score, levels = sort(unique(tour_rank$Score), decreasing = T))
    temp_score_table <- table(temp_score)
    temp <- NULL; temp_counter = 1
    for (i in 1:length(temp_score_table)) {
        temp <- c(temp, rep(temp_counter, temp_score_table[i]))
        temp_counter <- temp_counter + temp_score_table[i]
    }

    tour_rank$Rank <- temp
    tour_rank$PlotRank <- 1:nrow(tour_rank)
    return(tour_rank)
}
compute_hierarchy2 <- function(pairs_mock) {
    pairs_temp <- pairs_mock %>%
        select(Community, Species1, Species2, InteractionType, From, To)
    isolates_tournament <- tournament_rank(pairs_temp) %>% select(Species, Score)

    pairs_temp %>%
        left_join(rename_with(isolates_tournament, ~ paste0(., "1")), by = "Species1") %>%
        left_join(rename_with(isolates_tournament, ~ paste0(., "2")), by = "Species2") %>%
        filter(InteractionType == "exclusion") %>%
        mutate(WinnerScore = ifelse(From == Species1, Score1, Score2),
               LoserScore = ifelse(From == Species1, Score2, Score1)) %>%
        mutate(FollowRank = (WinnerScore > LoserScore) %>% factor(c(T,F))) %>%
        count(FollowRank, .drop = F, name = "Count") %>%
        mutate(FractionFollowRank = Count / sum(Count)) %>%
        filter(FollowRank == T) %>%
        pull(FractionFollowRank) %>%
        return()
}
randomize_pairs2 <- function(x) {
    # Shuffle pairs
    rng <- order(runif(nrow(x),0,1))
    x$InteractionType <- x$InteractionType[rng]
    # Shuffle dominance
    rng <- sample(1:nrow(x), size = nrow(x)/2, replace = F)
    temp <- x$From[rng]
    x$From[rng] <- x$To[rng]
    x$To[rng] <- temp

    return(x)
}


# 1. Pool pairs ----
poolPairs_N_outcome <- read_csv(paste0(folder_simulation, "12-aggregated_pairs/poolPairs_N_outcome.csv"), col_types = cols()) %>%
    mutate(Community = factor(Community, paste0("W", 0:19))) %>%
    select(Community, Species1, Species2, InteractionType, From, To)
monocultureSets_species <- read_csv(paste0(folder_simulation, "11-aggregated/monocultureSets_species.csv"), col_types = cols())
monocultureSets_richness <- read_csv(paste0(folder_simulation, "11-aggregated/monocultureSets_richness.csv"), col_types = cols())

poolPairs_N_outcome %>%
    filter(Community == "W0") %>%
    tournament_rank

monocultureSets_network <- monocultureSets_richness %>%
    rename(comm = Community) %>%
    rowwise() %>%
    mutate(Species = monocultureSets_species %>% filter(Community == comm) %>% list) %>%
    mutate(Pairs = poolPairs_N_outcome %>% filter(Community == comm) %>% list) %>%
    # Append the competitive ranks
    mutate(Species = list(left_join(Species, tournament_rank(Pairs), by = join_by(Species)))) %>%
    # Make network object from species and pair
    mutate(Network = make_network(Species, Pairs) %>% list) %>%
    rename(Community = comm) %>%
    # Calculate network hierarchy
    mutate(HierarchyScore = compute_hierarchy2(Pairs) %>% round(4)) %>%
    select(Community, Richness, HierarchyScore, Species, Pairs, Network) %>%
    ungroup()

save(monocultureSets_network, file = paste0(folder_simulation, "11-aggregated/monocultureSets_network.Rdata"))

monocultureSets_network %>%
    select(Community, HierarchyScore) %>%
    write_csv(paste0(folder_simulation, "11-aggregated/monocultureSets_hierarchy.csv"))


# 2. Within community pairs ----
communities_species <- read_csv(paste0(folder_simulation, "11-aggregated/communities_species.csv"), col_types = cols())
withinCommunityPairs_N_outcome <- read_csv(paste0(folder_simulation, "12-aggregated_pairs/withinCommunityPairs_N_outcome.csv"), col_types = cols()) %>%
    mutate(Community = factor(Community, paste0("W", 0:19))) %>%
    select(Community, Species1, Species2, InteractionType, From, To)
communities_species <- read_csv(paste0(folder_simulation, "11-aggregated/communities_species.csv"), col_types = cols())
communities_richness <- read_csv(paste0(folder_simulation, "11-aggregated/communities_richness.csv"), col_types = cols())

communities_network <- communities_richness %>%
    rename(comm = Community) %>%
    rowwise() %>%
    mutate(Species = communities_species %>% filter(Community == comm) %>% list) %>%
    mutate(Pairs = withinCommunityPairs_N_outcome %>% filter(Community == comm) %>% list) %>%
    # Append the competitive ranks
    mutate(Species = ifelse(nrow(Species) == 1, list(Species), list(left_join(Species, tournament_rank(Pairs), by = join_by(Species))))) %>%
    # Make network object from species and pair
    mutate(Network = make_network(Species, Pairs) %>% list) %>%
    rename(Community = comm) %>%
    # Calculate network hierarchy
    mutate(HierarchyScore = ifelse(nrow(Pairs) == 0, NA, compute_hierarchy2(Pairs) %>% round(4))) %>%
    select(Community, Richness, HierarchyScore, Species, Pairs, Network) %>%
    ungroup()

save(communities_network, file = paste0(folder_simulation, "11-aggregated/communities_network.Rdata"))

communities_network %>%
    select(Community, HierarchyScore) %>%
    write_csv(paste0(folder_simulation, "11-aggregated/communities_hierarchy.csv"))


