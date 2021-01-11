# Pairs from top-down communities

suppressWarnings(suppressMessages(library(tidyverse)))
suppressWarnings(suppressMessages(library(data.table)))
options(dplyr.summarise.inform = FALSE)
suppressWarnings(suppressMessages(library(cowplot)))
suppressWarnings(suppressMessages(library(tidygraph)))
suppressWarnings(suppressMessages(library(ggraph)))
suppressWarnings(suppressMessages(library(igraph)))
source("network_functions.R")

#
input_independent <- fread("../data/raw/simulation/mapping_files/input_independent.csv")
input_independent_community <- input_independent %>% filter(grepl("community-top_down", exp_id))

#
read_community_list <- function (community_composition) {
    community_composition %>% 
        filter(Transfer == max(Transfer), Type == "consumer") %>% 
        mutate(Community = paste0("Community", sub("W", "", Well))) %>% 
        select(Community, Type, ID, Abundance)
}
read_pair_from_commmunity_list <- function(pair_from_community_list) { 
    pair_from_community_list %>% 
        filter(Type == "consumer") %>% 
        select(Well, Community, Pair, InitialFrequency, ID) %>% 
        group_by(Well) %>% 
        arrange(Well) %>% 
        mutate(Isolate = 1:2) %>% 
        pivot_wider(names_from = Isolate, names_prefix = "Isolate", values_from = ID) %>% 
        ungroup()
}
read_pair_from_community_competition <- function(pair_from_community_data, pair_from_community_list) {
    pair_from_community_data %>% 
        filter(Type == "consumer") %>%
        left_join(select(pair_from_community_list, Well, Community, Pair, InitialFrequency), by = "Well") %>% 
        group_by(Community, Pair, InitialFrequency, Transfer) %>% 
        mutate(ID = factor(ID)) %>% 
        mutate(TotalAbundance = sum(Abundance), RelativeAbundance = Abundance/TotalAbundance) %>% 
        select(Community, Pair, InitialFrequency, Transfer, ID, RelativeAbundance) %>% 
        ungroup()
}
determine_pair_from_community_outcome <- function(pairs_consumer, pair_list) {
    frequency_changes <- pairs_consumer %>%
        left_join(pair_list, by = c("Community", "Pair", "InitialFrequency")) %>%
        filter(ID == Isolate1) %>% select(-ID) %>% 
        group_by(Community, Pair) %>% 
        pivot_wider(names_from = c(Transfer), names_prefix = "T", values_from = RelativeAbundance) %>% 
        mutate(FrequencyChange = ifelse((T5 - T0)>0, "T", "F")) %>%  # TRUE = Isolate1 increases, FALSE = Isolate1 decreases
        select(Community, Pair, Isolate1, Isolate2, FrequencyChange) %>% 
        group_by(Community, Pair, Isolate1, Isolate2) %>% 
        summarise(FrequencyChangePattern = paste0(FrequencyChange, collapse = "-"))
    
    map_frequency_pattern <- tibble(FrequencyChangePattern = c("T-T-T", "F-F-F", "T-F-F", "T-T-F", "T-F-T", "F-T-F", "F-F-T", "F-T-T"),
        InteractionType = c("exclusion", "exclusion", "coexistence", "coexistence", "coexistence", "coexistence", "mutual exclusion", "mutual exclusion"))
    
    frequency_changes$From = NA
    frequency_changes$To = NA
    frequency_changes$InteractionType = NA
    
    for (i in 1:nrow(frequency_changes)) {
        if (frequency_changes$FrequencyChangePattern[i] == "T-T-T") {
            from <- frequency_changes$Isolate1[i]
            to <- frequency_changes$Isolate2[i]
            interaction_type <- "exclusion"
        } else if (frequency_changes$FrequencyChangePattern[i] == "F-F-F") {
            from <- frequency_changes$Isolate2[i]
            to <- frequency_changes$Isolate1[i]
            interaction_type <- "exclusion"
        } else if (frequency_changes$FrequencyChangePattern[i] %in% c("T-F-F", "T-T-F", "T-F-T", "F-T-F")) {
            from <- frequency_changes$Isolate1[i]
            to <- frequency_changes$Isolate2[i]
            interaction_type <- "coexistence"
        } else {
            from <- NA
            to <- NA
            interaction_type <- NA
        }
        
        frequency_changes$From[i] <- from
        frequency_changes$To[i] <- to
        frequency_changes$InteractionType[i] <- interaction_type
    }
    
    frequency_changes %>%
        select(Community, Pair, Isolate1, Isolate2, InteractionType, From, To) %>% 
        ungroup() %>% 
        return()
}
determine_community_motif <- function(pairs_from_community) {
    # Make network
    isolate_list <- sort(unique(c(pairs_from_community$Isolate1, pairs_from_community$Isolate2)))
    nodes <- tibble(Isolate = 1:length(isolate_list), ID = isolate_list)
    edges <- pairs_from_community %>% 
        mutate(from = match(From, nodes$ID), to = match(To, nodes$ID)) %>% 
        select(InteractionType, from, to)
    edges_coext <- edges[edges$InteractionType == "coexistence",]
    edges_coext[,c("from", "to")] <- edges_coext[,c("to", "from")] # Add the mutual edges for coexistence links
    edges <- rbind(edges, edges_coext)
    
    # 
    #if (any(is.na(pairs_from_community$InteractionType))) {
    #    motif_counts <- tibble(Motif = 1:7, Count = rep(NA, 7))
    #} else {
        graph <- tbl_graph(nodes = nodes, edges = filter(edges, !is.na(InteractionType)), directed = T)
        motif_counts <- tibble(Motif = 1:7, Count = count_motif(graph))
    #}
    
    return(motif_counts)
}


