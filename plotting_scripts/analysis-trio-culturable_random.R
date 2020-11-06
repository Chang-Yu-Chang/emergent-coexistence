# Random culturable isolates

suppressWarnings(suppressMessages(library(tidyverse)))
suppressWarnings(suppressMessages(library(data.table)))
options(dplyr.summarise.inform = FALSE)
suppressWarnings(suppressMessages(library(cowplot)))
suppressWarnings(suppressMessages(library(tidygraph)))
suppressWarnings(suppressMessages(library(ggraph)))
suppressWarnings(suppressMessages(library(igraph)))

read_trio_list <- function (trio_culturable_list) {
    trio_culturable_list %>% 
        filter(Type == "consumer") %>% 
        select(Well, Trio, InitialFrequency, ID) %>% 
        group_by(Well) %>% 
        mutate(Isolate = 1:3) %>% 
        pivot_wider(names_from = Isolate, names_prefix = "Isolate", values_from = ID) %>% 
        ungroup()
}
read_trio_competition <- function(trio_culturable_data, trio_list) {
    trio_culturable_data %>% 
        filter(Type == "consumer") %>%
        left_join(select(trio_list, Well, Trio, InitialFrequency), by = "Well") %>% 
        group_by(Trio, InitialFrequency, Transfer) %>% 
        mutate(ID = factor(ID)) %>% 
        mutate(TotalAbundance = sum(Abundance), RelativeAbundance = Abundance/TotalAbundance) %>% 
        select(Trio, InitialFrequency, Transfer, ID, RelativeAbundance) %>% 
        ungroup()
}
determine_trio_outcome <- function(trio_outcome) {
    trio_outcome %>%
        filter(Transfer == max(Transfer)) %>% 
        filter(InitialFrequency == "33-33-33") %>% 
        group_by(Trio, InitialFrequency) %>%
        summarise(Richness = n())
}
plot_trio_temporal <- function(trio_competition) {
    trio_competition %>%
        filter(Trio %in% paste0("Trio", 1:5)) %>%
        ggplot() +
        geom_bar(aes(x = Transfer, y = RelativeAbundance, fill = ID), stat = "identity", color = 1) +
        facet_grid(InitialFrequency~Trio) +
        scale_x_continuous(expand = c(0,0)) +
        scale_y_continuous(expand = c(0,0)) +
        theme_bw() +
        theme(legend.position = "none")
}
read_pair_from_trio_list <- function(pair_from_trio_culturable_list) { 
    pair_from_trio_culturable_list %>% 
        filter(Type == "consumer") %>% 
        select(Well, Trio, Pair, InitialFrequency, ID) %>% 
        group_by(Well) %>% 
        mutate(Isolate = 1:2) %>% 
        pivot_wider(names_from = Isolate, names_prefix = "Isolate", values_from = ID) %>% 
        ungroup()
}
read_pair_from_trio_competition <- function(pair_from_trio_culturable_data, pair_from_trio_list) {
    pair_from_trio_culturable_data %>% 
        filter(Type == "consumer") %>%
        left_join(select(pair_from_trio_list, Well, Trio, Pair, InitialFrequency), by = "Well") %>% 
        group_by(Trio, Pair, InitialFrequency, Transfer) %>% 
        mutate(ID = factor(ID)) %>% 
        mutate(TotalAbundance = sum(Abundance), RelativeAbundance = Abundance/TotalAbundance) %>% 
        select(Trio, Pair, InitialFrequency, Transfer, ID, RelativeAbundance) %>% 
        ungroup()
}
#pairs_consumer <- df_pair_from_trio_competition
#pair_list<-df_pair_from_trio_list
determine_pair_from_trio_outcome <- function(pairs_consumer, pair_list) {
    n_trio_pair_freq <- nrow(pair_list)
    n_transfer <- max(pairs_consumer$Transfer)

    frequency_changes <- pairs_consumer %>%
        group_by(Trio, Pair, InitialFrequency) %>% 
        pivot_wider(names_from = c(Transfer), names_prefix = "T", values_from = RelativeAbundance) %>%
        replace_na(rep(list(0), n_transfer+1) %>% setNames(paste0("T", 0:n_transfer))) %>% 
        mutate(FrequencyChange = ifelse((T5 - T0)>0, "T", "F")) %>%  # TRUE = Isolate1 increases, FALSE = Isolate1 decreases
        left_join(pair_list, by = c("Trio", "Pair", "InitialFrequency")) %>%
        filter(ID == Isolate1) %>% 
        select(Trio, Pair, Isolate1, Isolate2, FrequencyChange) %>% 
        group_by(Trio, Pair, Isolate1, Isolate2) %>% 
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
        select(Trio, Pair, Isolate1, Isolate2, InteractionType, From, To) %>% 
        ungroup() %>% 
        return()
}
determine_trio_motif <- function(pairs_from_trio) {
    # Make network
    isolate_list <- sort(unique(c(pairs_from_trio$Isolate1, pairs_from_trio$Isolate2)))
    nodes <- tibble(Isolate = 1:length(isolate_list), ID = isolate_list)
    edges <- pairs_from_trio %>% 
        mutate(from = match(From, nodes$ID), to = match(To, nodes$ID)) %>% 
        select(InteractionType, from, to)
    edges_coext <- edges[edges$InteractionType == "coexistence",]
    edges_coext[,c("from", "to")] <- edges_coext[,c("to", "from")] # Add the mutual edges for coexistence links
    edges <- rbind(edges, edges_coext)
    
    # 
    if (any(is.na(pairs_from_trio$InteractionType))) {
        motif_counts <- tibble(Motif = 1:7, Count = rep(NA, 7))
    } else {
        graph <- tbl_graph(nodes = nodes, edges = edges, directed = T)
        motif_counts <- tibble(Motif = 1:7, Count = count_motif(graph))
    }
    
    return(motif_counts)
}
#source("network_functions.R")
# 
# 
# #args = commandArgs(trailingOnly = T)
# #input_independent <- fread(args[[1]])
# #input_independent <- fread("../data/raw/simulation/mapping_files/input_independent.csv")
# input_independent_trios <- input_independent %>% filter(grepl("trio-culturable_isolates", exp_id))
# #input_independent_trios <- input_independent_trios %>% filter(seed %in% 1:2)
# 
# temp_list <- rep(list(NA), nrow(input_independent_trios))
# 
# for (i in 1:nrow(input_independent_trios)) {
#     cat("\nexp_id = ", input_independent_trios$exp_id[i])
#     cat(",\tseed = ", input_independent_trios$seed[i])
#     
#     # Trio
#     df_trio_list <- fread(paste0("../data/raw/simulation/trio-culturable-", i, ".txt")) %>% 
#         read_trio_list()
#     df_trio_competition <- fread(paste0("../data/raw/simulation/trio-culturable_isolates-", i, "_composition.txt")) %>% 
#         read_trio_competition(df_trio_list)
#     df_trio_outcome <- determine_trio_outcome(df_trio_competition)
#     
#     # Pairs from trios
#     df_pair_from_trio_list <- fread(paste0("../data/raw/simulation/pair-culturable_from_trio-", i, ".txt")) %>% 
#         read_pair_from_trio_list()
#     df_pair_from_trio_competition <- fread(paste0("../data/raw/simulation/pair-culturable_from_trio-", i, "_composition.txt")) %>% 
#         read_pair_from_trio_competition(df_pair_from_trio_list)
#     df_pair_from_trio_outcome <- determine_pair_from_trio_outcome(df_pair_from_trio_competition, df_pair_from_trio_list)
#     df_trio_motif <- df_pair_from_trio_outcome %>% 
#         split.data.frame(f=.$Trio) %>% 
#         lapply(determine_trio_motif) %>% 
#         bind_rows(.id = "Trio") %>% 
#         filter(Count != 0) %>% 
#         select(Trio, Motif)
#     
#     temp_list[[i]] <- df_trio_motif
#     
# }
# 
# df_trio_motif_aggregate <- bind_rows(temp_list, .id = "Seed")
# df_trio_motif_aggregate <- mutate(Seed = 1, df_trio_motif)
# 
# #mapping_seed <- tibble(Seed = factor(1:4), l = c(0.2, 0, 0, 0.2), q = c(0.8, 0, 0.8, 0))
# 
# p1 <- df_trio_motif_aggregate %>% 
# #    left_join(mapping_seed) %>% 
#     #mutate(Motif = factor(Motif)) %>% 
#     left_join(df_trio_outcome) %>% 
#     ggplot() +
#     #geom_boxplot(aes(x = Motif, y = Richness, group = Motif)) +
#     geom_jitter(aes(x = Motif, y = Richness), shape = 21, width = 0.2, height = 0.2, size = 3) +
#     scale_x_continuous(limits = c(1,7), breaks = 1:7) +
#     scale_y_continuous(limits = c(0.5,3.5), breaks = 1:3) +
#     #facet_grid(Seed~., scales = "free_y") + 
# #    facet_wrap(l~q, labeller = label_both, nrow = 1) + 
#     theme_bw() +
#     theme(panel.grid.minor = element_blank())
#     #ggtitle(paste0(length(unique(df_trio_motif$Trio)), " trios"))
# p1
# 
# #cat("\nPlot the trio vs. richness")
# 
# ggsave("../plots/Fig2B.png", plot = p1, width = 10, height = 4)
# 






