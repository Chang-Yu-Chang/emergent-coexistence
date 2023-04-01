#' This script read the isolate and pair data and generate network objects
#' 1. generate network objects
#' 2. randomize the networks
#' 3. community network hierarchy
library(tidyverse)
library(tidygraph)
source(here::here("analysis/00-metadata.R"))

communities <- read_csv(paste0(folder_data, "output/communities_remained.csv"), show_col_types = F)
isolates <- read_csv(paste0(folder_data, "output/isolates_remained.csv"), show_col_types = F)
pairs <- read_csv(paste0(folder_data, "output/pairs_remained.csv"), show_col_types = F)

make_network <- function(isolates, pairs) {
    # Nodes
    nodes <- isolates

    # Edges
    ## Remove no-growth
    edges <- pairs %>%
        #filter(outcome %in% c("coexistence", "exclusion", "inconclusive")) %>%
        filter(outcome %in% c("1-exclusion", "2-exclusion", "3-coexistence", "4-coexistence", "5-inconclusive")) %>%
        mutate(from = From, to = To) %>% select(from, to, outcome)
    edges_coext <- edges[edges$outcome %in% c("3-coexistence","4-coexistence"),]
    edges_coext[,c("from", "to")] <- edges_coext[,c("to", "from")] # Add the mutual edges for coexistence links
    edges <- rbind(edges, edges_coext)

    # Network
    graph <- tbl_graph(nodes = nodes, edges = edges, directed = T)

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
    select(Community, CommunityLabel, CommunitySize, CommunityPairSize, Network) %>%
    ungroup

save(communities_network, file = paste0(folder_data, "temp/95-communities_network.Rdata"))


# # 2. Make a R list of randomized  networks ----
# randomize_network <- function(graph){
#     # Step1: remove the bidirection of coexistence
#     graph1 <- graph %>%
#         activate(edges) %>%
#         reroute(from = to, to = from, subset = from > to) %>%
#         distinct()
#
#     # Step2: shuffle interaction types
#     graph2 <- graph1 %>%
#         activate(edges) %>%
#         mutate(outcome = outcome[order(runif(n=n()))])
#
#     # Step3: for pairs that coexist, add back the reverse direction links
#     pairs_coexistence_reversed <- graph2 %>%
#         activate(edges) %>%
#         filter(outcome == "coexistence") %>%
#         reroute(from = to, to = from) %>%
#         activate(nodes) %>%
#         select(Isolate)
#
#     graph_coexistence <- graph_join(graph2, pairs_coexistence_reversed, by = "Isolate") %>%
#         filter(outcome == "coexistence")
#
#     # Step4: for exclusionary pairs, 50% of those has a flip direction
#     n_exclusion_pairs <- graph2 %>% filter(outcome == "exclusion") %>% pull(outcome) %>% length()
#     graph_exclusion <- graph2 %>%
#         activate(edges) %>%
#         filter(outcome == "exclusion") %>%
#         reroute(from = to, to = from, subset = sample(c(T,F), n_exclusion_pairs, replace = T, prob = c(0.5,0.5))) %>%
#         activate(nodes) %>%
#         select(Isolate)
#
#     graph3 <- graph_join(graph_coexistence, graph_exclusion, by = "Isolate") %>%
#         arrange(from, to)
#
#     return(graph3)
# }
# b <- 1000 # Number of bootstrapping/randomization
#
# # Make an empty two-layer R list
# net_randomized_list <- rep(list(rep(list(NA), b)), nrow(communities_network))
# names(net_randomized_list) <- communities$Community
#
# time_start <- proc.time()
# for (i in 1:nrow(communities_network)) {
#     time_t <- proc.time()
#     for (b_loop_index in 1:b) {
#         set.seed(b_loop_index)
#         net_randomized_list[[i]][[b_loop_index]] <- randomize_network(communities_network$Network[[i]])
#         if (b_loop_index %% 100 == 0) cat("\nboostrap =", b_loop_index)
#     }
#     # Print
#     cat("\n\n", communities$Community[i])
#     cat("\n", (proc.time() - time_t)[3], "seconds")
#
# }
# cat("\n\n total time:", (proc.time() - time_start)[3], "seconds\n\n")
#
# communities_network_randomized <- communities_network %>%
#     ungroup() %>%
#     mutate(NetworkRandomized = net_randomized_list)
#
# # Save the data file
# save(communities_network_randomized, file = paste0(folder_data, "temp/94-communities_network_randomized.Rdata"))



# 3. Calculate network hierarchy ----
#' Violation of ranks ----
## Compute fraction of pairs that follow the ranks (computed by the number of wins)
tournament_rank <- function(pairs) {
    #' Compute the rank of an isolate based on its number of wins and lose in pairwise competition
    if (igraph::is.igraph(pairs)) {
        net <- pairs
        pairs <- as_tibble(get.edgelist(net)) %>%
            setNames(c("From", "To")) %>%
            mutate(outcome = get.edge.attribute(net)$outcome) %>%
            rowwise() %>%
            mutate(Isolate1 = min(From, To), Isolate2 = max(From, To))
    }
    isolate_name <- pairs %>% select(Isolate1, Isolate2) %>% unlist %>% unique %>% sort()
    # Isolates' ranks in the tournament
    tour_rank <- data.frame(
        Isolate = isolate_name,
        # Win
        Win = filter(pairs, outcome == "exclusion") %>%
            select(From) %>% unlist() %>% factor(isolate_name) %>% table() %>% as.vector(),
        # Lose
        Lose = filter(pairs, outcome == "exclusion") %>%
            select(To) %>% unlist() %>% factor(isolate_name) %>% table() %>% as.vector(),
        # Draw; Note that I consider neturality and bistability as draw in the tournament
        Draw = filter(pairs, outcome %in% c("coexistence", "neutrality", "bistability")) %>%
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
        select(Community, Isolate1, Isolate2, outcome, From, To)
    isolates_tournament <- tournament_rank(pairs_temp) %>% select(Isolate, Score)

    pairs_temp %>%
        left_join(rename_with(isolates_tournament, ~ paste0(., "1")), by = "Isolate1") %>%
        left_join(rename_with(isolates_tournament, ~ paste0(., "2")), by = "Isolate2") %>%
        filter(outcome == "exclusion") %>%
        mutate(WinnerScore = ifelse(From == Isolate1, Score1, Score2),
               LoserScore = ifelse(From == Isolate1, Score2, Score1)) %>%
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
    x$outcome <- x$outcome[rng]
    # Shuffle dominance
    rng <- sample(1:nrow(x), size = nrow(x)/2, replace = F)
    temp <- x$From[rng]
    x$From[rng] <- x$To[rng]
    x$To[rng] <- temp

    return(x)
}

communities_hierarchy <- communities %>%
    select(comm = Community) %>%
    rowwise() %>%
    mutate(pairs_comm = pairs %>% filter(Community == comm) %>% list()) %>%
    mutate(HierarchyScore = compute_hierarchy2(pairs_comm)) %>%
    select(-pairs_comm) %>%
    select(Community = comm, HierarchyScore)

write_csv(communities_hierarchy, paste0(folder_data, "temp/95-communities_hierarchy.csv"))






