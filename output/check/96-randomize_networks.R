#' This script read the isolate and pair data and generate network objects
#' 1. generate network objects
#' 2. randomize the networks
#' 3. community network hierarchy
library(tidyverse)
library(tidygraph)
library(ggraph)
#source(here::here("plotting_scripts/misc.R"))

folder_script <- "~/Desktop/lab/emergent-coexistence/output/check/"
folder_main <- "~/Dropbox/lab/emergent-coexistence/data/raw/plate_scan/emergent_coexistence_plate_scan_check/"

communities <- read_csv("~/Dropbox/lab/emergent-coexistence/data/output/communities.csv", col_types = cols())
pairs_freq_ID <- read_csv(paste0(folder_main, "meta/00-pairs_freq_ID.csv"), show_col_types = F)
pairs_ID <- distinct(pairs_freq_ID, Batch, Community, Isolate1, Isolate2) %>% mutate(PairID = 1:n())
pairs_interaction <- read_csv(paste0(folder_main, "meta/95-pairs_interaction.csv"), show_col_types = F)
isolates_tournament <- read_csv(paste0(folder_main, "meta/95-isolates_tournament.csv"), show_col_types = F)
pairs_freq <- read_csv(paste0(folder_main, "meta/95-pairs_freq.csv"), show_col_types = F)


#paint_white_background <- function(x) theme(plot.background = element_rect(color = NA, fill = "white"))
make_network <- function(isolates, pairs) {
    # Nodes
    nodes <- isolates %>% select(Isolate, Rank, PlotRank)

    # Edges
    ## Remove no-growth
    edges <- pairs %>%
        filter(InteractionType %in% c("coexistence", "exclusion")) %>%
        mutate(from=From, to=To) %>% select(from, to, InteractionType)
    edges_coext <- edges[edges$InteractionType == "coexistence",]
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
    mutate(Isolates = filter(isolates_tournament, Community == comm) %>% list) %>%
    mutate(Pairs = filter(pairs_interaction, Community == comm) %>% list) %>%
    mutate(Network = make_network(Isolates, Pairs) %>% list) %>%
    rename(Community = comm) %>%
    select(Community, CommunityLabel, CommunitySize, CommunityPairSize, Network)

save(communities_network, file = paste0(folder_main, "meta/96-communities_network.Rdata"))


# 2. Make a R list of randomized  networks ----
randomize_network <- function(graph){
    # Step1: remove the bidirection of coexistence
    graph1 <- graph %>%
        activate(edges) %>%
        reroute(from = to, to = from, subset = from > to) %>%
        distinct()

    # Step2: shuffle interaction types
    graph2 <- graph1 %>%
        activate(edges) %>%
        mutate(InteractionType = InteractionType[order(runif(n=n()))])

    # Step3: for pairs that coexist, add back the reverse direction links
    pairs_coexistence_reversed <- graph2 %>%
        activate(edges) %>%
        filter(InteractionType == "coexistence") %>%
        reroute(from = to, to = from) %>%
        activate(nodes) %>%
        select(Isolate)

    graph_coexistence <- graph_join(graph2, pairs_coexistence_reversed, by = "Isolate") %>%
        filter(InteractionType == "coexistence")

    # Step4: for exclusionary pairs, 50% of those has a flip direction
    n_exclusion_pairs <- graph2 %>% filter(InteractionType == "exclusion") %>% pull(InteractionType) %>% length()
    graph_exclusion <- graph2 %>%
        activate(edges) %>%
        filter(InteractionType == "exclusion") %>%
        reroute(from = to, to = from, subset = sample(c(T,F), n_exclusion_pairs, replace = T, prob = c(0.5,0.5))) %>%
        activate(nodes) %>%
        select(Isolate)

    graph3 <- graph_join(graph_coexistence, graph_exclusion, by = "Isolate") %>%
        arrange(from, to)

    return(graph3)
}
b <- 1000 # Number of bootstrapping/randomization

# Make an empty two-layer R list
net_randomized_list <- rep(list(rep(list(NA), b)), nrow(communities_network))
names(net_randomized_list) <- communities$Community

time_start <- proc.time()
for (i in 1:nrow(communities_network)) {
    time_t <- proc.time()
    for (b_loop_index in 1:b) {
        set.seed(b_loop_index)
        net_randomized_list[[i]][[b_loop_index]] <- randomize_network(communities_network$Network[[i]])
        if (b_loop_index %% 100 == 0) cat("\nboostrap =", b_loop_index)
    }
    # Print
    cat("\n\n", communities$Community[i])
    cat("\n", (proc.time() - time_t)[3], "seconds")

}
cat("\n\n total time:", (proc.time() - time_start)[3], "seconds\n\n")

communities_network_randomized <- communities_network %>%
    ungroup() %>%
    mutate(NetworkRandomized = net_randomized_list)

# Save the data file
save(communities_network_randomized, file = paste0(folder_main, "meta/96-communities_network_randomized.Rdata"))








# 3. Calculate network hierarchy ----

#' Violation of ranks ----
## Compute fraction of pairs that follow the ranks (computed by the number of wins)
tournament_rank <- function(pairs) {
    #' Compute the rank of an isolate based on its number of wins and lose in pairwise competition
    if (igraph::is.igraph(pairs)) {
        net <- pairs
        pairs <- as_tibble(get.edgelist(net)) %>%
            setNames(c("From", "To")) %>%
            mutate(InteractionType = get.edge.attribute(net)$InteractionType) %>%
            rowwise() %>%
            mutate(Isolate1 = min(From, To), Isolate2 = max(From, To))
    }
    isolate_name <- pairs %>% select(Isolate1, Isolate2) %>% unlist %>% unique %>% sort()
    # Isolates' ranks in the tournament
    tour_rank <- data.frame(
        Isolate = isolate_name,
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
        select(Community, Isolate1, Isolate2, InteractionType, From, To)
    isolates_tournament <- tournament_rank(pairs_temp) %>% select(Isolate, Score)

    pairs_temp %>%
        left_join(rename_with(isolates_tournament, ~ paste0(., "1")), by = "Isolate1") %>%
        left_join(rename_with(isolates_tournament, ~ paste0(., "2")), by = "Isolate2") %>%
        filter(InteractionType == "exclusion") %>%
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
    x$InteractionType <- x$InteractionType[rng]
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
    mutate(pairs_comm = pairs_interaction %>% filter(Community == comm) %>% list()) %>%
    mutate(HierarchyScore = compute_hierarchy2(pairs_comm)) %>%
    select(-pairs_comm) %>%
    select(Community = comm, HierarchyScore)

write_csv(communities_hierarchy, paste0(folder_main, "meta/96-communities_hierarchy.csv"))






