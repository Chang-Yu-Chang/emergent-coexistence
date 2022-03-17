#' Count the number of pairwise coexistence as a function of distance to diagonal (|i-j|)
library(tidyverse)
library(data.table)
library(tidygraph)
library(igraph)
source(here::here("plotting_scripts/network_functions.R"))
isolates <- read_csv(here::here("data/output/isolates.csv"))
pairs <- read_csv(here::here("data/output/pairs.csv"))
pairs_meta <- read_csv(here::here("data/output/pairs_meta.csv"))
communities <- read_csv(here::here("data/output/communities.csv"))
load("~/Dropbox/lab/invasion-network/data/output/network_community.Rdata")
load("~/Dropbox/lab/invasion-network/data/output/network_randomized.Rdata")
#load(here::here("data/output/network_community.Rdata"))
#load(here::here("data/output/network_randomized.Rdata"))
# 20211215 remove random assembly network as it's not completed yet
net_list <- net_list[names(net_list) %>% str_detect("C\\d") %>% which]
net_randomized_list <- net_randomized_list[names(net_randomized_list) %>% str_detect("C\\d") %>% which]


# Observation
networks_diag <- communities %>%
    filter(str_detect(Community, "C\\d")) %>%
    mutate(Network = net_list) %>%
    rowwise() %>%
    mutate(Diagonal = count_diag_coexistence(Network) %>% mutate(Community) %>% list()) %>%
    pull(Diagonal) %>%
    bind_rows()

# Permutation
networks_diag_randomized_list <- rep(list(NA), length(net_list))
names(networks_diag_randomized_list) <- names(net_list)
tt <- proc.time()
for (i in 1:length(net_list)) {
    temp_tt <- proc.time()
    networks_diag_randomized_list[[i]] <- lapply(net_randomized_list[[i]], count_diag_coexistence) %>%
        bind_rows(.id = "Replicate")
    # Print
    cat("\n\n", communities$Community[i])
    cat("\n", (proc.time() - temp_tt)[3], "seconds")
    if (i == length(net_list)) cat("\n\n total time:", (proc.time() - tt)[3], "seconds")
}
networks_diag_randomized <- bind_rows(networks_diag_randomized_list, .id = "Community")

#
write_csv(networks_diag, file = here::here("data/output/networks_diag.csv"))
write_csv(networks_diag_randomized, file = here::here("data/output/networks_diag_randomized.csv"))

if (FALSE) {
net <- net_randomized_list$C11R2[[12]]
net_rank <- as_tibble(vertex_attr(net)) %>% arrange(PlotRank) %>% pull(Rank)
n_nodes <- length(net_rank)
net %>%
    activate(nodes) %>%
    select(Isolate, PlotRank) %>%
    activate(edges) %>%
    mutate(fromRank = .N()$PlotRank[match(from, .N()$Isolate)], toRank = .N()$PlotRank[match(to, .N()$Isolate)]) %>%
    #filter(fromRank <= toRank) %>%
    bind_edges(tibble(from = 1:n_nodes, to = 1:n_nodes, fromRank = 1:n_nodes, toRank = 1:n_nodes, InteractionType = "self")) %>%
    mutate(fromRank = factor(fromRank, n_nodes:1), toRank = factor(toRank)) %>%
    as_tibble() %>%
    ggplot() +
    geom_tile(aes(x = toRank, y = fromRank, fill = InteractionType), width = 0.9, height = 0.9) +
    #scale_x_discrete(position = "top", expand = c(0,0), labels = paste0("rank ", net_rank)) +
    #scale_y_discrete(position = "right", expand = c(0,0), labels = paste0("rank ", rev(net_rank))) +
    scale_fill_manual(values = c(assign_interaction_color("matrix"), "self" = "black"), breaks = c("exclusion", "coexistence")) +
    theme_classic() +
    labs()


net <- net_randomized_list$C11R2[[2]]
count_diag_coexistence(net)
pairs_comm <- as_edgelist(net) %>%
    magrittr::set_colnames(c("Isolate1", "Isolate2")) %>%
    as_tibble() %>%
    mutate(InteractionType = edge.attributes(net)$InteractionType) %>%
    filter(!(InteractionType == "coexistence" & Isolate1 > Isolate2))


get.adjacency(net, attr = "InteractionType", sparse = F)
pairs_comm %>%
    filter(Isolate1 == 1, Isolate2 == 4)

isolates_comm <- vertex.attributes(net) %>% as_tibble() %>% select(Isolate, Rank)
isolates_comm %>% arrange(Rank)

pairs_comm %>%
    filter(Isolate1 == 4 | Isolate2 == 4)

pairs_comm %>%
    filter(Isolate1 == 9 | Isolate2 == 9)

#
pairs_comm %>%
    # Join the isolate information
    left_join(rename_with(isolates_comm, ~paste0(., "1"), everything()), by = "Isolate1") %>%
    left_join(rename_with(isolates_comm, ~paste0(., "2"), everything()), by = "Isolate2") %>%
    # Calculate rank difference
    mutate(RankDifference = abs(Rank1 - Rank2)) %>%
    #arrange(desc(RankDifference))
    # Calculate the fraction of coexistence at each distance
    group_by(RankDifference, InteractionType) %>%
    summarize(Count = n(), .groups = "keep") %>%
    group_by(RankDifference) %>%
    mutate(TotalCount = sum(Count)) %>%
    mutate(Fraction = Count/TotalCount) %>%
    filter(InteractionType == "coexistence")



networks_diag_randomized %>%
    filter(Community == "C11R2") %>%
    filter(Replicate == 12)
filter(RankDifference == 11)


networks_diag_randomized %>%
    filter(Community == "C11R2") %>%
    group_by(RankDifference) %>%
    summarize(Count = n())

}
if (FALSE) {

    networks_diag %>%
        group_by(RankDifference) %>%
        summarize(Count = sum(Count), TotalCount = sum(TotalCount)) %>%
        mutate(Fraction = Count / TotalCount) %>%
        ggplot() +
        geom_point(aes(x = RankDifference, y = Fraction)) +
        scale_fill_manual(values = interaction_color) +
        scale_x_continuous(breaks = 0:12) +
        theme_classic()
    networks_diag_randomized_list[[1]] %>%
        ggplot() +
        geom_jitter(aes(x = RankDifference, y = Fraction), width = .1, height = .01, shape = 21) +
        scale_fill_manual(values = interaction_color) +
        scale_x_continuous(breaks = 0:12) +
        theme_classic()


    pairs_meta %>%
        filter(Assembly == "self_assembly") %>%
        select(Community, Isolate1, Isolate2, Rank1, Rank2, InteractionType) %>%
        mutate(RankDifference = abs(Rank1 - Rank2)) %>%
        group_by(RankDifference, InteractionType) %>%
        summarize(Count = n()) %>%
        group_by(RankDifference) %>%
        mutate(TotalCount = sum(Count)) %>%
        mutate(Fraction = Count/TotalCount) %>%
        filter(InteractionType == "coexistence") %>%
        ggplot() +
        #geom_col(aes(x = RankDifference, y = Fraction, fill = InteractionType)) +
        geom_point(aes(x = RankDifference, y = Fraction)) +
        scale_fill_manual(values = interaction_color) +
        scale_x_continuous(breaks = 0:12) +
        theme_classic()


    net <- net_list$C11R2
    net <- net_randomized_list$C1R2[[1]]
    pairs_comm <- as_edgelist(net) %>%
        magrittr::set_colnames(c("Isolate1", "Isolate2")) %>%
        as_tibble() %>%
        mutate(InteractionType = edge.attributes(net)$InteractionType)
    isolates_comm <- vertex.attributes(net) %>% as_tibble()


    pairs_comm %>%
        # Join the isolate information
        left_join(rename_with(isolates_comm, ~paste0(., "1"), everything()), by = "Isolate1") %>%
        left_join(rename_with(isolates_comm, ~paste0(., "2"), everything()), by = "Isolate2") %>%
        # Calculate rank difference
        mutate(RankDifference = abs(Rank1 - Rank2)) %>%
        # Calculate the fraction of coexistence at each distance
        group_by(RankDifference, InteractionType) %>%
        summarize(Count = n(), .groups = "keep") %>%
        group_by(RankDifference) %>%
        mutate(TotalCount = sum(Count)) %>%
        mutate(Fraction = Count/TotalCount) %>%
        filter(InteractionType == "coexistence")
}
if (FALSE) {
    # Calculating number of coexistence as the function of distance to diagonal
    adj_from_net <- function(net) {
        # Get adjacent matrix
        net_m <- get.adjacency(net, attr = "InteractionType", sparse = F)

        # Exclusion: win or lose
        temp_index <- which(net_m=="exclusion", arr.ind = T) %>% as.data.frame()
        for(i in 1:nrow(temp_index)) net_m[temp_index$col[i], temp_index$row[i]] <- "lose"
        net_m[net_m=="exclusion"] <- "win"

        # Bistability
        temp_index2 <- which(net_m=="bistability", arr.ind = T) %>% as.data.frame()
        for(i in 1:nrow(temp_index2)) net_m[temp_index2$col[i], temp_index2$row[i]] <- "bistability"
        net_m[net_m=="bistability"] <- "bistability"

        # Neutrality
        temp_index2 <- which(net_m=="neutrality", arr.ind = T) %>% as.data.frame()
        for(i in 1:nrow(temp_index2)) net_m[temp_index2$col[i], temp_index2$row[i]] <- "neutrality"
        net_m[net_m=="neutrality"] <- "neutrality"

        # Diagonal
        diag(net_m) <- "self"

        # Undefined
        temp_index <- which(net_m=="undefined", arr.ind = T) %>% as.data.frame()
        for(i in 1:nrow(temp_index)) net_m[temp_index$col[i], temp_index$row[i]] <- "undefined"

        # NA
        net_m[net_m == "" | is.na(net_m)] <- NA

        return(net_m)
    }
    # Find the fraction of coexistence as a function of distance to diagonal
    diag_distance <- function (net, observation = F) {
        # Convert network to matrix
        temp_matrix <- adj_from_net(net)

        # Re-order the matrix axis by the isolates' competitive rank
        temp_rank <- tournament_rank(net)$Isolate
        temp_matrix <- temp_matrix[temp_rank, temp_rank]

        # Order the matrix axis by competitive score
        t1 <- which(temp_matrix == "coexistence", arr.ind = T) %>%
            as_tibble() %>%
            filter(col > row) %>%
            mutate(DistanceToDiagonal = abs(row - col)) %>%
            group_by(DistanceToDiagonal) %>%
            summarize(CountCoexistence = n())

        # Total count for each distance
        t2 <- which(temp_matrix == "win" | temp_matrix == "lose" | temp_matrix == "coexistence", arr.ind = T) %>%
            as_tibble() %>%
            filter(col > row) %>%
            mutate(DistanceToDiagonal = abs(row - col)) %>%
            group_by(DistanceToDiagonal) %>%
            summarize(TotalCount = n())

        #
        if (observation) {
            full_join(t1, t2, by = "DistanceToDiagonal") %>%
                replace_na(list(CountCoexistence = 0)) %>%
                return()
        } else {
            return(t1)
        }
    }

    # Count the distance to diagonal in observed networks
    networks_diag <- lapply(net_list, function(x) diag_distance(x, observation = T)) %>% bind_rows(.id = "Community")

    # Count the distance to diagonal in randomized networks
    networks_diag_randomized_list <- rep(list(NA), length(net_list))
    names(networks_diag_randomized_list) <- names(net_list)

    tt <- proc.time()
    for (i in 1:length(net_list)) {
        temp_tt <- proc.time()
        networks_diag_randomized_list[[i]] <- lapply(net_randomized_list[[i]], diag_distance) %>%
            bind_rows(.id = "Replicate")
        # Print
        cat("\n\n", communities$Community[i])
        cat("\n", (proc.time() - temp_tt)[3], "seconds")
        if (i == length(net_list)) cat("\n\n total time:", (proc.time() - tt)[3], "seconds")
    }

    networks_diag_randomized <- bind_rows(networks_diag_randomized_list, .id = "Community")

}





















