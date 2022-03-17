# Clean up the AcrAss and RanAss community data
library(tidyverse)
library(tidygraph)
library(igraph)
library(ggraph)
library(cowplot)
source(here::here("plotting_scripts/network_functions.R"))

communities <- read_csv(here::here("data/output/communities.csv"), col_types = cols())
isolates <- read_csv(here::here("data/output/isolates.csv"), col_types = cols())
pairs <- read_csv(here::here("data/output/pairs.csv"), col_types = cols()) %>% mutate(InteractionType = ifelse(InteractionType == "neutrality", "coexistence", InteractionType))
pairs_freq <- read_csv(here::here("data/output/pairs_freq.csv"), col_types = cols())
#load("~/Dropbox/lab/invasion-network/data/output/network_community.Rdata")
community_factor <- c(communities %>% filter(str_detect(Community, "C\\d")) %>% arrange(CommunitySize) %>% pull(Community),
                      communities %>% filter(str_detect(Community, "Ass")) %>% pull(Community))
communities_size <- communities %>% mutate(Community = factor(Community, community_factor)) %>% arrange(Community) %>% pull(CommunitySize)

b = 10

# Hierarchy ----
compute_hierarchy1 <- function(pairs_mock) {
    pairs_temp <- pairs_mock %>%
        select(Isolate1, Isolate2, InteractionType, From, To)
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
randomize_pairs1 <- function(x) {
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

### Permutation
communities_randomized_list1 <- rep(list(NA), nrow(communities))
tt <- proc.time()
for (i in 1:length(communities_randomized_list1)) {
    temp_tt <- proc.time()
    pairs_temp <- pairs %>% filter(Community == communities$Community[i])
    isolates_temp <- tibble(Community = unique(pairs_temp$Community), Isolate = sort(unique(c(pairs_temp$Isolate1, pairs_temp$Isolate2))))

    communities_randomized_list1[[i]] <- tibble(Community = communities$Community[i], Replicate = 1:b) %>%
        mutate(pairs_comm = list(pairs_temp)) %>%
        mutate(pairs_randomized = map(pairs_comm, randomize_pairs1)) %>%
        rowwise() %>%
        mutate(h1 = compute_hierarchy1(pairs_randomized))
    cat("\n", communities$Community[i] %>% as.character())
    cat("\n", (proc.time() - temp_tt)[3], "seconds")
}

communities_hierarchy_randomized1 <- bind_rows(communities_randomized_list1) %>%
    select(Community, Replicate, h1) %>%
    mutate(Community = factor(Community))
### obv
communities_hierarchy1 <- pairs %>%
    nest(pairs_comm = -Community) %>%
    rowwise() %>%
    mutate(h1 = compute_hierarchy1(pairs_comm)) %>%
    select(Community, h1)

## Higgins 2017
compute_comp_score2 <- function(pairs_comm, target_isolate) {
    pairs_comm %>%
        filter(Isolate1 == target_isolate | Isolate2 == target_isolate) %>%
        mutate(Freq = ifelse(Isolate1 == target_isolate, Isolate1MeasuredFreq, 1-Isolate1MeasuredFreq)) %>%
        pull(Freq) %>%
        mean(na.rm = T)
}
compute_hierarchy2 <- function(isolates_mock, pairs_mock) {
    ranking <- isolates_mock %>%
        rowwise() %>%
        mutate(CompetitiveScore = compute_comp_score2(pairs_mock, Isolate)) %>%
        arrange(desc(CompetitiveScore)) %>%
        # Ranking by competitive score
        pull(Isolate)

    pairs_mock %>%
        mutate(Freq1 = Isolate1MeasuredFreq, Freq2 = 1-Freq1) %>%
        select(Isolate1, Isolate2, Freq1, Freq2) %>%
        mutate(Pair = 1:n()) %>%
        pivot_longer(cols = c(-Pair), names_to = ".value", names_pattern = "(.+)[12]") %>%
        mutate(Isolate = ordered(Isolate, ranking)) %>%
        group_by(Pair) %>%
        arrange(Pair, Isolate) %>%
        slice(1) %>%
        pull(Freq) %>%
        mean(na.rm = T) %>%
        return()

}
randomize_pairs2 <- function(x) {
    # Shuffle pairs
    rng <- order(runif(nrow(x),0,1))
    x$Isolate1MeasuredFreq <- x$Isolate1MeasuredFreq[rng]
    # Shuffle dominance
    rng <- sample(1:nrow(x), size = nrow(x)/2, replace = F)
    x$Isolate1MeasuredFreq[rng] <- 1 - x$Isolate1MeasuredFreq[rng]

    return(x)
}

### Permutation
communities_randomized_list2 <- rep(list(NA), nrow(communities))
tt <- proc.time()
for (i in 1:length(communities_randomized_list2)) {
    temp_tt <- proc.time()
    pairs_temp <- pairs_freq %>%
        filter(Community == communities$Community[i]) %>%
        select(Community, Isolate1, Isolate2, Isolate1MeasuredFreq)
    isolates_temp <- tibble(Community = unique(pairs_temp$Community), Isolate = sort(unique(c(pairs_temp$Isolate1, pairs_temp$Isolate2))))

    communities_randomized_list2[[i]] <- tibble(Community = communities$Community[i], Replicate = 1:b) %>%
        mutate(pairs_comm = list(pairs_temp), isolates_comm = list(isolates_temp)) %>%
        mutate(pairs_randomized = map(pairs_comm, randomize_pairs2)) %>%
        rowwise() %>%
        mutate(h2 = compute_hierarchy2(isolates_comm, pairs_randomized))
    cat("\n", communities$Community[i] %>% as.character())
    cat("\n", (proc.time() - temp_tt)[3], "seconds")
}

communities_hierarchy_randomized2 <- bind_rows(communities_randomized_list2) %>%
    select(Community, Replicate, h2) %>%
    mutate(Community = factor(Community))

### obv
communities_hierarchy2 <- pairs_freq %>%
    rename(comm = Community) %>%
    nest(pairs_comm = -comm) %>%
    rowwise() %>%
    mutate(isolates_comm = list(filter(isolates, Community == comm))) %>%
    mutate(h2 = compute_hierarchy2(isolates_comm, pairs_comm)) %>%
    select(Community = comm, h2)


# Join data
communities_hierarchy <- left_join(communities_hierarchy1, communities_hierarchy2) %>%
    mutate(Community = factor(Community, community_factor)) %>%
    pivot_longer(-Community, names_to = "Metric", values_to = "HierarchyScore")
communities_hierarchy_randomized <- left_join(communities_hierarchy_randomized1, communities_hierarchy_randomized2) %>%
    pivot_longer(-c(Community, Replicate), names_to = "Metric", values_to = "HierarchyScore")
write_csv(communities_hierarchy, here::here("data/output/communities_hierarchy.csv"))
write_csv(communities_hierarchy_randomized, here::here("data/output/communities_hierarchy_randomized.csv"))




# Make networks ----
## obv
temp <- communities %>%
    rename(comm = Community) %>%
    rowwise() %>%
    mutate(isolate_comm = filter(isolates, Community == comm) %>% list()) %>%
    mutate(pairs_comm = filter(pairs, Community == comm) %>% list()) %>%
    mutate(Network = make_network(isolate_comm, pairs_comm) %>% list())
net_list <- temp$Network %>% set_names(temp$comm)

## permutation
net_randomized_list <- rep(list(rep(list(NA), b)), length(net_list))
names(net_randomized_list) <- communities$Community

tt <- proc.time()
for (i in 1:length(net_list)) {
    temp_tt <- proc.time()
    for (b_loop_index in 1:b) {
        set.seed(b_loop_index)
        net_randomized_list[[i]][[b_loop_index]] <- randomize_network(net_list[[i]])
        if (b_loop_index %% 100 == 0) cat("\n boostrap =", b_loop_index)
    }
    # Print
    cat("\n\n", communities$Community[i])
    cat("\n", (proc.time() - temp_tt)[3], "seconds")
    if (i == length(net_list)) cat("\n\n total time:", (proc.time() - tt)[3], "seconds\n\n")
}

save(net_list, file = "~/Dropbox/lab/invasion-network/data/output/network.Rdata")
save(net_randomized_list, file = "~/Dropbox/lab/invasion-network/data/output/network_randomized.Rdata")



# Motif ----
load("~/Dropbox/lab/invasion-network/data/output/network.Rdata")
load("~/Dropbox/lab/invasion-network/data/output/network_randomized.Rdata")

## Obv
networks_motif <- communities %>%
    mutate(Network = net_list) %>%
    rowwise() %>%
    mutate(temp = list(tibble(Motif = 1:7, Count = count_motif(Network))), .keep = "unused") %>%
    unnest(cols = c(temp)) %>%
    group_by(Community) %>%
    mutate(Fraction = Count / sum(Count))

## Permutation
networks_motif_randomized <- communities %>%
    mutate(NetworkList = net_randomized_list) %>%
    rowwise() %>%
    mutate(temp = list(lapply(NetworkList, function(x) {
        tibble(Motif = 1:7, Count = count_motif(x))
    }) %>% bind_rows(.id = "Replicate")), .keep = "unused") %>%
    unnest(cols = c(temp)) %>%
    group_by(Replicate, Community) %>%
    mutate(Fraction = Count / sum(Count))

## Find percentiles
networks_motif_percentile <- networks_motif_randomized %>%
    group_by(Community, Motif) %>%
    arrange(desc(Count)) %>%
    slice(ceiling(b * 0.05), floor(b * 0.95)) %>%
    mutate(Percentile = c("p5", "p95")) %>%
    select(Community, Motif, Percentile, Count) %>%
    pivot_wider(names_from = Percentile, values_from = Count)

networks_motif <- networks_motif %>%
    left_join(networks_motif_percentile) %>%
    mutate(Sign = case_when(Count > p95 ~ "top",
                            Count < p5 ~ "bottom",
                            Count < p95 & Count > p5 ~ "n.s."))


write_csv(networks_motif, here::here("data/output/networks_motif.csv"))
write_csv(networks_motif_randomized, here::here("data/output/networks_motif_randomized.csv"))


# Diagonal analysis ----
## Observation
networks_diag <- communities %>%
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































