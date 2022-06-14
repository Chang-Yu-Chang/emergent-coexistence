#' Competitive hierarchy
#' Three measures
#' 1. Higgins et al 2017, based on relative abundances in pairs
#' 2. The fraction of pairs that follow the ranks

library(tidyverse)
source(here::here("plotting_scripts/misc.R"))
isolates_tournament <- read_csv("~/Dropbox/lab/emergent-coexistence/data/temp/isolates_tournament.csv", col_types = cols())
pairs_interaction <- read_csv("~/Dropbox/lab/emergent-coexistence/data/output/pairs_interaction.csv", col_types = cols()) %>%
    mutate(InteractionType = ifelse(InteractionType == "neutrality", "coexistence", InteractionType)) %>%
    filter(Set == "CFUandCASEU")
pairs_freq <- read_csv("~/Dropbox/lab/emergent-coexistence/data/output/pairs_freq.csv", col_types = cols())
communities <- read_csv("~/Dropbox/lab/emergent-coexistence/data/output/communities.csv", col_types = cols())
load("~/Dropbox/lab/emergent-coexistence/data/output/communities_network.RData")
#load("~/Dropbox/lab/emergent-coexistence/data/output/communities_network_randomized.RData")

b = 1000

# 1. Higgins et al 2017 ----
# Functions for computing hierarchy scores
## Compute competitive score for each isolates in a community
compute_comp_score1 <- function(pairs_comm, target_isolate) {
    pairs_comm %>%
        #select(-Isolate1InitialODFreq, -Isolate2InitialODFreq, -RawDataType, -Contamination) %>%
        filter(Isolate1 == target_isolate | Isolate2 == target_isolate) %>%
        mutate(Freq = ifelse(Isolate1 == target_isolate, Isolate1MeasuredFreq, 1-Isolate1MeasuredFreq)) %>%
        pull(Freq) %>%
        mean(na.rm = T)
}

## Compute hierarchy score for a community
compute_hierarchy1 <- function(isolates_mock, pairs_mock) {
    ranking <- isolates_mock %>%
        rowwise() %>%
        # Subset the pairs for each community
        mutate(CompetitiveScore = compute_comp_score1(pairs_mock, Isolate)) %>%
        select(Community, Isolate, CompetitiveScore) %>%
        arrange(Community, desc(CompetitiveScore)) %>%
        # Ranking by competitive score
        pull(Isolate)

    pairs_mock %>%
        mutate(Freq1 = Isolate1MeasuredFreq, Freq2 = 1-Freq1) %>%
        select(Community, Isolate1, Isolate2, Freq1, Freq2) %>%
        mutate(Pair = 1:n()) %>%
        pivot_longer(cols = c(-Community, -Pair), names_to = ".value", names_pattern = "(.+)[12]") %>%
        mutate(Isolate = ordered(Isolate, ranking)) %>%
        group_by(Community, Pair) %>%
        arrange(Community, Pair, Isolate) %>%
        slice(1) %>%
        pull(Freq) %>%
        mean() %>%
        return()

}
## Data: isolate and pairs
communities_hierarchy_obv1 <- communities %>%
    select(comm = Community) %>%
    rowwise() %>%
    mutate(isolates_comm = isolates_tournament %>% filter(Community == comm) %>% list()) %>%
    mutate(pairs_comm = pairs_freq %>% filter(Isolate1InitialODFreq == 50, Time == "Tend", Set == "CFUandCASEU") %>% filter(Community == comm) %>% list()) %>%
    mutate(HierarchyScore = compute_hierarchy1(isolates_comm, pairs_comm)) %>%
    select(-isolates_comm, -pairs_comm) %>%
    select(Community = comm, HierarchyScore)

if (FALSE) {

## Randomized pairs with shuffled relative abundances within a community
randomize_pairs1 <- function(x) {
    # Shuffle pairs
    rng <- order(runif(nrow(x),0,1))
    x$Isolate1MeasuredFreq <- x$Isolate1MeasuredFreq[rng]
    # Shuffle dominance
    rng <- sample(1:nrow(x), size = nrow(x)/2, replace = F)
    x$Isolate1MeasuredFreq[rng] <- 1 - x$Isolate1MeasuredFreq[rng]

    return(x)
}

# Randomized networks
communities_randomized_list <- rep(list(NA), nrow(communities))
tt <- proc.time()
for (i in 1:length(communities_randomized_list)) {
    temp_tt <- proc.time()
    pairs_temp <- pairs_freq %>%
        filter(Set == "CFUandCASEU") %>%
        filter(Isolate1InitialODFreq == 50, Time == "Tend") %>%
        filter(Community == communities$Community[i]) %>%
        select(Community, Isolate1, Isolate2, Isolate1MeasuredFreq)

    pairs_temp <- randomize_pairs1(pairs_temp)
    isolates_temp <- communities_network_randomized$Isolates[i][[1]] %>% select(Community, Isolate)

    communities_randomized_list[[i]] <- tibble(Community = communities$Community[i], Replicate = 1:b) %>%
        mutate(pairs_comm = list(pairs_temp), isolates_comm = list(isolates_temp)) %>%
        mutate(pairs_randomized = map(pairs_comm, randomize_pairs1)) %>%
        rowwise() %>%
        mutate(HierarchyScore = compute_hierarchy1(isolates_comm, pairs_randomized))
    cat("\n", communities$Community[i] %>% as.character())
    cat("\n", (proc.time() - temp_tt)[3], "seconds")
}

communities_hierarchy_randomized1 <- bind_rows(communities_randomized_list) %>%
    select(Community, Replicate, HierarchyScore) %>%
    mutate(Community = factor(Community))

## Compute p value
communities_hierarchy_pvalue <- communities_hierarchy_randomized1 %>%
    group_by(Community) %>%
    # Add a bottom row to each community
    bind_rows(communities %>% select(Community) %>% mutate(Replicate = b+1, HierarchyScore = 0)) %>%
    arrange(Community, desc(HierarchyScore)) %>%
    left_join(select(communities_hierarchy_obv, Community, HierarchyScoreObserved = HierarchyScore)) %>%
    mutate(Percentile = 0:b / b) %>%
    filter(HierarchyScoreObserved > HierarchyScore) %>%
    slice(1) %>%
    mutate(Significance = Percentile < 0.05) %>%
    select(Community, Percentile, Significance)
}


# 2. Violation of ranks ----
## Compute fraction of pairs that follow the ranks (computed by the number of wins)
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

## Data: isolate and pairs
communities_hierarchy_obv2 <- communities %>%
    select(comm = Community) %>%
    rowwise() %>%
    mutate(pairs_comm = pairs_interaction %>% filter(Community == comm) %>% list()) %>%
    mutate(HierarchyScore = compute_hierarchy2(pairs_comm)) %>%
    select(-pairs_comm) %>%
    select(Community = comm, HierarchyScore)

if (FALSE) {
# Randomized networks
communities_randomized_list <- rep(list(NA), nrow(communities))
tt <- proc.time()
for (i in 1:length(communities_randomized_list)) {
    temp_tt <- proc.time()
    pairs_temp <- pairs_interaction %>%
        filter(Community == communities$Community[i]) %>%
        randomize_pairs2()
    isolates_temp <- communities_network_randomized$Isolates[i][[1]] %>% select(Community, Isolate)
    communities_randomized_list[[i]] <- tibble(Community = communities$Community[i], Replicate = 1:b) %>%
        mutate(pairs_comm = list(pairs_temp)) %>%
        mutate(pairs_randomized = map(pairs_comm, randomize_pairs2)) %>%
        rowwise() %>%
        mutate(HierarchyScore = compute_hierarchy2(pairs_randomized))
    cat("\n", communities$Community[i] %>% as.character())
    cat("\n", (proc.time() - temp_tt)[3], "seconds")
}

communities_hierarchy_randomized2 <- bind_rows(communities_randomized_list) %>%
    select(Community, Replicate, HierarchyScore) %>%
    mutate(Community = factor(Community))

## Compute p value
communities_hierarchy_pvalue <- communities_hierarchy_randomized2 %>%
    group_by(Community) %>%
    # Add a bottom row to each community
    bind_rows(communities %>% select(Community) %>% mutate(Replicate = b+1, HierarchyScore = 0)) %>%
    arrange(Community, desc(HierarchyScore)) %>%
    left_join(select(communities_hierarchy_obv, Community, HierarchyScoreObserved = HierarchyScore)) %>%
    mutate(Percentile = 0:b / b) %>%
    filter(HierarchyScoreObserved > HierarchyScore) %>%
    slice(1) %>%
    mutate(Significance = Percentile < 0.05) %>%
    select(Community, Percentile, Significance)
communities_hierarchy2 <- communities_hierarchy_obv %>%
    left_join(communities_hierarchy_pvalue)

}




# Join data
communities_hierarchy <- communities_hierarchy_obv1 %>%
    rename_with(~ paste0(., "1"), !contains("Community")) %>%
    left_join(rename_with(communities_hierarchy_obv2, ~ paste0(., "2"), !contains("Community")))
if (FALSE) {
communities_hierarchy_randomized <- communities_hierarchy_randomized1 %>%
    rename_with(~ paste0(., "1"), !contains("Community") & !contains("Replicate")) %>%
    left_join(rename_with(communities_hierarchy_randomized2, ~ paste0(., "2"), !contains("Community") & !contains("Replicate")))

}


write_csv(communities_hierarchy, "~/Dropbox/lab/emergent-coexistence/data/output/communities_hierarchy.csv")
#write_csv(communities_hierarchy_randomized, "~/Dropbox/lab/emergent-coexistence/data/output/communities_hierarchy_randomized.csv")








