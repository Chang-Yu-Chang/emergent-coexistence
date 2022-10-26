#' This scripts reads the T0 vs. T8 frequencies and determine the competition outcomes
#' for each unique species pair by computing the fitness functions
#' 1. Combine T0 and T8 results into a pairs_freq table
#' 2. Determine the competition outcomes: match the significance of frequency changes to the fitness function table
#' 3. Calculate isolate tournament

library(tidyverse)
library(cowplot)
source(here::here("analysis/00-metadata.R"))

communities <- read_csv(paste0(folder_data, "temp/00c-communities.csv"), show_col_types = F)
pairs_ID <- read_csv(paste0(folder_data, "temp/00c-pairs_ID.csv"), show_col_types = F)
pairs_T0_boots <- read_csv(paste0(folder_data, "temp/06-pairs_T0_boots.csv"), show_col_types = F) # bootstraps using T0 mean and sd
pairs_T8_boots <- read_csv(paste0(folder_data, "temp/06-pairs_T8_boots.csv"), show_col_types = F) # bootstraps using random forest object probabilities=

# 1. Combine T0 and T8 frequencies ----
pairs_freq_T0_boots <- pairs_T0_boots %>%
    group_by(Community, Isolate1, Isolate2, Isolate1InitialODFreq, Time) %>%
    summarize(Isolate1CFUFreqMean = mean(Isolate1CFUFreq),
              Isolate1CFUFreqSd = sd(Isolate1CFUFreq)) %>%
    mutate(Community = factor(Community, communities$Community)) %>%
    arrange(Community, Isolate1, Isolate2, Isolate1InitialODFreq) %>%
    ungroup() %>%
    # PairFreqID
    mutate(PairFreqID = 1:n())
pairs_freq_T8_boots <- pairs_T8_boots %>%
    group_by(Community, Isolate1, Isolate2, Isolate1InitialODFreq, Time) %>%
    summarize(Isolate1CFUFreqMean = mean(Isolate1CFUFreq),
              Isolate1CFUFreqSd = sd(Isolate1CFUFreq)) %>%
    mutate(Community = factor(Community, communities$Community)) %>%
    arrange(Community, Isolate1, Isolate2, Isolate1InitialODFreq) %>%
    left_join(select(pairs_freq_T0_boots, PairFreqID, Community, Isolate1, Isolate2, Isolate1InitialODFreq))
pairs_freq <- bind_rows(pairs_freq_T0_boots, pairs_freq_T8_boots) %>%
    select(PairFreqID, everything())

write_csv(pairs_freq, paste0(folder_data, "temp/93-pairs_freq.csv"))

# 2. Determine pairwise competition ----
read_pairs_boots_table <- function () {
    #' This function batchly reads the table of boostrapped frequency change per pair
    temp <- rep(list(NA), nrow(pairs_ID))
    for (i in 1:nrow(pairs_ID)) {
        pair_name <- paste0(pairs_ID$Community[i], "_", pairs_ID$Isolate1[i], "_", pairs_ID$Isolate2[i])
        if (pair_name %in% pairs_no_colony) {cat("\nT8 has no colony, skip pair\t", pair_name); next}
        temp[[i]] <- read_csv(paste0(folder_pipeline, "images/", pairs_ID$Batch[i],"-09-bootstrap/", pair_name, ".csv"), show_col_types = F)
        cat("\n", i, "/", nrow(pairs_ID))
    }
    pairs_boots_table <- bind_rows(temp[which(!is.na(temp))])
}
compute_pairs_fitness <- function (pairs_boots_table) {
    #' This script takes the bootstrap samples and compute the significance of frequency changes
    pairs_boots_table %>%
        group_by(Community) %>%
        mutate(Community = factor(Community, paste0("C", rep(1:12, each = 8), "R", rep(1:8, 12)))) %>%
        arrange(Community, Isolate1, Isolate2, Isolate1InitialODFreq) %>%
        group_by(Community, Isolate1, Isolate2, Isolate1InitialODFreq, .drop = F) %>%
        mutate(Fraction = Count / sum(Count)) %>%
        # C11R1 isolate 6
        filter(!is.na(FreqChangeSign)) %>%
        pivot_wider(id_cols = c(Community, Isolate1, Isolate2, Isolate1InitialODFreq), names_from = FreqChangeSign, values_from = Fraction) %>%
        mutate(FitnessChange = case_when(
            increase > 0.95 ~ 1,
            decrease > 0.95 ~ -1,
            increase < 0.95 & decrease < 0.95 ~ 0
        ))
}
make_interaction_type <- function () {
    #' This function generates the fitness function table.
    #' There are a total of 27 possibilities
    interaction_type <- tibble(
        FromRare = rep(c(1, -1, 0), each = 9),
        FromMedium = rep(rep(c(1, -1, 0), each = 3), 3),
        FromAbundant = rep(c(1, -1, 0), 9),
        InteractionType = NA,
        InteractionTypeFiner = NA
    )
    ## Assign interaction types to combinations of frequency changes signs
    interaction_type$InteractionType[c(1,10,13,14)] <- "exclusion"
    interaction_type$InteractionType[c(2:6,8,11,20,23, 9,18,21,24,25,26,27)] <- "coexistence"

    ## Assign finer interaction types to combinations of frequency changes signs
    interaction_type$InteractionTypeFiner[c(1,14)] <- "competitive exclusion"
    interaction_type$InteractionTypeFiner[c(10,13)] <- "mutual exclusion"
    interaction_type$InteractionTypeFiner[c(2,5,8)] <- "stable coexistence"
    interaction_type$InteractionTypeFiner[c(4,6,11,20)] <- "frequency-dependent coexistence"
    interaction_type$InteractionTypeFiner[c(3)] <- "coexistence at 95%"
    interaction_type$InteractionTypeFiner[c(23)] <- "coexistence at 5%"
    interaction_type$InteractionTypeFiner[c(9,18,21,24:26)] <- "neutrality"
    interaction_type$InteractionTypeFiner[c(9,18,21,24:26)] <- "2-freq neutrality"
    interaction_type$InteractionTypeFiner[c(27)] <- "3-freq neutrality"
    interaction_type <- interaction_type %>%  mutate(FitnessFunction = paste(FromRare, FromMedium, FromAbundant, sep = "_"))
}
append_pairs_outocme <- function (pairs_fitness) {
    # This function appends the bootstrap result to the fitness function table
    pairs_fitness %>%
        pivot_wider(id_cols = c(Community, Isolate1, Isolate2), names_from = Isolate1InitialODFreq, values_from = FitnessChange) %>%
        rename(FromRare = `5`, FromMedium = `50`, FromAbundant = `95`) %>%
        left_join(interaction_type) %>%
        select(Community, Isolate1, Isolate2, InteractionType, InteractionTypeFiner, FitnessFunction)

}
interaction_type <- make_interaction_type()

pairs_interaction <- read_pairs_boots_table() %>%
    compute_pairs_fitness() %>%
    append_pairs_outocme() %>%
    mutate(From = case_when(
        FitnessFunction %in% c("1_1_1") ~ Isolate1, # Isolate1 wins
        FitnessFunction %in% c("-1_-1_-1") ~ Isolate2, # Isolate2 wins
        TRUE ~ Isolate1
    )) %>%
    mutate(To = case_when(
        FitnessFunction %in% c("1_1_1") ~ Isolate2,
        FitnessFunction %in% c("-1_-1_-1") ~ Isolate1,
        TRUE ~ Isolate2
    )) %>%
    ungroup()


# Flip the sign so that Isolate1 is the dominant and Isolate2 is subdominant
#' 50-50 frequency changes. To determine the dominance in coexistence pairs.
#' The isolate that increases frequency in the 50-50 competition is the dominant strain
pairs_coexist_dominant <- pairs_interaction %>%
    filter(InteractionType == "coexistence") %>%
    select(Community, Isolate1, Isolate2) %>%
    left_join(pairs_freq) %>%
    filter(Isolate1InitialODFreq == 50) %>%
    select(Community, Isolate1, Isolate2, Time, Isolate1CFUFreqMean) %>%
    pivot_wider(names_from = "Time", values_from = "Isolate1CFUFreqMean") %>%
    mutate(Isolate1Dominant = T8 > T0) %>%
    select(Community, Isolate1, Isolate2, Isolate1Dominant)

pairs_interaction <- left_join(pairs_interaction, pairs_coexist_dominant)

write_csv(pairs_interaction, paste0(folder_data, "temp/93-pairs_interaction.csv"))


# 3. isolate tournament ----
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

isolates_tournament <- communities %>%
    select(comm = Community, everything()) %>%
    rowwise() %>%
    mutate(pairs_comm = pairs_interaction %>% filter(Community == comm) %>% list()) %>%
    mutate(tournaments_comm = pairs_comm %>% tournament_rank() %>% list()) %>%
    select(Community = comm, tournaments_comm) %>%
    unnest(cols = tournaments_comm)
write_csv(isolates_tournament, paste0(folder_data, "temp/93-isolates_tournament.csv"))









