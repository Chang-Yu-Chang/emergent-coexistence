#' This scripts reads the T0 vs. T8 frequencies and determine the competition outcomes
#' for each unique species pair by computing the fitness functions
#' 1. Determine the competition outcomes: match the significance of frequency changes to the fitness function table
#' 2. Combine T0 and T8 results into a pairs_freq table
#' 3. Calculate isolate tournament

library(tidyverse)
library(cowplot)
folder_script <- "~/Desktop/lab/emergent-coexistence/analysis/"
folder_main <- "~/Dropbox/lab/emergent-coexistence/plate_scan_pipeline/"

pairs_freq_ID <- read_csv(paste0(folder_main, "meta/00-pairs_freq_ID.csv"), show_col_types = F)
pairs_T0_boots <- read_csv(paste0(folder_main, "meta/93-pairs_T0_boots.csv"), show_col_types = F) # bootstraps using T0 mean and sd
pairs_T8_boots <- read_csv(paste0(folder_main, "meta/93-pairs_T8_boots.csv"), show_col_types = F) # bootstraps using random forest object probabilities
pairs_T8 <- read_csv(paste0(folder_main, "meta/93-pairs_T8.csv"), show_col_types = F) # random forest classification
pairs_ID <- distinct(pairs_freq_ID, Batch, Community, Isolate1, Isolate2)

pairs_no_colony <- c(
    "C11R1_2_8",
    "C11R1_2_9",
    "C11R1_8_9",
    "C11R2_2_10"
)
plates_no_colony <- c(
    "B2_T8_C11R1_5-95_2_8",
    "B2_T8_C11R1_5-95_2_9",
    "B2_T8_C11R1_5-95_8_2",
    "B2_T8_C11R1_5-95_9_8",
    "B2_T8_C11R1_50-50_2_8",
    "B2_T8_C11R1_50-50_2_9",
    "C2_T8_C11R2_50-50_2_10",
    "C2_T8_C11R2_50-50_9_13",
    "C_T8_C11R1_50-50_1_2", # no plate
    "C_T8_C11R1_50-50_1_3" # no plate
)


# 1. Determine pairwise competition ----
interaction_type_finer <- c("competitive exclusion", "stable coexistence",
                            "mutual exclusion", "frequency-dependent coexistence",
                            "coexistence at 5%", "coexistence at 95%",
                            "2-freq neutrality", "3-freq neutrality")

read_pairs_boots_table <- function (T8_freq_type = "bootstrapped") {
    #' This function batchly reads the table of boostrapped frequency change per pair
    temp <- rep(list(NA), nrow(pairs_ID))
    for (i in 1:nrow(pairs_ID)) {
        pair_name <- paste0(pairs_ID$Community[i], "_", pairs_ID$Isolate1[i], "_", pairs_ID$Isolate2[i])
        if (pair_name %in% pairs_no_colony) {cat("\nT8 has no colony, skip pair\t", pair_name); next}
        temp[[i]] <- read_csv(paste0(folder_main, "images/meta-93-T0_T8_frequencies/", T8_freq_type, "/", pair_name, ".csv"), show_col_types = F)
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

pairs_interaction <- read_pairs_boots_table("bootstrapped") %>%
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

#write_csv(pairs_outcome_classifiedT8, paste0(folder_main, "meta/95-pairs_outcome_classifiedT8.csv"))
write_csv(pairs_interaction, paste0(folder_main, "meta/95-pairs_interaction.csv"))

pairs_ID %>%
    left_join(pairs_interaction) %>%
    mutate(InteractionTypeFiner = factor(InteractionTypeFiner, interaction_type_finer)) %>%
    group_by(InteractionType, InteractionTypeFiner, FitnessFunction) %>%
    count(name = "Count") %>%
    arrange(InteractionTypeFiner)


# 2. Combine T0 and T8 frequencies ----
pairs_freq_T0_boots <- pairs_T0_boots %>%
    group_by(Community, Isolate1, Isolate2, Isolate1InitialODFreq, Time) %>%
    summarize(Isolate1CFUFreqMean = mean(Isolate1CFUFreq),
              Isolate1CFUFreqSd = sd(Isolate1CFUFreq))
pairs_freq_T8_boots <- pairs_T8_boots %>%
    group_by(Community, Isolate1, Isolate2, Isolate1InitialODFreq, Time) %>%
    summarize(Isolate1CFUFreqMean = mean(Isolate1CFUFreq),
              Isolate1CFUFreqSd = sd(Isolate1CFUFreq))
# pairs_freq_T8 <- pairs_T8 %>%
#     select(Community, Isolate1, Isolate2, Isolate1InitialODFreq, Time,
#            Isolate1CFUFreqMean = Isolate1CFUFreq) %>%
#     mutate(Isolate1CFUFreqSd = NA)


# T8 frequency bootstrapping so there is SD at T8
pairs_freq <- bind_rows(pairs_freq_T0_boots, pairs_freq_T8_boots)

# No bootstrapping, the random forest predicted classification, so there is no SD at T8
#pairs_freq <- bind_rows(pairs_freq_T0_boots, pairs_freq_T8)

write_csv(pairs_freq, paste0(folder_main, "meta/95-pairs_freq.csv"))
#write_csv(pairs_freq, paste0(folder_main, "meta/95-pairs_freq.csv"))


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
communities <- read_csv("~/Dropbox/lab/emergent-coexistence/data/output/communities.csv", show_col_types = F)
pairs_interaction <- read_csv(paste0(folder_main, "meta/95-pairs_interaction.csv"), show_col_types = F)

isolates_tournament <- communities %>%
    select(comm = Community, everything()) %>%
    rowwise() %>%
    mutate(pairs_comm = pairs_interaction %>% filter(Community == comm) %>% list()) %>%
    mutate(tournaments_comm = pairs_comm %>% tournament_rank() %>% list()) %>%
    select(Community = comm, tournaments_comm) %>%
    unnest(cols = tournaments_comm)
write_csv(isolates_tournament, paste0(folder_main, "meta/95-isolates_tournament.csv"))









