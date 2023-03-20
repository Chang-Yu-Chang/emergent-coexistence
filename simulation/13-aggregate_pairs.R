#' This script aggregates the pairs data from 03a-poolPairs and 03b-withinCommunityPairs

library(tidyverse)
library(cowplot)
source(here::here("analysis/00-metadata.R"))
source(here::here("simulation/01-generate_input.R"))

# 0. parameters ----
input_parameters <- read_csv(here::here("simulation/01-input_parameters.csv"), col_types = cols())
input_poolPairs <- read_csv(here::here("simulation/03a-input_poolPairs.csv"), col_types = cols())
input_withinCommunityPairs <- read_csv(here::here("simulation/03b-input_withinCommunityPairs.csv"), col_types = cols())

#temp_end_time = paste0("T", input_poolPairs$n_timepoints[1])
#temp_end_time = input_poolPairs$n_timesteps[1]
temp_end_time <- input_communities$n_timepoints[1]

# Functions
## For reading and formating data
read_wide_file <- function (x, type = "N") {
    temp <- read_csv(x, col_types = cols(), name_repair = "unique_quiet") %>%
        pivot_longer(cols = starts_with("W"), names_to = "Well", values_to = "Abundance")
    if ("...1" %in% colnames(temp)) {
        if (type == "N") temp <- temp %>% rename(Family = ...1, Species = ...2)
        if (type == "R") temp <- temp %>% rename(Class = ...1, Resource = ...2)
    }

    return(temp)
}
format_columns <- function (x, type = "N") {
    # This function formats the Well, Family, Species, Class, Resource columns according to the designated order
    if (type == "N") {
    temp <- x %>%
        mutate(Well = factor(Well, paste0("W", 0:999)), .keep = "unused") %>%
        mutate(Family = factor(Family, unique(sal$Family))) %>%
        mutate(Species = factor(Species, sal$Species)) %>%
        mutate(Time = factor(Time, c("init", paste0("T", 1:20), "end")))
    }

    if (type == "R") {
        temp <- x %>%
            mutate(Well = factor(Well, paste0("W", 0:999)), .keep = "unused") %>%
            mutate(Class = factor(Class, unique(mal$Class))) %>%
            mutate(Resource = factor(Resource, mal$Resource)) %>%
            mutate(Time = factor(Time, c("init", paste0("T", 1:temp_end_time), "end")))
    }

    return(temp)
}
add_pairID <- function (x) {
    x %>%
        arrange(Well, Species) %>%
        mutate(Pair = paste0("P", rep(1:(n()/6), each = 6)))
}
read_init_composition <- function (input_mapping, treatment, comm, t = "init") {
    # This is a wrapper function for reading initial composition function
    # input_mapping <- input_withinCommunityPairs
    # treatment <- "withinCommunityPairs"
    # input_mapping <- input_poolPairs
    # treatment <- "poolPairs"
    # comm <- "W0"
    # t <- "init"
    paste0(input_mapping$output_dir[1], treatment, "_", comm, "-1-N_", t,".csv") %>%
        read_wide_file() %>%
        filter(Abundance != 0) %>%
        mutate(Community = comm, Time = t) %>%
        format_columns() %>%
        add_pairID() %>%
        arrange(Community, Well, Species) %>%
        select(Community, Pair, Species, Well, Time, Abundance)

}
read_later_composition <- function (input_mapping, treatment, comm, t = "end", pairs_N_sp) {
    #' This is a wrapper function
    #' The reason it's separated from the function read_init_compotition() is that
    #' the later time point would have dropped data if one species went extinct.
    #' This function deals with that by joining the incomplete data with a table
    #' of pairs
    # input_mapping <- input_withinCommunityPairs
    # treatment <- "withinCommunityPairs"
    # comm <- "W1"
    # t <- "init"
    # pair_N_sp <- withinCommunityPairs_N_sp
    paste0(input_mapping$output_dir[1], treatment, "_", comm, "-1-N_", t,".csv") %>%
        read_wide_file() %>%
        filter(Abundance != 0) %>%
        right_join(pairs_N_sp, by = join_by(Species, Well)) %>%
        replace_na(list(Abundance = 0, Time = t)) %>%
        mutate(Community = comm, Time = t) %>%
        format_columns() %>%
        arrange(Community, Well, Species) %>%
        select(Community, Pair, Species, Well, Time, Abundance)
}
bind_pair_frequency <- function (ini, end) {
    #bind_rows(poolPairs_N_init, poolPairs_N_end) %>%
    bind_rows(ini, end) %>%
        mutate(InitialFrequency = rep(rep(c(5, 50, 95), each = 2), n()/6)) %>%
        select(-Well) %>%
        select(Community, Pair, Species, InitialFrequency, Time, Abundance) %>%
        mutate(Pair = factor(Pair, paste0("P", 1:1000))) %>%
        # Relative abundance
        arrange(Community, Pair, InitialFrequency, Time, Species) %>%
        group_by(Community, Pair, InitialFrequency, Time) %>%
        mutate(Frequency = round(Abundance / sum(Abundance), 4)) %>%
        #
        select(Community, Pair, Species, InitialFrequency, Time, Frequency) %>%
        mutate(Isolate = c(1,2)) %>%
        group_by(Community, Pair, InitialFrequency) %>%
        pivot_wider(names_from = Isolate, values_from = c(Species, Frequency), names_sep = "") %>%
        ungroup() %>%
        select(Community, Pair, Species1, Species2, InitialFrequency, Time, Frequency1, Frequency2)
}

## For determining competition outcome
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
    interaction_type$InteractionType[c(7,12,15:17,19,22)] <- "unknown"

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
    interaction_type$InteractionTypeFiner[c(7,12,15:17,19,22)] <- "unknown"

    interaction_type <- interaction_type %>%  mutate(FitnessFunction = paste(FromRare, FromMedium, FromAbundant, sep = "_"))
}
append_pairs_outocme <- function (pairs_fitness) {
    # This function appends the bootstrap result to the fitness function table
    # Modified from 93-determine_competition.R
    pairs_fitness %>%
        pivot_wider(id_cols = c(Community, Species1, Species2), names_from = InitialFrequency, values_from = FitnessChange) %>%
        rename(FromRare = `5`, FromMedium = `50`, FromAbundant = `95`) %>%
        left_join(interaction_type, by = join_by(FromRare, FromMedium, FromAbundant)) %>%
        select(Community, Species1, Species2, InteractionType, InteractionTypeFiner, FitnessFunction)

}
interaction_type <- make_interaction_type()
compute_pairwise_outcome <- function (pairs_frequency) {
    #poolPairs_N_freq %>%
    pairs_frequency %>%
        # Compute fitness table
        select(Community, Pair, Species1, Species2, InitialFrequency, Time, Frequency1) %>%
        group_by(Community, Pair, Time) %>%
        pivot_wider(names_from = Time, values_from = Frequency1) %>%
        mutate(FitnessChange = case_when(
            end > init ~ 1,
            end < init ~ -1
        )) %>%
        # Compute fitness function
        append_pairs_outocme() %>%
        mutate(From = case_when(
            FitnessFunction %in% c("1_1_1") ~ Species1, # Isolate1 wins
            FitnessFunction %in% c("-1_-1_-1") ~ Species2, # Isolate2 wins
            TRUE ~ Species1
        )) %>%
        mutate(To = case_when(
            FitnessFunction %in% c("1_1_1") ~ Species2,
            FitnessFunction %in% c("-1_-1_-1") ~ Species1,
            TRUE ~ Species2
        )) %>%
        ungroup()
}


# 1. Pool pairs ----
communities_names <- paste0("W", 0:19)
poolPairs_N_freq <- rep(list(NA), length(communities_names))
poolPairs_N_outcome <- rep(list(NA), length(communities_names))

for (i in 1:20) {
    cat("\nMonocultureSet ", communities_names[i])
    # Read initial composition
    poolPairs_N_init <- read_init_composition(input_poolPairs, "poolPairs", comm = communities_names[i], t = "init")

    # Mapping file of pairs
    poolPairs_N_sp <- poolPairs_N_init %>% distinct(Community, Pair, Species, Well)

    # Read end composition
    poolPairs_N_end <- read_later_composition(input_poolPairs, "poolPairs", comm = communities_names[i],
                                              t = paste0("T", 1:temp_end_time), pairs_N_sp = poolPairs_N_sp) %>%
        mutate(Time = ifelse(Time == paste0("T", 1:temp_end_time), "end", Time))

    # Check if the init and end has the same number of rows
    stopifnot(nrow(poolPairs_N_init) == nrow(poolPairs_N_end))

    # Pair frequency
    poolPairs_N_freq[[i]] <- bind_pair_frequency(poolPairs_N_init, poolPairs_N_end)

    # Determine pairwise outcome
    poolPairs_N_outcome[[i]] <- compute_pairwise_outcome(poolPairs_N_freq[[i]])
}

poolPairs_N_freq <- bind_rows(poolPairs_N_freq[!is.na(poolPairs_N_freq)])
poolPairs_N_outcome <- bind_rows(poolPairs_N_outcome[!is.na(poolPairs_N_outcome)])

write_csv(poolPairs_N_freq, paste0(folder_simulation, "aggregated/13-poolPairs_N_freq.csv"))
write_csv(poolPairs_N_outcome, paste0(folder_simulation, "aggregated/13-poolPairs_N_outcome.csv"))


# 2. Community pairs ----
communities_richness <- read_csv(paste0(folder_simulation, "aggregated/12-communities_richness.csv"), col_types = cols())
communities_names <- paste0("W", 0:19)
withinCommunityPairs_N_freq <- rep(list(NA), length(communities_names))
withinCommunityPairs_N_outcome <- rep(list(NA), length(communities_names))

for (i in 1:20) {
    richness <- communities_richness %>%
        filter(Community == communities_names[i]) %>%
        pull(Richness)
    cat("\nCommunity", communities_names[i], " richness = ", richness)
    if (richness > 1) {
        # Read initial composition
        withinCommunityPairs_N_init <- read_init_composition(input_withinCommunityPairs, "withinCommunityPairs", comm = communities_names[i], t = "init")

        # Mapping file of pairs
        withinCommunityPairs_N_sp <- withinCommunityPairs_N_init %>% distinct(Community, Pair, Species, Well)

        # Read end composition
        withinCommunityPairs_N_end <- read_later_composition(input_withinCommunityPairs, "withinCommunityPairs", comm = communities_names[i],
                                                             t = paste0("T", 1:temp_end_time), pairs_N_sp = withinCommunityPairs_N_sp) %>%
            mutate(Time = ifelse(Time == paste0("T", 1:temp_end_time), "end", Time))

        # Check if the init and end has the same number of rows
        stopifnot(nrow(withinCommunityPairs_N_init) == nrow(withinCommunityPairs_N_end))

        # Pair frequency
        withinCommunityPairs_N_freq[[i]] <- bind_pair_frequency(withinCommunityPairs_N_init, withinCommunityPairs_N_end)

        # Determine pairwise outcome
        withinCommunityPairs_N_outcome[[i]] <- compute_pairwise_outcome(withinCommunityPairs_N_freq[[i]])
    }
}

withinCommunityPairs_N_freq <- bind_rows(withinCommunityPairs_N_freq[!is.na(withinCommunityPairs_N_freq)])
withinCommunityPairs_N_outcome <- bind_rows(withinCommunityPairs_N_outcome[!is.na(withinCommunityPairs_N_outcome)])

write_csv(withinCommunityPairs_N_freq, paste0(folder_simulation, "aggregated/13-withinCommunityPairs_N_freq.csv"))
write_csv(withinCommunityPairs_N_outcome, paste0(folder_simulation, "aggregated/13-withinCommunityPairs_N_outcome.csv"))





