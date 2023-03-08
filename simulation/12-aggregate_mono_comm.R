#' This script aggregates the data from 02a-monocultures and 02b-communities

library(tidyverse)
library(cowplot)
source(here::here("analysis/00-metadata.R"))

# 0. parameters ----
input_parameters <- read_csv(here::here("simulation/01-input_parameters.csv"), col_types = cols())
input_monocultures <- read_csv(here::here("simulation/02a-input_monocultures.csv"), col_types = cols())
input_communities <- read_csv(here::here("simulation/02b-input_communities.csv"), col_types = cols())
input_communitiesWithoutCrossfeeding <- read_csv(here::here("simulation/02c-input_communitiesWithoutCrossfeeding.csv"), col_types = cols())

# Generate family-species and class-resource table for matching
sa <- input_parameters$sa[1]
ma <- input_parameters$ma[1]
sal <- tibble(Family = paste0("F", c(rep(0, sa), rep(1, sa))), Species = paste0("S", 0:(sa * 2 - 1)))
mal <- tibble(Class = paste0("T", c(rep(0, ma), rep(1, ma))), Resource = paste0("R", 0:(ma * 2 - 1)))

read_wide_file <- function(x, type = "N") {
    temp <- read_csv(x, col_types = cols(), name_repair = "unique_quiet") %>%
        pivot_longer(cols = starts_with("W"), names_to = "Well", values_to = "Abundance")
    if ("...1" %in% colnames(temp)) {
        if (type == "N") temp <- temp %>% rename(Family = ...1, Species = ...2)
        if (type == "R") temp <- temp %>% rename(Class = ...1, Resource = ...2)
    }

    return(temp)
}

# 1. Monoculture ----
monocultures_abundance <- list.files(input_monocultures$output_dir[1], pattern = "monoculture-1-N") %>%
    map(c("N_T\\d+\\.csv", "init"), str_subset, string = .) %>%
    reduce(union) %>%
    paste0(input_monocultures$output_dir[1], .) %>%
    lapply(function(x) {
        time_point <- str_split_i(x, "-", -1) %>%
            str_replace("N_", "") %>%
            str_replace(".csv", "")
        read_wide_file(x) %>%
            mutate(Time = time_point) %>%
            return()
    }) %>%
    bind_rows %>%
    mutate(Well = ordered(Well, paste0("W", 0:(input_communities$n_wells[1]-1)))) %>%
    mutate(Time = ordered(Time, c("init", paste0("T", 1:20)))) %>%
    arrange(Well, Time) %>%
    filter(Abundance > 0) %>%
    select(Well, Time, Family, Species, Abundance)

write_csv(monocultures_abundance, paste0(folder_simulation, "11-aggregated/monocultures_abundance.csv"))


# 2. Self-assembled community composition ----
communities_abundance <- list.files(input_communities$output_dir[1], pattern = "communities-1-N") %>%
    map(c("N_T\\d+\\.csv", "init", "end"), str_subset, string = .) %>%
    reduce(union) %>%
    paste0(input_communities$output_dir[1], .) %>%
    lapply(function(x) {
        time_point <- str_split_i(x, "-", -1) %>%
            str_replace("N_", "") %>%
            str_replace(".csv", "")
        read_wide_file(x) %>%
            mutate(Time = time_point) %>%
            return()
    }) %>%
    bind_rows %>%
    mutate(Community = ordered(Well, paste0("W", 0:(input_communities$n_wells[1]-1)))) %>%
    select(-Well) %>%
    mutate(Time = ordered(Time, c("init", paste0("T", 1:20), "end"))) %>%
    arrange(Community, Time) %>%
    filter(Abundance > 0) %>%
    select(Community, Time, Family, Species, Abundance)

write_csv(communities_abundance, paste0(folder_simulation, "11-aggregated/communities_abundance.csv"))


# 3. Communities without crossfeeding ----
communitiesWithoutCrossfeeding_abundance <- list.files(input_communitiesWithoutCrossfeeding$output_dir[1], pattern = "communities-1-N") %>%
    map(c("N_T\\d+\\.csv", "init", "end"), str_subset, string = .) %>%
    reduce(union) %>%
    paste0(input_communitiesWithoutCrossfeeding$output_dir[1], .) %>%
    lapply(function(x) {
        time_point <- str_split_i(x, "-", -1) %>%
            str_replace("N_", "") %>%
            str_replace(".csv", "")
        read_wide_file(x) %>%
            mutate(Time = time_point) %>%
            return()
    }) %>%
    bind_rows %>%
    mutate(Community = ordered(Well, paste0("W", 0:(input_communitiesWithoutCrossfeeding$n_wells[1]-1)))) %>%
    select(-Well) %>%
    mutate(Time = ordered(Time, c("init", paste0("T", 1:20), "end"))) %>%
    arrange(Community, Time) %>%
    filter(Abundance > 0) %>%
    select(Community, Time, Family, Species, Abundance)

write_csv(communitiesWithoutCrossfeeding_abundance, paste0(folder_simulation, "11-aggregated/communitiesWithoutCrossfeeding_abundance.csv"))



















