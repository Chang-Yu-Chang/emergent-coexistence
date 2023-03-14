#' This script aggregates the data from 02a-monocultures and 02b-communities

library(tidyverse)
library(cowplot)
source(here::here("analysis/00-metadata.R"))
source(here::here("simulation/01-generate_input.R"))

# 0. parameters ----
input_parameters <- read_csv(here::here("simulation/01-input_parameters.csv"), col_types = cols())
input_monocultures <- read_csv(here::here("simulation/02a-input_monocultures.csv"), col_types = cols())
input_communities <- read_csv(here::here("simulation/02b-input_communities.csv"), col_types = cols())
input_communitiesWithoutCrossfeeding <- read_csv(here::here("simulation/02c-input_communitiesWithoutCrossfeeding.csv"), col_types = cols())

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
    #filter(Abundance > 0) %>%
    select(Well, Time, Family, Species, Abundance)

write_csv(monocultures_abundance, paste0(folder_simulation, "aggregated/12-monocultures_abundance.csv"))

#write_csv(monocultures_abundance_richness, paste0(folder_simulation, "aggregated/12-monocultureSets_abundance_richness.csv"))

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
    #filter(Abundance > 0) %>%
    select(Community, Time, Family, Species, Abundance)

write_csv(communities_abundance, paste0(folder_simulation, "aggregated/12-communities_abundance.csv"))

communities_abundance_richness <- communities_abundance %>%
    group_by(Community, .drop = F) %>%
    filter(Time == max(Time)) %>%
    group_by(Community) %>%
    filter(Abundance > 0.01*sum(Abundance)) %>%
    summarize(Richness = n())

write_csv(communities_abundance_richness, paste0(folder_simulation, "aggregated/12-communities_richness.csv"))


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

write_csv(communitiesWithoutCrossfeeding_abundance, paste0(folder_simulation, "aggregated/12-communitiesWithoutCrossfeeding_abundance.csv"))

communitiesWithoutCrossfeeding_abundance <- communitiesWithoutCrossfeeding_abundance %>%
    group_by(Community, .drop = F) %>%
    filter(Time == max(Time)) %>%
    group_by(Community) %>%
    filter(Abundance > 0.01*sum(Abundance)) %>%
    summarize(Richness = n())

write_csv(communitiesWithoutCrossfeeding_abundance, paste0(folder_simulation, "aggregated/12-communitiesWithoutCrossfeeding_richness.csv"))


















