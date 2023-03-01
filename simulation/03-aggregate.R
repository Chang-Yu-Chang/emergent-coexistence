#' This script aggregates the data from communities and monocultures

library(tidyverse)
library(cowplot)
source(here::here("analysis/00-metadata.R"))

# 0. parameters ----
# Read parameters
input_parameters <- read_csv(here::here("simulation/01-input_parameters.csv"), col_types = cols())
input_monocultures <- read_csv(here::here("simulation/01a-input_monocultures.csv"), col_types = cols())
input_communities <- read_csv(here::here("simulation/01b-input_communities.csv"), col_types = cols())
input_poolPairs <- read_csv(here::here("simulation/01c-input_poolPairs.csv"), col_types = cols())
input_withinCommunityPairs <- read_csv(here::here("simulation/01d-input_withinCommunityPairs.csv"), col_types = cols())
category_colors <- c(sugar = "#ED6A5A", acid = "#03CEA4", waste = "#51513D", fermenter = "#8A89C0", respirator = "#FFCB77")
input_row <- input_parameters[1,]

# Generate family-species and class-resource tibble for matching
sa <- input_parameters$sa[1]
ma <- input_parameters$ma[1]
sal <- tibble(Family = paste0("F", c(rep(0, sa), rep(1, sa))), Species = paste0("S", 0:(sa * 2 - 1)))
mal <- tibble(Class = paste0("T", c(rep(0, ma), rep(1, ma))), Resource = paste0("R", 0:(ma * 2 - 1)))

# 1. Self-assembled community composition ----
# Read simulation output files
read_wide_file <- function(x, type = "N") {
    temp <- read_csv(x, col_types = cols(), name_repair = "unique_quiet") %>%
        pivot_longer(cols = starts_with("W"), names_to = "Well", values_to = "Abundance")
    if ("...1" %in% colnames(temp)) {
        if (type == "N") temp <- temp %>% rename(Family = ...1, Species = ...2)
        if (type == "R") temp <- temp %>% rename(Class = ...1, Resource = ...2)
    }

    return(temp)
}

df_communities_N <- list.files(input_communities$output_dir[1], pattern = "selfAssembly-1-N") %>%
    map(c("N_T\\d+\\.csv", "init"), str_subset, string = .) %>%
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
    mutate(Well = ordered(Well, paste0("W", 0:(input_row$n_wells-1)))) %>%
    mutate(Time = ordered(Time, c("init", paste0("T", 1:20)))) %>%
    arrange(Well, Time)

write_csv(df_communities_N, paste0(folder_simulation, "11-aggregated/df_communities_N.csv"))

# 2. Monoculture ----
df_monocultures_N <- list.files(input_monocultures$output_dir[1], pattern = "monoculture-1-N") %>%
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
    mutate(Well = ordered(Well, paste0("W", 0:(input_row$n_wells-1)))) %>%
    mutate(Time = ordered(Time, c("init", paste0("T", 1:20)))) %>%
    arrange(Well, Time) %>%
    filter(Abundance > 0)

write_csv(df_monocultures_N, paste0(folder_simulation, "11-aggregated/df_monocultures_N.csv"))

# 3. Pool pairs ----
# Read pair data
df_poolPairs_N <- list.files(input_poolPairs$output_dir[1], pattern = "poolPairs_W0-1-N") %>%
    map(c("N_T\\d+\\.csv", "init"), str_subset, string = .) %>%
    reduce(union) %>%
    paste0(input_poolPairs$output_dir[1], .) %>%
    lapply(function(x) {
        time_point <- str_split_i(x, "-", -1) %>%
            str_replace("N_", "") %>%
            str_replace(".csv", "")
        read_wide_file(x) %>%
            mutate(Time = time_point) %>%
            return()
    }) %>%
    bind_rows %>%
    mutate(Well = ordered(Well, paste0("W", 0:(input_row$n_wells-1)))) %>%
    mutate(Time = ordered(Time, c("init", paste0("T", 1:20)))) %>%
    arrange(Well, Time)





