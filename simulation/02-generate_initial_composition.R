#' This script generates the initial composition for monocultures and communities
#'
#' Because these two simulations are independent of other simulations, this script can be run
#' along with the 01-generate_input.R

library(tidyverse)
source(here::here("analysis/00-metadata.R"))

# Read the parameters
input_parameters <- read_csv(here::here("simulation/01-input_parameters.csv"), col_types = cols())
input_monocultures <- read_csv(here::here("simulation/01a-input_monocultures.csv"), col_types = cols())
input_communities <- read_csv(here::here("simulation/01b-input_communities.csv"), col_types = cols())

# Generate family-species and class-resource matching tibble
# Note that input_independent has to use the same sa and ma throughout
sa <- input_parameters$sa[1]
ma <- input_parameters$ma[1]
sal <- tibble(Family = paste0("F", c(rep(0, sa), rep(1, sa))), Species = paste0("S", 0:(sa * 2 - 1)))
mal <- tibble(Class = paste0("T", c(rep(0, ma), rep(1, ma))), Resource = paste0("R", 0:(ma * 2 - 1)))

# Monoculture
draw_monoculture <- function(input_monocultures) {
    n_wells <- input_monocultures$n_wells
    tibble(
        Well = paste0("W", 0:(n_wells-1)),
        Species = paste0("S", sample(0:(sa*2-1), n_wells, replace = F)), # Sample up to 200 from the global pool
        Abundance = 1
    ) %>%
        full_join(sal, by = "Species") %>%
        # Fill the non-chosen species
        replace_na(list(Well = "W0", Abundance = 0)) %>%
        pivot_wider(id_cols = c(Family, Species), names_from = Well, values_from = Abundance, values_fill = 0) %>%
        mutate(Species = factor(Species, sal$Species)) %>%
        arrange(Species)

}
set_monoculture_resource <- function (input_monocultures) {
    n_wells <- input_monocultures$n_wells
    tibble(
        Well = paste0("W", 0:(n_wells-1)),
        Resource = "R0",
        Abundance = input_monocultures$R0_food
    ) %>%
        full_join(mal, by = "Resource") %>%
        pivot_wider(id_cols = c(Class, Resource), names_from = Well, values_from = Abundance, values_fill = 0) %>%
        select(-`NA`) %>%
        mutate(Resource = factor(Resource, mal$Resource)) %>%
        arrange(Resource)
}
draw_monoculture(input_monocultures) %>%
    write_csv(paste0(input_monocultures$output_dir, "monoculture-1-N_init.csv"))
set_monoculture_resource(input_monocultures) %>%
    write_csv(paste0(input_monocultures$output_dir, "monoculture-1-R_init.csv"))

# Self-assembly. Draw a fix number of S
draw_community <- function(input_communities, candidate_species = NA) {
    S <- input_communities$S
    n_communities <- input_communities$n_communities
    if (all(is.na(candidate_species))) {
        # Draw species. Default using a fix number of S. Identical to community-simulator
        species <- rep(list(NA), n_communities)
        for (i in 1:n_communities) species[[i]] <- sample(0:(2*sa-1), size = S, replace = F) %>% sort
        species <- unlist(species)
    }
    if (!all(is.na(candidate_species)) & is.vector(candidate_species)) {
        candidate_species <- str_replace(candidate_species, "S", "") %>% as.numeric()
        # Draw species. Default using a fix number of S. Identical to community-simulator
        species <- rep(list(NA), n_communities)
        for (i in 1:n_communities) species[[i]] <- sample(candidate_species, size = S, replace = F) %>% sort
        species <- unlist(species)
    }

    # Make a long list
    N <- tibble(Well = rep(paste0("W", 0:(n_communities-1)), each = S),
                Species = paste0("S", species),
                Abundance = rep(1, n_communities * S)) %>%
        full_join(sal, by = "Species") %>%
        # Fill the non-chosen species
        replace_na(list(Well = "W0", Abundance = 0)) %>%
        pivot_wider(id_cols = c(Family, Species), names_from = Well, values_from = Abundance, values_fill = 0) %>%
        mutate(Species = factor(Species, sal$Species)) %>%
        arrange(Species)

    return(N)
}
set_community_resource <- function (input_communities) {
    n_wells <- max(input_communities$n_wells, 1)
    tibble(
        Well = paste0("W", 0:(n_wells-1)),
        Resource = "R0",
        Abundance = input_communities$R0_food
    ) %>%
        right_join(mal, by = "Resource") %>%
        pivot_wider(id_cols = c(Class, Resource), names_from = Well, values_from = Abundance, values_fill = 0) %>%
        select(-`NA`) %>%
        mutate(Resource = factor(Resource, mal$Resource)) %>%
        arrange(Resource)
}

set.seed(1)
draw_community(input_communities) %>%
    write_csv(paste0(input_communities$output_dir, "selfAssembly-1-N_init.csv"))

set_community_resource(input_communities) %>%
    write_csv(paste0(input_communities$output_dir, "selfAssembly-1-R_init.csv"))

