#' This script generates the essential input for simulations
#'
#' 1. Generates the mapping files for simulations
#'     - 03a-input_poolPairs.csv
#'     - 03b-input_withinCommunityPairs.csv
#' 2. Generates the initial composition for these two dependent simulations

library(tidyverse)
source(here::here("analysis/00-metadata.R"))
source(here::here("simulation/01-generate_input.R"))
source(here::here("simulation/02-generate_input_mono_comm.R"))

# 1. generate mapping files ----
input_parameters <- read_csv(here::here("simulation/01-input_parameters.csv"), col_types = cols())

# Pairs from within communities
input_parameters %>%
    slice(rep(1, input_parameters$n_comm)) %>%
    mutate(save_timepoint = F) %>%
    mutate(output_dir = paste0(folder_simulation, "04a-withinCommunityTrios/")) %>%
    mutate(init_N0 = paste0("withinCommunityTrios_W", 0:(input_parameters$n_comm-1), "-1-N_init.csv")) %>%
    mutate(init_R0 = paste0("withinCommunityTrios_W", 0:(input_parameters$n_comm-1), "-1-R_init.csv")) %>%
    write_csv(here::here("simulation/04a-input_withinCommunityTrios.csv"))


# 2. generate initial composition files ----
# Note that n_wells in the mapping files will be modified depending on the diversity of the monocultureSets/community
input_monocultures <- read_csv(here::here("simulation/02a-input_monocultures.csv"), col_types = cols())
input_communities <- read_csv(here::here("simulation/02b-input_communities.csv"), col_types = cols())
input_withinCommunityTrios <- read_csv(here::here("simulation/04a-input_withinCommunityTrios.csv"), col_types = cols())

# 2.1. Generate within community pairs ----
draw_trios_from_community <- function(N_community_long) {
    # Communities with no or only one or two species
    if (nrow(N_community_long) <= 2) N_trios <- sal %>% mutate(W0 = 0)

    # Other communities with more one species
    if (nrow(N_community_long) > 2) {
        N_trios <- expand_grid(sp1 = N_community_long$Species, sp2 = N_community_long$Species, sp3 = N_community_long$Species) %>%
            mutate(across(everything(), ~ ordered(.x, sal$Species))) %>%
            filter(sp1 < sp2, sp2 < sp3, sp1 < sp3) %>%
            mutate(Trio = paste0("K", 0:(n()-1))) %>%
            # Make duplicate for four freq
            slice(rep(1:n(), each = 4)) %>%
            mutate(Well = paste0("W", 0:(n()-1))) %>%
            group_by(Well) %>%
            pivot_longer(cols = starts_with("sp"), names_to = "temp", values_to = "Species") %>%
            ungroup() %>%
            mutate(Abundance = rep(c(0.05, 0.05, 0.95,
                                     0.05, 0.95, 0.05,
                                     0.95, 0.05, 0.05,
                                     0.33, 0.33, 0.33), n()/12)) %>%
            full_join(sal, by = "Species") %>%
            replace_na(list(Well = "W0", Abundance = 0)) %>%
            select(Trio, Well, Family, Species, Abundance) %>%
            pivot_wider(id_cols = c(Family, Species), names_from = Well, values_from = Abundance, values_fill = 0) %>%
            mutate(Species = factor(Species, sal$Species)) %>%
            arrange(Species)
    }
    return(N_trios)
}

N_community <- read_csv(paste0(input_communities$output_dir[1], "communities-1-N_end.csv"), col_types = cols()) %>% rename(Family = ...1, Species = ...2)
N_community_split <- N_community %>%
    mutate_all(~replace(., .==0, NA)) %>%
    pivot_longer(cols = -c(Family, Species), names_to = "Community", values_to = "Abundance", values_drop_na = T) %>%
    arrange(Community) %>%
    # Remove rare species (relative abundance <0.01)
    group_by(Community) %>%
    mutate(RelativeAbundance = Abundance/sum(Abundance)) %>%
    filter(RelativeAbundance > 0.01) %>%
    # Order species and community
    mutate(Species = ordered(Species, sal$Species)) %>%
    mutate(Community = ordered(Community, colnames(N_community))) %>%
    arrange(Community, Species) %>%
    group_by(Community) %>%
    group_split()

# Some communities only have two species. Skip those communities
communities_richness <- N_community_split %>%
    bind_rows() %>%
    group_by(Community) %>%
    summarize(Richness = n())

communities_species <- N_community_split %>%
    bind_rows() %>%
    select(Community, Family, Species) %>%
    arrange(Community, Family, Species)

for (i in 1:length(N_community_split)) {
    cat("\n", paste0("Community W", i-1, " Richness=", communities_richness$Richness[i]))
    if (communities_richness$Richness[i] <= 2) next
    # init_N0
    draw_trios_from_community(N_community_split[[i]]) %>%
        write_csv(paste0(input_withinCommunityTrios$output_dir[i], "withinCommunityTrios_W", i-1, "-1-N_init.csv"))
    # Update n_wells in the input files
    input_withinCommunityTrios$n_wells[i] <- choose(nrow(N_community_split[[i]]), 3) * 4
    # init_R0
    set_community_resource(input_withinCommunityTrios[i,]) %>%
        write_csv(paste0(input_withinCommunityTrios$output_dir[i], "withinCommunityTrios_W", i-1, "-1-R_init.csv"))

}

# Update the mapping file
write_csv(input_withinCommunityTrios, here::here("simulation/04a-input_withinCommunityTrios.csv"))


