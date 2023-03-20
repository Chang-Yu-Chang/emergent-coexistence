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

# Pairs from the pool
n_comm <- input_parameters$n_communities
input_parameters %>%
    slice(rep(1, n_comm)) %>%
    mutate(save_timepoint = T) %>%
    mutate(output_dir = paste0(folder_simulation, "03a-poolPairs/")) %>%
    mutate(init_N0 = paste0("poolPairs_W", 0:(n_comm-1), "-1-N_init.csv"), exp_id = 1, S = 10) %>%
    mutate(init_R0 = paste0("poolPairs_W", 0:(n_comm-1), "-1-R_init.csv")) %>%
    write_csv(here::here("simulation/03a-input_poolPairs.csv"))

# Pairs from within self-assembled communities
input_parameters %>%
    slice(rep(1, n_comm)) %>%
    mutate(save_timepoint = T) %>%
    mutate(output_dir = paste0(folder_simulation, "03b-withinCommunityPairs/")) %>%
    mutate(init_N0 = paste0("withinCommunityPairs_W", 0:(n_comm-1), "-1-N_init.csv")) %>%
    mutate(init_R0 = paste0("withinCommunityPairs_W", 0:(n_comm-1), "-1-R_init.csv")) %>%
    write_csv(here::here("simulation/03b-input_withinCommunityPairs.csv"))

# Pairs for fitting LV models
n_comm <- input_parameters$n_communities
input_parameters %>%
    slice(rep(1, n_comm)) %>%
    mutate(save_timepoint = T, n_timepoint = 50, n_pass = 2) %>%
    mutate(output_dir = paste0(folder_simulation, "03c-LVPairs/")) %>%
    mutate(init_N0 = paste0("LVPairs_W", 0:(n_comm-1), "-1-N_init.csv"), exp_id = 1, S = 10) %>%
    mutate(init_R0 = paste0("LVPairs_W", 0:(n_comm-1), "-1-R_init.csv")) %>%
    write_csv(here::here("simulation/03c-input_LVPairs.csv"))


"Execute the chunks below after monocultures and communities are done"

# 2. generate initial composition files ----
# Note that n_wells in the mapping files will be modified depending on the diversity of the monocultureSets/community
input_monocultures <- read_csv(here::here("simulation/02a-input_monocultures.csv"), col_types = cols())
input_communities <- read_csv(here::here("simulation/02b-input_communities.csv"), col_types = cols())
input_poolPairs <- read_csv(here::here("simulation/03a-input_poolPairs.csv"), col_types = cols())
input_withinCommunityPairs <- read_csv(here::here("simulation/03b-input_withinCommunityPairs.csv"), col_types = cols())
input_LVPairs <- read_csv(here::here("simulation/03c-input_LVPairs.csv"), col_types = cols())

# 2.1. Generate competing pairs from culturable monocultures ----
draw_pairs_from_community <- function(N_community_long) {
    # Communities with no or only one species
    if (nrow(N_community_long) <= 1) N_pairs <- sal %>% mutate(W0 = 0)

    # Other communities with more one species
    if (nrow(N_community_long) > 1) {
        N_pairs <- expand_grid(sp1 = N_community_long$Species, sp2 = N_community_long$Species) %>%
            mutate(across(everything(), ~ ordered(.x, sal$Species))) %>%
            filter(sp1 < sp2) %>%
            mutate(Pair = paste0("P", 0:(n()-1))) %>%
            # Make duplicate for three freq
            slice(rep(1:n(), each = 3)) %>%
            mutate(Well = paste0("W", 0:(n()-1))) %>%
            group_by(Well) %>%
            pivot_longer(cols = starts_with("sp"), names_to = "temp", values_to = "Species") %>%
            ungroup() %>%
            mutate(Abundance = rep(c(0.05, 0.95, 0.5, 0.5, 0.95, 0.05), n()/6)) %>%
            full_join(sal, by = "Species") %>%
            replace_na(list(Well = "W0", Abundance = 0)) %>%
            select(Pair, Well, Family, Species, Abundance) %>%
            pivot_wider(id_cols = c(Family, Species), names_from = Well, values_from = Abundance, values_fill = 0) %>%
            mutate(Species = factor(Species, sal$Species)) %>%
            arrange(Species)
    }
    return(N_pairs)
}
set.seed(1)


species_mono <- read_csv(paste0(input_monocultures$output_dir[1], "monoculture-1-N_end.csv"), col_types = cols()) %>%
    rename(Family = ...1, Species = ...2) %>%
    mutate_all(~replace(., .==0, NA)) %>%
    pivot_longer(cols = starts_with("W"), names_to = "Well", values_to = "Abundance", values_drop_na = T) %>%
    pull(Species) %>%
    str_replace("S", "") %>% as.numeric()

N_monocultureSets <- draw_community(input_poolPairs[1,], species_mono)

N_monocultureSets_split <- N_monocultureSets %>%
    mutate_all(~replace(., .==0, NA)) %>%
    pivot_longer(cols = -c(Family, Species), names_to = "Community", values_to = "Abundance", values_drop_na = T) %>%
    # Order species and community
    mutate(Species = ordered(Species, sal$Species)) %>%
    mutate(Community = ordered(Community, colnames(N_monocultureSets))) %>%
    arrange(Community, Species) %>%
    group_by(Community) %>%
    group_split()


monocultureSets_richness <- N_monocultureSets_split %>%
    bind_rows() %>%
    group_by(Community) %>%
    summarize(Richness = n())

monocultureSets_species <- N_monocultureSets_split %>%
    bind_rows() %>%
    select(Community, Family, Species) %>%
    arrange(Community, Family, Species)

for (i in 1:length(N_monocultureSets_split)) {
    cat("\n", paste0("Pool set W", i-1, " Richness=",  monocultureSets_richness$Richness[i]))
    # init_N0
    draw_pairs_from_community(N_monocultureSets_split[[i]]) %>%
        write_csv(paste0(input_poolPairs$output_dir[i], "poolPairs_W", i-1, "-1-N_init.csv"))
    # Update n_wells in the input files
    input_poolPairs$n_wells[i] <- choose(nrow(N_monocultureSets_split[[i]]), 2) * 3
    # init_R0
    set_community_resource(input_poolPairs[i,]) %>%
        write_csv(paste0(input_poolPairs$output_dir[i], "poolPairs_W", i-1, "-1-R_init.csv"))
}

write_csv(input_poolPairs, here::here("simulation/03a-input_poolPairs.csv"))
write_csv(monocultureSets_richness, paste0(folder_simulation, "aggregated/03-monocultureSets_richness.csv"))
write_csv(monocultureSets_species, paste0(folder_simulation, "aggregated/03-monocultureSets_species.csv"))


# 2.2. Generate within community pairs ----
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
    # init_N0
    draw_pairs_from_community(N_community_split[[i]]) %>%
        write_csv(paste0(input_withinCommunityPairs$output_dir[i], "withinCommunityPairs_W", i-1, "-1-N_init.csv"))
    # Update n_wells in the input files
    input_withinCommunityPairs$n_wells[i] <- choose(nrow(N_community_split[[i]]), 2) * 3
    # init_R0
    set_community_resource(input_withinCommunityPairs[i,]) %>%
        write_csv(paste0(input_withinCommunityPairs$output_dir[i], "withinCommunityPairs_W", i-1, "-1-R_init.csv"))

}

write_csv(input_withinCommunityPairs, here::here("simulation/03b-input_withinCommunityPairs.csv"))
#write_csv(communities_richness, paste0(folder_simulation, "11-aggregated/communities_richness.csv"))
write_csv(communities_species, paste0(folder_simulation, "aggregated/03-communities_species.csv"))



# 2.3. Generate LV pair. Identical to pool pairs ----

set.seed(1)
species_mono <- read_csv(paste0(input_monocultures$output_dir[1], "monoculture-1-N_T5.csv"), col_types = cols()) %>%
    rename(Family = ...1, Species = ...2) %>%
    mutate_all(~replace(., .==0, NA)) %>%
    pivot_longer(cols = starts_with("W"), names_to = "Well", values_to = "Abundance", values_drop_na = T) %>%
    pull(Species) %>%
    str_replace("S", "") %>% as.numeric()

N_monocultureSets <- draw_community(input_poolPairs[1,], species_mono)

N_monocultureSets_split <- N_monocultureSets %>%
    mutate_all(~replace(., .==0, NA)) %>%
    pivot_longer(cols = -c(Family, Species), names_to = "Community", values_to = "Abundance", values_drop_na = T) %>%
    # Order species and community
    mutate(Species = ordered(Species, sal$Species)) %>%
    mutate(Community = ordered(Community, colnames(N_monocultureSets))) %>%
    arrange(Community, Species) %>%
    group_by(Community) %>%
    group_split()


monocultureSets_richness <- N_monocultureSets_split %>%
    bind_rows() %>%
    group_by(Community) %>%
    summarize(Richness = n())

monocultureSets_species <- N_monocultureSets_split %>%
    bind_rows() %>%
    select(Community, Family, Species) %>%
    arrange(Community, Family, Species)

for (i in 1) {
    cat("\n", paste0("Pool set W", i-1, " Richness=",  monocultureSets_richness$Richness[i]))
    # init_N0
    draw_pairs_from_community(N_monocultureSets_split[[i]]) %>%
        write_csv(paste0(input_LVPairs$output_dir[i], "LVPairs_W", i-1, "-1-N_init.csv"))
    # Update n_wells in the input files
    input_LVPairs$n_wells[i] <- choose(nrow(N_monocultureSets_split[[i]]), 2) * 3
    # init_R0
    set_community_resource(input_LVPairs[i,]) %>%
        write_csv(paste0(input_LVPairs$output_dir[i], "LVPairs_W", i-1, "-1-R_init.csv"))
}

if (FALSE) {
for (i in 1:length(N_monocultureSets_split)) {
    cat("\n", paste0("LV set W", i-1, " Richness=",  monocultureSets_richness$Richness[i]))
    # init_N0
    read_csv(paste0(input_poolPairs$output_dir[i], "poolPairs_W", i-1, "-1-N_init.csv"), show_col_types = F) %>%
        write_csv(paste0(input_LVPairs$output_dir[i], "LVPairs_W", i-1, "-1-N_init.csv"))
    # Update n_wells in the input files
    input_LVPairs$n_wells[i] <- choose(10, 2) * 3
    # init_R0
    read_csv(paste0(input_poolPairs$output_dir[i], "poolPairs_W", i-1, "-1-R_init.csv"), show_col_types = F) %>%
        write_csv(paste0(input_LVPairs$output_dir[i], "LVPairs_W", i-1, "-1-R_init.csv"))
}

}

write_csv(input_LVPairs, here::here("simulation/03c-input_LVPairs.csv"))
#write_csv(monocultureSets_richness, paste0(folder_simulation, "aggregated/03-monocultureSets_richness.csv"))
#write_csv(monocultureSets_species, paste0(folder_simulation, "aggregated/03-monocultureSets_species.csv"))



