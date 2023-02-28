#' This script generates the entry composition for pairwise competition
#'
#' Because these simulations depends on other simulations (monoculture and pairs),
#' this script can only be run after those simulations are done
#'
#' Note that the n_wells in input_withinCommunityPairs and input_poolPairs is modified after
#' knowing the number of species in the communities and pool

library(tidyverse)
source(here::here("analysis/00-metadata.R"))
source(here::here("simulation/02-generate_initial_composition.R"))

# 0. Read parameters ----
input_parameters <- read_csv(here::here("simulation/01-input_parameters.csv"), col_types = cols())
input_monocultures <- read_csv(here::here("simulation/01a-input_monocultures.csv"), col_types = cols())
input_communities <- read_csv(here::here("simulation/01b-input_communities.csv"), col_types = cols())
input_poolPairs <- read_csv(here::here("simulation/01c-input_poolPairs.csv"), col_types = cols())
input_withinCommunityPairs <- read_csv(here::here("simulation/01d-input_withinCommunityPairs.csv"), col_types = cols())


# 1. Generate competing pairs from culturable monocultures ----
# Execute this chunk when monoculture is done
draw_pairs_from_community <- function(N_community_long) {
    # Communities with no or only one species
    if (nrow(N_community_long) <= 1) N_pairs <- sal %>% mutate(W0 = 0)

    # Other communities with more one species
    if (nrow(N_community_long) > 1) {
        N_pairs <- expand_grid(sp1 = N_community_long$Species, sp2 = N_community_long$Species) %>%
            mutate(across(everything(), ~ ordered(.x, sal$Species))) %>%
            filter(sp1 < sp2) %>%
            mutate(Pair = paste0("P", 0:(n()-1))) %>%
            # Make duplicate for two freq
            slice(rep(1:n(), each = 2)) %>%
            mutate(Well = paste0("W", 0:(n()-1))) %>%
            group_by(Well) %>%
            pivot_longer(cols = starts_with("sp"), names_to = "temp", values_to = "Species") %>%
            ungroup() %>%
            mutate(Abundance = rep(c(0.05, 0.95, 0.95, 0.05), n()/4)) %>%
            full_join(sal, by = "Species") %>%
            replace_na(list(Well = "W0", Abundance = 0)) %>%
            select(Pair, Well, Family, Species, Abundance) %>%
            pivot_wider(id_cols = c(Family, Species), names_from = Well, values_from = Abundance, values_fill = 0) %>%
            mutate(Species = factor(Species, sal$Species)) %>%
            arrange(Species)
    }
    return(N_pairs)
}
# Pool pairs. Use isolates that can grow in monoculture
set.seed(1)
species_mono <- read_csv(paste0(input_monocultures$output_dir[1], "monoculture-1-N_end.csv"), col_types = cols()) %>%
    rename(Family = ...1, Species = ...2) %>%
    mutate_all(~replace(., .==0, NA)) %>%
    pivot_longer(cols = starts_with("W"), names_to = "Well", values_to = "Abundance", values_drop_na = T) %>%
    pull(Species) %>%
    str_replace("S", "") %>% as.numeric()

N_network <- draw_community(input_poolPairs[1,], species_mono)

N_network_split <- N_network %>%
    mutate_all(~replace(., .==0, NA)) %>%
    pivot_longer(cols = -c(Family, Species), names_to = "Community", values_to = "Abundance", values_drop_na = T) %>%
    # Order species and community
    mutate(Species = ordered(Species, sal$Species)) %>%
    mutate(Community = ordered(Community, colnames(N_network))) %>%
    arrange(Community, Species) %>%
    group_by(Community) %>%
    group_split()



for (i in 1:length(N_network_split)) {
    print(paste0("Pool set W", i-1, " Richness=", nrow(N_network_split[[i]])))
    # init_N0
    draw_pairs_from_community(N_network_split[[i]]) %>%
        write_csv(paste0(input_poolPairs$output_dir[i], "poolPairs_W", i-1, "-1-N_init.csv"))
    # Update n_wells in the input files
    input_poolPairs$n_wells[i] <- choose(nrow(N_network_split[[i]]), 2) * 2
    # init_R0
    set_community_resource(input_poolPairs[i,]) %>%
        write_csv(paste0(input_poolPairs$output_dir[i], "poolPairs_W", i-1, "-1-R_init.csv"))
}


write_csv(input_poolPairs, here::here("simulation/01c-input_poolPairs.csv"))


# 2. Generate within community pairs ----
# Execute this chunk when community assembly is done
N_community <- read_csv(paste0(input_communities$output_dir[1], "selfAssembly-1-N_end.csv"), col_types = cols()) %>% rename(Family = ...1, Species = ...2)
N_community_split <- N_community %>%
    mutate_all(~replace(., .==0, NA)) %>%
    pivot_longer(cols = -c(Family, Species), names_to = "Community", values_to = "Abundance", values_drop_na = T) %>%
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

for (i in 1:length(N_community_split)) {
    print(paste0("Community W", i-1, " Richness=", nrow(N_community_split[[i]])))
    # init_N0
    draw_pairs_from_community(N_community_split[[i]]) %>%
        write_csv(paste0(input_withinCommunityPairs$output_dir[i], "withinCommunityPairs_W", i-1, "-1-N_init.csv"))
    # Update n_wells in the input files
    input_withinCommunityPairs$n_wells[i] <- choose(nrow(N_community_split[[i]]), 2) * 2
    # init_R0
    set_community_resource(input_withinCommunityPairs[i,]) %>%
        write_csv(paste0(input_withinCommunityPairs$output_dir[i], "withinCommunityPairs_W", i-1, "-1-R_init.csv"))

}
write_csv(input_withinCommunityPairs, here::here("simulation/01d-input_withinCommunityPairs.csv"))



if (FALSE) {
    # Execute this chunk when monoculture is done
    draw_pairs_from_community <- function(N_community_long) {
        # Communities with no or only one species
        if (nrow(N_community_long) <= 1) N_pairs <- sal %>% mutate(W0 = 0)

        # Other communities with more one species
        if (nrow(N_community_long) > 1) {
            N_pairs <- expand_grid(sp1 = N_community_long$Species, sp2 = N_community_long$Species) %>%
                mutate(across(everything(), ~ ordered(.x, sal$Species))) %>%
                filter(sp1 < sp2) %>%
                mutate(Pair = paste0("P", 0:(n()-1))) %>%
                # Make duplicate for two freq
                slice(rep(1:n(), each = 2)) %>%
                mutate(Well = paste0("W", 0:(n()-1))) %>%
                group_by(Well) %>%
                pivot_longer(cols = starts_with("sp"), names_to = "temp", values_to = "Species") %>%
                ungroup() %>%
                mutate(Abundance = rep(c(0.05, 0.95, 0.95, 0.05), n()/4)) %>%
                full_join(sal, by = "Species") %>%
                replace_na(list(Well = "W0", Abundance = 0)) %>%
                select(Pair, Well, Family, Species, Abundance) %>%
                pivot_wider(id_cols = c(Family, Species), names_from = Well, values_from = Abundance, values_fill = 0) %>%
                mutate(Species = factor(Species, sal$Species)) %>%
                arrange(Species)
        }
        return(N_pairs)
    }
    # Pool pairs. Use isolates that can grow in monoculture
    ## Use isolates that can grow in monoculture
    set.seed(1)
    input_row <- input_pairs %>% filter(str_detect(init_N0, "poolPairs")) %>% slice(1)
    sp_mono <- read_csv(paste0(input_parameters$output_dir, "monoculture-1_end.csv")) %>%
        rename(Family = ...1, Species = ...2) %>%
        mutate_all(~replace(., .==0, NA)) %>%
        pivot_longer(cols = starts_with("W"), names_to = "Well", values_to = "Abundance", values_drop_na = T) %>%
        pull(Species) %>%
        str_replace("S", "") %>% as.numeric()
    N_network <- draw_community(input_row, sp_mono)
    temp <- N_network %>%
        mutate_all(~replace(., .==0, NA)) %>%
        pivot_longer(cols = -c(Family, Species), names_to = "Community", values_to = "Abundance", values_drop_na = T) %>%
        # Order species and community
        mutate(Species = ordered(Species, sal$Species)) %>%
        mutate(Community = ordered(Community, colnames(N_network))) %>%
        arrange(Community, Species) %>%
        group_by(Community) %>%
        group_split()
    ## Save the random network composition
    temp %>%
        bind_rows() %>%
        write_csv(paste0(input_parameters$output_dir, "poolNetwork-1_init.csv"))
    for (i in 1:length(temp)) {
        print(paste0("Pool set W", i-1, " Richness=", nrow(temp[[i]])))
        draw_pairs_from_community(temp[[i]]) %>%
            write_csv(paste0(input_parameters$output_dir, "poolPairs_W", i-1, "-1_init.csv"))
    }








    # Execute this chunk when community assembly is done
    input_row <- input_pairs %>% filter(str_detect(init_N0, "communityPairs")) %>% slice(1)

    # Community pairs
    N_community_end <- read_csv(paste0(input_parameters$output_dir, "selfAssembly-1_end.csv"), col_types = cols()) %>% rename(Family = ...1, Species = ...2)
    temp <- N_community_end %>%
        mutate_all(~replace(., .==0, NA)) %>%
        pivot_longer(cols = -c(Family, Species), names_to = "Community", values_to = "Abundance", values_drop_na = T) %>%
        # Remove rare species (relative abundance <0.01)
        group_by(Community) %>%
        mutate(RelativeAbundance = Abundance/sum(Abundance)) %>%
        filter(RelativeAbundance > 0.01) %>%
        # Order species and community
        mutate(Species = ordered(Species, sal$Species)) %>%
        mutate(Community = ordered(Community, colnames(N_community_end))) %>%
        arrange(Community, Species) %>%
        group_by(Community) %>%
        group_split()

    for (i in 1:length(temp)) {
        print(paste0("Community W", i-1, " Richness=", nrow(temp[[i]])))
        draw_pairs_from_community(temp[[i]]) %>% write_csv(paste0(input_parameters$output_dir, "communityPairs_W", i-1, "-1_init.csv"))
    }


}








