library(tidyverse)

# Generate the input_csv files ----
output_dir = "~/Dropbox/lab/invasion-network/simulation/data/raw12/"

# Example parameters
input_parameters <- tibble(
    output_dir = output_dir,
    save_timepoint = T,
    init_N0 = "init_pairs.csv",
    init_R0 = NA,
    exp_id = 1,
    seed = 1,
    sa = 500,
    ma = 20,
    S = 50,
    # c matrix
    c_symmetry = "empirical", # "symmetry", "asymmetry", "empirical"
    c_fs = 0.0415,  # mean uptake rates of fermenter on sugar
    sigc_fs = 0.0162, # standard deviation of uptake rates of fermenter on sugar
    c_fa = 0.0106,  # mean uptake rates of fermenter on acid
    sigc_fa = 0.0146, # standard deviation of uptake rates of fermenter on acid
    c_rs = 0.00911,  # mean uptake rates of respirator on sugar
    sigc_rs = 0.0102, # standard deviation of uptake rates of respirator on sugar
    c_ra = 0.0201,  # mean uptake rates of respirator on acid
    sigc_ra = 0.0146, # standard deviation of uptake rates of respirator on acid
    # D matrix
    metabolism = "empirical", # "common", "two-families", "empirical", "specific"
    ffss = 0, # fraction of flux from sugar to sugar in fermenter
    ffsa = 1, # fraction of flux from sugar to acid in fermenter
    ffas = 0, # fraction of flux from acid to sugar in fermenter
    ffaa = 1, # fraction of flux from acid to acid in fermenter
    frss = 0.49, # fraction of flux from sugar to sugar in respirator
    frsa = 0.51, # fraction of flux from sugar to acid in respirator
    fras = 0, # fraction of flux from acid to sugar in respirator
    fraa = 1, # fraction of flux from acid to acid in respirator
    # Leakiness
    l1 = 0.432,
    l2 = 0.00297,
    l1_sd = 0.105,
    l2_sd = 0.00252,
    n_communities = 20,
    n_wells = 100, # Note that the well number (column number) of init_N0 has to match n_wells
    rs = 0,
)

write_csv(input_parameters, file = here::here("simulation/input_parameters.csv"))

# Single species
temp1 <- input_parameters %>%
    slice(rep(1, 1)) %>%
    mutate(init_N0 = paste0("monoculture-1_init.csv"), exp_id = 1, n_wells = sa*2) %>%
    mutate(save_timepoint = F)
# Self-assembly
temp2 <- input_parameters %>%
    slice(rep(1, 1)) %>%
    mutate(save_timepoint = T) %>%
    mutate(init_N0 = paste0("selfAssembly-1_init.csv"), exp_id = 1, n_wells = 20, S = 50)

#input_independent <- bind_rows(temp1, temp2)
input_independent <- bind_rows(temp2)
write_csv(input_independent, here::here("simulation/input_independent.csv"))


# Pairs from the pool
n_comm <- input_parameters$n_communities
temp3 <- input_parameters %>%
    slice(rep(1, n_comm)) %>%
    mutate(save_timepoint = F) %>%
    mutate(init_N0 = paste0("poolPairs_W", 0:(n_comm-1), "-1_init.csv"), exp_id = 1, S = 10)
# Pairs from self-assembled communities
temp4 <- input_parameters %>%
    slice(rep(1, n_comm)) %>%
    mutate(save_timepoint = F) %>%
    mutate(init_N0 = paste0("communityPairs_W", 0:(n_comm-1), "-1_init.csv"))

#input_pairs <- bind_rows(temp3, temp4)
input_pairs <- bind_rows(temp4)
write_csv(input_pairs, here::here("simulation/input_pairs.csv"))


# Generate initial composition ----
# Generate family-species and class-resource matching tibble
input_independent <- read_csv(here::here("simulation/input_independent.csv"))
# Note that input_independent has to use the same sa and ma throughout
sa <- input_independent$sa[1]
ma <- input_independent$ma[1]
sal <- tibble(Family = paste0("F", c(rep(0, sa), rep(1, sa))), Species = paste0("S", 0:(sa * 2 - 1)))
mal <- tibble(Class = paste0("T", c(rep(0, ma), rep(1, ma), rep(2, ma))), Resource = paste0("R", 0:(ma * 3 - 1)))

# Monoculture
input_row <- input_independent %>% filter(str_detect(init_N0, "monoculture"))
draw_monoculture <- function(input_row) {
    n_wells <- input_row$n_wells
    tibble(Well = paste0("W", 0:(n_wells-1)),
           Species = paste0("S", 0:(sa*2-1)),
           Abundance = rep(1, n_wells)) %>%
        full_join(sal, by = "Species") %>%
        # Fill the non-chosen species
        replace_na(list(Well = "W0", Abundance = 0)) %>%
        pivot_wider(id_cols = c(Family, Species), names_from = Well, values_from = Abundance, values_fill = 0) %>%
        mutate(Species = factor(Species, sal$Species)) %>%
        arrange(Species)

}
N_mono <- draw_monoculture(input_row)
write_csv(N_mono, file = paste0(output_dir, "monoculture-1_init.csv"))

# Self-assembly. Draw a fix number of S
input_row <- input_independent %>% filter(str_detect(init_N0, "selfAssembly"))
draw_community <- function(input_row, candidate_species = NA) {
    S <- input_row$S
    n_communities <- input_row$n_communities
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
set.seed(1)
N_community <- draw_community(input_row)
write_csv(N_community, file = paste0(output_dir, "selfAssembly-1_init.csv"))












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
sp_mono <- read_csv(paste0(output_dir, "monoculture-1_end.csv")) %>%
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
    write_csv(paste0(output_dir, "poolNetwork-1_init.csv"))
for (i in 1:length(temp)) {
    print(paste0("Pool set W", i-1, " Richness=", nrow(temp[[i]])))
    draw_pairs_from_community(temp[[i]]) %>%
        write_csv(paste0(output_dir, "poolPairs_W", i-1, "-1_init.csv"))
}








# Execute this chunk when community assembly is done
input_row <- input_pairs %>% filter(str_detect(init_N0, "communityPairs")) %>% slice(1)

# Community pairs
N_community_end <- read_csv(paste0(output_dir, "selfAssembly-1_end.csv"), col_types = cols()) %>% rename(Family = ...1, Species = ...2)
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
    draw_pairs_from_community(temp[[i]]) %>% write_csv(paste0(output_dir, "communityPairs_W", i-1, "-1_init.csv"))
}










