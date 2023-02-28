#' This script generates required csv files for simlultion, including
#'
#' 1. the mapping files specifying the simulation parameters
#' 2. the initial composition of entry communities, pairs, and monocultures

library(tidyverse)
source(here::here("analysis/00-metadata.R"))

# Example parameters
input_parameters <- tibble(
    output_dir = NA,
    save_timepoint = T,
    init_N0 = "init_pairs.csv",
    init_R0 = "init_R0.csv",
    exp_id = 1,
    seed = 1,
    # Sampling pool
    sa = 500, # number of species in each specialist family
    ma = 20, # number of resources in each class
    S = 50, # number of species in the initial composition of each self-assembled community
    R0_food = 10000, # supplied R0 amount
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

# A universal parameter setting
write_csv(input_parameters, here::here("simulation/01-input_parameters.csv"))

# Single species, or monocultures
input_parameters %>%
    slice(rep(1, 1)) %>%
    mutate(save_timepoint = T) %>%
    mutate(output_dir = paste0(folder_simulation, "01-monocultures/")) %>%
    mutate(init_N0 = paste0("monoculture-1-N_init.csv"), exp_id = 1, n_wells = 200) %>%
    mutate(init_R0 = paste0("monoculture-1-R_init.csv")) %>%
    write_csv(here::here("simulation/01a-input_monocultures.csv"))

# Self-assembly, or communities
input_parameters %>%
    slice(rep(1, 1)) %>%
    mutate(save_timepoint = T) %>%
    mutate(output_dir = paste0(folder_simulation, "02-communities/")) %>%
    mutate(init_N0 = paste0("selfAssembly-1-N_init.csv"), exp_id = 1, n_wells = 20, S = 50) %>%
    mutate(init_R0 = paste0("selfAssembly-1-R_init.csv")) %>%
    write_csv(here::here("simulation/01b-input_communities.csv"))

# Pairs from the pool
n_comm <- input_parameters$n_communities
input_parameters %>%
    slice(rep(1, n_comm)) %>%
    mutate(save_timepoint = F) %>%
    mutate(output_dir = paste0(folder_simulation, "03-poolPairs/")) %>%
    mutate(init_N0 = paste0("poolPairs_W", 0:(n_comm-1), "-1-N_init.csv"), exp_id = 1, S = 10) %>%
    mutate(init_R0 = paste0("poolPairs_W", 0:(n_comm-1), "-1-R_init.csv")) %>%
    write_csv(here::here("simulation/01c-input_poolPairs.csv"))

# Pairs from within self-assembled communities
input_parameters %>%
    slice(rep(1, n_comm)) %>%
    mutate(save_timepoint = F) %>%
    mutate(output_dir = paste0(folder_simulation, "04-withinCommunityPairs/")) %>%
    mutate(init_N0 = paste0("withinCommunityPairs_W", 0:(n_comm-1), "-1-N_init.csv")) %>%
    mutate(init_R0 = paste0("withinCommunityPairs_W", 0:(n_comm-1), "-1-R_init.csv")) %>%
    write_csv(here::here("simulation/01d-input_withinCommunityPairs.csv"))











