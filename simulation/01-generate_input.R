#' This script generates the universal parameter mapping file 01-input_parameters.csv
#'
#' The 01-input_parameters.csv will later be read in
#' 02-generate_input_mono_comm.R and 03-generate_input_pairs.R
#' to generate specific input csv files for monoculture, communities,
#' and pairs simulation

library(tidyverse)
source(here::here("analysis/00-metadata.R"))

# Example parameters
input_parameters <- tibble(
    # File and folder names
    output_dir = folder_simulation,
    init_N0 = "init_pairs.csv",
    init_R0 = "init_R0.csv",
    exp_id = 1,
    seed = 1,
    # Species pool
    fa = 10,                       # Number of specialist families. Used to compute 'SA': 60*np.ones(3)
    sa = 100,                     # Number of species in each specialist family. Used to compute 'SA': 60*np.ones(3)
    ma = 1,                       # Number of resources in each resource type. Used to compute 'MA': 30*np.ones(3)
    Sgen = 0,                     # CM parameter. Number of species in the generalist family
    # Sampling parameters
    sampling_c = "Dirichlet",     # modified CM parameter. Sampling approach for consumer uptake rate: 'Binary', 'Gamma', 'Gaussian' (CM default), and 'Dirichlet' (Goldford et al 2018)
    sampling_D = "Uniform",       # Sampling approach for stoichiometric matrix: 'Dirichlet' (CM default) and 'Uniform' (Goldford et al 2018), and 'Zeros' to turn off cross-feeding
    mu_f = 0.4,                   # Goldford et al 2018 parameter. Mean family level concentration parameter for gaussian
    sigma_f = 0.01,               # Goldford et al 2018 parameter. SD family level concentration parameter for gaussian
    omega_f = 100,                # Goldford et al 2018 parameter. Degree of total variability with each family. A high omega ensure that species are very similar
    #
    muc = 10,                     # CM parameter. Mean sum of consumption rates, or consumption capacity (used in all models)
    sigc = 4,                     # CM parameter. Standard deviation of sum of consumption rates for Gaussian and Gamma models
    q = 0.5,                      # CM parameter. Preference strength of specialist families (0 for generalist and 1 for specialist)
    c0 = 0,                       # CM parameter. Sum of background consumption rates in binary model
    c1 = 1,                       # CM parameter. Specific consumption rate in binary model
    fs = 0.45,                    # CM parameter. Fraction of secretion flux with same resource type
    fw = 0.45,                    # CM parameter. Fraction of secretion flux to 'waste' resource
    sparsity = 0.2,               # CM parameter. Effective sparsity of metabolic matrix (between 0 and 1)
    food = 0,                     # CM parameter. Index of food source (when a single resource is supplied externally)
    R0_food = 10^3,               # CM parameter. Unperturbed fixed point for supplied food
    regulation = "independent",   # CM parameter. Metabolic regulation (see dRdt)
    response = "type I",          # CM parameter. Functional response (see dRdt)
    supply = "external",          # CM parameter. Resource supply (see dRdt) 'off' for batch culture. 'external' and 'self-renewing' for constant supply
    # assumptions = dictionary of metaparameters
    # 'sampling' = {'Gaussian','Binary','Gamma'} specifies choice of sampling algorithm
    # 'SA' = number of species in each family
    # 'MA' = number of resources of each type
    # 'Sgen' = number of generalist species
    # 'muc' = mean sum of consumption rates
    # 'sigc' = standard deviation for Gaussian sampling of consumer matrix
    # 'q' = family preference strength (from 0 to 1)
    # 'c0' = row sum of background consumption rates for Binary sampling
    # 'c1' = specific consumption rate for Binary sampling
    # 'fs' = fraction of secretion flux into same resource type
    # 'fw' = fraction of secretion flux into waste resource type
    # 'sparsity' = effective sparsity of metabolic matrix (from 0 to 1)
    # 'waste_type' = index of resource type to designate as "waste"
    # CM implicit parameters
    g = 1,                        # CM parameter. Conversion factor from energy uptake to growth rate (1/energy)
    w = 1,                        # CM parameter. Energy content of resource (energy/mass)
    l = 0,                        # CM parameter. Leakage fraction
    m = 1,                        # CM parameter. Minimal energy uptake for maintenance of species (energy/time). Set m=0 to turn off maintenance cost. i.e., no cell dies
    tau = 1,                      # CM parameter. Resource replenishment rate. Useful when 'supply' is set to 'external
    r = 1,                        # CM parameter. Rate of resource self-renewal (volume/mass/time). Useful when 'supply' is set to 'self-renewing'
    n = 2,                        # CM parameter. Hill coefficient for functional response (unitless). Useful when 'response' is 'type III'
    sigma_max = 1,                # CM parameter. Maximum input flux (mass/time). Useful when 'response' is 'type II' or 'type III'
    nreg = 10,                    # CM parameter. Hill coefficient for metabolic regulation (unitless). Useful when 'regulation' is 'energy' or 'mass'
    # Experiments
    n_timesteps = 10^6,           # Goldford et al 2018 parameter. Number of timesteps for ODE
    n_timepoints = 10,            # Goldford et al 2018 parameter. Number of timepoints for outputing the abundance for plotting
    save_timepoint = FALSE,
    #
    n_pass = 20,                  # Number of transfers
    t_propagation = 1,            # Length of propagation in one transfer
    dilution_factor = 1/1000,      # Dilution factor for passage
    n_wells = 50,                 # CM parameter. Number of independent wells
    #n_wells = 20,                # number of monocultures tested
    n_communities = 20,           # number of communities used
    S = 100                       # CM parameter. Number of species per well (randomly sampled from the pool of size Stot = sum(SA) + Sgen)
)

write_csv(input_parameters, here::here("simulation/01-input_parameters.csv"))

# Generate family-species and class-resource table for matching
fa <- input_parameters$fa[1] # Number of families
sa <- input_parameters$sa[1] # Number of species in each species family
ma <- input_parameters$ma[1] # Number of resources in each resource type
sal <- tibble(Family = paste0("F", rep(c(0:(fa-1)), each = sa)), Species = paste0("S", 0:(sa * fa - 1))) %>% mutate(SpeciesID = 1:n())
mal <- tibble(Class = paste0("T", rep(c(0:(fa-1)), each = ma)), Resource = paste0("R", 0:(ma * fa - 1))) %>% mutate(ResourceID = 1:n())



# Test
input_parameters %>%
    slice(rep(1, 1)) %>%
    mutate(save_timepoint = T) %>%
    mutate(output_dir = paste0(folder_simulation, "test/")) %>%
    mutate(exp_id = 1, n_wells = 20) %>%
    mutate(init_N0 = paste0("test-1-N_init.csv")) %>%
    mutate(init_R0 = paste0("test-1-R_init.csv")) %>%
    write_csv(here::here("simulation/01-input_test.csv"))


input_test <- read_csv(here::here("simulation/01-input_test.csv"), col_types = cols())

draw_community <- function(input_communities, candidate_species = NA) {
    S <- input_communities$S
    n_communities <- input_communities$n_communities
    if (all(is.na(candidate_species))) {
        # Draw species. Default using a fix number of S. Identical to community-simulator
        species <- rep(list(NA), n_communities)
        for (i in 1:n_communities) species[[i]] <- sample(0:(fa*sa-1), size = S, replace = FALSE) %>% sort
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
                Abundance = rep(1/S, n_communities * S)) %>%
        full_join(sal, by = "Species") %>%
        mutate(Species = factor(Species, sal$Species)) %>%
        # Fill the non-chosen species
        replace_na(list(Well = "W0", Abundance = 0)) %>%
        pivot_wider(id_cols = c(Family, Species), names_from = Well, values_from = Abundance, values_fill = 0) %>%
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
draw_community(input_test) %>%
    write_csv(paste0(input_test$output_dir, "test-1-N_init.csv"))

set_community_resource(input_test) %>%
    write_csv(paste0(input_test$output_dir, "test-1-R_init.csv"))

