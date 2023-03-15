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
    fa = 3,                       # Number of specialist families. Used to compute 'SA': 60*np.ones(3)
    sa = 500,                     # Number of species in each specialist family. Used to compute 'SA': 60*np.ones(3)
    ma = 10,                      # Number of resources in each class. Used to compute 'MA': 30*np.ones(3)
    Sgen = 0,                     # CM parameter. Number of species in the generalist family
    # Sampling parameters
    sampling = "Gamma",           # CM parameter. Sampling approach for consumer parameters: "Binary", "Gamma", "Gaussian"
    muc = 10,                     # CM parameter. Mean sum of consumption rates (used in all models)
    sigc = 4,                     # CM parameter. Standard deviation of sum of consumption rates for Gaussian and Gamma models
    q = 0.5,                      # CM parameter. Preference strength of specialist families (0 for generalist and 1 for specialist)
    c0 = 0,                       # CM parameter. Sum of background consumption rates in binary model
    c1 = 1,                       # CM parameter. Specific consumption rate in binary model
    fs = 0.45,                    # CM parameter. Fraction of secretion flux with same resource type
    fw = 0.45,                    # CM parameter. Fraction of secretion flux to 'waste' resource
    sparsity = 0.2,               # CM parameter. Effective sparsity of metabolic matrix (between 0 and 1)
    food = 0,                     # CM parameter. Index of food source (when a single resource is supplied externally)
    R0_food = 1000,               # CM parameter. Unperturbed fixed point for supplied food
    regulation = "independent",   # CM parameter. Metabolic regulation (see dRdt)
    response = "type I",          # CM parameter. Functional response (see dRdt)
    supply = "off",               # CM parameter. Resource supply (see dRdt) 'off' for batch culture. 'external' and 'self-renewing' for constant supply
    # CM internal parameters
    m = 0,                        # Set m=0 to turn off maintenance cost. i.e., no cell dies
    w = 1,                        # CM parameter. The resource use efficiency
    g = 1,
    l = 0.5,                      # CM parameter. Leakage fraction
    tau = 1,
    r = 1,
    sigma_max = 1,
    nreg = 10,
    n = 2,
    # Experiments
    n_pass = 20,                  # Number of transfers
    t_propagation = 1,            # Length of propagation in one transfer
    save_timepoint = FALSE,
    n_timepoint = 50,
    dilution_factor = 1/1000,      # Dilution factor for passage
    n_wells = 50,                 # CM parameter. Number of independent wells
    #n_wells = 20,                # number of monocultures tested
    n_communities = 20,           # number of communities used
    S = 100                       # CM parameter. Number of species per well (randomly sampled from the pool of size Stot = sum(SA) + Sgen)
)

write_csv(input_parameters, here::here("simulation/01-input_parameters.csv"))


# Generate family-species and class-resource table for matching
fa <- input_parameters$fa[1]
sa <- input_parameters$sa[1]
ma <- input_parameters$ma[1]
sal <- tibble(Family = paste0("F", rep(c(0:(fa-1)), each = sa)), Species = paste0("S", 0:(sa * fa - 1))) %>% mutate(SpeciesID = 1:n())
mal <- tibble(Class = paste0("T", rep(c(0:(fa-1)), each = ma)), Resource = paste0("R", 0:(ma * fa - 1))) %>% mutate(ResourceID = 1:n())


