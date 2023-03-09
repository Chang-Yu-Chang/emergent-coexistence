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
    output_dir = folder_simulation,
    save_timepoint = T,
    init_N0 = "init_pairs.csv",
    init_R0 = "init_R0.csv",
    exp_id = 1,
    seed = 1,
    # Testing community-simulator parameters
    w = 1,
    l = 0.8,
    # Passaging
    n_pass = 20,            # number of transfer or passages
    t_propagation = 10,     # time length of propagation, or length of growth cycle
    # Sampling pool
    sa = 500,               # number of species in each specialist family
    ma = 20,                # number of resources in each class
    S = 50,                 # number of species in the initial composition of each self-assembled community
    R0_food = 1000,         # supplied R0 amount
    n_communities = 20,     # number of communities used
    n_wells = 200,          # number of monocultures tested
    sampling = "empirical", # sampling approach for consumer parameters "empirical", "binary", "gamma"
    # c matrix
    c_fs = 0.183,           # mean uptake rates of fermenter on sugar
    sigc_fs = 0.0478,       # standard deviation of uptake rates of fermenter on sugar
    c_fa = 0.192,           # mean uptake rates of fermenter on acid
    sigc_fa = 0.0407,       # standard deviation of uptake rates of fermenter on acid
    c_rs = 0.0268,          # mean uptake rates of respirator on sugar
    sigc_rs = 0.0314,       # standard deviation of uptake rates of respirator on sugar
    c_ra = 0.236,           # mean uptake rates of respirator on acid
    sigc_ra = 0.0828,       # standard deviation of uptake rates of respirator on acid
    # D matrix
    ffss = 0,               # fraction of flux from sugar to sugar in fermenter
    ffsa = 1,               # fraction of flux from sugar to acid in fermenter
    ffas = 0,               # fraction of flux from acid to sugar in fermenter
    ffaa = 1,               # fraction of flux from acid to acid in fermenter
    frss = 0.487,           # fraction of flux from sugar to sugar in respirator
    frsa = 0.513,           # fraction of flux from sugar to acid in respirator
    fras = 0,               # fraction of flux from acid to sugar in respirator
    fraa = 1,               # fraction of flux from acid to acid in respirator
    # l matrix
    l1 = 0.442,             # mean leakiness of fermenter
    l1_sd = 0.108,          # sd leakiness of fermenter
    l2 = 0.00241,           # mean leakiness of respirator
    l2_sd = 0.00224         # sd leakiness of fermenter
)
    # Test
    # mutate(
    #     # c matrix
    #     c_fs = 1,  # mean uptake rates of fermenter on sugar
    #     sigc_fs = 0, # standard deviation of uptake rates of fermenter on sugar
    #     c_fa = 1,  # mean uptake rates of fermenter on acid
    #     sigc_fa = 0, # standard deviation of uptake rates of fermenter on acid
    #     c_rs = 1,  # mean uptake rates of respirator on sugar
    #     sigc_rs = 0, # standard deviation of uptake rates of respirator on sugar
    #     c_ra = 1,  # mean uptake rates of respirator on acid
    #     sigc_ra = 0, # standard deviation of uptake rates of respirator on acid
    #     # D matrix
    #     metabolism = "empirical", # "common", "two-families", "empirical", "specific"
    #     ffss = 0, # fraction of flux from sugar to sugar in fermenter
    #     ffsa = 0, # fraction of flux from sugar to acid in fermenter
    #     ffas = 0, # fraction of flux from acid to sugar in fermenter
    #     ffaa = 0, # fraction of flux from acid to acid in fermenter
    #     frss = 0, # fraction of flux from sugar to sugar in respirator
    #     frsa = 0, # fraction of flux from sugar to acid in respirator
    #     fras = 0, # fraction of flux from acid to sugar in respirator
    #     fraa = 0, # fraction of flux from acid to acid in respirator
    #     # l matrix
    #     l1 = 0,
    #     l1_sd = 0,
    #     l2 = 0,
    #     l2_sd = 0
    # )
    #
write_csv(input_parameters, here::here("simulation/01-input_parameters.csv"))


# Generate family-species and class-resource table for matching
sa <- input_parameters$sa[1]
ma <- input_parameters$ma[1]
sal <- tibble(Family = paste0("F", c(rep(0, sa), rep(1, sa))), Species = paste0("S", 0:(sa * 2 - 1)))
mal <- tibble(Class = paste0("T", c(rep(0, ma), rep(1, ma))), Resource = paste0("R", 0:(ma * 2 - 1)))









# c_fs = 0.0415,  # mean uptake rates of fermenter on sugar
# sigc_fs = 0.0162, # standard deviation of uptake rates of fermenter on sugar
# c_fa = 0.0106,  # mean uptake rates of fermenter on acid
# sigc_fa = 0.0146, # standard deviation of uptake rates of fermenter on acid
# c_rs = 0.00911,  # mean uptake rates of respirator on sugar
# sigc_rs = 0.0102, # standard deviation of uptake rates of respirator on sugar
# c_ra = 0.0201,  # mean uptake rates of respirator on acid
# sigc_ra = 0.0146, # standard deviation of uptake rates of respirator on acid
# l1 = 0.432,
# l1_sd = 0.105,
# l2 = 0.00297,
# l2_sd = 0.00252,
