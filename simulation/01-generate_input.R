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
    output_dir = NA,
    save_timepoint = T,
    init_N0 = "init_pairs.csv",
    init_R0 = "init_R0.csv",
    exp_id = 1,
    seed = 1,
    # Sampling pool
    sa = 500, # number of species in each specialist family
    ma = 10, # number of resources in each class
    S = 100, # number of species in the initial composition of each self-assembled community
    R0_food = 1000, # supplied R0 amount
    # c matrix
    c_symmetry = "empirical", # "symmetry", "asymmetry", "empirical"
    # c_fs = 0.0415,  # mean uptake rates of fermenter on sugar
    # sigc_fs = 0.0162, # standard deviation of uptake rates of fermenter on sugar
    # c_fa = 0.0106,  # mean uptake rates of fermenter on acid
    # sigc_fa = 0.0146, # standard deviation of uptake rates of fermenter on acid
    # c_rs = 0.00911,  # mean uptake rates of respirator on sugar
    # sigc_rs = 0.0102, # standard deviation of uptake rates of respirator on sugar
    # c_ra = 0.0201,  # mean uptake rates of respirator on acid
    # sigc_ra = 0.0146, # standard deviation of uptake rates of respirator on acid
    c_fs = 0.218,  # mean uptake rates of fermenter on sugar
    sigc_fs = 0.0444, # standard deviation of uptake rates of fermenter on sugar
    c_fa = 0.21,  # mean uptake rates of fermenter on acid
    sigc_fa = 0.0357, # standard deviation of uptake rates of fermenter on acid
    c_rs = 0.07,  # mean uptake rates of respirator on sugar
    sigc_rs = 0.0341, # standard deviation of uptake rates of respirator on sugar
    c_ra = 0.31,  # mean uptake rates of respirator on acid
    sigc_ra = 0.115, # standard deviation of uptake rates of respirator on acid
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
    # l1 = 0.432,
    # l1_sd = 0.105,
    # l2 = 0.00297,
    # l2_sd = 0.00252,
    l1 = 0.442,
    l1_sd = 0.108,
    l2 = 0.00241,
    l2_sd = 0.00224,
    n_communities = 20,
    n_wells = 100, # Note that the well number (column number) of init_N0 has to match n_wells
    rs = 0,
)

# A universal parameter setting
write_csv(input_parameters, here::here("simulation/01-input_parameters.csv"))










