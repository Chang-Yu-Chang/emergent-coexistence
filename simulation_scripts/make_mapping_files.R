#' Make the input mapping files

suppressWarnings(suppressMessages(library(tidyverse)))
suppressWarnings(suppressMessages(library(data.table)))

seeds = 1:20 # Random seed. Default 1:100
cat("\nSeeds = ", seeds)
cat("\nTotal seeds are ", seeds, "\n")
# data_directory = "../data/raw/simulation/"
data_directory = "/home/cc2553/project/invasion-network/data/"
mapping_file_directory = "../data/raw/simulation/mapping_files/"
make_input_csv <- function(...){
    args = list(...)

    # List of parameters
    df_default <-
        data.frame(
            stringsAsFactors = FALSE,

            selected_function = "f1_additive", #Function that is under selection
            protocol = "simple_screening", #protocol to implement
            seed = 1, #Seed for species poo l
            exp_id = paste("f1_additive", "simple_screening", 1, sep = "-"), # ID for simulation (will determine filenames
            overwrite_plate = NA, # If not NA, then must be a text file. Overwrite the initial plate with this composititon saved in this text file
            passage_overwrite_plate = F, # If overwrite_plate != NA, set TRUE if the overwrite_plate is at equilibrium and need an addititonal transfer
            output_dir = "data/", # Output directory. Default is a data subfolder
            save_function = T, # Save Function data
            save_composition = T, # Save Composition Data
            save_plate = F, #Save initial plate
            function_lograte = 1, #How often do you save the function in transfers
            composition_lograte = 20, #How often do you save the composition in transfers

            #Experiment Parameters (applies to for all protocols)

            scale = 1000000, # Number of cells when N_i = 1
            n_inoc = 1000000, # Number of cells sampled from the regional species at start
            rich_medium = T, #Whether to generate a rich medium sampled from a a random distribution or a minimal media with only a single resource
            monoculture = F, # whether to run simple screening with monoculture
            dilution = 0.001, #Dilution factor at every transfers
            n_wells = 96, # Number of welss on a plate
            n_propagation = 1, # Incubation time
            n_transfer = 40, #Number of Transfers total number of transfers
            n_transfer_selection = 20, #Number of tranfers implementing selection regime
            metacommunity_sampling = "Power", # {"Power", "Lognormal", "Default"} Sampling method for initial metacommunity
            power_alpha = NA, # Default = 0.01
            lognormal_mean = NA, # Default = 8
            lognormal_sd = NA, # Default = 8


            #Parameters for community function, #parameters that determine properties of function
            phi_distribution = "Norm", # {"Norm", "Uniform"}
            phi_mean = 0, #
            phi_sd = 1, # Standard deviation for drawing specifc speices/interaction function
            phi_lower = 0,
            phi_upper = 1,
            ruggedness = 0.8, # (1-ruggedness) percent of function are set to 0
            function_ratio = 1, # Scaling factor between species- and interaction-specific function variances
            binary_threshold = 1, #Threshold for binary functions
            g0 = 1, # The baseline conversion factor of biomass per energy
            cost_distribution = "Norm", # {"Norm", "Uniform"}
            cost_mean = 0, # Mean fraction of cost feeded into a gamma distribution. Suggested up to 0.05
            cost_sd = 0, # Sd fraction of cost feeded into a gamma distribution. cost_sd = 0 if cost_mean = 0, cost_sd= 0.01 if cost_mean >0
            cost_lower = 0, # Lower bound for cost if cost_distribution="Uniform"
            cost_upper = 1, # Upper bound for cost if cost_distribution="Uniform"
            invader_index =  2,
            invader_sampling = "Gamma",
            invader_strength = 10,
            target_resource = NA, # Target resource production when selected_function=f6_target_resourece

            #Paramaters for Directed Selection (for directed selection protocols that can't be coded up in experiment paramaters)

            directed_selection = F, # If true whenever select_top is selected the highest performing Community is propagated asexually and some kind of pertubations can ne applied
            knock_out = F, #If True performs knock out pertubations
            knock_in = F, #If True performs knock in pertubation
            knock_in_threshold = NA, # value determines threshold for isolates to knock in, #If NA isolates are chosen at random.
            bottleneck = F, #If True perform bottleneck pertubations
            bottleneck_size = NA, #Magnitude of bottleneck. If not set it default to dilution
            migration = F, #If true perform migration pertubations
            n_migration = 1e6, # Number of cells to migration in the directed selection
            s_migration = NA, # Number of species to migrate. If s_migration is NA defaults to power law migration (so this is normal).
            coalescence = F, #If true perform coalescence pertubation
            frac_coalescence = NA, # fraction of coalesced community that is champion. Defaults to 0.5 if NA
            resource_shift = F, #If true performs resource pertubations
            r_type = NA, # Type of resource pertubation. rescale_add, rescale_remove, add, remove, old. if NA defaults to resource swap
            r_percent = NA, # Tunes the magnitude of resource pertubation if NA does not perform resource pertubation

            #Paramaters for community simulator package, note that we have split up a couple of paramaters that are inputed as list (SA and SGen). In the mapping file
            #if paramater is set as NA it takes the default value in community_simulator package. Also some paramaters could actually be inputed as lists but this is beyond the scope of this structure of mapping file i.e m, w,g r

            sampling = "Binary_Gamma", #{'Gaussian','Binary','Gamma', 'Binary_Gamma'} specifies choice of sampling algorithm
            sn = 2100, #number of species per specialist family
            sf = 1, #number of specialist families, # note SA = sn *np.ones(sf)
            Sgen = 0, #number of generalist species
            rn = 90, #number of resources per resource clas
            rf = 1, #number of resource classes, #Note RA = rn*np.ones(rf)
            R0_food = 1000, #Total amount of supplied food
            food = NA, #index of food source (when a single resource is supplied externally). ONly work for single-resource medium
            supply = NA, #resource supply (see dRdt)
            muc = NA, #mean sum of comnsumption rate
            sigc = NA, #Standard deviation of sum of consumption rates for Gaussian and Gamma models
            c0 = NA, #Sum of background consumption rates in binary model
            c1 = NA, #Specific consumption rate in binary model
            q = NA, #preference strength
            sparsity = NA, #Effective sparsity of metabolic matrix (between 0 and 1)
            fs = NA, #Fraction of secretion flux with same resource type
            fw = NA , #Fraction of secretion flux to 'waste' resource
            g = NA, #energy to biomass conversion
            w = NA, #resource energy value
            l = 0, # Leakage fraction
            m = 0, # Minimal resurce uptake; mortality
            n = NA, #hill coefficient when n= 1 response is type 2
            response = "type III", #functional response (see dRdt)
            sigma_max = NA, # default = 1. sigma max for functional response
            regulation = NA, #metabolic regulation (see dRdt)
            nreg = NA, #Hill coefficient that tunes steepness. of metabolic regulation.
            tau = NA, # external resource supply  rate (for chemostat)
            r = NA, #renewal rate for self renewing resources
            S = 100 # number of species in the initial community, legacy of the community-simulator

        )


    argument_names <- colnames(df_default)
    # If there is no change, output a line of default
    if (length(args) == 0) {
        output_row <- df_default
    } else { # If there is any variable specified, change that
        output_row <- df_default
        # The variable has to be in the argument list
        to_changed_args <- names(args)
        if (!all(to_changed_args %in% argument_names)) stop("Errors: arguments do not exist in the list")
        for (i in 1:length(to_changed_args)) output_row[,to_changed_args[i]] <- args[[i]][1]

        # Dependency
        ## Set exp_id names by the seed, selected function, and protocol (and directed selection type if protocol is directed selection)
        output_row$exp_id = paste(output_row$selected_function, output_row$protocol, output_row$seed, sep = "-")

        ## Monoculture
        if (output_row$monoculture == TRUE) {
            if (output_row$protocol != "simple_screening") stop("Errors: monoculture plate has to be simple_screening")
            #output_row$protocol = "monoculture"
            output_row$exp_id = paste(output_row$selected_function, "monoculture", output_row$seed, sep = "-")
        }

        ## Cost per function. Fixed sd of cost is specified
        if (output_row$cost_mean != 0) output_row$cost_sd <- 0.01

        ## Check on the dependency of arguments on directed selection
        list_directed_selections <- c("knock_out", "knock_in", "bottleneck", "migration", "coalescence", "resource_shift")

        if (any(unlist(output_row[,list_directed_selections]))) {
            # Set the flag TRUE
            #output_row$protocol <- "directed_selection"
            output_row$directed_selection <- TRUE

            # Check on the dependency of arguments. For example, bottleneck_size is not NA when bottleneck is TRUE
            if (all(output_row[list_directed_selections] == FALSE)) stop("Errors: A directed selection approach must be speicified")

            # exp_id
            if (!("exp_id" %in% names(args))) output_row$exp_id = paste(output_row$selected_function, output_row$protocol, list_directed_selections[unlist(output_row[,list_directed_selections])], output_row$seed, sep = "-")

            if (output_row$knock_in == TRUE) {
                temp1 <- output_row$knock_in_threshold
                if (is.na(temp1)) temp1 <- 0.95
                output_row$exp_id = paste(output_row$selected_function, output_row$protocol, "knock_in", paste0("p", temp1*100), output_row$seed, sep = "-")
            } else if (output_row$bottleneck == TRUE) {
                temp1 <- output_row$bottleneck_size
                if (is.na(temp1)) temp1 <- output_row$dilution
                output_row$exp_id = paste(output_row$selected_function, output_row$protocol, "bottleneck", 1/temp1, output_row$seed, sep = "-")
            } else if(output_row$migration == TRUE) {
                temp1 <- output_row$s_migration
                temp2 <- output_row$n_migration
                if (is.na(temp1)) temp1 <- "log"
                if (is.na(temp2)) temp2 <- 1000000
                output_row$exp_id = paste(output_row$selected_function, output_row$protocol, "migration", temp1, temp2, output_row$seed, sep = "-")
            } else if(output_row$coalescence == TRUE) {
                temp1 <- output_row$frac_coalescence
                if (is.na(temp1)) temp1 <- 0.5
                output_row$exp_id = paste(output_row$selected_function, output_row$protocol, "coalescence", paste0("p", temp1*100), output_row$seed, sep = "-")
            } else if(output_row$resource_shift == TRUE) {
                temp1 <- output_row$r_type
                temp2 <- output_row$r_percent
                if (is.na(temp1)) temp1 <- "add"
                if (is.na(temp2)) temp2 <- 0.01
                output_row$exp_id = paste(output_row$selected_function, output_row$protocol, "resource_shift", temp1, paste0("p", temp2*100), output_row$seed, sep = "-")
            }
        }

        # Check on the possible bug. For example, n_transfer must be larger than n_transfer_selection
        if (output_row$n_transfer < output_row$n_transfer_selection) stop("Errors: n_transfer must be greater than n_transfer_selection")
        if (output_row$n_transfer < output_row$function_lograte) stop("Errors: n_transfer must be greater than the function_lograte")
        if (output_row$n_transfer < output_row$composition_lograte) stop("Errors: n_transfer must be greater than the composition_lograte")
        #if (output_row$protocol != "directed_selection" & output_row$directed_selection == T) stop("protocol name needs to be changed to directed selection")

        # Placeholder for checking whether the protocol is available

    }

    # Modify the parameter format so it's readible by python
    output_row[sapply(output_row, isTRUE)] <- "True"
    output_row[sapply(output_row, isFALSE)] <- "False"

    # If exp_id speficied in the arguments, use it
    if ("exp_id" %in% names(args)) output_row$exp_id <- args$exp_id

    # Turn off scientific notation printout for large / small fractions numeric paramters
    for (i in 1:length(output_row[1,])) if (is.numeric(output_row[1,i])) output_row[1,i] <- format(output_row[1,i], scientific = FALSE)

    #
    output_row$seed <- as.numeric(output_row$seed)
    output_row$composition_lograte <- as.numeric(output_row$composition_lograte)
    output_row$n_transfer <- as.numeric(output_row$n_transfer)
    output_row$n_transfer_selection  <- as.numeric(output_row$n_transfer_selection)

    return(output_row)
}

list_treatments <- tibble(
    folder_id = c(rep("simple_medium", 8), rep("rich_medium", 8)),
    exp_id = c(paste0("simple_medium", 1:8), paste0("rich_medium", 1:8)),
    rich_medium = c(rep(F, 8), rep(T, 8)),
    n_inoc = c(rep(10^3, 8), rep(10^6, 8)),
    n_propagation = c(rep(20, 8), rep(1,8)),
    dilution = c(rep(0.001, 16)),
    l = rep(c(0, 0, 0, 0, 0.5, 0.5, 0.5, 0.5), 2),
    q = rep(c(0, 0, 0.8, 0.8, 0, 0, 0.8, 0.8), 2),
    sn = c(rep(c(2100, 700), 8)),
    sf = c(rep(c(1, 3), 8)),
    Sgen = 0,
    rn = c(rep(c(90, 30), 8)),
    rf = c(rep(c(1, 3), 8)),
    n_wells = 10,
    power_alpha = 0.01,
    sampling = "Binary_Gamma",
    n_transfer = 10,
    n_transfer_selection = 10,
    save_function = F,
    composition_lograte = 1
)
# 4 and 8 does no work
list_treatments <- filter(list_treatments, !grepl("medium4", exp_id) & !grepl("medium8", exp_id))


input_independent_wrapper <- function (i, treatment) {
    df <- make_input_csv()
    df <- df[rep(1,22),]; row.names(df) <- 1:nrow(df)
    experiments <- c("monoculture", "top_down_community", paste0("pair_from_random_species_", 1:10), paste0("pair_from_top_down_community_", 1:10))
    monoculture <- c(T, rep(F, 21))
    df$exp_id <- paste0(treatment$exp_id, "-", experiments, "-", i)
    df$monoculture <- monoculture
    df$overwrite_plate <- c(NA, NA, paste0(treatment$output_dir, treatment$exp_id, "-", experiments[-c(1,2)], "-", i, ".txt"))
    df$n_wells = c(NA, 10, rep(NA, 20))

    for (j in 2:ncol(treatment)) df[,names(treatment)[j]] = treatment[,names(treatment)[j]]
    df[is.na(df)] <- "NA"
    return(df)
}

input_set_wrapper <- function (i, treatment) {
    df <- make_input_csv()
    df <- df[rep(1,2),]; row.names(df) <- 1:nrow(df)
    experiments <- c("monoculture", "top_down_community")
    monoculture <- c(T, F)
    df$exp_id <- paste0(treatment$exp_id, "-", experiments, "-", i)
    df$monoculture <- monoculture
    df$overwrite_plate <- c(NA, NA)
    df$n_wells = c(NA, 10)
    df$output_dir <- paste0(data_directory, "set_", treatment$folder_id, "/")
    for (j in 3:ncol(treatment)) df[,names(treatment)[j]] = treatment[,names(treatment)[j]]
    df[is.na(df)] <- "NA"
    return(df)
}

input_synthetic_wrapper <- function (i, treatment) {
    df <- make_input_csv()
    df <- df[rep(1,20),]; row.names(df) <- 1:nrow(df)
    experiments <- c(paste0("pair_from_random_species_", 1:10), paste0("pair_from_top_down_community_", 1:10))
    monoculture <- c(rep(F, 20))
    df$exp_id <- paste0(treatment$exp_id, "-", experiments, "-", i)
    df$monoculture <- monoculture
    df$n_wells = c(rep(NA, 20))
    df$output_dir <- paste0(data_directory, "synthetic_", treatment$folder_id, "/")
    df$overwrite_plate <- c(paste0(df$output_dir, treatment$exp_id, "-", experiments, "-", i, ".txt"))
    for (j in 3:ncol(treatment)) df[,names(treatment)[j]] = treatment[,names(treatment)[j]]
    df[is.na(df)] <- "NA"
    return(df)

}

#input_independent_list <- rep(list(rep(list(NA), length(seeds))), nrow(list_treatments))
input_set_list <- rep(list(rep(list(NA), length(seeds))), nrow(list_treatments))
input_synthetic_list <- rep(list(rep(list(NA), length(seeds))), nrow(list_treatments))

for (k in 1:nrow(list_treatments)) {
    #input_independent_list[[k]][[1]] <- input_independent_wrapper(i = 1, treatment = list_treatments[k,])
    input_set_list[[k]][[1]] <- input_set_wrapper(i = 1, treatment = list_treatments[k,])
    input_synthetic_list[[k]][[1]] <- input_synthetic_wrapper(i = 1, treatment = list_treatments[k,])
    for (i in seeds) {
        # input_independent_list[[k]][[i]] <- input_independent_list[[k]][[1]] %>%
        #     mutate(seed = i, exp_id = sub("-\\d$", paste0("-", i), exp_id)) %>%
        #     mutate(seed = i, overwrite_plate = sub("-\\d.txt", paste0("-", i, ".txt"), overwrite_plate))
        input_set_list[[k]][[i]] <- input_set_list[[k]][[1]] %>%
            mutate(seed = i, exp_id = sub("-\\d$", paste0("-", i), exp_id)) %>%
            mutate(seed = i, overwrite_plate = sub("-\\d.txt", paste0("-", i, ".txt"), overwrite_plate))
        input_synthetic_list[[k]][[i]] <- input_synthetic_list[[k]][[1]] %>%
            mutate(seed = i, exp_id = sub("-\\d$", paste0("-", i), exp_id)) %>%
            mutate(seed = i, overwrite_plate = sub("-\\d.txt", paste0("-", i, ".txt"), overwrite_plate))

    }
    #input_independent <- rbindlist(input_independent_list[[k]])
    input_set <- rbindlist(input_set_list[[k]])
    input_synthetic <- rbindlist(input_synthetic_list[[k]])
    #fwrite(input_independent, paste0(mapping_file_directory, "input_independent.csv"))
    #fwrite(input_independent, paste0(mapping_file_directory, paste0("input_independent_", list_treatments$exp_id[k],".csv")))
    #fwrite(input_set, paste0(mapping_file_directory, paste0("input_set_", list_treatments$exp_id[k],".csv")))
    #fwrite(input_synthetic, paste0(mapping_file_directory, paste0("input_synthetic_", list_treatments$exp_id[k],".csv")))
}

# cat("\nMaking pooled input_independent.csv\n")
# input_independent <- input_independent_list %>% lapply(rbindlist) %>% rbindlist()
# fwrite(input_independent, paste0(mapping_file_directory, "input_independent.csv"))

# cat("\nMaking pooled input_set.csv\n")
# input_set <- input_set_list %>% lapply(rbindlist) %>% rbindlist()
# fwrite(input_set, paste0(mapping_file_directory, "input_set.csv"))
#
# cat("\nMaking pooled input_synthetic.csv\n")
# input_synthetic <- input_synthetic_list %>% lapply(rbindlist) %>% rbindlist()
# fwrite(input_synthetic, paste0(mapping_file_directory, "input_synthetic.csv"))

bind_rows(input_set_list) %>%
    filter(grepl("simple_medium", exp_id)) %>%
    fwrite(paste0(mapping_file_directory, "input_set_simple_medium.csv"))

bind_rows(input_synthetic_list) %>%
    filter(grepl("simple_medium", exp_id)) %>%
    fwrite(paste0(mapping_file_directory, "input_synthetic_simple_medium.csv"))


bind_rows(input_set_list) %>%
    filter(grepl("rich_medium", exp_id)) %>%
    fwrite(paste0(mapping_file_directory, "input_set_rich_medium.csv"))

bind_rows(input_synthetic_list) %>%
    filter(grepl("rich_medium", exp_id)) %>%
    fwrite(paste0(mapping_file_directory, "input_synthetic_rich_medium.csv"))

