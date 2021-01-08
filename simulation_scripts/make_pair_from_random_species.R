#' Make synthetic pairs from culturable isolates

suppressWarnings(suppressMessages(library(tidyverse)))
suppressWarnings(suppressMessages(library(data.table)))

args = commandArgs(trailingOnly = T)
input_independent <- fread(args[[1]])
#input_independent <- fread("../data/raw/simulation/mapping_files/input_independent_simple_medium.csv")
temp <- args[[1]] %>% strsplit("/") %>% `[[`(1)
treatment <- sub("input_synthetic_", "", temp[length(temp)]) %>% sub(".csv", "", .) # simple_medium
input_independent_random_species <- input_independent %>% filter(grepl("pair_from_random_species", exp_id))

# Pair parameters
frequency1 <- c(0.05, 0.5, 0.95)
frequency2 <- c(0.95, 0.5, 0.05)
stopifnot(length(frequency1) == length(frequency2))
n_initial_frequencies <- length(frequency1)
n_species = 10
n_pairs = choose(10, 2)
list_species <- rep(list(NA), nrow(input_independent_random_species))

cat("\nMaking random pairs from culturable isolates")
for (i in 1:nrow(input_independent_random_species)) {
    seed1 <- input_independent_random_species$seed[i]
    set1 <- sub(paste0(treatment, "-pair_from_random_species_"), "", input_independent_random_species$exp_id[i]) %>% strsplit("-") %>% unlist %>% `[`(1) %>% as.numeric()
    scale <- input_independent_random_species$scale[i]
    output_dir <- input_independent_random_species$output_dir[i]
    cat("\nexp_id = ", input_independent_random_species$exp_id[i])

    #
    df_monoculture_culturable <- fread(paste0("../data/raw/simulation/", treatment, "-monoculture_culturable-", seed1, ".txt"))
    culturable_isolates = unique(df_monoculture_culturable$ID)

    # Consumer for one well
    set.seed(seed1 * set1)
    sampled_species <- sort(sample(culturable_isolates, n_species)) # Species ID is 0-indexed
    sampled_pairs <- t(combn(sampled_species, 2))
    list_species[[i]] <- tibble(exp_id = input_independent_random_species$exp_id[i], Seed = seed1, Set = set1, ID = sampled_species)

    df_culturable_pairs <- tibble(
        Well = paste0("W", 0:((n_pairs*n_initial_frequencies)-1)),
        Pair = paste0("Pair", rep(1:n_pairs, n_initial_frequencies)),
        Isolate1 = rep(sampled_pairs[,1], n_initial_frequencies),
        Isolate2 = rep(sampled_pairs[,2], n_initial_frequencies),
        InitialFrequency = rep(paste0(100*frequency1, "-", 100*frequency2), each = n_pairs),
        Abundance1 = rep(frequency1, each = n_pairs),
        Abundance2 = rep(frequency2, each = n_pairs)
    ) %>%
        pivot_longer(cols = starts_with("Isolate"), names_to = "Isolate", values_to = "ID") %>%
        mutate(Abundance = ifelse(Isolate == "Isolate1", Abundance1, Abundance2), Type = "consumer") %>%
        select(Well, Pair, InitialFrequency, Type, ID, Abundance)

    # Resource for one well
    df_resource_one_well <- fread(paste0(output_dir, "synthetic_", sub("\\d+$", "", treatment), "/", treatment, "-monoculture-", seed1, "_composition.txt")) %>%
        filter(Transfer == 0, Well == "W0", Type %in% c("resource", "R0")) %>%
        select(Type, ID, Abundance)
    df_resource <- df_culturable_pairs %>% select(Well, Pair, InitialFrequency) %>% distinct() %>%
        right_join(df_resource_one_well, by = character())
    df_culturable_pairs_resource <- bind_rows(df_culturable_pairs, df_resource) %>% mutate(Transfer = 0)

    #
    fwrite(df_culturable_pairs_resource, file = paste0(output_dir, "synthetic_", sub("\\d+$", "", treatment), "/", input_independent_random_species$exp_id[i], ".txt"))
    cat("\t ", set1)

}

# list_species %>%
#     rbindlist() %>%
#     fwrite(file = paste0(output_dir, treatment, "-random_species.txt"))
