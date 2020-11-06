#' Make synthetic pairs from culturable isolates

suppressWarnings(suppressMessages(library(tidyverse)))
suppressWarnings(suppressMessages(library(data.table)))

args = commandArgs(trailingOnly = T)
input_independent <- fread(args[[1]])
#input_independent <- fread("../data/raw/simulation/mapping_files/input_independent.csv")
input_independent_pairs <- input_independent %>% filter(grepl("pair-culturable_isolates", exp_id))

# Pair parameters
frequency1 <- c(0.05, 0.5, 0.95)
frequency2 <- c(0.95, 0.5, 0.05)
stopifnot(length(frequency1) == length(frequency2))
n_initial_frequencies <- length(frequency1)

cat("\nMaking random pairs from culturable isolates")
for (i in 1:nrow(input_independent_pairs)) {
    n_pairs <- input_independent_pairs$n_pairs[i]
    seed <- input_independent_pairs$seed[i]
    scale <- input_independent_pairs$scale[i] 
    output_dir <- input_independent_pairs$output_dir[i]
    cat("\nexp_id = ", input_independent_pairs$exp_id[i])
    
    df_monoculture_culturable <- fread(paste0("../data/raw/simulation/monoculture-culturable-", seed, ".txt"))
    n_culturable_isolates = nrow(df_monoculture_culturable)
    pairs_pool <- t(combn(df_monoculture_culturable$ID, 2)) # Species ID is 0-indexed
    n_pairs_pool <- choose(n_culturable_isolates, 2)
    set.seed(seed)
    pairs_sampled <- pairs_pool[sample(n_pairs_pool, size = n_pairs, replace = F), ]
    
    # Make the initial plate for pairs of culturable isolates
    df_culturable_pairs <- tibble(
        Well = paste0("W", 0:((n_pairs*n_initial_frequencies)-1)),
        Pair = paste0("Pair", rep(1:n_pairs, n_initial_frequencies)),
        Isolate1 = rep(pairs_sampled[,1], n_initial_frequencies),
        Isolate2 = rep(pairs_sampled[,2], n_initial_frequencies),
        InitialFrequency = rep(paste0(100*frequency1, "-", 100*frequency2), each = n_pairs),
        Abundance1 = rep(frequency1, each = n_pairs),
        Abundance2 = rep(frequency2, each = n_pairs)
    ) %>% 
        pivot_longer(cols = starts_with("Isolate"), names_to = "Isolate", values_to = "ID") %>% 
        mutate(Abundance = ifelse(Isolate == "Isolate1", Abundance1, Abundance2), Type = "consumer") %>% 
        select(Well, Pair, InitialFrequency, Type, ID, Abundance)
    
    df_resource_one_well <- fread(paste0("../data/raw/simulation/monoculture-", seed, "_composition.txt")) %>% 
        filter(Transfer == 0, Well == "W0", Type %in% c("resource", "R0")) %>% 
        select(Type, ID, Abundance)
    df_resource <- df_culturable_pairs %>% select(Well, Pair, InitialFrequency) %>% distinct() %>%
        right_join(df_resource_one_well, by = character())
    
    df_culturable_pairs_resource <- bind_rows(df_culturable_pairs, df_resource) %>% mutate(Transfer = 0)
    fwrite(df_culturable_pairs_resource, file = paste0(output_dir, "pair-culturable-", seed, ".txt"))
    cat("\tmade ", n_pairs, " random pairs of culturable isolates")
}



