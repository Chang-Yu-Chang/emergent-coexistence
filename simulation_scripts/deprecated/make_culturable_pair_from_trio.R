#' Make every possible synthetic pairs from a given trio

suppressWarnings(suppressMessages(library(tidyverse)))
suppressWarnings(suppressMessages(library(data.table)))

args = commandArgs(trailingOnly = T)
input_independent <- fread(args[[1]])
#input_independent <- fread("../data/raw/simulation/mapping_files/input_independent.csv")
input_independent_trios <- input_independent %>% filter(grepl("trio-culturable_isolates", exp_id))

#
read_trio_list <- function (trio_culturable_list) {
    trio_culturable_list %>% 
        filter(Type == "consumer") %>% 
        select(Well, Trio, InitialFrequency, ID) %>% 
        group_by(Well) %>% 
        mutate(Isolate = 1:3) %>% 
        pivot_wider(names_from = Isolate, names_prefix = "Isolate", values_from = ID) %>% 
        ungroup()
}

# Trio parameters
n_pairs_per_trio = choose(3,2) 
frequency1 <- c(0.05, 0.5, 0.95)
frequency2 <- c(0.95, 0.5, 0.05)
stopifnot(length(frequency1) == length(frequency2))
n_initial_frequencies <- length(frequency1)

cat("\nMaking pairs from the trios")
for (i in 1:nrow(input_independent_trios)) {
    n_trios <- input_independent_trios$n_trios[i]
    seed <- input_independent_trios$seed[i]
    scale <- input_independent_trios$scale[i] 
    output_dir <- input_independent_trios$output_dir[i]
    cat("\nexp_id = ", input_independent_trios$exp_id[i])
    
    df_monoculture_culturable <- fread(paste0("../data/raw/simulation/monoculture-culturable-", seed, ".txt"))
    df_trio_list <- fread(paste0("../data/raw/simulation/trio-culturable-", i, ".txt")) %>% 
        read_trio_list() %>% select(Trio, Isolate1, Isolate2, Isolate3) %>% 
        distinct()
    df_trio_pair_list <- df_trio_list %>% 
        split.data.frame(f = .$Trio) %>% 
        lapply(function(x) {
            temp <- t(combn(c(x$Isolate1, x$Isolate2, x$Isolate3), 2)) 
            tibble(Isolate1 = temp[,1], Isolate2 = temp[,2]) %>% 
                mutate(Pair = paste0("Pair", 1:n_pairs_per_trio)) %>% 
                select(Pair, everything())
        }) %>% 
        rbindlist(idcol = "Trio")

    # Make the initial plate for trios of culturable isolates
    df_culturable_pair_from_trio <- tibble(
        Well = paste0("W", 0:((n_trios*n_pairs_per_trio*n_initial_frequencies)-1)),
        Trio = paste0("Trio", rep(1:n_trios, each = n_pairs_per_trio*n_initial_frequencies)),
        Pair = paste0("Pair", rep(1:n_pairs_per_trio, n_initial_frequencies * n_trios)),
        # Isolate1 = rep(trios_sampled[,1], n_initial_frequencies),
        # Isolate2 = rep(trios_sampled[,2], n_initial_frequencies),
        # Isolate3 = rep(trios_sampled[,3], n_initial_frequencies),
        InitialFrequency = rep(rep(paste0(100*frequency1, "-", 100*frequency2), each = n_pairs_per_trio), n_trios),
        Abundance1 = rep(rep(frequency1, each = n_pairs_per_trio), n_trios),
        Abundance2 = rep(rep(frequency2, each = n_pairs_per_trio), n_trios),
    ) %>% 
        left_join(df_trio_pair_list, by = c("Trio", "Pair")) %>% 
        group_by(Well) %>% 
        pivot_longer(cols = starts_with("Isolate"), names_to = "Isolate", values_to = "ID") %>% 
        mutate(Abundance = ifelse(Isolate == "Isolate1", Abundance1,  Abundance2), Type = "consumer") %>% 
        select(Well, Trio, Pair, InitialFrequency, Type, ID, Abundance)
    
    df_resource_one_well <- fread(paste0("../data/raw/simulation/monoculture-", seed, "_composition.txt")) %>% 
        filter(Transfer == 0, Well == "W0", Type %in% c("resource", "R0")) %>% 
        select(Type, ID, Abundance)
    df_resource <- df_culturable_pair_from_trio %>% select(Well, Trio, Pair, InitialFrequency) %>% distinct() %>%
        right_join(df_resource_one_well, by = character())
    
    df_culturable_trios_resource <- bind_rows(df_culturable_pair_from_trio, df_resource) %>% mutate(Transfer = 0)
    fwrite(df_culturable_trios_resource, file = paste0(output_dir, "pair-culturable_from_trio-", seed, ".txt"))
    cat("\tmade a total of ", n_trios*3, " pairwise combination from ", n_trios," trios")
}
