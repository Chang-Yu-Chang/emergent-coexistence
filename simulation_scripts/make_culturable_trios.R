#' Make synthetic trios from culturable isolates

suppressWarnings(suppressMessages(library(tidyverse)))
suppressWarnings(suppressMessages(library(data.table)))

args = commandArgs(trailingOnly = T)
input_independent <- fread(args[[1]])
#input_independent <- fread("../data/raw/simulation/mapping_files/input_independent.csv")
input_independent_trios <- input_independent %>% filter(grepl("trio-culturable_isolates", exp_id))

# Trio parameters
frequency1 <- c(0.33)
frequency2 <- c(0.33)
frequency3 <- c(0.33)
stopifnot(length(frequency1) == length(frequency2) & length(frequency1) == length(frequency3))
n_initial_frequencies <- length(frequency1)

cat("\nMaking trios from culturable isolates")
for (i in 1:nrow(input_independent_trios)) {
    n_trios <- input_independent_trios$n_trios[i]
    seed <- input_independent_trios$seed[i]
    scale <- input_independent_trios$scale[i] 
    output_dir <- input_independent_trios$output_dir[i]
    cat("\nexp_id = ", input_independent_trios$exp_id[i])
    
    df_monoculture_culturable <- fread(paste0("../data/raw/simulation/monoculture-culturable-", seed, ".txt"))
    n_culturable_isolates = nrow(df_monoculture_culturable)
    set.seed(seed)
    trios_sampled <- matrix(NA, n_trios, 3)
    for (j in 1:n_trios) trios_sampled[j,] <- sort(sample(df_monoculture_culturable$ID, size = 3, replace = F))

    # Make the initial plate for trios of culturable isolates
    df_culturable_trios <- tibble(
        Well = paste0("W", 0:((n_trios*n_initial_frequencies)-1)),
        Trio = paste0("Trio", rep(1:n_trios, n_initial_frequencies)),
        Isolate1 = rep(trios_sampled[,1], n_initial_frequencies),
        Isolate2 = rep(trios_sampled[,2], n_initial_frequencies),
        Isolate3 = rep(trios_sampled[,3], n_initial_frequencies),
        InitialFrequency = rep(paste0(100*frequency1, "-", 100*frequency2, "-", 100*frequency3), each = n_trios),
        Abundance1 = rep(frequency1, each = n_trios),
        Abundance2 = rep(frequency2, each = n_trios),
        Abundance3 = rep(frequency3, each = n_trios)
    ) %>% 
        pivot_longer(cols = starts_with("Isolate"), names_to = "Isolate", values_to = "ID") %>% 
        mutate(Abundance = ifelse(Isolate == "Isolate1", Abundance1, ifelse(Isolate == "Isolate2", Abundance2, Abundance3)), Type = "consumer") %>% 
        select(Well, Trio, InitialFrequency, Type, ID, Abundance)
    
    df_resource_one_well <- fread(paste0("../data/raw/simulation/monoculture-", seed, "_composition.txt")) %>% 
        filter(Transfer == 0, Well == "W0", Type %in% c("resource", "R0")) %>% 
        select(Type, ID, Abundance)
    df_resource <- df_culturable_trios %>% select(Well, Trio, InitialFrequency) %>% distinct() %>%
        right_join(df_resource_one_well, by = character())
    
    df_culturable_trios_resource <- bind_rows(df_culturable_trios, df_resource) %>% mutate(Transfer = 0)
    fwrite(df_culturable_trios_resource, file = paste0(output_dir, "trio-culturable-", seed, ".txt"))
    cat("\tmade ", n_trios, " random trios of culturable isolates")
}
