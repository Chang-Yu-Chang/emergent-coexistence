#' Make pairwise combinations from top-down assembled communities

suppressWarnings(suppressMessages(library(tidyverse)))
suppressWarnings(suppressMessages(library(data.table)))

args = commandArgs(trailingOnly = T)
input_independent <- fread(args[[1]])
#input_independent <- fread("../data/raw/simulation/mapping_files/input_independent.csv")
input_independent_comm <- input_independent %>% filter(grepl("community-top_down", exp_id))

# Read communities 
read_community_composition <- function (community_composition) {
    community_composition %>% 
        filter(Transfer == max(Transfer), Type == "consumer") %>% 
        mutate(Community = paste0("Community", sub("W", "", Well))) %>% 
        select(Community, Type, ID, Abundance)
}

# Pair parameters
frequency1 <- c(0.05, 0.5, 0.95)
frequency2 <- c(0.95, 0.5, 0.05)
stopifnot(length(frequency1) == length(frequency2))
n_initial_frequencies <- length(frequency1)


for (i in 1:nrow(input_independent_comm)) {
    comm_seed <- input_independent_comm$seed[i]
    scale <- input_independent_comm$scale[i] 
    output_dir <- input_independent_comm$output_dir[i]
    cat("\nexp_id = ", input_independent_comm$exp_id[i])
    
    # Number of communities to be examined
    input_independent_pair_from_comm <- input_independent %>%  filter(grepl("pair-from_top_down", exp_id), seed == comm_seed)
    n_comms <- nrow(input_independent_pair_from_comm)
    comms <- input_independent_pair_from_comm$exp_id %>% gsub(paste0("pair-from_top_down_community-", i, "-community"), "", .)
    df_community_composition <- fread(paste0("../data/raw/simulation/community-top_down-", comm_seed, "_composition.txt")) %>% 
        read_community_composition() %>% 
        filter(Community %in% paste0("Community", comms))
    
    # Resource
    df_resource_one_well <- fread(paste0("../data/raw/simulation/monoculture-", comm_seed, "_composition.txt")) %>% 
        filter(Transfer == 0, Well == "W0", Type %in% c("resource", "R0")) %>% 
        select(Type, ID, Abundance)
    
    temp_list <- rep(list(NA), n_comms)
    names(temp_list) <- comms
    
    for (j in 1:n_comms) {
        
        df_community_composition_subset <- df_community_composition %>% filter(Community == paste0("Community", comms[j]))
        df_isolates <- df_community_composition_subset %>%
            mutate(RelativeAbundance = Abundance/sum(Abundance)) %>%
            filter(RelativeAbundance >= 0.01) %>% 
            pull(ID)
        n_isolates <- length(df_isolates)
        n_pairs <- choose(n_isolates, 2)
        cat("\nCommunity ", comms[j], " has ", n_isolates, " isolates, ", n_pairs, " possible pairs, and ", n_pairs * n_initial_frequencies, " unique pair-frequency")
        
        # Make pair list
        temp <- t(combn(df_isolates, 2))
        df_comm_pair_list <- tibble(Community = comms[j], Pair = paste0("Pair", 1:n_pairs), Isolate1 = temp[,1], Isolate2 = temp[,2])
        
        # Make the initial plate of pairs from the communities
        df_pair_from_community <- tibble(
            Well = paste0("W", 0:(n_pairs*n_initial_frequencies-1)),
            Community = rep(comms[j], n_pairs*n_initial_frequencies),
            Pair = paste0("Pair", rep(1:n_pairs, n_initial_frequencies)),
            InitialFrequency = rep(paste0(100*frequency1, "-", 100*frequency2), each = n_pairs),
            Abundance1 = rep(frequency1, each = n_pairs),
            Abundance2 = rep(frequency2, each = n_pairs),
        ) %>% 
            left_join(df_comm_pair_list, by = c("Community", "Pair")) %>% 
            group_by(Well) %>% 
            pivot_longer(cols = starts_with("Isolate"), names_to = "Isolate", values_to = "ID") %>% 
            mutate(Abundance = ifelse(Isolate == "Isolate1", Abundance1,  Abundance2), Type = "consumer") %>% 
            select(Well, Community, Pair, InitialFrequency, Type, ID, Abundance)
        
        df_resource <- df_pair_from_community %>% select(Well, Community, Pair, InitialFrequency) %>% distinct() %>%
            right_join(df_resource_one_well, by = character())
        
        bind_rows(df_pair_from_community, df_resource) %>% mutate(Transfer = 0) %>% 
            fwrite(file = paste0(output_dir, "pair-from_top_down_community-", comm_seed, "-community", comms[j], ".txt"))
    }
    
}


