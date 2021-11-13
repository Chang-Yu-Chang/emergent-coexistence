library(tidyverse)
source("misc.R")

# Read data
input_randomNetworks <- read_csv("input_randomNetworks.csv")


S = 300 # number of species in the pool; change this line when the pool size changes
n_networks = 10
size_network = 6

species_list <- paste0("S", sample(0:(S-1), size = n_networks * size_network, replace = F))
well_list <- rep(paste0("W", 0:(n_networks-1)), each = size_network)

df <- input_randomNetworks[rep(1:nrow(input_randomNetworks), each = n_networks*size_network),] %>%
    mutate(Well = rep(well_list, nrow(input_randomNetworks)), Species = rep(species_list, nrow(input_randomNetworks)))

# Generate pairs
generate_pairs <- function(x) {
    if (length(x) <= 1) return(tibble(Isolate1 = character(), Isolate2 = character()))
    if (length(x) >= 2) {
        x %>%
            sort() %>%
            combn(2) %>%
            t() %>%
            as_tibble %>%
            setNames(c("Isolate1", "Isolate2")) %>%
            return()
    }
}

df_pairs <- df %>%
    group_by(exp_id, seed, vamp, q2, Well) %>%
    summarize(Pairs = list(generate_pairs(Species))) %>%
    group_by(exp_id, seed, vamp, q2) %>%
    summarize(Pairs = list(bind_rows(Pairs) %>% distinct))


for (i in 1:nrow(df_pairs)) {
    n_pairs <- nrow(df_pairs$Pairs[i][[1]])
    N0 <- as_tibble(matrix(0, nrow = S, ncol = n_pairs)) %>% setNames(paste0("W", 0:(n_pairs-1)))
    species_list <- paste0("S", 0:(S-1))
    pair_list <- df_pairs$Pairs[i][[1]]
    n_pairs <- nrow(pair_list)
    N0 <- as_tibble(matrix(0, nrow = S, ncol = n_pairs*2)) %>% setNames(paste0("W", 0:(2*n_pairs-1)))
    for (j in 1:n_pairs) {
        index_pair <- c(match(pair_list$Isolate1[j], species_list), match(pair_list$Isolate2[j], species_list))
        N0[index_pair,2*j-1] <- c(0.1, 0.9)
        N0[index_pair,2*j] <- c(0.9, 0.1)
    }
    write_csv(N0, file = paste0(input_randomNetworks$output_dir[1], 'randomNetworks_', df_pairs$exp_id[i], '_init.csv'))
    cat(paste0("exp_id=", df_pairs$exp_id[i], "\n"))
}




