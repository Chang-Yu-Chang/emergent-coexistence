# Run simulation of the with given prability of pairwise coexistence

library(tidygraph)
library(tidyverse)
library(data.table)
source("network_functions.R")

#
n_sizes <- c(3, 4, 8, 12, 20) # Number of species
p_range <- seq(0, 1, by = 0.1) # Probability of pairwise coexistence
b = 100 # Number of seeds

temp_list <- rep(list(NA), length(n_sizes) * length(p_range) * b)
counter = 1

for (i in 1:length(n_sizes)) {
    cat ("\n\nCommunity size = ", n_sizes[i], "\n")
    for (j in 1:length(p_range)) {
        for (k in 1:b) {
            set.seed(k)
            temp_list[[counter]] <- make_random_network(n = n_sizes[i], p = p_range[j]) %>%
                count_motif() %>% 
                tibble(Motif = factor(1:7), Count = .) %>% 
                mutate(CommunitySize = n_sizes[i], ProbPairCoexistence = p_range[j], Seed = k)
            counter = counter + 1
            if (k%%100 == 0) cat(k, " ")
        }
    }
}

simulated_motif_counts <- temp_list %>% rbindlist() %>% 
    select(CommunitySize, ProbPairCoexistence, Seed,  Motif, Count)


fwrite(simulated_motif_counts, file = "../data/temp/simulated_motif_counts.txt")