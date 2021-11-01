#' Make a R list of randomized  networks
library(tidyverse)
library(igraph)
library(tidygraph)
source(here::here("plotting_scripts/network_functions.R"))

isolates <- read_csv(here::here("data/output/isolates.csv"))
pairs <- read_csv(here::here("data/output/pairs.csv"))
communities <- read_csv(here::here("data/output/communities.csv"))
load(here::here("data/output/network_community.Rdata"))

# Make random networks and save them in a list ----
b <- 1000 # Number of bootstrapping/randomization

# Make an empty two-layer R list
net_randomized_list <- rep(list(rep(list(NA), b)), length(net_list))
names(net_randomized_list) <- communities$Community

tt <- proc.time()
for (i in 1:length(net_list)) {
  temp_tt <- proc.time()
  for (b_loop_index in 1:b) {
    set.seed(b_loop_index)
    net_randomized_list[[i]][[b_loop_index]] <- randomize_network(net_list[[i]])
    if (b_loop_index %% 1000 == 0) cat("\n boostrap =", b_loop_index)
  }
  # Print
  cat("\n\n", communities$Community[i])
  cat("\n", (proc.time() - temp_tt)[3], "seconds")
  if (i == length(net_list)) cat("\n\n total time:", (proc.time() - tt)[3], "seconds\n\n")

}

# Save the data file
save(net_randomized_list, file = here::here("data/output/network_randomized.Rdata"))
