#' Make a R list of randomized  networks
library(tidyverse)
library(igraph)
library(tidygraph)
source(here::here("plotting_scripts/misc.R"))

load("~/Dropbox/lab/emergent-coexistence/data/output/communities_network.Rdata")
communities <- read_csv("~/Dropbox/lab/emergent-coexistence/data/output/communities.csv", col_types = cols())


# Randomize networks
b <- 1000 # Number of bootstrapping/randomization

# Make an empty two-layer R list
net_randomized_list <- rep(list(rep(list(NA), b)), nrow(communities_network))
names(net_randomized_list) <- communities$Community

time_start <- proc.time()
for (i in 1:nrow(communities_network)) {
    time_t <- proc.time()
    for (b_loop_index in 1:b) {
        set.seed(b_loop_index)
        net_randomized_list[[i]][[b_loop_index]] <- randomize_network(communities_network$Network[[i]])
        if (b_loop_index %% 1000 == 0) cat("\n boostrap =", b_loop_index)
    }
    # Print
    cat("\n\n", communities$Community[i])
    cat("\n", (proc.time() - time_t)[3], "seconds")

}
cat("\n\n total time:", (proc.time() - time_start)[3], "seconds\n\n")

communities_network_randomized <- communities_network %>%
    ungroup() %>%
    mutate(NetworkRandomized = net_randomized_list)

# Save the data file
save(communities_network_randomized, file = "~/Dropbox/lab/emergent-coexistence/data/output/communities_network_randomized.Rdata")
