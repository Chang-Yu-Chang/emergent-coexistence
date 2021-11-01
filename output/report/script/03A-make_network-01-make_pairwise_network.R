#' Make pairwise network from outcome of pairwise competitions
library(tidyverse)
library(data.table)
library(tidygraph)
library(ggraph)
source(here::here("plotting_scripts/network_functions.R"))

# Read data ----
communities <- read_csv(here::here("data/output/communities.csv"))
isolates <- read_csv(here::here("data/output/isolates.csv"))
pairs <- read_csv(here::here("data/output/pairs.csv"))

# Make network, which self-contain all isolates and pairs information
net_list <- rep(list(NA), length(communities$Community)) # tbl graph object
names(net_list) <- communities$Community
#names(net_tbl_list) <- communities$Community

for (i in 1:length(communities$Community)) {
    net_list[[i]] <- make_network(isolates = filter(isolates, Community == communities$Community[i]),
                                  pairs = filter(pairs, Community == communities$Community[i]))
}

# Plot
p_net_list <- rep(list(NA), length(net_list))
for (i in 1:length(net_list)) p_net_list[[i]] <- plot_competitive_network(net_list[[i]])

save(net_list, p_net_list, file = here::here("data/output/network_community.Rdata"))
