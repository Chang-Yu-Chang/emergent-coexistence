#' Make pairwise network from outcome of pairwise competitions
library(tidyverse)
library(data.table)
library(tidygraph)
library(ggraph)
source(here::here("plotting_scripts/network_functions.R"))

# Read data ----
communities <- fread(here::here("data/output/communities.csv"))
community_name <- communities$Community
community_size <- communities$CommunitySize
isolates <- fread(here::here("data/output/isolates.csv"))
pairs <- fread(here::here("data/output/pairs.csv"))

# Make network, which self-contain all isolates and pairs information
net_list <- rep(list(NA), length(community_name)) # tbl graph object
names(net_list) <- community_name
names(net_tbl_list) <- community_name

for (i in 1:length(community_name)) {
    net_list[[i]] <- make_network(isolates = filter(isolates, Community == community_name[i]),
                                  pairs = filter(pairs, Community == community_name[i]))
}

# Plot
p_net_list <- rep(list(NA), length(net_list))
for (i in 1:length(net_list)) p_net_list[[i]] <- plot_competitive_network(net_list[[i]])

save(net_list, p_net_list, file = here::here("data/output/network_community.Rdata"))
