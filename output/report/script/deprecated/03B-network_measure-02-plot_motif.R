#' Plot the observed motif counts and randomized network motif counts
library(tidyverse)
library(data.table)
library(igraph)
library(tidygraph)
library(ggraph)
library(cowplot)
source(here::here("plotting_scripts/network_functions.R"))

isolates <- read_csv(here::here("data/output/isolates.csv"))
pairs <- read_csv(here::here("data/output/pairs.csv"))
communities <- read_csv(here::here("data/output/communities.csv"))
networks_motif <- read_csv(here::here("data/output/networks_motif.csv")) %>% mutate(Community = ordered(Community, communities$Community))
networks_motif_randomized <- read_csv(here::here("data/output/networks_motif_randomized.csv")) %>% mutate(Community = ordered(Community, communities$Community))
networks_motif_randomized_percentile <- read_csv(here::here("data/temp/networks_motif_randomized_percentile.csv")) %>% mutate(Community = ordered(Community, communities$Community))
b <- max(networks_motif_randomized$Randomization) # Extract the number of bootstrapping


# Plot the motif examples ----
# Motif list
motif_list <- rep(list(NA), 7)
temp_id <- c(11, 7, 8, 12, 13, 14, 15)
g <- igraph::graph.isocreate(size = 3, 1)
layout <- ggraph::create_layout(g, layout = 'circle')

for (i in 1:7) {
  g <- igraph::graph.isocreate(size = 3, temp_id[i])
  E(g)$InteractionType <- ifelse(which_mutual(g), "coexistence", "exclusion")
  g <- tidygraph::as_tbl_graph(g) %>%
    mutate(x = layout$x, y = layout$y, graph = paste0("motif", i))
  motif_list[[i]] <- g
}

save(motif_list, file = here::here("data/output/motif_list.Rdata"))




if (FALSE) {
# Plot
#set_graph_style(plot_margin = margin(1,1,1,1))
interaction_type <- c("exclusion", "coexistence", "lose", "bistability", "neutrality", "self", "undefined")
myColor = c("#DB7469", "#557BAA", "#73C966", "#EECF6D", "#8650C4", "black", "grey80")
names(myColor) <- interaction_type
node_size <- 5

# Example motif plot
p_motif_example <- plot_grid(plotlist = lapply(motif_list, function(x) plot_competitive_network(x, node_size=3)), nrow = 1)

}


