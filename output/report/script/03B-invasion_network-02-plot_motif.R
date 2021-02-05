#' Plot the observed motif counts and randomized network motif counts ----
library(tidyverse)
library(data.table)
library(igraph)
library(tidygraph)
library(ggraph)
library(cowplot)
source(here::here("plotting_scripts/network_functions.R"))

# Data
isolates <- fread(here::here("data/output/isolates.csv"))
pairs <- fread(here::here("data/output/pairs.csv"))
communities <- fread(here::here("data/output/communities.csv"))
communities_name <- communities$Community
networks_motif <- fread(here::here("data/temp/networks_motif.csv")) %>% mutate(Community = ordered(Community, communities_name))
networks_motif_randomized <- fread(here::here("data/temp/networks_motif_randomized.csv")) %>% mutate(Community = ordered(Community, communities_name))
networks_motif_randomized_percentile <- fread(here::here("data/temp/networks_motif_randomized_percentile.csv")) %>% mutate(Community = ordered(Community, communities_name))
b <- max(networks_motif_randomized$Randomization ) # Extract the number of bootstrapping


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

# Plot
#set_graph_style(plot_margin = margin(1,1,1,1))
interaction_type <- c("exclusion", "coexistence", "lose", "bistability", "neutrality", "self", "undefined")
myColor = c("#DB7469", "#557BAA", "#73C966", "#EECF6D", "#8650C4", "black", "grey80")
names(myColor) <- interaction_type
node_size <- 5

# Example motif plot
p_motif_example <- plot_grid(plotlist = lapply(motif_list, function(x) plot_competitive_network(x, node_size=3)), nrow = 1)


# Plot observations and 5% and 95% percentiles ----
p_motif_randomized <-
  ggplot() +
  # 5% and 95% percentiles in randomized networks
  geom_point(data = networks_motif_randomized_percentile, aes(x = Motif, y = CountMotif, group = Motif), col = "black") +
  geom_segment(data = spread(networks_motif_randomized_percentile, Percentile, CountMotif),
               aes(x = Motif, xend = Motif, y = p5, yend = p95), col = "black") +
  # Observations
  geom_point(data = networks_motif, aes(x = Motif, y = CountMotif), col = "red") +
  #  scale_colour_manual(name="Error Bars",values=c("black", "red")) +
  facet_wrap(Community ~., scale = "free_y") +
  theme_bw() +
  ggtitle("Motif counts in the invasion networks") +
  labs(x = "motif", y = "Count")


# Plot histogram in motif count distribution ----
p_motif_randomized_histo <-
  networks_motif_randomized %>%
#  filter(Community == "C11R1") %>%
  ggplot() +
  geom_histogram(aes(CountMotif), fill = NA, col = 1) +
  geom_vline(data = networks_motif, aes(xintercept = CountMotif), col = "red") +
  facet_grid(Community ~ Motif, scale = "free") +
  coord_flip() +
  theme_bw() +
  theme(strip.text.x = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 90)) +
  labs(x = "motif counts", y = "counts in randomization")

