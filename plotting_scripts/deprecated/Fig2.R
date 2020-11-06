# Detect motif in the empirical network

library(cowplot)
library(tidygraph)
library(ggraph)
library(tidyverse)
library(data.table)
source("network_functions.R")

# Read data
simulated_motif_counts <- fread("../data/temp/simulated_motif_counts.txt")
load("../data/temp/graph_list.Rdata") # Load observed networks graph_list
load("../data/temp/example_motif_list.Rdata") # Load example motif graphs example_motifs
# communities_abundance <- read_csv("../data/temp/communities_abundance.csv")
# communities <- read_csv("../data/temp/communities.csv")
# community_names <- communities$Community
# community_sizes <- communities$CommunitySize
# isolates <- read_csv("../data/output/isolates.csv")
# isolates_melted <- read_csv("../data/output/isolates_melted.csv")
# pairs <- read_csv("../data/output/pairs.csv")
# community_names_ordered_by_size <- communities %>% arrange(CommunitySize) %>% pull(Community)


# Panel XX: motif distribution, compared to randomized network
b = 100

cat("\n Randomizing the empirical graphs")
temp_list <- rep(list(rep(list(NA), b)), length(graph_list))
names(temp_list) <- names(graph_list)
for (j in 1:length(graph_list)){
    cat("\ngraph:", names(graph_list)[j], "\n")
    for (i in 1:b) {
        temp_list[[j]][[i]] <- count_motif(randomize_network(graph_list[[j]]))
        if (i%%10 == 0) cat(i, " ")
    }
}
motif_counts <- temp_list %>%
    lapply(function(x) {
        lapply(x, function(y) {tibble(Motif = factor(1:7), Count = y)}) %>%
            rbindlist(idcol = "Seed")
    }) %>% 
    bind_rows(.id = "Community")

motif_counts_p95 <- motif_counts %>% 
    group_by(Community, Motif) %>% 
    filter(Count >= quantile(Count, 0.95)) %>% 
    distinct(Community, Motif, Count) %>% 
    arrange(Community, Motif, Count) %>% 
    slice_min(Count) %>% 
    mutate(Percentile = "p95")
motif_counts_p05 <- motif_counts %>% 
    group_by(Community, Motif) %>% 
    filter(Count <= quantile(Count, 0.05)) %>% 
    distinct(Community, Motif, Count) %>% 
    arrange(Community, Motif, Count) %>% 
    slice_max(Count) %>% 
    mutate(Percentile = "p05")
motif_counts_percentile <- bind_rows(motif_counts_p05, motif_counts_p95) %>% 
    mutate(Community = ordered(Community, levels = community_names_ordered_by_size))


colors <- c("observed" = "red", "random [5th and 95th percentiles]" = "black")

p3 <- summary_network_motifs %>% 
    mutate(Community = ordered(Community, levels = community_names_ordered_by_size)) %>% 
    ggplot(aes(x = Motif, y = Count)) +
    geom_point(data = motif_counts_percentile, aes(x = Motif, y = Count, color = "random [5th and 95th percentiles]")) +
    geom_segment(data = pivot_wider(motif_counts_percentile, names_from = Percentile, values_from = Count), 
        aes(x = Motif, xend = Motif, y = p05, yend = p95, color = "random [5th and 95th percentiles]")) +
    geom_point(aes(color = "observed")) +
    scale_color_manual(values = colors) +
    facet_wrap(Community~., scales = "free_y", nrow = 2) +
    theme_cowplot() + 
    theme(legend.position = "bottom") +
    panel_border(color = "black") +
    labs(color = "")
p3

ggsave("../plots/Fig2C.png", plot = p3, width = 14, height = 4)


# Panel XX: pooled networks
motif_counts_aggregated <- motif_counts %>% 
    group_by(Seed, Motif) %>% 
    summarize(Count = sum(Count))

motif_counts_aggregated_p95 <- motif_counts_aggregated %>% 
    group_by(Motif) %>% 
    filter(Count >= quantile(Count, 0.95)) %>% 
    distinct(Motif, Count) %>% 
    arrange(Motif, Count) %>% 
    slice_min(Count) %>% 
    mutate(Percentile = "p95")
motif_counts_aggregated_p05 <- motif_counts_aggregated %>% 
    group_by(Motif) %>% 
    filter(Count <= quantile(Count, 0.05)) %>% 
    distinct(Motif, Count) %>% 
    arrange(Motif, Count) %>% 
    slice_min(Count) %>% 
    mutate(Percentile = "p05")
motif_counts_aggregated_percentile <- bind_rows(motif_counts_aggregated_p05, motif_counts_aggregated_p95) 

summary_network_motifs_aggregated <- summary_network_motifs %>% 
    group_by(Motif) %>% 
    summarize(Count = sum(Count))


plot_example_motifs <- function(node_size=5) {
    temp_list <- rep(list(NA), 7)
    temp_id <- c(11, 7, 8, 12, 13, 14, 15)
    
    for (i in 1:7) {
        g <- as_tbl_graph(igraph::graph.isocreate(size = 3, temp_id[i]))
        layout <- create_layout(g, layout = 'circle')
        g <- activate(g, edges) %>% mutate(InteractionType = ifelse(edge_is_mutual(), "coexistence", "exclusion"))
        temp_list[[i]] <- g %>% activate(nodes) %>% mutate(x = layout$x, y = layout$y, graph = paste0(i))
    }
    
    merged_graph <- bind_graphs(temp_list)
    
    plot_competitive_network(merged_graph, layout = "example_motif", node_size = node_size) + 
        facet_nodes(~graph, nrow = 1)
}
p_example <- plot_example_motifs(node_size = 3)


colors <- c("observed" = "red", "random [5th and 95th percentiles]" = "black")
p4 <- summary_network_motifs_aggregated %>% 
    ggplot(aes(x = Motif, y = Count)) +
    geom_point(data = motif_counts_aggregated_percentile, aes(x = Motif, y = Count, color = "random [5th and 95th percentiles]"), size = 3) +
    geom_segment(data = pivot_wider(motif_counts_aggregated_percentile, names_from = Percentile, values_from = Count), 
        aes(x = Motif, xend = Motif, y = p05, yend = p95, color = "random [5th and 95th percentiles]")) +
    geom_point(aes(color = "observed"), size = 3) +
    scale_color_manual(values = colors) +
    facet_grid(.~Motif, scales = "free_x") +
    theme_cowplot() + 
    theme(legend.position = "bottom", strip.text = element_blank(), strip.background = element_blank(),
        axis.text.x = element_blank()) +
    labs(color = "")

p <- plot_grid(p_example, p4, ncol = 1, axis = "rl", align = "vh", rel_heights = c(2, 7))
p
ggsave("../plots/Fig2D.png", plot = p, width = 8, height = 5)



# Combining the plots
p <- plot_grid(p2, p3, ncol = 1, align = "v", rel_heights = c(1, 2))

ggsave("../plots/Fig2.png", plot = p, width = 10, height = 5)









if (FALSE) {
    
    graph_list[[11]] %>%
        plot_competitive_network()
    
    # BArplot
    motif_counts %>%
        mutate(Community = ordered(Community, level = community_names)) %>% 
        group_by(Community, Seed) %>% 
        mutate(TotalMotifCount = sum(Count)) %>% 
        group_by(Community, Seed, Motif) %>% 
        summarize(RelativeMotifCount = Count/TotalMotifCount) %>% 
        ggplot(aes(x = Seed, y = RelativeMotifCount, fill = Motif)) +
        geom_bar(position = "stack", stat = "identity") +
        scale_x_continuous(expand = c(0,0)) +
        scale_y_continuous(expand = c(0,0)) +
        facet_grid(Community~.) + 
        theme_bw()
    
    
}
#










