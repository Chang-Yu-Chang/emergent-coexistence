# Detect motif in the empirical network

library(cowplot)
library(tidygraph)
library(ggraph)
library(tidyverse)
library(data.table)
source("network_functions.R")

# Read data
communities_abundance <- read_csv("../data/temp/communities_abundance.csv")
communities <- read_csv("../data/temp/communities.csv")
community_names <- communities$Community
community_sizes <- communities$CommunitySize
isolates <- read_csv("../data/output/isolates.csv")
isolates_melted <- read_csv("../data/output/isolates_melted.csv")
pairs <- read_csv("../data/output/pairs.csv")


# Panel XX: Make networks
graph_list <- rep(list(NA), length(community_names))
names(graph_list) <- community_names

for (i in 1:length(community_names)) {
    isolates_subset <- isolates %>% filter(Community == community_names[i])
    pairs_subset <- pairs %>% filter(Community == community_names[i])
    graph_list[[i]] <- make_network(isolates = isolates_subset, pairs = pairs_subset)
}

# Panel XX: fraction of pairwise coexistence as a function of community size
summary_network_pairs <- graph_list %>% 
    lapply(summarize_network_pairs) %>%
    bind_rows()

p1 <- summary_network_pairs %>%
    ggplot(aes(x = NumberNodes, y = FractionCoexistence)) +
    geom_jitter(size = 3, shape = 21) +
    geom_smooth(method = "lm", formula = y ~ x) +
    scale_x_continuous(breaks = 1:13) +
    theme_cowplot() +
    panel_border(color = "black") +
    labs(x = "Community size", y = "Fraction of pairwise coexistence")

ggsave("../plots/Fig2B.png", p1, width = 4, height = 4)

# Panel XX: example of one netowkr and adjacent matrix 
p_net <- plot_competitive_network(graph_list$C11R2, node_size = 4, layout = "circle")
p_mat <- plot_adjacent_matrix(graph_list$C11R2)

p_merged <- plot_grid(p_net, p_mat)
ggsave("../plots/Fig_example.png", p_merged, width = 8, height = 4)


# Panel XX: motif count as a function of community size
summary_network_motifs <- graph_list %>%
    lapply(summarize_network_motif) %>%
    bind_rows(.id = "Community")

p2 <- summary_network_motifs %>%
    ggplot(aes(x = NumberNodes, y = RelativeMotifCount)) +
    geom_jitter(size = 3, shape = 21, width = 0.1) +
    geom_smooth(method = "lm", formula = y ~ x) +
    scale_x_continuous(limits = c(2, 13), breaks = c(2, 7, 12)) +
    scale_y_continuous(limits = c(-0.001,1.001)) +
    facet_wrap(.~Motif, nrow = 1) +
    theme_cowplot() +
    theme(strip.background = element_blank(), strip.text = element_blank()) +
    panel_border(color = "black") +
    labs(x = "Community size", y = "Fraction of motif")


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
motif_counts_percentile <- bind_rows(motif_counts_p05, motif_counts_p95) 

p5 <- summary_network_motifs %>% 
    mutate(Community = ordered(Community, level = community_names)) %>% 
    ggplot(aes(x = Motif, y = Count)) +
    geom_point(data = motif_counts_percentile, aes(x = Motif, y = Count)) +
    geom_segment(data = pivot_wider(motif_counts_percentile, names_from = Percentile, values_from = Count), 
        aes(x = Motif, xend = Motif, y = p05, yend = p95)) +
    geom_point(color = "red") +
    facet_wrap(Community~., scales = "free_y") +
    theme_cowplot() + 
    panel_border(color = "black")

ggsave("../plots/Fig2C.png", plot = p5, width = 10, height = 10)

# Combining the plots
p <- plot_grid(p2, p2, ncol = 1, align = "v", rel_heights = c(1, 2))

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










