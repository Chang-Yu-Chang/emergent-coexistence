# Simulate the network topography

library(tidygraph)
library(ggraph)
library(tidyverse)
library(data.table)
library(cowplot)
source("network_functions.R")

# Read data
communities <- fread("../data/temp/communities.csv")
community_names_ordered_by_size <- communities %>% arrange(CommunitySize) %>% pull(Community)
simulated_motif_counts <- fread("../data/temp/simulated_motif_counts.txt")
observed_motif_counts <- fread("../data/temp/observed_motif_counts.txt")
random_motif_counts <- fread("../data/temp/random_motif_counts.txt")
random_motif_counts_percentile <- fread("../data/temp/random_motif_counts_percentile.txt")
load("../data/temp/graph_list.Rdata") # Load observed networks graph_list
load("../data/temp/example_motif_list.Rdata") # Load example motif graphs example_motifs

# Panel A cartoon for experiment
# p_A <- ggdraw() + draw_image("../data/experimental_scheme/Fig1A.png")

# Panel B: motif distribution as a function of pairwise coexistence
## Motif demo
colors_grey <- grey(seq(1,0, length.out = length(example_motif_list)))
names(colors_grey) = 1:7
for (i in 1:length(example_motif_list)) {
    p_motif_list[[i]] <- example_motif_list[[i]] %>% 
        plot_competitive_network(node_size = 3) +
        theme(panel.background = element_rect(fill = colors_grey[i], color = NA))
    }
p_motifs <- plot_grid(plotlist = p_motif_list, nrow = 1)

## Motif count
motif_type <- c("others", "0-scored")
motif_color = c("#DB7469", "#557BAA")
names(motif_color) <- motif_type

p1 <- simulated_motif_counts %>% 
    mutate(Motif = factor(Motif)) %>% 
    group_by(CommunitySize, ProbPairCoexistence, Motif) %>% 
    summarize(MeanCount = mean(Count)) %>% 
    group_by(CommunitySize, ProbPairCoexistence) %>% 
    mutate(SumMeanCount = sum(MeanCount), RelativeMeanCount = MeanCount/SumMeanCount) %>% 
    filter(CommunitySize == 12) %>% 
    ggplot() +
    geom_area(aes(x = ProbPairCoexistence, y = RelativeMeanCount, fill = Motif), color = 1) +
    scale_x_continuous(expand = c(0,0), breaks = c(0, 0.5, 1)) +
    scale_y_continuous(expand = c(0,0), breaks = c(0, 0.5, 1)) +
    scale_fill_manual(values = colors_grey) +
    theme_cowplot() +
    panel_border(color = 1) +
    theme(legend.title = element_blank(),
        legend.position = "none", legend.direction = "horizontal") +
    labs(x = "Probability of pairwise coexistence", y = "Relative motif count")

p_B <- plot_grid(p_motifs, p1, ncol = 1, rel_heights = c(1,5))
ggsave("../plots/Fig1B.png", plot = p_B, width = 6, height = 6)

if (FALSE) {
    ## 0-scored motif vs others
    p2 <- simulated_motif_counts %>% 
        mutate(MotifType = ifelse(Motif %in% c(1, 7), "0-scored", "others")) %>% 
        group_by(CommunitySize, Seed, ProbPairCoexistence, MotifType) %>% 
        summarize(Count = sum(Count)) %>% 
        group_by(CommunitySize, ProbPairCoexistence, MotifType) %>% 
        summarize(MeanCount = mean(Count)) %>% 
        filter(CommunitySize == 12) %>% 
        ggplot() +
        geom_area(aes(x = ProbPairCoexistence, y = MeanCount, fill = MotifType), color = 1) +
        scale_x_continuous(expand = c(0,0), breaks = c(0, 0.5, 1)) +
        scale_y_continuous(expand = c(0,0), breaks = c(0, 100, 200)) +
        scale_fill_manual(values = motif_color) +
        theme_cowplot() +
        theme(legend.title = element_blank(),
            legend.position = "top") +
        labs(x = "Fraction of pairwise coexistence", y = "Mean motif count")
    
    ggsave("../plots/Fig1B_example.png", plot = p2, width = 4, height = 4)
    
    
    # Possible supp
    # Panel C
    ps1 <- simulated_motif_counts %>% 
        mutate(Motif = factor(Motif)) %>% 
        group_by(CommunitySize, ProbPairCoexistence, Motif) %>% 
        summarize(MeanCount = mean(Count)) %>% 
        group_by(CommunitySize, ProbPairCoexistence) %>% 
        mutate(SumMeanCount = sum(MeanCount), RelativeMeanCount = MeanCount/SumMeanCount) %>% 
        ggplot() +
        geom_area(aes(x = ProbPairCoexistence, y = RelativeMeanCount, fill = Motif), color = 1) +
        scale_x_continuous(expand = c(0,0), breaks = c(0, 0.5, 1)) +
        scale_y_continuous(expand = c(0,0), breaks = c(0, 0.5, 1)) +
        scale_fill_manual(values = colors_grey) +
        facet_grid(.~CommunitySize) +
        guides(fill = guide_legend(nrow = 1)) +
        theme_cowplot() +
        theme(legend.position = "top", 
            legend.direction = "horizontal") +
        panel_border(color = 1) +
        labs(x = "Probability of pairwise coexistence", y = "Relative motif count")
    
    ggsave("../plots/FigS1.png", plot = ps1, width = 10, height = 4)
    
    # Panel XX: motif count as a function of community size
    p_motifs <- example_motif_list %>% 
        lapply(function(x) plot_competitive_network(x, node_size = 3)) %>% 
        plot_grid(plotlist = ., nrow = 1)
    
    summary_network_motifs <- graph_list %>%
        lapply(summarize_network_motif) %>%
        bind_rows(.id = "Community")
    
    ps2 <- summary_network_motifs %>%
        ggplot(aes(x = CommunitySize, y = RelativeMotifCount)) +
        geom_jitter(size = 3, shape = 21, width = 0.1) +
        geom_smooth(method = "lm", formula = y ~ x) +
        scale_x_continuous(limits = c(2, 13), breaks = c(2, 7, 12)) +
        scale_y_continuous(limits = c(-0.001,1.001), breaks = c(0, 0.5, 1)) +
        facet_wrap(.~Motif, nrow = 1) +
        theme_cowplot() +
        theme(strip.background = element_blank(), strip.text = element_blank()) +
        panel_border(color = "black") +
        labs(x = "Community size", y = "Relative motif count")
    
    p <- plot_grid(p_motifs, ps2, ncol = 1, axis = "rl", align = "hv", rel_heights = c(2,5))
    ggsave("../plots/FigS2.png", plot = p, width = 10, height = 4)

 }

# Panel C: fraction of pairwise coexistence as a function of community size
summary_network_pairs <- graph_list %>% 
    lapply(summarize_network_pairs) %>%
    bind_rows()

p_C <- summary_network_pairs %>%
    ggplot(aes(x = NumberNodes, y = FractionCoexistence)) +
    geom_jitter(size = 3, shape = 21) +
    geom_smooth(method = "lm", formula = y ~ x) +
    scale_x_continuous(breaks = 1:13) +
    theme_cowplot() +
    panel_border(color = "black") +
    labs(x = "Community size", y = "Fraction of pairwise coexistence")

ggsave("../plots/Fig1C.png", p_C, width = 4, height = 4)

if (FALSE) {
    # # Panel XX: example of one netowkr and adjacent matrix 
    # p_net <- plot_competitive_network(graph_list$C11R2, node_size = 4, layout = "circle")
    # p_mat <- plot_adjacent_matrix(graph_list$C11R2)
    # p_merged <- plot_grid(p_net, p_mat)
    # ggsave("../plots/Fig_example.png", p_merged, width = 8, height = 4)
    
    # Panel: adjacent matrix
    graph_list_ordered_by_size <- rep(list(NA), length(graph_list))
    for (i in 1:length(graph_list)) graph_list_ordered_by_size[[i]] <- graph_list[[community_names_ordered_by_size[i]]]
    names(graph_list_ordered_by_size) <- community_names_ordered_by_size
    p_mats <- graph_list_ordered_by_size %>% lapply(plot_adjacent_matrix)
    ps3 <- plot_grid(plotlist = p_mats, nrow = 2, labels = community_names_ordered_by_size)
    
    ggsave("../plots/FigS3.png", plot = ps3, width = 20, height = 6)
    
}



# Panel D: the normalized motif counts, relative to the expectation at random
random_motif_counts_confidence_intervals <- random_motif_counts %>%
    group_by(Community, Motif) %>%
    summarize(MeanCount = mean(Count), SdCount = sd(Count))

p_D <- observed_motif_counts %>% 
    left_join(random_motif_counts_confidence_intervals) %>% 
    mutate(StandardizedCount = (Count - MeanCount)/SdCount) %>%
    mutate(CommunitySize = factor(CommunitySize)) %>% 
    ggplot() +
    # geom_rect(xmin = -Inf, xmax = Inf, ymin = 0.34, ymax = -0.34, fill = "grey", alpha = 0.5) +
    # geom_rect(xmin = -Inf, xmax = Inf, ymin = 0.475, ymax = -0.475, fill = NA, color = 1, linetype = 2) +
    geom_ribbon(aes(x = Motif,  ymin = 0.38, ymax = -0.38), fill = "grey", alpha = 0.5) +
    geom_ribbon(aes(x = Motif,  ymin = 0.475, ymax = -0.475), fill = NA, color = 1, linetype = 2) +
    geom_point(aes(x = Motif, y = StandardizedCount, color = Community, size = CommunitySize), shape = 21) +
    scale_x_continuous(breaks = 1:7) +
    # facet_grid(.~Motif, scales = "free_x") +
    theme_cowplot() +
    theme(legend.position = "top") + 
    guides(color = F, size = guide_legend(nrow = 1, title.position = "top")) +
    labs(x ="Motif", y = "Standardized count") 

#p <- plot_grid(p_motifs, p4, ncol = 1, axis = "lr", align = "v", rel_heights = c(1,5))
ggsave("../plots/Fig1D.png", p_D, width = 4, height = 4)

if (FALSE) {
    random_motif_counts_percentile <- random_motif_counts_percentile %>% 
        mutate(Motif = factor(Motif)) %>% 
        mutate(Community = ordered(Community, levels = community_names_ordered_by_size))
    colors <- c("observed" = "red", "random [5th and 95th percentiles]" = "black")
    
    ps4 <- summary_network_motifs %>% 
        mutate(Community = ordered(Community, levels = community_names_ordered_by_size)) %>% 
        ggplot() +
        geom_point(data = random_motif_counts_percentile, aes(x = Motif, y = Count, color = "random [5th and 95th percentiles]")) +
        geom_segment(data = pivot_wider(random_motif_counts_percentile, names_from = Percentile, values_from = Count), 
            aes(x = Motif, xend = Motif, y = p05, yend = p95, color = "random [5th and 95th percentiles]")) +
        geom_point(aes(color = "observed", x = Motif, y = Count)) +
        scale_color_manual(values = colors) +
        facet_wrap(Community~., scales = "free_y", nrow = 2) +
        theme_cowplot() + 
        theme(legend.position = "bottom") +
        panel_border(color = "black") +
        labs(color = "")
    
    ggsave("../plots/FigS4.png", ps4, width = 10, height = 4)
    
}


# Panel E: competitive hierarchy
#graph 
example_motifs[[1]] %>%
    


#

p <- plot_grid(p_B, p_C, p_D, ncol = 2, axis = "tblr", align = "hv")

ggsave("../plots/Fig1.png", p, width = 10, height = 10)




if (FALSE){
    simulated_motif_counts_mean <- simulated_motif_counts %>% 
        group_by(CommunitySize, ProbPairCoexistence, Motif) %>% 
        summarize(MeanCount = mean(Count))
    
    simulated_motif_counts_mean %>% 
        ggplot() +
        geom_area(aes(x = ProbPairCoexistence, y = MeanCount, fill = Motif), color = 1) +
        scale_x_continuous(expand = c(0,0)) +
        scale_y_continuous(expand = c(0,0)) +
        facet_grid(CommunitySize ~ ., scales = "free_y") +
        theme_cowplot()
    
    
    ggraph(graph, layout = "linear", circular = T) +
    geom_node_point() +
    geom_edge_arc(aes(color = interaction), show.legend = T) +
    theme_graph()
    
    
    temp_list <- rep(list(rep(list(NA), b)), length(p_range))
    names(temp_list) <- p_range
    
    for (j in 1:length(p_range)) {
        cat("\np =", p_range[j], "\n")
        for (i in 1:b) {
            temp_list[[j]][[i]] <- count_motif(make_random_network(n = n, p = p_range[j]))
            if (i%%10 == 0) cat(i, " ")
        }
    }
    
    motif_counts <- temp_list %>%
        lapply(function(x) {
            lapply(x, function(y) {tibble(Motif = factor(1:7), Count = y)}) %>%
                rbindlist(idcol = "Seed")
        }) %>% 
        rbindlist(idcol = "p")
    
    p1 <- motif_count_mean %>%  
        ggplot(aes(x = p, y = MeanCount, color = Motif, group = Motif)) +
        geom_point() + geom_line() +
        theme_bw()
    
}


