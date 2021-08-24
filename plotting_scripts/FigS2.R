# Fig S2

library(tidygraph)
library(ggraph)
library(tidyverse)
library(data.table)
library(cowplot)
source("network_functions.R")

communities <- fread("../data/output/communities.csv")
community_names_ordered_by_size <- communities %>% arrange(CommunitySize) %>% pull(Community)
community_abundance <- fread("../data/temp/communities_abundance.csv")
load("../data/output/network_community.Rdata") # Load observed networks net_list
load("../data/output/example_motif_list.Rdata") # Load example motif graphs example_motifs
load("../data/output/network_randomized.Rdata") # Load net_randomized_list
interaction_color <- assign_interaction_color()

# Panel: adjacent matrix
net_list_ordered_by_size <- rep(list(NA), length(net_list))
for (i in 1:length(net_list)) net_list_ordered_by_size[[i]] <- net_list[[community_names_ordered_by_size[i]]]
names(net_list_ordered_by_size) <- community_names_ordered_by_size
p_mats <- net_list_ordered_by_size %>% lapply(function(x) {plot_adjacent_matrix(x) + ggtitle("")})
p_mats[[14]] <- (plot_adjacent_matrix(net_list_ordered_by_size[[13]]) +
                     theme(legend.position = "right", legend.title = element_blank(), legend.text = element_text(size = 25)) +
                     guides(fill = guide_legend())) %>%
    ggpubr::get_legend() %>%
    ggpubr::as_ggplot()
p1 <- plot_grid(plotlist = p_mats, nrow = 2)#, labels = community_names_ordered_by_size)

#ggsave("../plots/FigS2.png", plot = p1, width = 18, height = 6)
net_randomized_list %>%
    lapply(count_motif)
"
generate random_motif_counts_percentile
"
count_motif(net_randomized_list$C1R2[[1]])
random_motif_counts_percentile

random_motif_counts_percentile <- random_motif_counts_percentile %>%
    mutate(Motif = factor(Motif)) %>%
    mutate(Community = ordered(Community, levels = community_names_ordered_by_size))
colors <- c("observed" = "red", "random [5th and 95th percentiles]" = "black")

p2 <- summary_network_motifs %>%
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

plot_grid(p1, p2, nrow = 2)

ggsave("../plots/FigS2.png", p, width = 10, height = 8)


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

    summary_network_motifs <- net_list %>%
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
