#' This scripts generate the network figures for pool pairs and community pairs
#'
#'  1. monoculture set networks
#'  2. community networks
#'  3. hierarchy scores

library(tidyverse)
library(cowplot)
library(tidygraph)
library(ggraph)
source(here::here("analysis/00-metadata.R"))

# 0. load network objects ----
load(paste0(folder_simulation, "11-aggregated/monocultureSets_network.Rdata"))
load(paste0(folder_simulation, "11-aggregated/communities_network.Rdata"))
# monocultureSets_network
# communities_network

# 0.1 plotting functions ----
set.seed(1)
node_size = 3
edge_width = .8
net <- monocultureSets_network$Network[[1]]
plot_network_hierarchy <- function(net, n_rank = 12, n_break = 12) {
    graph_ranked <- net %>%
        activate(nodes) %>%
        select(Species, Rank, PlotRank) %>%
        activate(edges) %>%
        mutate(fromRank = .N()$PlotRank[match(from, .N()$Species)],
               toRank = .N()$PlotRank[match(to, .N()$Species)])

    n_nodes <- igraph::vcount(graph_ranked)
    n_exclusion_violation <- graph_ranked %>%
        activate(edges) %>%
        filter(fromRank > toRank, InteractionType == "exclusion") %>%
        igraph::ecount()

    if (n_exclusion_violation != 0) {
        pairs_exclusion_violation <- graph_ranked %>%
            activate(edges) %>%
            filter(fromRank > toRank, InteractionType == "exclusion") %>%
            mutate(temp = fromRank, fromRank = toRank, toRank = temp) %>%
            select(-temp) %>%
            mutate(InteractionType = "exclusion violating rank") %>%
            as_tibble

        graph_ranked <- graph_ranked %>%
            activate(edges) %>%
            filter(fromRank <= toRank) %>%
            bind_edges(pairs_exclusion_violation)
    }

    temp <- graph_ranked %>%
        # Node position
        activate(nodes) %>%
        mutate(y = -Rank) %>%
        group_by(Rank) %>%
        mutate(x = {seq(0, 1, length.out = n() + 2) %>% `[`(c(-1, -length(.)))}) %>% # + rnorm(n(), 0, .5)) %>%
        ungroup() %>%
        # Filter out coexistence edges
        activate(edges) %>%
        filter(!(InteractionType == "coexistence" & fromRank > toRank)) %>%
        mutate(Temp = sample(c(-1, 1), size = n(), replace = T))
    strength_angle <- as_tibble(temp)$Temp * 0.06

    temp %>%
        ggraph(layout = "nicely") +
        geom_hline(yintercept = c(-n_rank:-1), color = "grey90") +
        #geom_segment(aes(x = 0.5, xend = 0.5, y = -Inf, yend = -1), color = "grey90") +
        geom_node_point(size = node_size, shape = 21, fill = "grey", stroke = node_size/5, color = "black") +
        geom_edge_arc(strength = strength_angle, alpha = 0.5,
                      aes(color = InteractionType), width = edge_width,
                      arrow = arrow(length = unit(edge_width, "mm"), type = "closed", angle = 30, ends = "last"),
                      start_cap = circle(node_size/2, "mm"),
                      end_cap = circle(node_size/2, "mm")) +
        scale_edge_color_manual(values = assign_interaction_color(level = "hierarchy"),
                                breaks = c("exclusion", "exclusion violating rank", "coexistence", "unknown"),
                                labels = c("exclusion pair that follows rank", "exclusion pair that violates rank", "coexistence", "unknown")) +
        scale_x_continuous(limits = c(0.1, 0.9), expand = c(0,0)) +
        scale_y_continuous(limits = c(-n_break-1, 0), breaks = -n_break:-1, labels = n_break:1) +
        theme_void() +
        theme(
            legend.position = "none",
            legend.title = element_blank(),
            axis.title = element_blank(),
            strip.text = element_blank(),
            plot.margin=unit(c(0,0,0,0),"mm")
        ) +
        labs(y = "Rank")

}

# 1. Monoculture set network ----
monocultureSets_network_hierarchy <- monocultureSets_network %>%
    mutate(Community = factor(Community, Community)) %>%
    rowwise() %>%
    mutate(NetworkHierarchyPlot = plot_network_hierarchy(Network, n_rank = Richness) %>% list())

p <- plot_grid(plotlist = monocultureSets_network_hierarchy$NetworkHierarchyPlot) +
    paint_white_background()
ggsave(here::here("simulation/plots/24-monocultureSets_network.png"), p, width = 15, height = 10)


# 2. Community networks ----
communities_network_hierarchy <- communities_network %>%
    mutate(Community = factor(Community, Community)) %>%
    rowwise() %>%
    mutate(NetworkHierarchyPlot = ifelse(nrow(Species) == 1, list(NA), list(plot_network_hierarchy(Network, n_rank = Richness))))

p <- plot_grid(plotlist = communities_network_hierarchy$NetworkHierarchyPlot) +
    paint_white_background()
ggsave(here::here("simulation/plots/24-communities_network.png"), p, width = 15, height = 10)


# 3. Hierarchy score ----
p <- bind_rows(monocultureSets_network, communities_network) %>%
    mutate(Treatment = factor(rep(c("monocultureSet", "community"), each = 20), c("monocultureSet", "community"))) %>%
    ggplot(aes(x = Treatment, y = HierarchyScore)) +
    geom_boxplot(linewidth = .8, outlier.color = NA) +
    geom_jitter(shape = 1, size = 2, stroke = .8, height = 0, width = .1) +
    scale_y_continuous(limits = c(0,1.01), breaks = c(0, .25, .5, .75, 1)) +
    theme_classic() +
    theme(panel.grid.major.y = element_line(color = "grey", linetype = 2),
          panel.spacing = unit(0, "mm"),
          panel.border = element_rect(fill = NA, color = 1, size = 1.5),
          axis.text = element_text(size = 10, color = 1),
          #axis.text.x = element_blank(),
          axis.title = element_text(size = 10, color = 1),
          legend.title = element_blank(),
          plot.title = element_text(size = 10, color = 1)
    ) +
    labs(x = "", y = "Hierarchy score") +
    guides(color = "none")

ggsave(here::here("simulation/plots/24-hierarchy_scores.png"), p, width = 4, height = 3)


# p_net_hierarchy_list[[13]] <- p_net_hierarchy_list[[13]] +
#     scale_y_continuous(limits = c(-12-1, 0), breaks = -12:-1, labels = 12:1, position = "right") +
#     theme(axis.title.y = element_text(color = 1, size = 10, angle = 270, margin = margin(l = 2, unit = "mm")),
#           axis.text.y = element_text(color = 1, size = 10, margin = margin(l = 1, unit = "mm")))
# pB_axistitle <- ggdraw() + draw_label("Community", fontface = 'plain', x = .5, hjust = .5) + theme(plot.margin = margin(5, 0, 5, 7))
# p_temp <- plot_grid(plotlist = p_net_hierarchy_list,
#                     rel_widths = c(monocultureSets_network_hierarchy$CommunitySize / max(monocultureSets_network_hierarchy$CommunitySize))^1.5,
#                     labels = 1:13, label_fontface = "plain", label_x = c(rep(0.5, 12), 0.45), hjust = c(rep(.5, 12), 1),
#                     nrow = 1, axis = "tb", align = "h") + paint_white_background()
# pB <- plot_grid(pB_axistitle, p_temp, ncol = 1, rel_heights = c(.1, 1)) + paint_white_background()
#
# ## legend
# edge_width = 0.5
# temp <- p_net_hierarchy_list[[13]] +
#     geom_edge_arc(strength = 10,
#                   aes(color = InteractionType), width = edge_width*1.3,
#                   arrow = arrow(length = unit(edge_width*5, "mm"), type = "closed", angle = 30, ends = "last"),
#                   start_cap = circle(node_size/2, "mm"),
#                   end_cap = circle(node_size/2, "mm")) +
#     #scale_alpha_manual(values = 1) +
#     theme(legend.key.size = unit(2, "line"),
#           legend.key.height = unit(6, "mm"),
#           legend.position = "right",
#           legend.direction = "vertical",
#           legend.text = element_text(size = 10),
#           legend.background = element_rect(fill = NA, color = NA))
# p_legend <- get_legend(temp)
# pB <- ggdraw(pB) + draw_plot(p_legend,.2,.2,.1,.1)






# Matrix form of networks ----


plot_matrix <- function (species, pairs) {
    # species <- monocultureSets_network$Species[[i]]$Species
    # pairs <- monocultureSets_network$Pairs[[i]]
    pairs %>%
        bind_rows(tibble(Species1 = species$Species, Species2 = Species1)) %>%
        mutate(Species1 = factor(Species1, rev(species$Species))) %>%
        mutate(Species2 = factor(Species2, species$Species)) %>%
        ggplot() +
        geom_tile(aes(x = Species2, y = Species1, fill = InteractionType), color = 1, linewidth = .5) +
        scale_fill_manual(values = interaction_color) +
        scale_x_discrete(drop = FALSE, expand = c(0, 0), position = "top") +
        scale_y_discrete(drop = FALSE, expand = c(0, 0), position = "right") +
        theme_classic() +
        theme(
            legend.position = c(0.3, 0.4),
            legend.background = element_blank(),
            axis.text = element_blank(),
            axis.title = element_blank(),
            plot.margin = unit(c(1,0.2,0.2,0.2), "cm")
        ) +
        guides(fill = guide_legend(title = "")) +
        labs()
}

monocultureSets_network <- monocultureSets_network %>%
    rowwise() %>%
    mutate(MatrixPlot = plot_matrix(Species, Pairs) %>% list())

p <- plot_grid(plotlist = monocultureSets_network$MatrixPlot, labels = monocultureSets_network$Community,
               hjust = 0, label_x = 0.01) + paint_white_background()

ggsave(here::here("simulation/plots/24-monocultureSets_matrix.png"), p, width = 15, height = 10)

communities_network <- communities_network %>%
    rowwise() %>%
    mutate(MatrixPlot = plot_matrix(Species, Pairs) %>% list())

p <- plot_grid(plotlist = communities_network$MatrixPlot, labels = communities_network$Community,
               hjust = 0, label_x = 0.01) + paint_white_background()

ggsave(here::here("simulation/plots/24-communities_matrix.png"), p, width = 15, height = 10)
























