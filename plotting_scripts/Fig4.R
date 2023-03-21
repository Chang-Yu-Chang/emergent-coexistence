library(tidyverse)
library(cowplot)
library(tidygraph)
library(ggraph)
library(igraph)
source(here::here("analysis/00-metadata.R"))

isolates <- read_csv(paste0(folder_data, "output/isolates.csv"), show_col_types = F)
pairs <- read_csv(paste0(folder_data, "output/pairs.csv"), show_col_types = F)
communities <- read_csv(paste0(folder_data, "temp/00c-communities.csv"), show_col_types = F)
pairs_freq <- read_csv(paste0(folder_data, "temp/93-pairs_freq.csv"), show_col_types = F)
load(paste0(folder_data, "temp/95-communities_network.Rdata"))
communities_hierarchy <- read_csv(paste0(folder_data, "temp/95-communities_hierarchy.csv"), show_col_types = F)


# Figure 4A One example network ----
p_net <- communities_network %>%
    filter(Community == "C2R6") %>%
    pull(Network) %>% `[[`(1) %>%
    activate(nodes) %>%
    mutate(y = -Rank) %>%
    arrange(Rank) %>%
    mutate(x = c(.3, .7, .3, .7)) %>%
    ungroup() %>%
    ggraph(layout = "nicely") +
    geom_edge_link(aes(color = InteractionType), width = 1,
                   arrow = arrow(length = unit(2, "mm"), type = "closed", angle = 30, ends = "last"),
                   start_cap = circle(3, "mm"),
                   end_cap = circle(3, "mm")) +
    scale_edge_color_manual(values = assign_interaction_color(),
                            breaks = c("exclusion", "exclusion violating rank"),
                            labels = c("exclusion following rank", "exclusion violating rank")) +
    scale_x_continuous(limits = c(.2, .8), expand = c(0,0)) +
    scale_y_continuous(limits = c(-4, 0.5), breaks = -4:-1, labels = 4:1) +
    theme_void() +
    theme(
        legend.position = "none",
        panel.grid.major.y = element_line(color = "grey90"),
        plot.margin=unit(c(0,0,0,0),"mm"),
        axis.text.y = element_text(color = 1, size = 10, margin = margin(r = 2, unit = "mm")),
        axis.title.y = element_text(color = 1, size = 10, angle = 90, margin = margin(r = 5, unit = "mm"))
    ) +
    guides(color = "none") +
    labs(y = "Rank")

pA <- p_net +
    draw_image(here::here("plots/cartoons/Fig2B_1.png"), x = -.2, y = -1, vjust = 0.5, hjust = 0, clip = "on", scale = .8) +
    draw_image(here::here("plots/cartoons/Fig2B_3.png"), x = .2, y = -1, vjust = 0.5, hjust = 0, clip = "on", scale = .8) +
    draw_image(here::here("plots/cartoons/Fig2B_4.png"), x = -.2, y = -3, vjust = 0.5, hjust = 0, clip = "on", scale = .8) +
    draw_image(here::here("plots/cartoons/Fig2B_2.png"), x = .2, y = -3, vjust = .6, hjust = 0, clip = "on", scale = .8) +
    paint_white_background()


# Figure 4B. Hierarchy network plot ----
set.seed(1)
node_size = 3
edge_width = .8
plot_network_hierarchy <- function(net, n_rank = 12, n_break = 12) {
    graph_ranked <- net %>%
        activate(nodes) %>%
        select(Isolate, Rank, PlotRank) %>%
        activate(edges) %>%
        mutate(fromRank = .N()$PlotRank[match(from, .N()$Isolate)],
               toRank = .N()$PlotRank[match(to, .N()$Isolate)])

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
        geom_node_point(size = node_size, shape = 21, fill = "grey", stroke = node_size/5, color = "black") +
        geom_edge_arc(strength = strength_angle, alpha = 0.5,
                      aes(color = InteractionType), linewidth = edge_width,
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
communities_network_hierarchy <- communities_network %>%
    mutate(Community = factor(Community, Community)) %>%
    arrange(CommunitySize) %>%
    rowwise() %>%
    mutate(NetworkHierarchyPlot = plot_network_hierarchy(Network, n_rank = CommunitySize) %>% list())

p_net_hierarchy_list <- communities_network_hierarchy$NetworkHierarchyPlot
p_net_hierarchy_list[[13]] <- p_net_hierarchy_list[[13]] +
    scale_y_continuous(limits = c(-12-1, 0), breaks = -12:-1, labels = 12:1, position = "right") +
    theme(axis.title.y = element_text(color = 1, size = 10, angle = 270, margin = margin(l = 2, unit = "mm")),
          axis.text.y = element_text(color = 1, size = 10, margin = margin(l = 1, unit = "mm")))
pB_axistitle <- ggdraw() + draw_label("Community", fontface = 'plain', x = .5, hjust = .5) + theme(plot.margin = margin(5, 0, 5, 7))
p_temp <- plot_grid(plotlist = p_net_hierarchy_list,
                    rel_widths = c(communities_network_hierarchy$CommunitySize / max(communities_network_hierarchy$CommunitySize))^1.5,
                    labels = 1:13, label_fontface = "plain", label_x = c(rep(0.5, 12), 0.45), hjust = c(rep(.5, 12), 1),
                    nrow = 1, axis = "tb", align = "h") + paint_white_background()
pB <- plot_grid(pB_axistitle, p_temp, ncol = 1, rel_heights = c(.1, 1)) + paint_white_background()

# line legend
edge_width = 0.5
p_legend <- get_legend({
    p_net_hierarchy_list[[13]] +
    geom_edge_arc(strength = 10,
                  aes(color = InteractionType), width = edge_width*1.3,
                  arrow = arrow(length = unit(edge_width*5, "mm"), type = "closed", angle = 30, ends = "last"),
                  start_cap = circle(node_size/2, "mm"),
                  end_cap = circle(node_size/2, "mm")) +
    #scale_alpha_manual(values = 1) +
    theme(legend.key.size = unit(2, "line"),
          legend.key.height = unit(6, "mm"),
          legend.position = "right",
          legend.direction = "vertical",
          legend.text = element_text(size = 10),
          legend.background = element_rect(fill = NA, color = NA))
 })

# # Figure 4C: Hierarchy ----
# pC <- mutate(communities_hierarchy, Treatment = "experiment") %>%
#     ggplot(aes(x = Treatment, y = HierarchyScore)) +
#     geom_boxplot(width = .5, lwd = .8, outlier.color = NA) +
#     geom_jitter(shape = 1, size = 2, stroke = .8, height = 0, width = .1) +
#     scale_y_continuous(limits = c(0,1.01), breaks = c(0, .25, .5, .75, 1)) +
#     theme_classic() +
#     theme(panel.grid.major.y = element_line(color = "grey", linetype = 2),
#           panel.spacing = unit(0, "mm"),
#           panel.border = element_rect(fill = NA, color = 1, linewidth = 1.5),
#           axis.text = element_text(size = 10, color = 1),
#           axis.text.x = element_blank(),
#           axis.title = element_text(size = 10, color = 1),
#           legend.title = element_blank(),
#           plot.title = element_text(size = 10, color = 1)
#     ) +
#     labs(x = "", y = "Score") +
#     guides(color = "none") +
#     ggtitle("Hierarchy")



# Figure 4C: transitivity
triad_possible <- tibble(CommunitySize = 3:12, TotalTriads = choose(CommunitySize, 3))

# Randomize the pairs by shuffling the network
load(paste0(folder_data, "temp/94-communities_network_randomized.Rdata"))
communities_network_randomized_rps <- communities_network_randomized %>%
    rowwise() %>%
    mutate(CountRPS_randomized = list(sapply(NetworkRandomized, function (x) triad_census(x)[10]))) %>%
    select(Community, CommunityLabel, CommunitySize, CountRPS_randomized) %>%
    unnest(CountRPS_randomized)


if (FALSE) {

# Boxplot
communities_network_randomized_rps %>%
    mutate(Community = factor(Community, communities$Community)) %>%
    ggplot(aes(x = Community, y = CountRPS_randomized)) +
    geom_boxplot(outlier.fill = NA) +
    geom_jitter(height = 0, shape = 21) +
    theme_classic() +
    theme() +
    labs()

# Violin
communities_network_randomized_rps %>%
    mutate(Community = factor(Community, communities$Community)) %>%
    ggplot(aes(x = Community, y = CountRPS_randomized)) +
    geom_violin(draw_quantiles = c(0.25, 0.5, 0.75), scale = "width") +
    geom_jitter(height = 0, shape = 21) +
    theme_classic() +
    theme() +
    labs()


# Histogram violin
communities_network_randomized_rps %>%
    mutate(Community = factor(Community, communities$Community)) %>%
    group_by(Community, CountRPS_randomized) %>%
    #group_by(CountRPS_randomized) %>%
    summarize(Count = n()) %>%
    ggplot() +
    geom_vline(xintercept = 0, linetype = 1, color = grey(0.9)) +
    #geom_hline(yintercept = 0:15, linetype = 1, color = grey(0.9)) +
    geom_tile(aes(x = 0, y = CountRPS_randomized, width = Count, height = 1), color = NA, fill = grey(0.8), linewidth = .5) +
    #geom_histogram(color = 1, fill = NA, binwidth = 1) +
    #geom_jitter(height = 0, shape = 21) +
    scale_x_continuous(breaks = seq(-500, 500, 500)) +
    scale_y_continuous(breaks = 0:15) +
    facet_grid(.~Community, labeller = labeller(setNames(communities$CommunityLabel, communities$Community))) +
    theme_classic() +
    theme(panel.spacing.x = unit(0,"mm"),
          strip.background = element_rect(color = NA)) +
    labs(x = "boostrap count", y = "number of RPS")

}

communities_network_rps <- communities_network %>%
    rowwise() %>%
    #mutate(ExclusionNetwork = extract_exclusion_network(Network) %>% list()) %>%
    mutate(CountRPS = triad_census(Network)[10]) # count rock-paper-scissor triads

#
pC <- communities_network_randomized_rps %>%
    mutate(Community = factor(Community, communities$Community)) %>%
    group_by(CommunityLabel, CountRPS_randomized) %>%
    summarize(Count = n()) %>%
    group_by(CommunityLabel) %>%
    mutate(Fraction = Count / sum(Count)) %>%
    ggplot() +
    geom_hline(yintercept = 0, linetype = 1, color = grey(0.9)) +
    geom_col(aes(x = CountRPS_randomized, y = Count, fill = "bootstrap"), color = NA, width = .9) +
    geom_point(data = communities_network_rps, y = 0, aes(x = CountRPS, shape = "observed RPS"), stroke = .5) +
    scale_x_continuous(breaks = 0:15) +
    scale_y_continuous(breaks = c(0, 500), label = c(0,500), expand = c(0.15,1)) +
    scale_shape_manual(values = c("observed RPS" = 21))+
    scale_fill_manual(values = c("bootstrap" = grey(0.8)))+
    coord_flip() +
    facet_grid(.~CommunityLabel) +
    theme_classic() +
    theme(
        legend.position = c(0.2, 0.85),
        legend.background = element_rect(fill = "white"),
        legend.margin = margin(0,0,0,0, "cm"),
        legend.spacing.y = unit(0, "cm"),
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.4, "cm"),
        panel.spacing.x = unit(0,"mm"),
        #panel.border = element_rect(fill = NA, color = 1),
        strip.background = element_rect(color = NA),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 5)
    ) +
    guides(fill = guide_legend(title = ""), shape = guide_legend(title = "")) +
    labs(x = "number of RPS", y = "bootstrap count")



if (FALSE) {
pC <- communities_network_rps %>%
    ggplot() +
    #geom_col(data = triad_possible, aes(x = CommunitySize, y = TotalTriads, fill = "possible triads"), color = NA, position = position_identity()) +
    geom_point(aes(x = CommunitySize, y = CountRPS, shape = "observed RPS"), position = position_dodge2(width = 0.3), stroke = .5) +
    scale_x_continuous(breaks = 1:12) +
    scale_shape_manual(values = c("observed RPS" = 21))+
    scale_fill_manual(values = c("possible triads" = grey(0.8)))+
    theme_classic() +
    theme(
        legend.position = c(0.3, 0.8),
        legend.background = element_blank(),
        legend.margin = margin(0,0,0,0, "cm"),
        legend.spacing.y = unit(-.1, "cm"),
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.4, "cm")
    ) +
    guides(fill = guide_legend(title = ""), shape = guide_legend(title = "")) +
    labs(x = "# of species", y = "# of triads")

}

#
p_left <- plot_grid(pA, NULL, ncol = 1, rel_heights = c(1,1.5), scale = c(.8, .9), labels = c("A", "C"), axis = "lr", align = "v")
p <- plot_grid(p_left, pB, nrow = 1, rel_widths = c(1,4), scale = c(1, .9), labels = c("", "B")) + paint_white_background()
p <- ggdraw(p) +
    draw_plot(pC, x = 0, y = 0, width = 0.4, height = 0.5,  hjust = 0, vjust = 0) +
    draw_plot(p_legend, x = 0.5, y = 0.15, width = 0.1, height = 0.1,  hjust = 0, vjust = 0)
ggsave(here::here("plots/Fig4.png"), p, width = 10, height = 5)


















