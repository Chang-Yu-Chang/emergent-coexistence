# Scripts for figures
library(tidyverse)
library(cowplot)
library(officer)
library(tidygraph)
library(ggraph)
library(gridExtra)
library(officer)
library(flextable)
source(here::here("plotting_scripts/network_functions.R"))

sequences_abundance <- read_csv(here::here("data/output/sequences_abundance.csv"), col_types = cols())
isolates <- read_csv(here::here("data/output/isolates.csv"), col_types = cols())
pairs <- read_csv(here::here("data/output/pairs.csv"), col_types = cols()) %>% mutate(InteractionType = ifelse(InteractionType == "neutrality", "coexistence", InteractionType))
pairs_freq <- read_csv(here::here("data/output/pairs_freq.csv"), col_types = cols())
pairs_example_outcomes_finer <- read_csv(here::here("data/output/pairs_example_outcomes_finer.csv"), col_types = cols())
communities <- read_csv(here::here("data/output/communities.csv"), col_types = cols()) %>%
    filter(Assembly == "self_assembly") %>%
    arrange(CommunitySize) %>%
    mutate(CommunityLabel = 1:13) %>%
    select(Community, CommunityLabel, CommunitySize, CommunityPairSize)
community_factor <- communities$Community
communities_size <- communities$CommunitySize
communities_hierarchy <- read_csv(here::here("data/output/communities_hierarchy.csv"), col_types = cols()) %>% right_join(select(communities, Community))
load("~/Dropbox/lab/invasion-network/data/output/network.Rdata")
net_list <- net_list %>% `[`(as.character(community_factor)) # Re-order the net_list according to the communities


# Figure 1 ----
p <- ggdraw() + draw_image(here::here("plots/cartoons/Fig1.png")) + theme(plot.background = element_rect(fill = "white", color = NA))
ggsave(here::here("plots/Fig1.png"), p, width = 27, height = 15)


# Figure 2 ----
# Figure 2A: networks of one community. C7R1
## Main network
net <- net_list$C7R1 %>%
    activate(nodes) %>%
    mutate(x = c(0, 0, 1, 1), y = c(1, 0, 1, 0))
node_size = 15
p1 <- net %>%
    ggraph(layout = "nicely") +
    geom_edge_link(aes(color = InteractionType), width = 2,
                   arrow = arrow(length = unit(node_size/2-1, "mm"), type = "closed", angle = 30, ends = "last"),
                   start_cap = circle(node_size/2+1, "mm"),
                   end_cap = circle(node_size/2+1, "mm")) +
    scale_edge_color_manual(values = interaction_color) +
    scale_x_continuous(limits = c(-0.4, 1.4), breaks = c(0, .5, 1)) +
    scale_y_continuous(limits = c(-0.4, 1.4), breaks = c(0, .5, 1)) +
    theme_void() +
    theme(legend.position = "none",
          panel.background = element_blank(),
          plot.margin=unit(c(3,3,3,3),"mm"),
          plot.background = element_rect(fill = NA, color = NA)) +
    labs() +
    draw_image(here::here("plots/cartoons/Fig1B_1.png"), x = -0.5, y = 0.75, vjust = 0.25, hjust = 0, clip = "on", scale = .3) +
    draw_image(here::here("plots/cartoons/Fig1B_2.png"), x = -0.5, y = -0.25, vjust = 0.25, hjust = 0, clip = "on", scale = .3) +
    draw_image(here::here("plots/cartoons/Fig1B_3.png"), x = 0.5, y = 0.75, vjust = 0.25, hjust = 0, clip = "on", scale = .3) +
    draw_image(here::here("plots/cartoons/Fig1B_4.png"), x = 0.5, y = -0.25, vjust = 0.25, hjust = 0, clip = "on", scale = .3)

## network legend
plot_network_legend <- function(net) {
    net %>%
    ggraph(layout = "nicely") +
    geom_edge_link(aes(color = InteractionType), width = 2,
                   arrow = arrow(length = unit(2, "mm"), type = "closed", angle = 30, ends = "last")) +
    scale_edge_color_manual(values = interaction_color) +
    theme_void() +
    theme(legend.key.size = unit(3,"line"),
          legend.text = element_text(size = 12),
          legend.position = c(0.5, 0.5),
          legend.title = element_blank(),
          legend.direction = "vertical",
          legend.background = element_blank())
}
p_temp1 <- net %>%
    activate(edges) %>%
    filter(InteractionType == "exclusion") %>%
    plot_network_legend()
p_temp2 <- net %>%
    activate(edges) %>%
    filter(InteractionType == "coexistence") %>%
    plot_network_legend()
p_legend_network <- plot_grid(get_legend(p_temp1), get_legend(p_temp2), nrow = 2, align = "v", axis = "l")

## frequency plots
pairs_example_freq <- pairs %>%
    filter(Community == "C7R1") %>%
    select(Community, starts_with("Isolate"), starts_with("Interaction")) %>%
    mutate(Pair = 1:n()) %>%
    left_join(pairs_freq, by = c("Community", "Isolate1", "Isolate2")) %>%
    mutate(Isolate1InitialODFreq = factor(Isolate1InitialODFreq)) %>%
    mutate(InteractionType = factor(InteractionType, c("exclusion", "coexistence")))
plot_example_freq <- function(pairs_freq) {
    pairs_freq %>%
        mutate(Isolate1InitialODFreq = factor(Isolate1InitialODFreq, c(95,50,5))) %>%
        ggplot(aes(x = Time, y = Isolate1MeasuredFreq, color = Isolate1InitialODFreq, group = Isolate1InitialODFreq)) +
        geom_point(size = 2) +
        geom_line(size = 1) +
        scale_y_continuous(breaks = c(0, .5, 1), limits = c(0,1)) +
        scale_x_discrete(labels = c(0,8)) +
        scale_color_manual(values = frequency_color, label = c("95%", "50%", "5%")) +
        facet_wrap(.~Pair) +
        theme_bw() +
        theme(panel.spacing = unit(2, "mm"), strip.text.x = element_blank(),
              panel.border = element_rect(color = 1, fill = NA, size = 1),
              panel.grid.minor.y = element_blank(),
              #axis.title = element_text(size = 10), axis.text = element_text(color = 1, size = 8),
              axis.title = element_blank(), axis.text = element_blank(),
              axis.ticks = element_blank(),
              panel.background = element_rect(fill = "white"),
              plot.background = element_blank()) +
        guides(color = "none") +
        labs(x = "Time (days)", y = "Frequency")
}
p_pairs_example_freq_list <- pairs_example_freq %>% group_split(Pair) %>% lapply(plot_example_freq)
p_pairs_example_freq_list[[1]] <- p_pairs_example_freq_list[[1]] + theme(axis.text = element_text(size = 10), axis.title = element_text(size = 10))
p_legend_color <- {p_pairs_example_freq_list[[1]] +
        theme(legend.background = element_blank(), legend.text = element_text(size = 12), legend.title = element_text(size = 12)) +
        guides(color = guide_legend(title = "Initial frequency"))} %>%
    cowplot::get_legend()
ss <- .2
pA <- ggdraw(p1) +
    draw_plot(p_pairs_example_freq_list[[1]], x = .05, y = .5, width = ss*1.5, height = ss*1.5, hjust = .5, vjust = .5) +
    draw_plot(p_pairs_example_freq_list[[2]], x = .5, y = .85, width = ss, height = ss, hjust = .5, vjust = .5) +
    draw_plot(p_pairs_example_freq_list[[3]], x = .4, y = .6, width = ss, height = ss, hjust = .5, vjust = .5) +
    draw_plot(p_pairs_example_freq_list[[4]], x = .6, y = .6, width = ss, height = ss, hjust = .5, vjust = .5) +
    draw_plot(p_pairs_example_freq_list[[5]], x = .5, y = .15, width = ss, height = ss, hjust = .5, vjust = .5) +
    draw_plot(p_pairs_example_freq_list[[6]], x = .85, y = .5, width = ss, height = ss, hjust = .5, vjust = .5) +
    draw_plot(p_legend_color, x = .15, y = .9, width = ss/2, height = ss/2, hjust = 0.5, vjust = .5) +
    draw_plot(p_legend_network, x = .15, y = .15, width = ss/2, height = ss/2, hjust = 0.5, vjust = .5) +
    theme(panel.background = element_blank(), plot.background = element_rect(color = NA, fill = "white"),
          plot.margin = unit(c(10,0,0,10), "mm"))

#ggsave(here::here("plots/Fig2A-example_network.png"), pA, width = 5, height = 5)


# Figure 2B: All 13 self-assembled community graphs
plot_competitive_network_grey <- function(x, node_size, edge_width){
    plot_competitive_network(x, node_size = node_size, edge_width = edge_width) +
        theme(plot.background = element_rect(fill = "grey90", color = NA),
              panel.background = element_rect(fill = "grey90", color = NA))

}
p_net_list <- communities %>%
    mutate(Network = net_list[1:13]) %>%
    mutate(NetworkPlotSize = max(CommunitySize) / CommunitySize / 4) %>%
    rowwise() %>%
    mutate(p_net = plot_competitive_network_grey(Network, 0, NetworkPlotSize) %>% list()) %>%
    pull(p_net)
p1 <- plot_grid(plotlist = p_net_list, nrow = 1, scale = 1.3) + paint_white_background()

## pairwise outcomes per community
p2 <- pairs %>%
    filter(Assembly == "self_assembly") %>%
    group_by(Community, InteractionType) %>%
    count(name = "Count") %>%
    group_by(Community) %>% mutate(Fraction = Count / sum(Count)) %>% ungroup() %>%
    mutate(Community = factor(Community, community_factor)) %>%
    arrange(Community) %>%
    mutate(CommunityLabel = rep(1:13, each = 2) %>% factor()) %>%
    ggplot() +
    geom_col(aes(x = CommunityLabel, fill = InteractionType, y = Fraction), color = 1, width = .8, size = .5) +
    geom_text(data = communities, aes(x = CommunityLabel, y = .1, label = paste0("n=", CommunityPairSize)), vjust = -.5, size = 3) +
    scale_fill_manual(values = assign_interaction_color()) +
    scale_y_continuous(breaks = c(0,.5,1), limit = c(0, 1), expand = c(0,0)) +
    theme_classic() +
    theme(legend.text = element_text(size = 12),
          axis.text = element_text(color = 1, size = 12),
          axis.title = element_text(color = 1, size = 12),
          legend.title = element_blank(),
          legend.position = "top") +
    guides(fill = guide_legend(reverse = T)) +
    labs(x = "Community", y = "Fraction")


pB <- plot_grid(p1, p2, ncol = 1, scale = .9, rel_heights = c(1, 4), axis = "lr", align = "v") + paint_white_background()
#ggsave(here::here("plots/Fig2B-all_networks.png"), pB, width = 8, height = 3)

#
p <- plot_grid(pA, pB, nrow = 1, labels = c("A", "B"), rel_widths = c(1, 2), axis = "tr", align = "h") + paint_white_background()
ggsave(here::here("plots/Fig2.png"), p, width = 12, height = 4)



# Figure 3 ----
pairs_interaction_finer <- pairs %>%
    filter(Assembly == "self_assembly") %>%
    mutate(InteractionType = factor(InteractionType, c("exclusion", "coexistence"))) %>%
    mutate(InteractionTypeFiner = factor(InteractionTypeFiner, c("competitive exclusion", "mutual exclusion", "stable coexistence", "frequency-dependent coexistence", "neutrality"))) %>%
    group_by(InteractionType, InteractionTypeFiner) %>%
    count(name = "Count") %>% ungroup() %>%
    mutate(Fraction = Count / sum(Count)) %>%
    mutate(Label = str_replace(InteractionTypeFiner, " ", "\n"))

## Legends
temp <- pairs %>%
    ggplot() +
    geom_tile(aes(x = Isolate1, y = Isolate2, fill = InteractionTypeFiner), height = .8, width = .8, alpha = .9) +
    scale_fill_manual(values = assign_interaction_color(level = "finer"),
                      breaks = c("competitive exclusion", "mutual exclusion", "stable coexistence", "frequency-dependent coexistence", "neutrality"),
                      labels = paste0(pairs_interaction_finer$InteractionTypeFiner, " (", round(pairs_interaction_finer$Fraction, 3) * 100,"%)")) +
    theme(legend.title = element_blank(),
          legend.position = "right",
          legend.spacing.y = unit("2", "mm"),
          legend.text = element_text(size = 12)
    ) +
    guides(fill = guide_legend(byrow = T)) +
    paint_white_background()
p_legend_fill <- get_legend(temp)


## Plot the waffle
pairs_example_freq <- pairs %>%
    filter(Assembly == "self_assembly") %>%
    mutate(InteractionType = factor(InteractionType, c("exclusion", "coexistence"))) %>%
    mutate(InteractionTypeFiner = ordered(InteractionTypeFiner, c("competitive exclusion", "mutual exclusion", "stable coexistence", "frequency-dependent coexistence", "neutrality"))) %>%
    arrange(InteractionType, InteractionTypeFiner) %>%
    # Add the coordinate in the grid
    bind_rows(tibble(InteractionType = rep(NA, 4), InteractionTypeFiner = rep(NA, 4))) %>%
    #mutate(x = rep(1:19, each = 10), y = rep(1:10, 19), PlotOrder = paste0(x, "_", y)) %>%
    filter(!is.na(tibble(InteractionType))) %>%
    mutate(PairID = factor(1:n())) %>%
    # Join the frequency data
    left_join(pairs_freq, by = c("Community", "Isolate1", "Isolate2")) %>%
    select(PairID, InteractionType, InteractionTypeFiner, Isolate1InitialODFreq, Time, Isolate1MeasuredFreq,
           Community, Isolate1, Isolate2) %>%
    mutate(Time = str_replace(Time, "T", "")) %>%
    mutate(Isolate1InitialODFreq = factor(Isolate1InitialODFreq))

plot_example_freq <- function(pairs_freq) {
    # Extract params
    comm <- unique(pairs_freq$Community)
    isolate1 <- unique(pairs_freq$Isolate1)
    isolate2 <- unique(pairs_freq$Isolate2)
    interaction_type_finer <- pairs %>%
        filter(Community == comm, Isolate1 == isolate1, Isolate2 == isolate2) %>%
        pull(InteractionTypeFiner)

    #
    pairs_freq %>%
        mutate(Isolate1InitialODFreq = factor(Isolate1InitialODFreq, c(95,50,5))) %>%
        ggplot() +
        geom_rect(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, aes(fill = InteractionTypeFiner), alpha = .1) +
        geom_point(size = 2, aes(x = Time, y = Isolate1MeasuredFreq, color = Isolate1InitialODFreq, group = Isolate1InitialODFreq)) +
        geom_line(size = 1, aes(x = Time, y = Isolate1MeasuredFreq, color = Isolate1InitialODFreq, group = Isolate1InitialODFreq)) +
        scale_y_continuous(breaks = c(0, .5, 1), limits = c(0,1)) +
        scale_color_manual(values = frequency_color, label = c("95%", "50%", "5%")) +
        scale_fill_manual(values = assign_interaction_color(level = "finer")) +
        theme_bw() +
        theme(panel.spacing = unit(2, "mm"),
              panel.border = element_rect(color = 1, fill = NA, size = 1),
              panel.grid = element_blank(),
              panel.grid.minor.y = element_blank(),
              axis.title = element_blank(), axis.text = element_blank(),
              axis.ticks = element_blank(),
              plot.background = element_blank(),
              plot.margin = margin(0,0,0,0, "mm")) +
        guides(color = "none", fill = "none") +
        labs(x = "Time", y = "Frequency")
}
temp_list <- pairs_example_freq %>%
    as_tibble() %>%
    arrange(PairID) %>%
    group_split(PairID) %>%
    lapply(plot_example_freq)

## Grid layout
m <- matrix(c(1:186, rep(NA, 4)), nrow = 10)
p_waffle <- arrangeGrob(grobs = temp_list, layout_matrix = m)
p_waffle <- plot_grid(p_waffle, NULL, rel_widths = c(3, 1), scale = c(.9, 1)) + paint_white_background()

#p_right <- plot_grid(p_legend_fill, p_legend_color, nrow = 2)
ss = .3
p <- ggdraw(p_waffle) +
    draw_plot(p_legend_fill, x = .86, y = .7, width = ss/2, height = ss/2, hjust = 0.5, vjust = .5) +
    draw_plot(p_legend_color, x = .77, y = .2, width = ss/2, height = ss/2, hjust = 0.5, vjust = .5) +
    theme(panel.background = element_blank(), plot.background = element_rect(color = NA, fill = "white"),
          plot.margin = unit(c(0,0,0,0), "mm"))

ggsave(here::here("plots/Fig3.png"), p, width = 13, height = 5)





# Figure 4 ----
# Figure 4A One example network
node_size = 5
edge_width = 1.5
p_net <- net_list$C7R1 %>%
    activate(nodes) %>%
    mutate(y = -Rank) %>%
    group_by(Rank) %>%
    mutate(x = {seq(0, 1, length.out = n() + 2) %>% `[`(c(-1, -length(.)))}) %>%
    ungroup() %>%
    ggraph(layout = "nicely") +
    geom_edge_link(aes(color = InteractionType), width = edge_width,
                  arrow = arrow(length = unit(edge_width, "mm"), type = "closed", angle = 30, ends = "last"),
                  start_cap = circle(node_size/2+1, "mm"),
                  end_cap = circle(node_size/2+1, "mm")) +
    scale_edge_color_manual(values = assign_interaction_color(level = "matrix"),
                            breaks = c("exclusion", "exclusion violating rank"),
                            labels = c("exclusion following rank", "exclusion violating rank")) +
    scale_x_continuous(limits = c(.2, .8), expand = c(0,0)) +
    scale_y_continuous(limits = c(-6, 0), breaks = -4:-1, labels = 4:1) +
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
    draw_image(here::here("plots/cartoons/Fig1B_1.png"), x = 0, y = -1, vjust = 0.5, hjust = 0, clip = "on", scale = .8) +
    draw_image(here::here("plots/cartoons/Fig1B_3.png"), x = -.17, y = -2, vjust = 0.5, hjust = 0, clip = "on", scale = .8) +
    draw_image(here::here("plots/cartoons/Fig1B_4.png"), x = .17, y = -2, vjust = 0.5, hjust = 0, clip = "on", scale = .8) +
    draw_image(here::here("plots/cartoons/Fig1B_2.png"), x = 0, y = -4, vjust = .6, hjust = 0, clip = "on", scale = .8) +
    draw_plot(p_legend_network, x = -0.1, y = -6.2, height = 1.5) +
    paint_white_background()

#ggsave(here::here("plots/Fig4A-example.png"), pA, width = 3, height = 3)


# Figure 4B. Hierarchy network plot
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
        filter(InteractionType != "coexistence") %>%
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
        scale_edge_color_manual(values = assign_interaction_color(level = "matrix"),
                                breaks = c("exclusion", "exclusion violating rank"),
                                labels = c("exclusion pair that follows rank", "exclusion pair that violates rank")) +
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
communities_net <- communities %>%
    mutate(Community = factor(Community, community_factor)) %>%
    arrange(Community) %>%
    left_join(as_tibble_col(net_list, column_name = "Network") %>% mutate(Community = names(net_list))) %>%
    select(Community, CommunitySize, Network) %>%
    rowwise() %>%
    mutate(p_net = plot_network_hierarchy(Network, n_rank = CommunitySize) %>% list())

p_net_hierarchy_list <- communities_net %>% pull(p_net)
p_net_hierarchy_list[[13]] <- p_net_hierarchy_list[[13]] +
    scale_y_continuous(limits = c(-12-1, 0), breaks = -12:-1, labels = 12:1, position = "right") +
    theme(axis.title.y = element_text(color = 1, size = 10, angle = 270, margin = margin(l = 2, unit = "mm")),
          axis.text.y = element_text(color = 1, size = 10, margin = margin(l = 1, unit = "mm")))
pB_axistitle <- ggdraw() + draw_label("Community", fontface = 'plain', x = .5, hjust = .5) + theme(plot.margin = margin(5, 0, 5, 7))
p_temp <- plot_grid(plotlist = p_net_hierarchy_list,
                    rel_widths = c(communities_net$CommunitySize / max(communities_net$CommunitySize))^1.5,
                    labels = 1:13, label_x = c(rep(0.5, 12), 0.45), hjust = c(rep(.5, 12), 1),
                    nrow = 1, axis = "tb", align = "h") + paint_white_background()
pB <- plot_grid(pB_axistitle, p_temp, ncol = 1, rel_heights = c(.1, 1)) + paint_white_background()

## legend
temp <- p_net_hierarchy_list[[13]] +
    geom_edge_arc(strength = 10,
                  aes(color = InteractionType), width = edge_width*2.5,
                  arrow = arrow(length = unit(edge_width*2.5, "mm"), type = "closed", angle = 30, ends = "last"),
                  start_cap = circle(node_size/2, "mm"),
                  end_cap = circle(node_size/2, "mm")) +
    theme(legend.key.size = unit(3, "line"),
          legend.key.height = unit(6, "mm"),
          legend.position = "right",
          legend.direction = "vertical",
          legend.text = element_text(size = 12),
          legend.background = element_rect(fill = NA, color = NA))
p_legend <- get_legend(temp)
pB <- ggdraw(pB) + draw_plot(p_legend,.2,.3,.1,.1)

#ggsave(here::here("plots/Fig4B-network_hierarchy_experiment.png"), pB, width = 10, height = 3)



# Figure 4C: Hierarchy
pC <- mutate(communities_hierarchy, Treatment = "experiment") %>%
    filter(Metric == "h1") %>%
    ggplot(aes(x = Treatment, y = HierarchyScore)) +
    geom_boxplot(width = .5, lwd = .8, outlier.color = NA) +
    geom_jitter(shape = 1, size = 2, stroke = .8, height = 0, width = .1) +
    scale_y_continuous(limits = c(0,1.01), breaks = c(0, .25, .5, .75, 1)) +
    theme_classic() +
    theme(panel.grid.major.y = element_line(color = "grey", linetype = 2),
          panel.spacing = unit(0, "mm"),
          panel.border = element_rect(fill = NA, color = 1, size = 1.5),
          axis.text = element_text(size = 10, color = 1),
          axis.text.x = element_blank(),
          axis.title = element_text(size = 10, color = 1),
          legend.title = element_blank(),
          plot.title = element_text(size = 10, color = 1)
    ) +
    labs(x = "", y = "Score") +
    guides(color = "none") +
    ggtitle("Hierarchy")

#ggsave(here::here("plots/Fig4C-hierarchy.png"), pC, width = 4, height = 4)

#
p_left <- plot_grid(pA, pC, ncol = 1, rel_heights = c(1,1), scale = c(.8, .9), labels = c("A", "C"), axis = "lr", align = "v")
p <- plot_grid(p_left, pB, nrow = 1, rel_widths = c(1,4), scale = c(1, .9), labels = c("", "B")) + paint_white_background()
ggsave(here::here("plots/Fig4.png"), p, width = 10, height = 5)








# Supplement
# Figure S1 ----
## Panel A. cartoon of self-assembly experiment and isolate abundance
p1 <- ggdraw() + draw_image(here::here("plots/cartoons/FigS1.png")) + theme(plot.background = element_rect(fill = "white", color = NA))
## Panel B. isolate abundance
temp <- sequences_abundance %>%
    filter(AlignmentType == "local") %>%
    filter(AllowMismatch == Inf) %>%
    filter(BasePairMismatch <= 4) %>%
    mutate(Community = ordered(Community,  communities$Community)) %>%
    group_by(Community) %>%
    mutate(RankRelativeAbundance = rank(-RelativeAbundance))

isolates_abundance <- isolates %>%
    filter(Assembly == "self_assembly") %>%
    left_join(temp)
color_sets <- tibble(Color = c("yellow", "deepskyblue3", "blue", "darkorchid2", "firebrick", "orange2", "grey"),
                     Family = c("Aeromonadaceae", "Enterobacteriaceae", "Moraxellaceae", "Pseudomonadaceae","Comamonadaceae","Alcaligenaceae", "Sphingobacteriaceae"))
p2 <- isolates_abundance %>%
    mutate(Community = factor(Community, community_factor)) %>%
    arrange(Community) %>%
    ggplot() +
    geom_bar(aes(x = Community, y = RelativeAbundance, fill = Family), size = .3, color = "grey30", position = "stack", stat = "identity") +
    theme_bw() +
    scale_fill_manual(values = setNames(color_sets$Color, color_sets$Family)) +
    scale_x_discrete(labels = 1:13) +
    scale_y_continuous(breaks = c(0, .5, 1), expand = c(0,0), limits = c(0, 1)) +
    theme(#axis.text.x = element_text(angle = 90, vjust = 0.5),
          panel.grid = element_blank(),
          axis.title = element_text(size = 12),
          panel.border = element_rect(color = 1, size = 1)) +
    labs(y = "Relative abundance")

## Stats
temp %>%
    group_by(Community) %>%
    summarize(Total = sum(RelativeAbundance)) %>%
    summarize(Mean = mean(Total))

p <- plot_grid(p1, p2, nrow = 1, scale = c(.8, .9), rel_widths = c(1, 1.5), labels = c("A", "B")) + paint_white_background()
ggsave(here::here("plots/FigS1.png"), p, width = 10, height = 3)



# Figure S2 ----
# Rank versus ranked abundance
p1 <- isolates_abundance %>%
    mutate(Fermenter = ifelse(Fermenter, "fermenter", "respirator")) %>%
    geom_smooth(aes(x = Rank, y = RelativeAbundance), method = "lm") +
    ggplot() +
    geom_point(aes(x = Rank, y = RelativeAbundance, color = Fermenter, fill = Fermenter),
               shape = 21, size = 3, stroke = 0, alpha = 0.7) +
    scale_color_manual(values = category_color, breaks = c("fermenter", "respirator")) +
    scale_fill_manual(values = category_color, breaks = c("fermenter", "respirator")) +
    scale_x_continuous(breaks = 1:12) +
    scale_y_continuous(limits = c(0,1)) +
    theme_classic() +
    theme(legend.position = "top", legend.title = element_blank()) +
    labs(x = "Rank", y = "Relative abundance")

# Rank versus ranked abundance
p2 <- isolates_abundance %>%
    mutate(Fermenter = ifelse(Fermenter, "fermenter", "respirator")) %>%
    #group_by(Rank, RankRelativeAbundance) %>%
    #count(name = "Count") %>%
    ggplot() +
    geom_smooth(aes(x = Rank, y = RankRelativeAbundance), method = "lm") +
    geom_point(aes(x = Rank, y = RankRelativeAbundance, color = Fermenter, fill = Fermenter),
               shape = 21, size = 3, stroke = 0, alpha = 0.7,
               position = position_jitter(width = .1, height = .1)) +
    scale_color_manual(values = category_color, breaks = c("fermenter", "respirator")) +
    scale_fill_manual(values = category_color, breaks = c("fermenter", "respirator")) +
    scale_x_continuous(breaks = 1:12) +
    scale_y_continuous(breaks = 1:12) +
    theme_classic() +
        theme(legend.position = "top", legend.title = element_blank()) +
    labs(x = "Rank", y = "Ranked relative abundance")

#
isolates_abundance %>%
    select(Rank, RelativeAbundance) %>%
    glm(RelativeAbundance ~ Rank, data = .) %>%
    broom::tidy()

isolates_abundance %>%
    select(Rank, RankRelativeAbundance) %>%
    glm(RankRelativeAbundance ~ Rank, data = .) %>%
    broom::tidy()

#
p <- plot_grid(p1, p2, nrow = 1, labels = c("A", "B"))
#ggsave(here::here("plots/FigS2.png"), p, width = 6, height = 3)

#











# Table S1. Pairwise interaction tables ----
ft1 <- read_csv(here::here("data/output/pairs_interaction_table.csv")) %>%
    mutate(InteractionType = ifelse(InteractionType == "neutrality", "coexistence", InteractionType)) %>%
    arrange(InteractionType, InteractionTypeFiner, FromRare, FromMedium, FromAbundant) %>%
    setNames(c("From 5%", "From 50%", "From 95%", "Outcome", "Finer outcome", "Count")) %>%
    flextable() %>%
    width(j = 1:3, width = 1) %>%
    width(j = 5, width = 2.5)
save_as_image(ft1, here::here("plots/TableS1.png"))


# Table S2. Isolates ----
t2 <- isolates_abundance %>%
    left_join(select(communities, Community, CommunityLabel)) %>%
    select(CommunityLabel, Fermenter, Community, Isolate, Family, Genus, RelativeAbundance) %>%
    select(-Community) %>%
    mutate(Fermenter = ifelse(Fermenter, "fermenter", "respirator")) %>%
    group_by(CommunityLabel) %>%
    mutate(TotalAbundance = sum(RelativeAbundance, na.rm = T)) %>%
    mutate(RelativeAbundance = ifelse(is.na(RelativeAbundance), RelativeAbundance, paste0(round(RelativeAbundance*100, 2), "%"))) %>%
    mutate(TotalAbundance = paste0(round(TotalAbundance*100, 2), "%")) %>%
    select(-Fermenter) %>%
    #mutate(Sequence = str_replace_all(Sequence, '(?=(?:.{50})+$)', "\n")) %>%
    arrange(CommunityLabel, Isolate)

t2_1 <- t2 %>%
    rename(Community = CommunityLabel, `Relative abundance` = RelativeAbundance, `Total abundance` = TotalAbundance) %>%
    filter(Community %in% c(1:9))
border_rows1 <- t2_1 %>%
    group_by(Community) %>%
    mutate(temp = Isolate == max(Isolate)) %>%
    pull(temp) %>%
    which()
t2_2 <- t2 %>%
    rename(Community = CommunityLabel, `Relative abundance` = RelativeAbundance, `Total abundance` = TotalAbundance) %>%
    filter(Community %in% c(10:13))
border_rows2 <- t2_2 %>%
    group_by(Community) %>%
    mutate(temp = Isolate == max(Isolate)) %>%
    pull(temp) %>%
    which()



ft2_1 <- t2_1 %>%
    flextable() %>%
    width(j = 1:6, width = 1) %>%
    merge_v(j = c("Community", "Total abundance")) %>%
    border(i = border_rows1, border.bottom = fp_border()) %>%
    border(j = 1:6, border = fp_border(width = 2), part = "header") %>%
    border(j = 1:6, border = fp_border(), part = "body")
ft2_2 <- t2_2 %>%
    flextable() %>%
    width(j = 1:6, width = 1) %>%
    merge_v(j = c("Community", "Total abundance")) %>%
    border(i = border_rows2, border.bottom = fp_border()) %>%
    border(j = 1:6, border = fp_border(width = 2), part = "header") %>%
    border(j = 1:6, border = fp_border(), part = "body")


save_as_image(ft2_1, here::here("plots/TableS2_1.png"))
save_as_image(ft2_2, here::here("plots/TableS2_2.png"))




















