#' Script for figures

library(tidyverse)
library(cowplot)
library(tidygraph)
library(ggraph)
library(gridExtra)
source(here::here("analysis/00-metadata.R"))

isolates <- read_csv(paste0(folder_data, "output/isolates.csv"), show_col_types = F)
pairs <- read_csv(paste0(folder_data, "output/pairs.csv"), show_col_types = F)
communities <- read_csv(paste0(folder_data, "temp/00c-communities.csv"), show_col_types = F)
pairs_freq <- read_csv(paste0(folder_data, "temp/93-pairs_freq.csv"), show_col_types = F)
load(paste0(folder_data, "temp/95-communities_network.Rdata"))
communities_hierarchy <- read_csv(paste0(folder_data, "temp/95-communities_hierarchy.csv"), show_col_types = F)

# 0. Clean up column factors ----
# Arrange communities by size
communities <- communities %>%
    arrange(CommunitySize) %>%
    mutate(Community = factor(Community, Community))

# Clean up the pairs data ----
pairs <- pairs %>%
    # Remove no-colony pairs
    unite(col = "Pair", Community, Isolate1, Isolate2, sep = "_", remove = F) %>%
    filter(!(Pair %in% pairs_no_colony)) %>%
    # Remove low-accuracy model pairs
    filter(AccuracyMean > 0.9)

# Count pairwise competition outcomes ----
pairs %>%
    group_by(InteractionType) %>%
    count(name = "Count") %>%
    ungroup() %>%
    mutate(Fraction = Count / sum(Count))
pairs %>%
    group_by(InteractionType, InteractionTypeFiner) %>%
    count(name = "Count") %>%
    ungroup() %>%
    mutate(Fraction = Count / sum(Count)) %>%
    mutate(InteractionTypeFiner = factor(InteractionTypeFiner, interaction_type_finer)) %>%
    arrange(InteractionTypeFiner)


# Figure 1 ----
#p <- ggdraw() + draw_image(here::here("plots/cartoons/Fig1.pdf")) + paint_white_background()
#ggsave(here::here("plots/Fig1.png"), p, width = 27, height = 15)


# Figure 2 ----
# Figure 2A: cartoon
pA <- ggdraw() + draw_image(here::here("plots/cartoons/Fig2A.png")) + paint_white_background()

# Figure 2B: example network C2R6
## Main network
node_size = 15
net <- communities_network %>%
    filter(Community == "C2R6") %>%
    pull(Network) %>% `[[`(1) %>%
    activate(nodes) %>%
    mutate(x = c(0, 0, 1, 1), y = c(1, 0, 1, 0))
p1 <- net %>%
    ggraph(layout = "nicely") +
    geom_edge_link(aes(color = InteractionType), width = 1,
                   arrow = arrow(length = unit(node_size/2-3, "mm"), type = "closed", angle = 30, ends = "last"),
                   start_cap = circle(node_size-7, "mm"),
                   end_cap = circle(node_size-7, "mm")) +
    scale_edge_color_manual(values = interaction_color, breaks = c("coexistence", "exclusion")) +
    scale_x_continuous(limits = c(-0.2, 1.2), breaks = c(0, .5, 1)) +
    scale_y_continuous(limits = c(-0.2, 1.2), breaks = c(0, .5, 1)) +
    theme_void() +
    theme(legend.position = "none",
          legend.key.size = unit(1.3, "line"),
          panel.background = element_blank(),
          plot.margin = unit(c(0,0,0,0), "mm"),
          plot.background = element_rect(fill = NA, color = NA)) +
    draw_image(here::here("plots/cartoons/Fig2B_1.png"), x = -0.5, y = 0.75, vjust = 0.25, hjust = 0, clip = "on", scale = .3) +
    draw_image(here::here("plots/cartoons/Fig2B_2.png"), x = 0.5, y = 0.75, vjust = 0.25, hjust = 0, clip = "on", scale = .3) +
    draw_image(here::here("plots/cartoons/Fig2B_3.png"), x = -0.5, y = -0.25, vjust = 0.25, hjust = 0, clip = "on", scale = .3) +
    draw_image(here::here("plots/cartoons/Fig2B_4.png"), x = 0.5, y = -0.25, vjust = 0.25, hjust = 0, clip = "on", scale = .3)

## frequency plots
pairs_example_freq <- pairs %>%
    filter(Community == "C2R6") %>%
    select(Community, starts_with("Isolate"), starts_with("Interaction")) %>%
    mutate(Pair = 1:n()) %>%
    left_join(pairs_freq, by = c("Community", "Isolate1", "Isolate2")) %>%
    mutate(Isolate1InitialODFreq = factor(Isolate1InitialODFreq)) %>%
    mutate(InteractionType = factor(InteractionType, c("exclusion", "coexistence")))

plot_example_freq <- function(pairs_freq) {
    pairs_freq %>%
        mutate(Isolate1InitialODFreq = factor(Isolate1InitialODFreq, c(95,50,5))) %>%
        ggplot(aes(x = Time, y = Isolate1CFUFreqMean, color = Isolate1InitialODFreq, group = Isolate1InitialODFreq)) +
        geom_line(size = 1) +
        geom_point(size = 1.5) +
        geom_segment(size = 1, aes(x = Time, xend = Time,
                                    y = Isolate1CFUFreqMean + 2*Isolate1CFUFreqSd,
                                    yend = Isolate1CFUFreqMean - 2*Isolate1CFUFreqSd,
                                    color = Isolate1InitialODFreq)) +
        scale_y_continuous(breaks = c(0, .5, 1), limits = c(-0.2, 1.2)) +
        scale_x_discrete(labels = c("start", "end")) +
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
        labs(x = "Time", y = "Frequency")
}
p_pairs_example_freq_list <- pairs_example_freq %>% group_split(Pair) %>% lapply(plot_example_freq)
p_pairs_example_freq_list[[1]] <- p_pairs_example_freq_list[[1]] + theme(axis.text = element_text(size = 8), axis.title = element_text(size = 10))
p_legend_color <- get_legend({p_pairs_example_freq_list[[1]] +
        theme(legend.background = element_blank(), legend.text = element_text(size = 10), legend.title = element_blank()) +
        guides(color = guide_legend(title = "Initial frequency"))})
ss <- .25
pB <- ggdraw(p1) +
    draw_plot(p_pairs_example_freq_list[[1]], x = -.02, y = .5, width = ss*1.5, height = ss*1.5, hjust = .5, vjust = .5) +
    draw_plot(p_pairs_example_freq_list[[2]], x = .5, y = .85, width = ss, height = ss*1.1, hjust = .5, vjust = .5) +
    draw_plot(p_pairs_example_freq_list[[3]], x = .35, y = .58, width = ss, height = ss*1.1, hjust = .5, vjust = .5) +
    draw_plot(p_pairs_example_freq_list[[4]], x = .65, y = .58, width = ss, height = ss*1.1, hjust = .5, vjust = .5) +
    draw_plot(p_pairs_example_freq_list[[5]], x = .5, y = .15, width = ss, height = ss*1.1, hjust = .5, vjust = .5) +
    draw_plot(p_pairs_example_freq_list[[6]], x = .95, y = .55, width = ss, height = ss*1.1, hjust = .5, vjust = .5) +
    #draw_plot(p_legend_color, x = .05, y = .9, width = ss/2, height = ss/2, hjust = 0.5, vjust = .5) +
    theme(panel.background = element_blank(), plot.background = element_rect(color = NA, fill = "white"),
          plot.margin = unit(c(0,3,3,5), "mm"))

# Figure 2C: All 13 self-assembled community graphs
plot_competitive_network_grey <- function(g, node_size = 10, edge_width = 1){
    # Layout
    graph_layout <- create_layout(g, "circle")
    mean_x_coord <- mean(graph_layout$x)
    mean_y_coord <- mean(graph_layout$x)
    g <- g %>% activate(nodes) %>% mutate(x = graph_layout$x - mean_x_coord, y = graph_layout$y - mean_y_coord)

    # Axis range
    nodes_axis_x <- activate(g, nodes) %>% pull(x) %>% range()
    nodes_axis_y <- activate(g, nodes) %>% pull(y) %>% range()

    # Graph
    g %>%
        ggraph(layout = "nicely") +
        geom_node_point(fill = "grey", size = node_size, shape = 21, colour = "black", stroke = node_size/5) +
        geom_edge_link(aes(color = InteractionType), width = edge_width) +
        scale_edge_color_manual(values = interaction_color) +
        scale_x_continuous(limits = nodes_axis_x*1.3) +
        scale_y_continuous(limits = nodes_axis_y*1.3) +
        theme_graph() +
        theme(
            legend.position = "none",
            legend.title = element_blank(),
            strip.text = element_blank(),
            plot.margin = unit(c(3,3,3,3), "mm"),
            plot.background = element_rect(fill = "grey90", color = NA),
            panel.background = element_rect(fill = "grey90", color = NA)
        )
}
p_net_list <- communities_network %>%
    mutate(Community = factor(Community, Community)) %>%
    arrange(CommunitySize) %>%
    mutate(NetworkPlotSize = max(CommunitySize) / CommunitySize / 4) %>%
    rowwise() %>%
    mutate(NetworkPlot = plot_competitive_network_grey(Network, 0, NetworkPlotSize) %>% list()) %>%
    pull(NetworkPlot)
p1 <- plot_grid(plotlist = p_net_list, nrow = 1, scale = 1.3) + paint_white_background()

## pairwise outcomes per community
p2 <- pairs %>%
    filter(!is.na(FitnessFunction)) %>%
    group_by(Community, InteractionType) %>%
    count(name = "Count") %>%
    group_by(Community) %>% mutate(Fraction = Count / sum(Count), TotalCount = sum(Count)) %>%
    ungroup() %>%
    mutate(Community = factor(Community, communities$Community)) %>%
    arrange(Community) %>%
    replace_na(list(InteractionType = "unknown")) %>%
    ggplot() +
    geom_col(aes(x = Community, fill = InteractionType, y = Fraction), color = 1, width = .8, size = .5) +
    geom_text(aes(x = Community, y = .9, label = paste0("n=", TotalCount))) +
    scale_fill_manual(values = assign_interaction_color(), breaks = c("coexistence", "exclusion", "unknown")) +
    scale_y_continuous(breaks = c(0,.5,1), limit = c(0, 1), expand = c(0,0)) +
    theme_classic() +
    theme(legend.text = element_text(size = 12),
          axis.text = element_text(color = 1, size = 12),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          axis.title = element_text(color = 1, size = 12),
          legend.title = element_blank(),
          legend.position = "top") +
    labs(x = "Community", y = "Fraction")

pC <- plot_grid(NULL, p1, p2, ncol = 1, scale = .9, rel_heights = c(0.3, 1, 5), axis = "lr", align = "v") + paint_white_background()

#
p_bottom <- plot_grid(pB, pC, nrow = 1, labels = c("B", "C"), scale = c(0.9, 1), rel_widths = c(1, 1.8), axis = "t", align = "h")
p <- plot_grid(pA, p_bottom, nrow = 2, labels = c("A", ""), scale = c(.95, .95), rel_heights = c(.8, 1)) + paint_white_background()
ggsave(here::here("plots/Fig2.png"), p, width = 12, height = 7)



# Figure 3 ----
get_interaction_legend <- function (pairs) {
    panel_colors <- c(interaction_color[1], rep(interaction_color[2], 6), interaction_color[3]) %>% setNames(pairs_interaction_finer$InteractionTypeFiner)
    panel_fills <- assign_interaction_color(level = "finer")[-3]
    temp <- pairs %>%
        mutate(InteractionType = factor(InteractionType)) %>%
        mutate(InteractionTypeFiner = factor(InteractionTypeFiner)) %>%
        ggplot() +
        geom_tile(aes(x = Isolate1, y = Isolate2, color = InteractionType, fill = InteractionTypeFiner), height = .8, width = .8, alpha = .9, size = 1) +
        scale_color_manual(values = panel_colors,
                           labels = paste0(pairs_interaction_finer$InteractionTypeFiner, " (", round(pairs_interaction_finer$Fraction, 3) * 100,"%)")) +
        scale_fill_manual(values = panel_fills,
                          labels = paste0(pairs_interaction_finer$InteractionTypeFiner, " (", round(pairs_interaction_finer$Fraction, 3) * 100,"%)")) +
        theme(legend.position = "right",
              legend.spacing.y = unit("2", "mm"),
              legend.text = element_text(size = 12),
              legend.background = element_blank()) +
        labs(color = "", fill = "") +
        paint_white_background()
    return(get_legend(temp))
}
plot_example_freq <- function(pairs_freq) {
    axis_lower <- -0.5
    axis_upper <- 1.5
    pairs_freq %>%
        mutate(Time = case_when(Time == "T0" ~ 0,Time == "T8" ~ 1)) %>%
        mutate(Isolate1InitialODFreq = factor(Isolate1InitialODFreq, c(95,50,5))) %>%
        ggplot() +
        geom_rect(xmin = axis_lower, xmax = axis_upper, ymin = axis_lower+0.2, ymax = axis_upper-0.2, aes(color = InteractionType, fill = InteractionTypeFiner), alpha = .08, size = .8) +
        geom_hline(size = .2, yintercept = c(0,1), linetype = 1, color = "grey90") +
        geom_line(size = .4, aes(x = Time, y = Isolate1CFUFreqMean, color = Isolate1InitialODFreq, group = Isolate1InitialODFreq)) +
        geom_point(size = .2, aes(x = Time, y = Isolate1CFUFreqMean, color = Isolate1InitialODFreq, group = Isolate1InitialODFreq)) +
        geom_segment(size = .4, aes(x = Time, xend = Time,
                                    y = Isolate1CFUFreqMean + 2*Isolate1CFUFreqSd,
                                    yend = Isolate1CFUFreqMean - 2*Isolate1CFUFreqSd,
                                    color = Isolate1InitialODFreq)) +
        scale_x_continuous(expand = c(0, 0), limits = c(axis_lower, axis_upper)) +
        scale_y_continuous(expand = c(0, 0), breaks = c(0, .5, 1), limits = c(axis_lower+0.2, axis_upper-0.2)) +
        scale_color_manual(values = c(frequency_color, interaction_color), label = c("95%", "50%", "5%", "exclusion", "coexistence")) +
        scale_fill_manual(values = assign_interaction_color(level = "finer")) +
        theme_classic() +
        theme(panel.spacing = unit(2, "mm"),
              panel.border = element_blank(),
              panel.grid.minor.y = element_blank(),
              axis.line = element_blank(),
              axis.title = element_blank(), axis.text = element_blank(),
              axis.ticks = element_blank(),
              plot.background = element_blank(),
              plot.title = element_text(size = 5, margin = margin(0,0,0,0)),
              plot.margin = margin(0,0,0,0, "mm")) +
        guides(color = "none", fill = "none") +
        labs(x = "Time", y = "Frequency") +
        #ggtitle(unique(pairs_freq$PairID)) +
        NULL
}
plot_example_freq2 <- function(pairs_freq) {
    axis_lower <- -0.5
    axis_upper <- 1.5
    pairs_freq %>%
        mutate(Time = case_when(Time == "T0" ~ 0,Time == "T8" ~ 1)) %>%
        mutate(Isolate1InitialODFreq = factor(Isolate1InitialODFreq, c(95,50,5))) %>%
        ggplot() +
        geom_rect(xmin = axis_lower, xmax = axis_upper, ymin = axis_lower+0.2, ymax = axis_upper-0.2, aes(color = InteractionType, fill = InteractionTypeFiner), alpha = .08, size = .8) +
        geom_hline(size = .5, yintercept = c(0,1), linetype = 1, color = "grey90") +
        geom_line(size = 1, aes(x = Time, y = Isolate1CFUFreqMean, color = Isolate1InitialODFreq, group = Isolate1InitialODFreq)) +
        geom_point(size = 1.5, aes(x = Time, y = Isolate1CFUFreqMean, color = Isolate1InitialODFreq, group = Isolate1InitialODFreq)) +
        geom_segment(size = 1, aes(x = Time, xend = Time,
                                    y = Isolate1CFUFreqMean + 2*Isolate1CFUFreqSd,
                                    yend = Isolate1CFUFreqMean - 2*Isolate1CFUFreqSd,
                                    color = Isolate1InitialODFreq)) +
        scale_x_continuous(expand = c(0, 0), limits = c(axis_lower, axis_upper)) +
        scale_y_continuous(expand = c(0, 0), breaks = c(0, .5, 1), limits = c(axis_lower+0.2, axis_upper-0.2)) +
        scale_color_manual(values = c(frequency_color, interaction_color), label = c("95%", "50%", "5%", "exclusion", "coexistence")) +
        scale_fill_manual(values = assign_interaction_color(level = "finer")) +
        theme_classic() +
        theme(panel.spacing = unit(2, "mm"),
              panel.border = element_blank(),
              panel.grid.minor.y = element_blank(),
              axis.line = element_blank(),
              axis.title = element_blank(), axis.text = element_blank(),
              axis.ticks = element_blank(),
              plot.background = element_blank(),
              plot.title = element_text(size = 5, margin = margin(0,0,0,0)),
              plot.margin = margin(0,0,0,0, "mm")) +
        guides(color = "none", fill = "none") +
        labs(x = "Time", y = "Frequency") +
        #ggtitle(unique(pairs_freq$PairID)) +
        NULL
}

# Legend for fill
pairs_interaction_finer <- pairs %>%
    group_by(InteractionType, InteractionTypeFiner) %>%
    count(name = "Count") %>%
    ungroup() %>%
    mutate(Fraction = Count / sum(Count)) %>%
    mutate(InteractionTypeFiner = factor(InteractionTypeFiner, interaction_type_finer)) %>%
    arrange(InteractionTypeFiner)
p_legend_fill <- pairs %>%
    get_interaction_legend()

# Append competition outcome to frequencies
pairs_example_freq <- pairs %>%
    filter(!is.na(FitnessFunction)) %>%
    mutate(InteractionType = factor(InteractionType, c("exclusion", "coexistence"))) %>%
    mutate(InteractionTypeFiner = factor(InteractionTypeFiner, interaction_type_finer)) %>%
    arrange(InteractionType, InteractionTypeFiner) %>%
    left_join(pairs_freq) %>%
    select(PairID, InteractionType, InteractionTypeFiner, Isolate1InitialODFreq, Time, Isolate1CFUFreqMean, Isolate1CFUFreqSd) %>%
    mutate(Isolate1InitialODFreq = factor(Isolate1InitialODFreq), Time = factor(Time, c("T0", "T8")))

temp_list <- pairs_example_freq %>%
    arrange(InteractionTypeFiner, PairID) %>%
    group_split(InteractionTypeFiner, PairID) %>%
    lapply(plot_example_freq)

p_example <- pairs_example_freq %>%
    filter(PairID == 3) %>%
    plot_example_freq2() +
    theme(axis.text = element_text(size = 10), axis.title = element_text(size = 10)) +
    scale_x_continuous(breaks = c(0,1), labels = c("start", "end"), limits = c(-0.5, 1.5))

## Grid layout
m <- matrix(c(1:nrow(pairs), rep(NA, 10-(nrow(pairs)%%10))), nrow = 10)
p_waffle <- arrangeGrob(grobs = temp_list, layout_matrix = m)
p_waffle <- plot_grid(NULL, p_waffle, rel_widths = c(1, 2.5), scale = c(1, 0.9)) + paint_white_background()

ss = .3
p <- ggdraw(p_waffle) +
    draw_plot(p_example, x = .1, y = .75, width = ss*0.4, height = ss*0.8, hjust = 0.5, vjust = .5) +
    draw_plot(p_legend_color, x = .2, y = .8, width = ss/2, height = ss/2, hjust = 0.5, vjust = .5) +
    draw_plot(p_legend_fill, x = .16, y = .4, width = ss, height = ss*2, hjust = 0.5, vjust = .5) +
    theme(panel.background = element_blank(), plot.background = element_rect(color = NA, fill = "white"),
          plot.margin = unit(c(0,0,0,0), "mm"))

ggsave(here::here("plots/Fig3.png"), p, width = 12, height = 5)



# Figure 4 ----
# Figure 4A One example network
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
    draw_image(here::here("plots/cartoons/Fig2B_1.png"), x = -.2, y = -1, vjust = 0.5, hjust = 0, clip = "on", scale = .8) +
    draw_image(here::here("plots/cartoons/Fig2B_3.png"), x = .2, y = -1, vjust = 0.5, hjust = 0, clip = "on", scale = .8) +
    draw_image(here::here("plots/cartoons/Fig2B_4.png"), x = -.2, y = -3, vjust = 0.5, hjust = 0, clip = "on", scale = .8) +
    draw_image(here::here("plots/cartoons/Fig2B_2.png"), x = .2, y = -3, vjust = .6, hjust = 0, clip = "on", scale = .8) +
    paint_white_background()


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

## legend
edge_width = 0.5
temp <- p_net_hierarchy_list[[13]] +
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
p_legend <- get_legend(temp)
pB <- ggdraw(pB) + draw_plot(p_legend,.2,.2,.1,.1)

# Figure 4C: Hierarchy
pC <- mutate(communities_hierarchy, Treatment = "experiment") %>%
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

#
p_left <- plot_grid(pA, pC, ncol = 1, rel_heights = c(1,1), scale = c(.8, .9), labels = c("A", "C"), axis = "lr", align = "v")
p <- plot_grid(p_left, pB, nrow = 1, rel_widths = c(1,4), scale = c(1, .9), labels = c("", "B")) + paint_white_background()
ggsave(here::here("plots/Fig4.png"), p, width = 10, height = 5)



