library(tidyverse)
library(cowplot)
library(tidygraph)
library(ggraph)
library(igraph)
source(here::here("analysis/00-metadata.R"))

isolates <- read_csv(paste0(folder_data, "output/isolates_remained.csv"), show_col_types = F)
pairs <- read_csv(paste0(folder_data, "output/pairs_remained.csv"), show_col_types = F)
communities <- read_csv(paste0(folder_data, "temp/00c-communities.csv"), show_col_types = F)
pairs_freq <- read_csv(paste0(folder_data, "temp_old/93a-pairs_freq.csv"), show_col_types = F)
load(paste0(folder_data, "temp_old/95-communities_network.Rdata"))
#communities_hierarchy <- read_csv(paste0(folder_data, "temp/95-communities_hierarchy.csv"), show_col_types = F)


# Check if for all exclusion pairs, it's the higher ranked species excluding the lower rank
pairs_exclusion <- pairs %>%
    filter(outcome %in% c("1-exclusion", "2-exclusion")) %>%  # 103 exclusion pairs
    select(Community, From, To, outcome)

isolates_rank <- isolates %>%
    select(Community, Isolate, Rank)


pairs_exclusion %>%
    left_join(rename(isolates_rank, From = Isolate, FromRank = Rank)) %>%
    left_join(rename(isolates_rank, To = Isolate, ToRank = Rank)) %>%
    filter(FromRank > ToRank)


isolates %>%
    select(Community, Isolate, Rank, Win, Draw, Lose, Game, Score) %>%
    filter(Community == "C11R1")

#
pairs_freq %>%
    filter(Community == "C11R1", Isolate1 == 3, Isolate2 == 9) %>%
    ggplot() +
    geom_line(aes(x = Time, y = Isolate1CFUFreqMean, color = factor(Isolate1InitialODFreq), group = Isolate1InitialODFreq)) +
    theme_classic() +
    theme() +
    guides() +
    labs()


# Plot individual networks
set.seed(1)
node_size = 3
edge_width = .8
plot_network_hierarchy <- function(net, n_rank = 12, n_break = 12) {
    #net <- communities_network$Network[[2]]
    graph_ranked <- net %>%
        activate(nodes) %>%
        select(Isolate, Rank, PlotRank) %>%
        activate(edges) %>%
        mutate(fromRank = .N()$PlotRank[match(from, .N()$Isolate)],
               toRank = .N()$PlotRank[match(to, .N()$Isolate)])

    n_nodes <- igraph::vcount(graph_ranked)
    n_exclusion_violation <- graph_ranked %>%
        activate(edges) %>%
        filter(fromRank > toRank, outcome == "1-exclusion" | outcome == "2-exclusion") %>%
        igraph::ecount()

    if (n_exclusion_violation != 0) {
        pairs_exclusion_violation <- graph_ranked %>%
            activate(edges) %>%
            filter(fromRank > toRank, outcome == "1-exclusion" | outcome == "2-exclusion") %>%
            mutate(temp = fromRank, fromRank = toRank, toRank = temp) %>%
            select(-temp) %>%
            mutate(outcome = "6-exclusion violating rank") %>%
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
        #filter(fromRank > toRank) %>%
        filter((outcome == "1-exclusion" | outcome == "2-exclusion")) %>%
        mutate(Temp = sample(c(-1, 1), size = n(), replace = T))
    strength_angle <- as_tibble(temp)$Temp * 0.06

    temp %>%
        ggraph(layout = "nicely") +
        geom_hline(yintercept = c(-n_rank:-1), color = "grey90") +
        geom_node_point(size = node_size, shape = 21, fill = "grey", stroke = node_size/5, color = "black") +
        geom_edge_arc(strength = strength_angle, alpha = 0.5,
                      aes(color = outcome), width = edge_width,
                      arrow = arrow(length = unit(edge_width, "mm"), type = "closed", angle = 30, ends = "last"),
                      start_cap = circle(node_size/2, "mm"),
                      end_cap = circle(node_size/2, "mm")) +
        scale_edge_color_manual(values = outcome_colors, labels = outcome_labels) +
        # scale_edge_color_manual(values = assign_interaction_color(level = "hierarchy"),
        #                         breaks = c("exclusion", "exclusion violating rank", "coexistence", "unknown"),
        #                         labels = c("exclusion pair that follows rank", "exclusion pair that violates rank", "coexistence", "unknown")) +
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

# Check the isolate identity in hierarchy ----
p_net_hierarchy_list <- communities_network_hierarchy$NetworkHierarchyPlot
#for (i in 1:13) p_net_hierarchy_list[[i]] <- p_net_hierarchy_list + theme(legend.position = "none")
p_net_hierarchy_list[[13]] <- p_net_hierarchy_list[[13]] +
    scale_y_continuous(limits = c(-12-1, 0), breaks = -12:-1, labels = 12:1, position = "right") +
    theme(axis.title.y = element_text(color = 1, size = 10, angle = 270, margin = margin(l = 2, unit = "mm")),
          axis.text.y = element_text(color = 1, size = 10, margin = margin(l = 1, unit = "mm")))
p_axistitle <- ggdraw() + draw_label("Community", fontface = 'plain', x = .5, hjust = .5) + theme(plot.margin = margin(5, 0, 5, 7))
p_temp <- plot_grid(plotlist = p_net_hierarchy_list,
                    rel_widths = c(communities_network_hierarchy$CommunitySize / max(communities_network_hierarchy$CommunitySize))^1.5,
                    labels = 1:13, label_fontface = "plain", label_x = c(rep(0.5, 12), 0.45), hjust = c(rep(.5, 12), 1),
                    nrow = 1, axis = "tb", align = "h") + paint_white_background()
p <- plot_grid(p_axistitle, p_temp, ncol = 1, rel_heights = c(.1, 1)) + paint_white_background()


# line legend
edge_width = 0.5
p_legend <- get_legend({
    p_net_hierarchy_list[[13]] +
        geom_edge_arc(strength = 10,
                      aes(color = outcome), width = edge_width*1.3,
                      arrow = arrow(length = unit(edge_width*5, "mm"), type = "closed", angle = 30, ends = "last"),
                      start_cap = circle(node_size/2, "mm"),
                      end_cap = circle(node_size/2, "mm")) +
        scale_alpha_manual(values = 1) +
        theme(
            legend.key.size = unit(2, "line"),
            legend.key.height = unit(6, "mm"),
            legend.position = "right",
            legend.direction = "vertical",
            legend.text = element_text(size = 10),
            legend.background = element_rect(fill = NA, color = NA),
        ) +
        guides(color = guide_legend(title = "pairwise competition outcome", override.aes = list(alpha = 1, linewidth = 2)))
})

p <- ggdraw(p) +
    draw_plot(p_legend, x = 0.5, y = 0.15, width = 0.1, height = 0.1,  hjust = 0, vjust = 0)

#ggsave(here::here("plots/Fig4.png"), p, width = 10, height = 5)

ggsave(paste0(folder_data, "temp/36-1-graph.png"), p, width = 10, height = 5)

# Check the isolate identity ----
plot_network_hierarchy_family <- function(net, n_rank = 12, n_break = 12) {
    #net <- communities_network$Network[[2]]
    graph_ranked <- net %>%
        activate(nodes) %>%
        #select(Isolate, Rank, PlotRank) %>%
        activate(edges) %>%
        mutate(fromRank = .N()$PlotRank[match(from, .N()$Isolate)],
               toRank = .N()$PlotRank[match(to, .N()$Isolate)])

    n_nodes <- igraph::vcount(graph_ranked)
    n_exclusion_violation <- graph_ranked %>%
        activate(edges) %>%
        filter(fromRank > toRank, outcome == "1-exclusion" | outcome == "2-exclusion") %>%
        igraph::ecount()

    if (n_exclusion_violation != 0) {
        pairs_exclusion_violation <- graph_ranked %>%
            activate(edges) %>%
            filter(fromRank > toRank, outcome == "1-exclusion" | outcome == "2-exclusion") %>%
            mutate(temp = fromRank, fromRank = toRank, toRank = temp) %>%
            select(-temp) %>%
            mutate(outcome = "6-exclusion violating rank") %>%
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
        #filter(fromRank > toRank) %>%
        filter((outcome == "1-exclusion" | outcome == "2-exclusion")) %>%
        mutate(Temp = sample(c(-1, 1), size = n(), replace = T))
    strength_angle <- as_tibble(temp)$Temp * 0.06

    temp %>%
        ggraph(layout = "nicely") +
        geom_hline(yintercept = c(-n_rank:-1), color = "grey90") +
        geom_node_point(aes(fill = Family), size = node_size, shape = 21, stroke = node_size/5, color = "black") +
        geom_edge_arc(strength = strength_angle, alpha = 0.5,
                      aes(color = outcome), width = edge_width,
                      arrow = arrow(length = unit(edge_width, "mm"), type = "closed", angle = 30, ends = "last"),
                      start_cap = circle(node_size/2, "mm"),
                      end_cap = circle(node_size/2, "mm")) +
        scale_edge_color_manual(values = outcome_colors, labels = outcome_labels) +
        # scale_edge_color_manual(values = assign_interaction_color(level = "hierarchy"),
        #                         breaks = c("exclusion", "exclusion violating rank", "coexistence", "unknown"),
        #                         labels = c("exclusion pair that follows rank", "exclusion pair that violates rank", "coexistence", "unknown")) +
        scale_x_continuous(limits = c(0.1, 0.9), expand = c(0,0)) +
        scale_y_continuous(limits = c(-n_break-1, 0), breaks = -n_break:-1, labels = n_break:1) +
        theme_void() +
        theme(
#            legend.position = "none",
            legend.title = element_blank(),
            axis.title = element_blank(),
            strip.text = element_blank(),
            plot.margin=unit(c(0,0,0,0),"mm")
        ) +
        labs(y = "Rank")

}

p <- plot_network_hierarchy_family(communities_network$Network[[12]]) + paint_white_background()

ggsave(paste0(folder_data, "temp/36-2-graph_family.png"), p, width = 5, height = 5)


# Bionomial distribution of tranistivity ----
#triad_possible <- tibble(CommunitySize = 3:12, TotalTriads = choose(CommunitySize, 3))

subset_exclusion_network <- function (net) {
    net %>%
        activate(edges) %>%
        filter(outcome %in% c("1-exclusion", "2-exclusion"))
}

communities_network <- communities_network %>%
    rowwise() %>%
    mutate(ExclusionNetwork = list(subset_exclusion_network(Network)))

networks <- union(
    communities_network$Network[[1]], communities_network$Network[[2]],
    communities_network$Network[[3]], communities_network$Network[[4]],
    communities_network$Network[[5]], communities_network$Network[[6]],
    communities_network$Network[[7]], communities_network$Network[[8]],
    communities_network$Network[[9]], communities_network$Network[[10]],
    communities_network$Network[[11]], communities_network$Network[[12]],
    communities_network$Network[[13]], byname = "auto")

triad_census(networks) %>% sum() # 220 possible triads with 0, 1, 2, or 3 links
triad_census(networks) %>% `[`(c(9, 10, 12:16)) %>% sum() # 182 fully connected nodes

exclusion_networks <- union(
    communities_network$ExclusionNetwork[[1]], communities_network$ExclusionNetwork[[2]],
    communities_network$ExclusionNetwork[[3]], communities_network$ExclusionNetwork[[4]],
    communities_network$ExclusionNetwork[[5]], communities_network$ExclusionNetwork[[6]],
    communities_network$ExclusionNetwork[[7]], communities_network$ExclusionNetwork[[8]],
    communities_network$ExclusionNetwork[[9]], communities_network$ExclusionNetwork[[10]],
    communities_network$ExclusionNetwork[[11]], communities_network$ExclusionNetwork[[12]],
    communities_network$ExclusionNetwork[[13]], byname = "auto")

triad_census(exclusion_networks) %>% `[`(c(9, 10, 12:16)) %>% sum() # 77 fully connected nodes


union(communities_network$ExclusionNetwork[[1]], communities_network$ExclusionNetwork[[2]],
      communities_network$ExclusionNetwork[[3]], communities_network$ExclusionNetwork[[4]]) %>%
    triad_census() %>% `[`(c(9, 10, 12:16)) %>% sum() # These small networks should only have total 2 fully connected triads

set.seed(2)
N = 10000
tb <- tibble(
    BootstrapID = 1:N,
    NumberExpectedRPS = rbinom(N, 77, prob = 1/4)
)

range(tb$NumberExpectedRPS)

p <- tb %>%
    ggplot() +
    geom_histogram(aes(x = NumberExpectedRPS), color = "black", fill = "white") +
    scale_x_continuous(limits = c(0, max(tb$NumberExpectedRPS)+2)) +
    annotate("text", x = Inf, y = Inf, label = "X ~ Binomial(n=77, p = 1/4)", vjust = 5, hjust = 2) +
    theme_classic() +
    theme() +
    guides() +
    labs(x = "number of expected nontransitivity", y = "count")

ggsave(paste0(folder_data, "temp/36-3-bionomial.png"), p, width = 5, height = 5)




# Test to form a directed network
tbl_graph(nodes = tibble(node = 1:3), edges = tibble(from = c(1,1,2), to = c(2,3,3))) %>%
    triad_census() %>% `[`(c(9, 10, 12:16)) %>% sum()













count_exclusion_RPS <- function (net) {
    temp <- net %>%
        activate(edges) %>%
        filter(outcome == "exclusion") %>%
        triad_census()
    temp[10]
}
count_exclusion_triads <- function (net) {
    temp <- net %>%
        activate(edges) %>%
        filter(outcome == "exclusion") %>%
        triad_census()
    sum(temp[9:10])
}

communities_network$Network[[6]] %>%
    count_exclusion_triads()


communities_network_rps <- communities_network %>%
    rowwise() %>%
    #mutate(ExclusionNetwork = extract_exclusion_network(Network) %>% list()) %>%
    mutate(CountRPS = count_exclusion_RPS(Network)) %>% # count rock-paper-scissor triads
    mutate(Xaxis = "haha")

expected_nT <- tibble(Xaxis = c("haha", "haha"), Fraction = c(3/4, 1/4), Motif = c("Transitive", "Nontransitive"))


pC <- communities_network_rps %>%
    ggplot() +
    geom_col(data = mutate(expected_nT, Motif = factor(Motif, c("Transitive", "Nontransitive"))), aes(x = Xaxis, Fraction, group = Motif), fill = "white", color = 1, width = .6) +
    geom_point(aes(x = Xaxis, y = CountRPS), position = position_jitter(width = .2, height = 0), shape = 21, size = 2, stroke = 0.8) +
    #scale_fill_manual(value = ) +
    #scale_x_continuous(limits = c(0.5, 1.5)) +
    scale_y_continuous(expand = c(0,0.05), breaks = seq(0,1, 0.25)) +
    theme_classic() +
    theme(
        panel.grid.major.y = element_line(color = "grey", linetype = 2),
        axis.text = element_text(size = 10, color = 1),
        axis.text.x = element_blank(),
        panel.border = element_rect(fill = NA, color = 1, linewidth = 1)
        #panel.border = element_rect(color = 1, fill = NA, linewidth = 0.5)
    ) +
    guides() +
    labs(x = "", y = "fraction")



# Test Hierarchy ----
p <- mutate(communities_hierarchy, Treatment = "experiment") %>%
    ggplot(aes(x = Treatment, y = HierarchyScore)) +
    geom_boxplot(width = .5, lwd = .8, outlier.color = NA) +
    geom_jitter(shape = 1, size = 2, stroke = .8, height = 0, width = .1) +
    #scale_y_continuous(limits = c(0,1.01), breaks = c(0, .25, .5, .75, 1)) +
    theme_classic() +
    theme(
        panel.grid.major.y = element_line(color = "grey", linetype = 2),
        panel.spacing = unit(0, "mm"),
        panel.border = element_rect(fill = NA, color = 1, linewidth = 1),
        axis.text = element_text(size = 10, color = 1),
        axis.text.x = element_blank(),
        axis.title = element_text(size = 10, color = 1)
    ) +
    labs(x = "", y = "hierarchy score") +
    guides(color = "none")


ggsave(paste0(folder_data, "temp/36-2-hierarchy.png"), p, width = 4, height = 4)

#
ss <- 0.13
p_left <- plot_grid(pA, NULL, ncol = 1, rel_heights = c(1,1.5), scale = c(.8, .9), labels = c("A", ""), axis = "lr", align = "v")
p_bottom <- plot_grid(pC, pD, nrow = 1, align = "h", axis = "tb", labels = c("C", "D"))
p <- plot_grid(p_left, pB, nrow = 1, rel_widths = c(1,4), scale = c(1, .9), labels = c("", "B")) + paint_white_background()
p <- ggdraw(p) +
    draw_plot(p_bottom, x = 0, y = 0, width = 0.4, height = 0.5,  hjust = 0, vjust = 0) +
    #draw_plot(pD, x = 0.3, y = 0, width = 0.2, height = 0.5,  hjust = 0, vjust = 0) +
    draw_plot(p_T, x = 0.12, y = .25, width = ss*0.6, height = ss, hjust = 0, vjust = 0) +
    draw_plot(p_nT, x = 0.12, y = .07, width = ss*0.6, height = ss, hjust = 0, vjust = 0) +
    draw_plot(p_legend, x = 0.5, y = 0.15, width = 0.1, height = 0.1,  hjust = 0, vjust = 0)
ggsave(here::here("plots/Fig4.png"), p, width = 10, height = 5)


















