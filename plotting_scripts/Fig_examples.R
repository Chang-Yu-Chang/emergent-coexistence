# Make example figures

library(tidyverse)
library(data.table)
library(tidygraph)
library(ggraph)
library(cowplot)
source("network_functions.R")

# Figure: barplots of abundance changes ----
t_coexist <- tibble(
    ID = rep(rep(c(1,2), each = 5), 3),
    Abundance = c(
        c(0.95, 0.75, 0.65, 0.40, 0.30, 0.05, 0.25, 0.35, 0.60, 0.70),
        c(0.50, 0.45, 0.42, 0.35, 0.30, 0.50, 0.55, 0.58, 0.65, 0.70),
        c(0.05, 0.15, 0.20, 0.25, 0.30, 0.95, 0.85, 0.80, 0.75, 0.70)),
    Transfer = rep(rep(c(1:5), 2), 3),
    InitialFrequency = rep(c("95-5", "50-50", "5-95"), each = 10),
    InteractionType = rep("coexistence", 30)
)


t_exclusion <- tibble(
    ID = rep(rep(c(1,2), each = 5), 3),
    Abundance = c(
        c(0.95, 0.50, 0.30, 0.15, 0.00, 0.05, 0.50, 0.70, 0.85, 1.00),
        c(0.50, 0.30, 0.15, 0.02, 0.00, 0.50, 0.70, 0.85, 0.98, 1.00),
        c(0.05, 0.03, 0.00, 0.00, 0.00, 0.95, 0.97, 1.00, 1.00, 1.00)),
    Transfer = rep(rep(c(1:5), 2), 3),
    InitialFrequency = rep(c("95-5", "50-50", "5-95"), each = 10),
    InteractionType = rep("exclusion", 30)
)

p_bar <- bind_rows(t_coexist, t_exclusion) %>%
    mutate(ID = factor(ID)) %>%
    ggplot() +
    geom_bar(aes(x = Transfer, y = Abundance, fill = ID), stat = "identity", color = 1) +
    scale_y_continuous(breaks = c(0, 0.5, 1)) +
    scale_fill_manual(values = c("white", "grey40")) +
    facet_grid(InitialFrequency~InteractionType) +
    theme_cowplot() +
    theme(legend.position = "none",
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank()) +
    labs(x = "Time")

p_line <- bind_rows(t_coexist, t_exclusion) %>%
    filter(ID == 2, Transfer %in% c(1,4)) %>%
    mutate(Transfer = ifelse(Transfer == 1, "T1", "T8")) %>%
    mutate(ID = factor(ID), Transfer = factor(Transfer)) %>%
    ggplot(aes(x = Transfer, y = Abundance, group = InitialFrequency)) +
    geom_line() +
    geom_point() +
    scale_y_continuous(limits = c(0,1), breaks = c(0, 0.5, 1)) +
    scale_fill_manual(values = c("white", "grey40")) +
    facet_grid(.~InteractionType) +
    theme_bw() +
    theme(legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_blank()) +
    labs(x = "Time", y = "Frequency")
p_line

isolates1 <- tibble(Isolate = 2:1, PlotRank = 2:1, Rank = 2:1)
edges1 <- tibble(From = 1, To = 2, InteractionType = "coexistence")
g1 <- make_network(isolates1, edges1) %>%
    activate(nodes) %>% mutate(graph = 1)

isolates2 <- tibble(Isolate = 2:1, PlotRank = 2:1, Rank = 2:1)
edges2 <- tibble(From = 1, To = 2, InteractionType = "exclusion")
g2 <- make_network(isolates2, edges2) %>%
    activate(nodes) %>% mutate(graph = 2)

p_coexistence <- plot_competitive_network(g1, g_layout = "linear", node_size = 10) + theme(legend.position = "none", panel.background = element_blank())
p_exclusion <- plot_competitive_network(g2, g_layout = "linear", node_size = 10) + theme(legend.position = "none", panel.background = element_blank())
p_links <- plot_grid(p_coexistence, p_exclusion, nrow = 1)


#p1 <- plot_grid(p_links, p_bar, ncol = 1, rel_heights = c(2, 8), scale = c(1, 1), align = "hv", axis = "lr")
p1 <- plot_grid(p_links, p_line, ncol = 1, rel_heights = c(2, 2), scale = c(1, 1), align = "hv", axis = "lr")
ggsave("../plots/Ex1_mutual_invasion.png", plot = p1, width = 5, height = 5)
ggsave("../plots/Ex1_mutual_invasion.pdf", plot = p1, width = 5, height = 5)



# Figure: diagram of individual motifs ----
make_example_motifs <- function() {
    temp_list <- rep(list(NA), 7)
    temp_id <- c(11, 7, 8, 12, 13, 14, 15)

    for (i in 1:7) {
        g <- as_tbl_graph(igraph::graph.isocreate(size = 3, temp_id[i]))
        layout <- create_layout(g, layout = 'circle')
        g <- activate(g, edges) %>% mutate(InteractionType = ifelse(edge_is_mutual(), "coexistence", "exclusion"))
        temp_list[[i]] <- g %>% activate(nodes) %>% mutate(x = layout$x, y = layout$y, graph = paste0("motif", i))
    }

    return(temp_list)
}

example_motif_list <- make_example_motifs()
p2 <- bind_graphs(example_motif_list)  %>%
    plot_competitive_network(g_layout = "example_motif", node_size = 5) +
    facet_nodes(~graph, nrow = 1)

ggsave("../plots/Ex2_motifs.png", plot = p2, width = 10, height = 1.5)
ggsave("../plots/Ex2_motifs_enlarged.png", plot = p2, width = 8, height = 1.2)
save(example_motif_list, file = "../data/output/example_motif_list.Rdata")


# Figure: example competitive networks ----
if (FALSE) {
    plot_competitive_network2 <- function(graph, node_size = 10, layout = "circle") {
        # Layout
        graph_layout <- create_layout(graph, layout)
        graph <- graph %>% activate(nodes) %>% mutate(x = graph_layout$x, y = graph_layout$y)

        # Nodes
        nodex_axis_x <- activate(graph, nodes) %>% pull(x) %>% range()
        nodex_axis_y <- activate(graph, nodes) %>% pull(y) %>% range()

        # Edges
        interaction_type <- c("exclusion", "coexistence", "lose", "bistability", "neutrality", "self", "undefined")
        interaction_color = c("#DB7469", "#557BAA", "#73C966", "#EECF6D", "#8650C4", "black", "grey80")
        names(interaction_color) <- interaction_type

        graph %>%
            ggraph(layout = "nicely") +
            geom_node_point(size = node_size, shape = 21, fill = "gray", colour = "black", stroke = node_size/5) +
            geom_edge_link(aes(color = InteractionType), width = node_size/10,
                arrow = arrow(length = unit(node_size/2, "mm"), type = "closed", angle = 30, ends = "last"),
                start_cap = circle(node_size/2+1, "mm"),
                end_cap = circle(node_size/2+1, "mm")) +
            scale_edge_color_manual(values = interaction_color) +
            scale_x_continuous(limits = nodex_axis_x*1.2) +
            scale_y_continuous(limits = nodex_axis_y*1.2) +
            theme_graph() +
            theme(
                legend.position = "top",
                legend.title = element_blank(),
                legend.text = element_text(size=20, aes(color = InteractionType)),
                plot.margin=unit(c(3,3,3,3),"mm")
            )
    }
}
set.seed(1)
n_nodes = 6
example_pairs <- as_tibble(t(combn(n_nodes,2))) %>%
    setNames(c("From", "To")) %>%
    mutate(InteractionType = sample(c("coexistence", "exclusion"), size = n(), replace = T, prob = c(.5, .5)))
example_isolates <- tibble(Isolate = 1:n_nodes, Rank = 1:n_nodes, PlotRank = 1:n_nodes)
example_graph <- make_network(isolates = example_isolates, pairs = example_pairs)
p3_network <- plot_competitive_network(example_graph, g_layout = "circle")
#ggsave("../plots/Ex3_graph.png", plot = p3, width = 5, height = 5)

## Motif distribution
p3_motif <- tibble(Motif = 1:7, Count = count_motif(example_graph)) %>%
    ggplot() +
    geom_bar(aes(x = Motif, y = Count), stat = "identity", fill = "grey80", color = 1) +
    scale_x_continuous(limits = c(0,8), breaks = 1:7, expand = c(0,0)) +
    scale_y_continuous(breaks = c(0, 5, 10), limits = c(0,10)) +
    theme_cowplot() +
    theme(legend.position = c(0.05, .5)) +
    ggtitle("")
p3 <- plot_grid(p3_network, p3_motif, nrow = 1, rel_widths = c(4, 8))
ggsave("../plots/Ex3_graph.png", plot = p3, width = 12, height = 4)

# Figure: example global competitive network ----
set.seed(1)
n_nodes = 20

## Full coexistence
example_pairs <- as_tibble(t(combn(n_nodes,2))) %>%
    setNames(c("From", "To")) %>%
    mutate(InteractionType = sample(c("coexistence", "exclusion"), size = n(), replace = T, prob = c(1,0)))
example_isolates <- tibble(Isolate = 1:n_nodes, Rank = 1:n_nodes, PlotRank = 1:n_nodes)
example_graph <- make_network(isolates = example_isolates, pairs = example_pairs)
p_coexistence <- plot_competitive_network(example_graph, g_layout = "linear", node_size = 3)


example_pairs <- as_tibble(t(combn(n_nodes,2))) %>%
    setNames(c("From", "To")) %>%
    mutate(InteractionType = sample(c("coexistence", "exclusion"), size = n(), replace = T, prob = c(0.1, 0.9)))
isolates_win <- example_pairs %>%
    group_by(From, InteractionType) %>%
    summarise(Count = n()) %>%
    ungroup() %>%
    pivot_wider(names_from = InteractionType, values_from = Count) %>%
    replace_na(list(coexistence = 0, exclusion = 0)) %>%
    mutate(Isolate = From, Win = exclusion, Draw = coexistence) %>%
    select(Isolate, Win, Draw)
isolates_lose <- example_pairs %>%
    group_by(To, InteractionType) %>%
    summarise(Count = n()) %>%
    ungroup() %>%
    pivot_wider(names_from = InteractionType, values_from = Count) %>%
    replace_na(list(coexistence = 0, exclusion = 0)) %>%
    mutate(Isolate = To, Lose = exclusion, Draw = coexistence) %>%
    select(Isolate, Lose)

example_isolates <- full_join(isolates_win, isolates_lose, by = "Isolate") %>%
    replace_na(list(Win = 0, Draw = 0, Lose = 0)) %>%
    mutate(Score = Win - Lose) %>%
    arrange(desc(Score)) %>%
    mutate(Rank = 1:n_nodes, PlotRank = 1:n_nodes) %>%
    arrange(Isolate)

example_graph <- make_network(isolates = example_isolates, pairs = example_pairs)
p_hierarchy <- plot_competitive_network(example_graph, g_layout = "linear", node_size = 3)

p4 <- plot_grid(p_coexistence, p_hierarchy, nrow = 1)
ggsave("../plots/Ex4_graph_global.png", plot = p4, width = 8, height = 4)



# Figure: example motif count ----
simulated_motif_counts <- fread("../data/temp/simulated_motif_counts.txt")
motif_color <- assign_motif_color()

p5 <- simulated_motif_counts %>%
    mutate(MotifType = ifelse(Motif == 7, "Motif7", "Others") %>% ordered(levels = c("Others", "Motif7"))) %>%
    filter(CommunitySize == 12) %>%
    #mutate(Motif = factor(Motif)) %>%
    group_by(CommunitySize, ProbPairCoexistence, MotifType) %>%
    summarize(MeanCount = mean(Count)) %>%
    group_by(CommunitySize, ProbPairCoexistence) %>%
    mutate(SumMeanCount = sum(MeanCount), RelativeMeanCount = MeanCount/SumMeanCount) %>%
    ggplot() +
    geom_area(aes(x = ProbPairCoexistence, y = RelativeMeanCount, fill = MotifType), color = 1) +
    scale_x_continuous(expand = c(0,0), breaks = c(0, 0.5, 1)) +
    scale_y_continuous(expand = c(0,0), breaks = c(0, 0.5, 1)) +
    scale_fill_manual(values = c("Others" = grey(1), `Motif7` = grey(0.5))) +
    theme_cowplot() +
    panel_border(color = 1) +
    theme(legend.title = element_blank(), legend.position = "top", legend.direction = "horizontal") +
    labs(x = "Pairwise coexistence in pool", y = "Expected motif count")

#p_B <- plot_grid(p_motifs, p1, ncol = 1, rel_heights = c(1,5))
ggsave("../plots/Ex_motif_random.png", plot = p5, width = 4, height = 4)
