# Make the example plots

library(tidyverse)
library(tidygraph)
library(ggraph)
library(cowplot)
source("network_functions.R")
interaction_type <- c("exclusion", "coexistence", "lose", "bistability", "neutrality", "self", "undefined")
interaction_color = c("#DB7469", "#557BAA", "#73C966", "#EECF6D", "#8650C4", "black", "grey80")
names(interaction_color) <- interaction_type


# Diagram of 7 motifs
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
p1 <- bind_graphs(example_motif_list)  %>%
    plot_competitive_network(g_layout = "example_motif", node_size = 5) +
    facet_nodes(~graph, nrow = 1)

ggsave("../plots/talk-exmple_motifs.png", plot = p1, width = 10, height = 1.5)
#save(example_motif_list, file = "../data/temp/example_motif_list.Rdata")

# Make the example network
isolates <- read_csv(here::here("data/output/isolates.csv"))
pairs <- read_csv(here::here("data/output/pairs.csv"))

n_species = 6
set.seed(1)

plot_competitive_network2 <- function(graph, interaction_color, node_size = 10, layout = "circle") {
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
        theme_void() +
        theme(legend.position = "none", plot.margin = margin(0,0,0,0)) +
        labs(x = "") %>%
        return()
}

generate_example_graphs <- function (n_species) {
    layout = "circle"
    example_pairs <- as_tibble(t(combn(n_species,2))) %>%
        setNames(c("From", "To")) %>%
        mutate(InteractionType = sample(c("coexistence", "exclusion"), size = n(), replace = T))
    example_isolates <- tibble(Isolate = 1:n_species, Rank = 1:n_species, PlotRank = 1:n_species)
    example_graph <- make_network(isolates = example_isolates, pairs = example_pairs)

    graph_layout <- create_layout(example_graph, layout)
    example_graph <- example_graph %>% activate(nodes) %>% mutate(x = graph_layout$x, y = graph_layout$y)
    nodex_axis_x <- activate(example_graph, nodes) %>% pull(x) %>% range()
    nodex_axis_y <- activate(example_graph, nodes) %>% pull(y) %>% range()

    ## grey links
    node_size = 10
    p2 <- example_graph %>%
        ggraph(layout = "nicely") +
        geom_node_point(size = node_size, shape = 21, fill = "gray", colour = "black", stroke = node_size/5) +
        geom_edge_link(aes(color = InteractionType), width = node_size/10,
                       start_cap = circle(node_size/2+1, "mm"),
                       end_cap = circle(node_size/2+1, "mm")) +
        scale_edge_color_manual(values = c("grey", "grey")) +
        scale_x_continuous(limits = nodex_axis_x*1.2) +
        scale_y_continuous(limits = nodex_axis_y*1.2) +
        theme_void() +
        theme(legend.position = "none", plot.margin = margin(0,0,0,0)) +
        labs(x = "")

    ## colored links
    p3 <- example_graph %>%
        ggraph(layout = "nicely") +
        geom_node_point(size = node_size, shape = 21, fill = "gray", colour = "black", stroke = node_size/5) +
        geom_edge_link(aes(color = InteractionType), width = node_size/10,
                       arrow = arrow(length = unit(node_size/2, "mm"), type = "closed", angle = 30, ends = "last"),
                       start_cap = circle(node_size/2+1, "mm"),
                       end_cap = circle(node_size/2+1, "mm")) +
        scale_edge_color_manual(values = interaction_color) +
        scale_x_continuous(limits = nodex_axis_x*1.2) +
        scale_y_continuous(limits = nodex_axis_y*1.2) +
        theme_void() +
        theme(legend.position = "none", plot.margin = margin(0,0,0,0)) +
        labs(x = "")

    ## Pairs barplot example
    p_pairs_interaction <- example_pairs %>%
        group_by(InteractionType) %>%
        summarize(Count = n()) %>%
        ggplot() +
        geom_col(aes(x = InteractionType, y = Count, fill = InteractionType), color = 1) +
        scale_fill_manual(values = interaction_color, breaks = c("exclusion", "coexistence")) +
        scale_y_continuous(limits = c(0, nrow(example_pairs)), expand = c(0,0)) +
        theme_classic() +
        theme(axis.title.x = element_blank(), legend.position = "none",
              axis.text.x = element_blank(), axis.text.y = element_text(color = "black"),
              axis.title.y = element_text(size = 10), plot.background = element_blank(), panel.background = element_blank()) +
        labs(x = "", y = "Number of pairs", fill = "")

    node_size =5
    p_a <- tbl_graph(nodes = tibble(Isolate = 1:2, G = "g1", x = c(1, -1), y = c(0, 0)), edges = tibble(G = "g1", from = 1:2, to = 2:1, InteractionType = c("coexistence", "coexistence"))) %>%
        ggraph(layout = "nicely") +
        geom_node_point(size = node_size, shape = 21, fill = "gray", colour = "black", stroke = node_size/5) +
        geom_edge_link(aes(color = InteractionType), width = node_size/10,
                       arrow = arrow(length = unit(node_size/2, "mm"), type = "closed", angle = 30, ends = "last"),
                       start_cap = circle(node_size/2+1, "mm"),
                       end_cap = circle(node_size/2+1, "mm")) +
        annotate("text", x = 0, y = 0, vjust = 2, label = "coexistence", color = "#557BAA", size = 3) +
        scale_edge_color_manual(values = interaction_color) +
        scale_x_continuous(limits = c(-1.2, 1.2)) +
        scale_y_continuous(limits = c(-1.2, 1.2)) +
        theme_void() +
        theme(legend.position = "none", plot.margin = margin(0,0,0,0), plot.background = element_blank()) +
        labs(x = "")

    p_b <- tbl_graph(nodes = tibble(Isolate = 1:2, G = "g2", x = c(1, -1), y = c(0, 0)), edges = tibble(G = "g2", from = 2, to = 1, InteractionType = "exclusion")) %>%
        ggraph(layout = "nicely") +
        geom_node_point(size = node_size, shape = 21, fill = "gray", colour = "black", stroke = node_size/5) +
        geom_edge_link(aes(color = InteractionType), width = node_size/10,
                       arrow = arrow(length = unit(node_size/2, "mm"), type = "closed", angle = 30, ends = "last"),
                       start_cap = circle(node_size/2+1, "mm"),
                       end_cap = circle(node_size/2+1, "mm")) +
        annotate("text", x = 0, y = 0, vjust = 2, label = "exclusion", color = "#DB7469", size = 3) +
        scale_edge_color_manual(values = interaction_color) +
        scale_x_continuous(limits = c(-1.2, 1.2)) +
        scale_y_continuous(limits = c(-1.2, 1.2)) +
        theme_void() +
        theme(legend.position = "none", plot.margin = margin(0,0,0,0), plot.background = element_blank()) +
        labs(x = "")
    p_pair_example <- plot_grid(p_a, p_b)
    p4 <- plot_grid(p_pairs_interaction, p_pair_example, nrow = 2, axis = "lrbt", align = "vh", rel_heights = c(4,1))
    #ggsave("../plots/talk-example_bar.png", plot = p3, width = 2.5, height = 3)

    ## Arrow
    p <- tibble(x1 = -1, x2 = 1, y1 = 0, y2 = 0) %>%
        ggplot() +
        geom_segment(aes(x = x1, xend = x2, y = y1, yend = y2), arrow = arrow(ends = "last", type = "closed", angle = 45), size = 3) +
        theme_void()

    return(list(p2, p, p4, p, p3))
}

p1 <- plot_grid(plotlist = generate_example_graphs(5), nrow = 1, axis = "lrbt", align = "vh", scale = c(0.8, 1, 1, 1, 0.8), rel_widths = c(3, .5, 2, .5, 3))
p2 <- plot_grid(plotlist = generate_example_graphs(6), nrow = 1, axis = "lrbt", align = "vh", scale = c(0.8, 1, 1, 1, 0.8), rel_widths = c(3, .5, 2, .5, 3))
p3 <- plot_grid(plotlist = generate_example_graphs(7), nrow = 1, axis = "lrbt", align = "vh", scale = c(0.8, 1, 1, 1, 0.8), rel_widths = c(3, .5, 2, .5, 3))
p <- plot_grid(p1, p2, p3, ncol = 1, labels = paste0("Community ", 1:3))
ggsave("../plots/talk-example.png", plot = p, width = 12, height = 10)


# Network with hiererachy
load(here::here("data/output/network_community.Rdata"))


node_size = 10
example_graph <- net_list$C11R1 %>%
    activate("nodes") %>%
    arrange(PlotRank)
graph_layout <- create_layout(example_graph, layout = "auto") %>%
    arrange(PlotRank) %>%
    mutate(x = c(-.5, .5, 0, -1, 1, -2, 2, -1, 1),
           y = c(5, 5, 4, 3, 3, 2, 2, 1, 1))
example_graph %>%
    ggraph(layout = graph_layout) +
    geom_node_point(aes(x = x, y = y), size = node_size, shape = 21, fill = "gray", colour = "black", stroke = node_size/5) +
    geom_edge_link(aes(color = InteractionType), width = node_size/10,
                   arrow = arrow(length = unit(node_size/2, "mm"), type = "closed", angle = 30, ends = "last"),
                   start_cap = circle(node_size/2+1, "mm"),
                   end_cap = circle(node_size/2+1, "mm")) +
    scale_edge_color_manual(values = interaction_color) +
    theme_void() +
    theme(legend.position = "none", plot.margin = margin(0,0,0,0)) +
    labs(x = "")


p <- example_graph %>%
    ggraph(layout = "linear") +
    geom_node_point(size = node_size, shape = 21, fill = "gray", colour = "black", stroke = node_size/5) +
    #geom_node_text(aes(label = Isolate), size = 10) +
    geom_edge_arc(aes(color = InteractionType, alpha = InteractionType), edge_width = node_size / 10,
                  arrow = arrow(length = unit(node_size/2, "mm"), type = "open", angle = 30, ends = "last"),
                  start_cap = circle(node_size/2+1, "mm"), end_cap = circle(node_size/2+1, "mm")) +
    coord_flip() +
    scale_edge_color_manual(values = interaction_color) +
    scale_edge_alpha_manual(values = c("coexistence" = 0.3, "exclusion" = 0.9)) +
    scale_y_continuous(limits = c(-4, 2)) +
    scale_x_reverse() +
    scale_color_manual(values = c(`TRUE` = "#DB7469", `FALSE` = "#557BAA")) +
    guides(color = "none") +
    theme_void() +
    theme(legend.position = "right", plot.background = element_rect(fill = "white", color = NA),
          legend.title = element_blank(), plot.margin = margin(10,30,10,10, unit = "pt"),
          legend.text = element_text(size = 20)) +
    labs(x = "")
p
ggsave("../plots/talk-example_hiererchy.png", p, width = 9, height = 5)

# Expectation plot for coexisence pairs
p <- tibble(InteractionType = c("coexistence", "exclusion"), Count = c(130, 0)) %>%
    ggplot() +
    geom_col(aes(x = InteractionType, y = Count, fill = InteractionType), color = 1) +
    scale_fill_manual(values = interaction_color, breaks = c("exclusion", "coexistence")) +
    scale_y_continuous(limits = c(0, 150), expand = c(0,0)) +
    theme_classic() +
    theme(axis.title.x = element_blank(), legend.position = "top",
          axis.text.x = element_text(size = 10, color = "black"), axis.text.y = element_text(color = "black"),
          axis.title.y = element_text(size = 10)) +
    labs(x = "", y = "Number of pairs", fill = "")
ggsave(here::here("plots/talk-pairs_expectation.png"), p, width = 3, height = 4)


# Example of pairwise mechanisms
plot_pairwise_mechanisms <- function(graph) {
    graph %>%
        ggraph(layout = "nicely") +
        geom_node_point(aes(color = Color, shape = Shape), size = node_size, fill = "gray", stroke = node_size/5) +
        geom_edge_link(aes(color = InteractionType), width = node_size/3,
                       arrow = arrow(length = unit(node_size/2, "mm"), type = "closed", angle = 45, ends = "last"),
                       start_cap = circle(node_size/2+3, "mm"),
                       end_cap = circle(node_size/2+3, "mm")) +
        scale_edge_color_manual(values = interaction_color) +
        scale_color_manual(values = c("black" = "black", "white" = NA)) +
        scale_shape_manual(values = c("real" = 21, "fake" = NA)) +
        scale_x_continuous(limits = c(-2, 2)) +
        scale_y_continuous(limits = c(-1.5, 1)) +
        theme_void() +
        theme(legend.position = "none", plot.background = element_blank()) +
        labs(x = "")
}
node_size = 10
p_a <- tbl_graph(nodes = tibble(Isolate = 1:2, G = "g2", x = c(1, -1), y = c(-0.5, -0.5), Color = rep("black", 2), Shape = rep("real", 2)),
                 edges = tibble(G = "g2", from = 2, to = 1, InteractionType = "exclusion")) %>%
    plot_pairwise_mechanisms()
p_b <- tbl_graph(nodes = tibble(Isolate = 1:3, G = "g2", x = c(-1, 1, 0), y = c(0, 0, -1), Color = rep("black", 3), Shape = rep("real", 3)),
                 edges = tibble(G = "g2", from = c(1, 2), to = c(2, 3), InteractionType = "exclusion")) %>%
    plot_pairwise_mechanisms()
p_c <- tbl_graph(nodes = tibble(Isolate = 1:4, G = "g2", x = c(-1, 1, 0, 0), y = c(0, 0, -1, 0), Color = c(rep("black", 3), "white"), Shape = c(rep("real", 3), "fake")),
                 edges = tibble(G = "g2", from = c(1, 2, 3), to = c(2, 1, 4), InteractionType = c("coexistence", "coexistence", "bistability"))) %>%
    plot_pairwise_mechanisms()


p <- plot_grid(p_a, p_b, p_c, nrow = 1,
               #labels = c("Pairwise interaction", "Interaction chain", "Higher-order interaction"),
               axis = "lrtb", align = "vh", label_size = 20) +
    theme(plot.background = element_rect(fill = NA, color = NA))
p
ggsave(here::here("plots/talk-pairwise_mechanisms.png"), p, width = 9, height = 3)


set.seed(2)
p <- tbl_graph(nodes = tibble(Isolate = 1:6),
          edges = as_tibble(t(combn(6, 2))) %>%
              setNames(c("From", "To")) %>%
              mutate(InteractionType = sample(c("coexistence", "exclusion"), size = n(), replace = T))) %>%
    ggraph(layout = "circle") +
    geom_node_point(shape = 21, size = node_size, fill = "gray", stroke = node_size/5) +
    geom_edge_link(aes(color = InteractionType), width = node_size/5,
                   arrow = arrow(length = unit(node_size/2, "mm"), type = "closed", angle = 30, ends = "last"),
                   start_cap = circle(node_size/2+3, "mm"),
                   end_cap = circle(node_size/2+3, "mm")) +
    scale_edge_color_manual(values = interaction_color) +
    scale_color_manual(values = c("black" = "black", "white" = NA)) +
    scale_shape_manual(values = c("real" = 21, "fake" = NA)) +
    scale_x_continuous(limits = c(-1.2, 1.2)) +
    scale_y_continuous(limits = c(-1.2, 1.2)) +
    theme_void() +
    theme(legend.position = "none", plot.background = element_rect(fill = NA, color = NA), plot.margin = margin(0,0,0,0, unit = "pt")) +
    labs(x = "")

ggsave(here::here("plots/talk-pairwise_mechanisms_network.png"), p, width = 3, height = 3)






