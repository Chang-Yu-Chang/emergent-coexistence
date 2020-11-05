# Make the example plots

library(tidyverse)
library(data.table)
library(tidygraph)
library(ggraph)
library(cowplot)
source("network_functions.R")

# Panel XX: diagram of motifs
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
    plot_competitive_network(layout = "example_motif", node_size = 5) + 
    facet_nodes(~graph, nrow = 1)

ggsave("../plots/Ex1_motifs.png", plot = p1, width = 10, height = 1.5)
save(example_motif_list, file = "../data/temp/example_motif_list.Rdata")

# Panel XX: Make networks
isolates <- read_csv("../data/output/isolates.csv")
isolates_melted <- read_csv("../data/output/isolates_melted.csv")
pairs <- read_csv("../data/output/pairs.csv")

set.seed(1)
example_pairs <- as_tibble(t(combn(6,2))) %>%
    setNames(c("From", "To")) %>%
    mutate(InteractionType = sample(c("coexistence", "exclusion"), size = n(), replace = T))
example_isolates <- tibble(Isolate = 1:6, Rank = 1:6, PlotRank = 1:6)
example_graph <- make_network(isolates = example_isolates, pairs = example_pairs)

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
p2 <- plot_competitive_network2(example_graph, layout = "circle")
p2
ggsave("../plots/Ex2_graph.png", plot = p2, width = 5, height = 5)



