library(tidyverse)
library(cowplot)
library(tidygraph)
library(ggraph)
source(here::here("analysis/00-metadata.R"))

communities <- read_csv(paste0(folder_data, "output/communities_remained.csv"), show_col_types = F)
isolates <- read_csv(paste0(folder_data, "output/isolates_remained.csv"), show_col_types = F)
pairs <- read_csv(paste0(folder_data, "output/pairs_remained.csv"), show_col_types = F)

make_network <- function(isolates_comm, pairs_comm) {
    # Nodes
    nodes <- isolates_comm

    # Edges
    edges <- pairs_comm %>%
        filter(outcome %in% c("1-exclusion", "2-exclusion", "3-coexistence", "4-coexistence", "5-inconclusive")) %>%
        #' tbl_graph cannot handle it when a link connects from 1 to 12 (In pairs.csv this means the Isolate)
        #' but there are only 10 isolates. For instance, C11R2 has 10 isolates from 1-12 but 9 and 11
        #' are removed because of bad Sanger-ESV alignments. This creates incontinuous numbered isolates.
        #' To account for this, the `from` and `to` in the edge tibble of network object means the `row` number in the node tibble
        mutate(from = match(From, isolates_comm$Isolate), to = match(To, isolates_comm$Isolate)) %>%
        select(from, to, outcome, PairID)
    edges_coext <- edges[edges$outcome %in% c("3-coexistence","4-coexistence"),]
    edges_coext[,c("from", "to")] <- edges_coext[,c("to", "from")] # Add the mutual edges for coexistence links
    edges <- bind_rows(edges, edges_coext)

    # Somehow it works
    nodes <- nodes %>%
        mutate(Isolate = 1:n())

    # Network
    graph <- tbl_graph(nodes = nodes, edges = edges, directed = T)


    return(graph)
}

isolates_comm <- filter(isolates, Community == "C11R2")
pairs_comm <- filter(pairs, Community == "C11R2")
# graph %>% activate(nodes) %>% filter(!is.na(Community))

# Check to make network object from isolate and pair ----
communities_network <- communities %>%
    rename(comm = Community) %>%
    rowwise() %>%
    mutate(Isolates = isolates %>% filter(Community == comm) %>% list) %>%
    mutate(Pairs = pairs %>% filter(Community == comm) %>% list) %>%
    mutate(Network = make_network(Isolates, Pairs) %>% list) %>%
    rename(Community = comm) %>%
    select(Community, CommunityLabel, CommunitySize, CommunityPairSize, Network) %>%
    ungroup

save(communities_network, file = paste0(folder_data, "temp/38-communities_network.Rdata"))

# Check the number of links is correct in plooted networks ----
plot_competitive_network_grey <- function(g, node_size = 10, edge_width = 1){
    # Layout
    graph_layout <- create_layout(g, "circle")
    mean_x_coord <- mean(graph_layout$x)
    mean_y_coord <- mean(graph_layout$x)
    g <- g %>% activate(nodes) %>% mutate(x = graph_layout$x - mean_x_coord, y = graph_layout$y - mean_y_coord)

    # Axis range
    nodes_axis_x <- (activate(g, nodes) %>% pull(x) %>% range()) * 1.1
    nodes_axis_y <- (activate(g, nodes) %>% pull(y) %>% range()) * 1.1

    # Graph
    g %>%
        ggraph(layout = "nicely") +
        geom_edge_link(aes(color = outcome), width = edge_width/2,
                       arrow = arrow(length = unit(2, "mm"), type = "closed", angle = 30, ends = "last"),
                       start_cap = circle(node_size, "mm"),
                       end_cap = circle(node_size, "mm")) +
        geom_node_point(fill = "white", size = node_size*1.2, shape = 21, colour = "black", stroke = node_size/3) +
        scale_edge_color_manual(values = outcome_colors, label = outcome_labels) +
        scale_x_continuous(limits = nodes_axis_x*1) +
        scale_y_continuous(limits = nodes_axis_y*1) +
        theme_graph() +
        theme(
            legend.position = "none",
            legend.title = element_blank(),
            strip.text = element_blank(),
            plot.margin = unit(c(3,3,3,3), "mm")
        ) +
        guides() +
        labs()
}
p_net_list <- communities_network %>%
    mutate(Community = factor(Community, Community)) %>%
    arrange(CommunityLabel) %>%
    mutate(NetworkPlotSize = max(CommunitySize) / CommunitySize) %>%
    rowwise() %>%
    mutate(NetworkPlot = plot_competitive_network_grey(Network, NetworkPlotSize, NetworkPlotSize) %>% list()) %>%
    pull(NetworkPlot)
p <- plot_grid(plotlist = p_net_list, nrow = 2, scale = .9, labels = 1:13) + paint_white_background()
ggsave(paste0(folder_data, "temp/38-01-graph.png"), p, width = 13, height = 4)
#ggsave(here::here("plots/FigS91-community_graph.png"), p, width = 13, height = 4)

# Check the number of links is correct in exclusion networks ----
subset_exclusion <- function (graph) {
    graph %>%
    activate(edges) %>%
    filter(outcome %in% c("1-exclusion", "2-exclusion"))
}
communities_network_exclusion <- communities_network %>%
    mutate(Community = factor(Community, Community)) %>%
    arrange(CommunityLabel) %>%
    mutate(NetworkPlotSize = max(CommunitySize) / CommunitySize) %>%
    rowwise() %>%
    mutate(ExclusionNetwork = subset_exclusion(Network) %>% list) %>%
    mutate(ExclusionNetworkPlot = plot_competitive_network_grey(ExclusionNetwork, NetworkPlotSize, NetworkPlotSize) %>% list())

p_net_list <- communities_network_exclusion$ExclusionNetworkPlot
p <- plot_grid(plotlist = p_net_list, nrow = 2, scale = .9, labels = 1:12) + paint_white_background()

ggsave(paste0(folder_data, "temp/38-02-graph_exclusion.png"), p, width = 13, height = 4)

# Number of the exclusion pairs per community. Manually checked
pairs %>%
    left_join(communities) %>%
    filter(outcome %in% c("1-exclusion", "2-exclusion")) %>%
    mutate(outcome = factor(outcome, c("1-exclusion", "2-exclusion"))) %>%
    group_by(Community, CommunityLabel, outcome, .drop = F) %>%
    count() %>%
    arrange(CommunityLabel, outcome) %>%
    pivot_wider(names_from = outcome, values_from = n)

# Check the isolate identity ----
plot_competitive_network_family <- function(g, node_size = 10, edge_width = 1){
    # Layout
    graph_layout <- create_layout(g, "circle")
    mean_x_coord <- mean(graph_layout$x)
    mean_y_coord <- mean(graph_layout$x)
    g <- g %>% activate(nodes) %>% mutate(x = graph_layout$x - mean_x_coord, y = graph_layout$y - mean_y_coord)

    # Axis range
    nodes_axis_x <- (activate(g, nodes) %>% pull(x) %>% range()) * 1.1
    nodes_axis_y <- (activate(g, nodes) %>% pull(y) %>% range()) * 1.1

    # Graph
    g %>%
        ggraph(layout = "nicely") +
        geom_edge_link(aes(color = outcome), width = edge_width/2,
                       arrow = arrow(length = unit(2, "mm"), type = "closed", angle = 30, ends = "last"),
                       start_cap = circle(node_size, "mm"),
                       end_cap = circle(node_size, "mm")) +
        geom_node_point(aes(fill = Family), size = node_size*1.2, shape = 21, colour = "black", stroke = node_size/3) +
        geom_node_text(aes(label = Isolate), size = node_size*1.2) +
        scale_edge_color_manual(values = outcome_colors, label = outcome_labels) +
        scale_fill_manual(values = family_colors) +
        scale_x_continuous(limits = nodes_axis_x*1) +
        scale_y_continuous(limits = nodes_axis_y*1) +
        theme_graph() +
        theme(
            legend.position = "none",
            legend.title = element_blank(),
            strip.text = element_blank(),
            plot.margin = unit(c(3,3,3,3), "mm")
        ) +
        guides() +
        labs()
}

communities_network_exclusion <- communities_network %>%
    mutate(Community = factor(Community, Community)) %>%
    arrange(CommunityLabel) %>%
    mutate(NetworkPlotSize = max(CommunitySize) / CommunitySize) %>%
    rowwise() %>%
    mutate(ExclusionNetwork = subset_exclusion(Network) %>% list) %>%
    mutate(ExclusionNetworkPlot = plot_competitive_network_family(ExclusionNetwork, NetworkPlotSize, NetworkPlotSize) %>% list())

p_net_list <- communities_network_exclusion$ExclusionNetworkPlot
p_net_list[[13]] <- get_legend(plot_competitive_network_family(bind_graphs(communities_network_exclusion$ExclusionNetwork), node_size = 2) + theme(legend.position = "right", legend.spacing.y = unit(2, "mm"), legend.key.size = unit(3, "mm")))
p <- plot_grid(plotlist = p_net_list, nrow = 2, scale = .9, labels = 1:12) + paint_white_background()
ggsave(paste0(folder_data, "temp/38-03-graph_exclusion_family.png"), p, width = 13, height = 4)

#
pairs %>%
    filter(Community == "C11R5") %>%
    select(Community, Isolate1, Isolate2, From, To, Genus1, Genus2, outcome)

isolates %>%
    filter(Community == "C11R5") %>%
    select(Community, Isolate, Genus) # C11R5 4 should be a entero not a pseudo
















