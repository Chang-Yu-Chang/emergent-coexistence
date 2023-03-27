#' Script for supplement figures

library(tidyverse)
library(cowplot)
library(ggraph)
library(tidygraph)
source(here::here("analysis/00-metadata.R"))

communities <- read_csv(paste0(folder_data, "temp/00c-communities.csv"), show_col_types = F)
load(paste0(folder_data, "temp/95-communities_network.Rdata"))



# Figure S10 All 13 self-assembled community graphs ----
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
        geom_edge_link(aes(color = InteractionType), width = edge_width/2, arrow = arrow()) +
        geom_node_point(fill = "white", size = node_size*1.2, shape = 21, colour = "black", stroke = node_size/3) +
        scale_edge_color_manual(values = interaction_color) +
        scale_x_continuous(limits = nodes_axis_x*1) +
        scale_y_continuous(limits = nodes_axis_y*1) +
        theme_graph() +
        theme(
            legend.position = "none",
            legend.title = element_blank(),
            strip.text = element_blank(),
            plot.margin = unit(c(3,3,3,3), "mm")
            # plot.background = element_rect(fill = "grey90", color = NA),
            # panel.background = element_rect(fill = "grey90", color = NA)
        )
}
p_net_list <- communities_network %>%
    mutate(Community = factor(Community, Community)) %>%
    arrange(CommunitySize) %>%
    mutate(NetworkPlotSize = max(CommunitySize) / CommunitySize) %>%
    rowwise() %>%
    mutate(NetworkPlot = plot_competitive_network_grey(Network, NetworkPlotSize, NetworkPlotSize) %>% list()) %>%
    pull(NetworkPlot)
p_network <- plot_grid(plotlist = p_net_list, nrow = 2, scale = .9, labels = 1:13) + paint_white_background()
p_legend <- get_legend({tibble(InteractionType = c("coexistence", "exclusion", "unknown"), x = 1:3, y = 1:3) %>%
        ggplot() +
        geom_line(aes(color = InteractionType, x = x, y = y, group = InteractionType), linewidth = 2) +
        scale_color_manual(values = interaction_color) +
        theme(legend.position = "right",
              legend.title = element_blank(),
              legend.key.size = unit(.8, "cm"),
              legend.text = element_text(size = 12))})
p <- ggdraw(p_network) +
    draw_plot(p_legend, x = .93, y = .25, hjust = 0.5, vjust = .5)

ggsave(here::here("plots/FigS91-community_graph.png"), p, width = 13, height = 4)


communities_network$Network[[9]] %>%
    activate(edges) %>%
    #filter(InteractionType == "exclusion") %>%
    triad_census()
    plot_competitive_network_grey(5, 5)

