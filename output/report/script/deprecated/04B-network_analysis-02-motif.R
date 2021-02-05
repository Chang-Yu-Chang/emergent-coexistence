#'
#'



# Plot the motif examples ----
# Motif list
motif_list <- rep(list(NA), 7)
temp_id <- c(11, 7, 8, 12, 13, 14, 15)
g <- igraph::graph.isocreate(size = 3, 1)
layout <- create_layout(g, layout = 'circle')


for (i in 1:7) {
  g <- igraph::graph.isocreate(size = 3, temp_id[i])
  E(g)$InteractionType <- ifelse(which_mutual(g), "coexistence", "exclusion")
  g <- tidygraph::as_tbl_graph(g) %>%
    mutate(x = layout$x, y = layout$y, graph = paste0("motif", i))
  motif_list[[i]] <- g
}

# Plot
set_graph_style(plot_margin = margin(1,1,1,1))
# Color setting
interaction_type <- c("exclusion", "coexistence", "lose", "bistability", "neutrality", "self", "undefined")
myColor = c("#DB7469", "#557BAA", "#73C966", "#EECF6D", "#8650C4", "black", "grey80")
names(myColor) <- interaction_type
node_size <- 5

p_motif_example <-
  tidygraph::bind_graphs(motif_list[[1]], motif_list[[2]], motif_list[[3]], motif_list[[4]], motif_list[[5]], motif_list[[6]], motif_list[[7]]) %>%
  ggraph(layout = "nicely") +
  geom_edge_link(aes(color = InteractionType), width = 1,
    arrow = arrow(length = unit(node_size, "mm"), type = "open", angle = 20),
    start_cap = circle(node_size, "mm"),
    end_cap = circle(node_size, "mm")) +
  geom_node_point(size = node_size + 2, color = "black") +
  scale_edge_color_manual(values = myColor) +
  facet_nodes(~graph, nrow = 1) +
  theme_graph() +
  theme(
    legend.position = "none",
    legend.direction = "none",
    legend.title = element_blank(),

    strip.text = element_blank()
  )
