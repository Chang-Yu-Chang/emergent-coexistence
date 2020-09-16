# Simulate the network topography

library(tidygraph)
library(ggraph)
library(tidyverse)
library(data.table)
source("network_functions.R")



n = 10 # Number of species
p_range <- seq(0, 1, by = 0.1)
#p = 0.5 # Probability of pairwise coexistence
b = 100

temp_list <- rep(list(rep(list(NA), b)), length(p_range))
names(temp_list) <- p_range

for (j in 1:length(p_range)) {
    cat("\np =", p_range[j], "\n")
    for (i in 1:b) {
        temp_list[[j]][[i]] <- count_motif(make_random_network(n = n, p = p_range[j]))
        if (i%%10 == 0) cat(i, " ")
    }
}

motif_counts <- temp_list %>%
    lapply(function(x) {
        lapply(x, function(y) {tibble(Motif = factor(1:7), Count = y)}) %>%
            rbindlist(idcol = "Seed")
    }) %>% 
    rbindlist(idcol = "p")


motif_counts %>%
    ggplot(aes(x = Seed, y = Count, fill = Motif)) +
    geom_bar(position = "stack", stat = "identity") +
    scale_x_continuous(expand = c(0,0)) +
    scale_y_continuous(expand = c(0,0)) +
    facet_grid(p~.) + 
    theme_bw()


# Panel XX: motif distribution as a function of coexistence fraction
p1 <- motif_counts %>% 
    group_by(p, Motif) %>% 
    summarize(MeanCount = mean(Count)) %>% 
    ggplot(aes(x = p, y = MeanCount, color = Motif, group = Motif)) +
    geom_point() + geom_line() +
    theme_bw()

ggsave("../plots/Fig1A.png", plot = p1, width = 4, height = 3)


# Panel XX. matrix representation 
plot_graph_list <- graph_list %>% lapply(plot_adjacent_matrix)
p2 <-  plot_grid(plotlist = plot_graph_list)

ggsave("../plots/Fig1B.png", plot = p2, width = 4, height = 3)



# Panel XX. network 
plot_competitive_network(graph_list[[11]], layout = "circle")
plot_competitive_network(graph_list[[11]], layout = "circle")

#







if (FALSE){
ggraph(graph, layout = "linear", circular = T) +
    geom_node_point() +
    geom_edge_arc(aes(color = interaction), show.legend = T) +
    theme_graph()
}


