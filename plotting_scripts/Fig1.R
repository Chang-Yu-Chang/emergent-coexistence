# Simulate the network topography

library(tidygraph)
library(ggraph)
library(tidyverse)
library(data.table)
library(cowplot)
source("network_functions.R")

n_sizes <- c(4, 8, 12) # Number of species
p_range <- seq(0, 1, by = 0.1) # Probability of pairwise coexistence
b = 100 # Number of seeds

temp_list <- rep(list(NA), length(n_sizes) * length(p_range) * b)
counter = 1

for (i in 1:length(n_sizes)) {
    for (j in 1:length(p_range)) {
        for (k in 1:b) {
            set.seed(k)
            temp_list[[counter]] <- make_random_network(n = n_sizes[i], p = p_range[j]) %>%
                count_motif() %>% 
                tibble(Motif = factor(1:7), Count = .) %>% 
                mutate(CommunitySize = n_sizes[i], FractionCoexistence = p_range[j], Seed = k)
            counter = counter + 1
            if (counter%%100 == 0) cat(counter, " ")
        }
    }
}
motif_counts <- temp_list %>% rbindlist()

motif_counts %>%
    ggplot(aes(x = Seed, y = Count, fill = Motif)) +
    geom_bar(position = "stack", stat = "identity") +
    scale_x_continuous(expand = c(0,0)) +
    scale_y_continuous(expand = c(0,0)) +
    facet_grid(p~.) + 
    theme_bw()

# Panel XX: motif distribution as a function of pairwise coexistence
motif_count_mean <- motif_counts %>% 
    group_by(CommunitySize, FractionCoexistence, Motif) %>% 
    summarize(MeanCount = mean(Count))

motif_count_mean %>% 
    ggplot() +
    geom_area(aes(x = FractionCoexistence, y = MeanCount, fill = Motif), color = 1) +
    scale_x_continuous(expand = c(0,0)) +
    scale_y_continuous(expand = c(0,0)) +
    facet_grid(CommunitySize ~ ., scales = "free_y") +
    theme_cowplot()

motif_type <- c("others", "0-scored")
motif_color = c("#DB7469", "#557BAA")
names(motif_color) <- motif_type

p1 <- motif_counts %>% 
    #mutate(MotifType = ifelse(Motif %in% c(1, 7), "0-scored", "others")) %>% 
    group_by(CommunitySize, FractionCoexistence, Motif) %>% 
    summarize(MeanCount = mean(Count)) %>% 
    filter(CommunitySize == 12) %>% 
    ggplot() +
    geom_area(aes(x = FractionCoexistence, y = MeanCount, alpha = Motif), color = 1) +
    scale_x_continuous(expand = c(0,0), breaks = c(0, 0.5, 1)) +
    scale_y_continuous(expand = c(0,0), breaks = c(0, 100, 200)) +
    #    facet_grid(CommunitySize ~ ., scales = "free_y") +
    theme_cowplot() +
    theme(legend.title = element_blank(),
        legend.position = "top", legend.direction = "horizontal") +
    labs(x = "Fraction of pairwise coexistence", y = "Mean motif count")


p2 <- motif_counts %>% 
    mutate(MotifType = ifelse(Motif %in% c(1, 7), "0-scored", "others")) %>% 
    group_by(CommunitySize, Seed, FractionCoexistence, MotifType) %>% 
    summarize(Count = sum(Count)) %>% 
    group_by(CommunitySize, FractionCoexistence, MotifType) %>% 
    summarize(MeanCount = mean(Count)) %>% 
    filter(CommunitySize == 12) %>% 
    ggplot() +
    geom_area(aes(x = FractionCoexistence, y = MeanCount, fill = MotifType), color = 1) +
    scale_x_continuous(expand = c(0,0), breaks = c(0, 0.5, 1)) +
    scale_y_continuous(expand = c(0,0), breaks = c(0, 100, 200)) +
    scale_fill_manual(values = motif_color) +
    theme_cowplot() +
    theme(legend.title = element_blank(),
        legend.position = "top") +
    labs(x = "Fraction of pairwise coexistence", y = "Mean motif count")


p <- plot_grid(p1, p2, nrow = 1, align = "hv")
ggsave("../plots/Fig1A.png", plot = p, width = 8, height = 4)
ggsave("../plots/Fig1A_example.png", plot = p2, width = 4, height = 4)




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
    
    p1 <- motif_count_mean %>%  
        ggplot(aes(x = p, y = MeanCount, color = Motif, group = Motif)) +
        geom_point() + geom_line() +
        theme_bw()
    
}


