#
library(tidyverse)
library(tidygraph)
library(tidymodels)
library(ggraph)
library(cowplot)
library(ggsci)
library(ggsignif)
source(here::here("plotting_scripts/network_functions.R"))
interaction_type <- c("exclusion", "coexistence", "lose", "bistability", "neutrality", "self", "undefined")
interaction_color = c("#DB7469", "#557BAA", "#73C966", "#EECF6D", "#8650C4", "black", "grey80")
names(interaction_color) <- interaction_type
add_interaction_direction <- function(pairs) {
    pairs %>%
        mutate(from = ifelse(InteractionType == "exclusion" & From == "species1", Species1, Species2),
               to = ifelse(InteractionType == "exclusion" & From == "species1", Species2, Species1),
               from = ifelse(InteractionType == "coexistence", Species1, from),
               to = ifelse(InteractionType == "coexistence", Species2, to)) %>%
        select(-From)
}


l = 0.5
q = 0.5
v = 1
S = 300
#
pairs_poolPairs <- read_csv("~/Dropbox/lab/invasion-network/simulation/data/temp/pairs_poolPairs.csv") %>%
    filter(l1 == l, q2 == q, vamp == v) %>% add_interaction_direction
pairs_randomNetworks <- read_csv("~/Dropbox/lab/invasion-network/simulation/data/temp/pairs_randomNetworks.csv") %>%
    filter(l1 == l, q2 == q, vamp == v) %>% add_interaction_direction
communities_selfAssembly <- read_csv("~/Dropbox/lab/invasion-network/simulation/data/temp/communities_selfAssembly.csv") %>%
    filter(l1 == l, q2 == q, vamp == v)
pairs_communityPairs <- read_csv("~/Dropbox/lab/invasion-network/simulation/data/temp/pairs_communityPairs.csv") %>%
    filter(l1 == l, q2 == q, vamp == v) %>% add_interaction_direction

# Subset community pairs that are "within" communities
communities_species <- communities_selfAssembly %>%
    arrange(exp_id, Well, Species) %>%
    select(exp_id, Community = Well, Species)
generate_pairs <- function(x) {
    if (length(x) <= 1) return(tibble(Species1 = character(), Species2 = character()))
    if (length(x) >= 2) {
        x %>%
            ordered(level = paste0("S", 0:(S-1))) %>%
            sort() %>%
            as.character() %>%
            combn(2) %>%
            t() %>%
            as_tibble %>%
            setNames(c("Species1", "Species2")) %>%
            return()
    }
}
pairs_communityPairs_within <- communities_species %>%
    group_by(Community) %>%
    summarize(Pairs = list(generate_pairs(Species))) %>%
    unnest(cols = c(Pairs)) %>%
    left_join(pairs_communityPairs)

# Combine all pairs data
pairs <- bind_rows(pairs_poolPairs,
                   mutate(pairs_randomNetworks, Experiment = "poolPair"),
                   pairs_communityPairs_within)


# Compare the pair interactions  across pools
p <- pairs %>%
    mutate(Experiment = ordered(Experiment, c("poolPair", "communityPairs"))) %>%
    group_by(Experiment, InteractionType) %>% summarize(Count = n()) %>%
    ggplot() +
    geom_col(aes(x = Experiment, y = Count, fill = InteractionType), position = "fill", color = 1) +
    geom_text(data = group_by(pairs, Experiment) %>% summarize(Count = n()), aes(x = Experiment, label = paste0("n = ", Count)), y = 1, vjust = 2) +
    #facet_grid(.~PairFermenter, scales = "free_y", labeller = label_both) +
    scale_fill_manual(values = interaction_color, breaks = c("coexistence", "exclusion")) +
    scale_y_continuous(breaks = c(0, 0.5, 1), expand = c(0,0)) +
    scale_x_discrete(labels = c("poolPair" = "Species Pool", "communityPairs" = "Communities"),
                     breaks = c("poolPair", "communityPairs")) +
    theme_classic() +
    theme(legend.position = "top", axis.text.x = element_text(size = 12),
          axis.title.y = element_text(size = 12), axis.text.y = element_text(size = 10),
          legend.text = element_text(size = 12)) +
    labs(x = "", y = "Fraction", fill = "")
p
ggsave("../plots/Fig_simulation-pairs.png", p, width = 4, height = 5)

## Stat
observed_stat <- pairs %>%
    filter(Experiment %in% c("poolPair", "communityPairs")) %>%
    chisq_test(InteractionType ~ Experiment) %>% pull(statistic)
null_stat <- pairs %>%
    filter(Experiment %in% c("poolPair", "communityPairs")) %>%
    specify(InteractionType ~ Experiment, success = "coexistence") %>%
    hypothesize(null = "independence") %>%
    generate(reps = 1000, type = "permute") %>%
    calculate(stat = "Chisq", order = c("poolPair", "communityPairs"))
null_stat %>%
    get_p_value(obs_stat = observed_stat, direction = "right")




# Random networks
make_network_from_pairs <- function(from, to, InteractionType ) {
    pairs <- tibble(from = from , to = to, InteractionType = InteractionType)
    # Nodes
    nodes <- pairs %>%
        select(from, to) %>%
        unlist() %>% unique %>% ordered(paste0("S", 0:(S-1))) %>% sort() %>%
        tibble(Isolate = 1:length(.), Species = .)

    # Edges
    edges <- pairs %>% select(from, to, InteractionType) %>%
        mutate(from = match(from, nodes$Species),
               to = match(to, nodes$Species))

    edges_coext <- edges[edges$InteractionType == "coexistence",]
    edges_coext[,c("from", "to")] <- edges_coext[,c("to", "from")] # Add the mutual edges for coexistence links
    edges <- rbind(edges, edges_coext)

    # Network
    graph <- tbl_graph(nodes = nodes, edges = edges, directed = T)

    return(graph)
}
df_graph_random <- pairs_randomNetworks  %>%
    mutate(Community = rep(1:10,  each = choose(6, 2)) %>% as.character()) %>%
    select(Experiment, Community, from, to, InteractionType, PairFermenter) %>%
    group_by(Experiment, Community) %>%
    summarize(Graph = list(make_network_from_pairs(from, to, InteractionType))) %>%
    rowwise() %>%
    mutate(Motif = list(count_motif(Graph)))

# Community networks
df_graph_comm <- pairs_communityPairs_within %>%
    select(Experiment, Community, from, to, InteractionType, PairFermenter) %>%
    group_by(Experiment, Community) %>%
    summarize(Graph = list(make_network_from_pairs(from, to, InteractionType))) %>%
    rowwise() %>%
    mutate(Motif = list(count_motif(Graph)))

# Combine the networks
df_graph <- bind_rows(df_graph_comm, df_graph_random)

df_motif <- df_graph %>%
    select(Experiment, Community, Motif) %>%
    unnest_longer(col = Motif, indices_to = "Motif", values = "Count") %>%
    group_by(Experiment, Community) %>%
    mutate(Fraction = Count / sum(Count))

# Plot motif
p1 <- df_motif %>%
    mutate(Community = factor(Community), Experiment = ordered(Experiment, c("randomNetworks", "communityPairs"))) %>%
    mutate(Motif = as.character(Motif)) %>%
    filter(Motif %in% c(2)) %>%
    ggplot(aes(x = Experiment, y = Fraction, color = Experiment)) +
    geom_boxplot(position = position_dodge(width = 1)) +
    geom_jitter(shape = 21, size = 2, width = 0.2, height = 0) +
    #geom_signif(y_position = 1.1, xmin = 1, xmax = 2, annotation = c("NS"), color = "black") +
    geom_signif(comparisons = list(c("randomNetworks", "communityPairs")), test = "wilcox.test", map_signif_level = T, color = "black", tip_length = 0.01) +
    scale_y_continuous(breaks = c(0, 0.5, 1), limits = c(0, 1.1)) +
    scale_color_manual(values = c("randomNetworks" = "orange", "communityPairs" = "cornflower blue")) +
    scale_x_discrete(labels = c("communityPairs" = "Communities", "randomNetworks" = "Random\nspecies")) +
    #facet_grid(.~Experiment, scales = "free_x") +
    theme_classic() +
    theme(legend.position = "none", legend.title = element_blank(), strip.text = element_blank(),
          panel.spacing = unit(0, "pt"),
          axis.text.x = element_text(size = 12),
          axis.title.y = element_text(size = 12), axis.text.y = element_text(size = 12),
          legend.text = element_text(size = 12)) +
    labs(x = "", y = "Hierarchy")
p1
ggsave("../plots/Fig_simulation-motif_pool.png", plot = p1, width = 3, height = 5)


# Stat

observed_stat <- df_motif %>%
    filter(Motif %in% c(2)) %>%
    t_test(Fraction ~ Experiment, order = c("communityPairs", "randomNetworks"))

null_stat <- df_motif %>%
    filter(Motif %in% c(2)) %>%
    specify(Fraction ~ Experiment) %>%
    hypothesize(null = "independence") %>%
    generate(reps = 1000, type = "permute") %>%
    calculate(stat = "t", order = c("communityPairs", "randomNetworks"))

null_stat %>%
    get_p_value(obs_stat = observed_stat, direction = "right")




if (FALSE) {
    node_size = 3
    load(here::here("data/output/motif_list.Rdata"))
    p_motif1 <- plot_competitive_network(motif_list[[1]], node_size = node_size) + theme(plot.margin = margin(0,0,0,0))
    p_motif2 <- plot_competitive_network(motif_list[[2]], node_size = node_size) + theme(plot.margin = margin(0,0,0,0))
    p2 <- plot_grid(p_motif1, p_motif2, nrow = 1) + theme(plot.background = element_rect(color = NA, fill = "white"))
    p <- plot_grid(p1, p2, ncol = 1, axis = "lrtb", align = "vh", rel_heights = c(3,1))
    ggsave("../plots/Fig_simulation-motif_pool.png", plot = p, width = 4, height = 5)

    library(igraph)
    g <- df_graph$Graph[1][[1]]


    igraph::cohesive_blocks(g)

    load(here::here("data/output/network_community.Rdata"))
    net_list[[1]] %>%
        as.undirected() %>%
        hierarchy()
    cohesive.blocks()

    transitivity(net_list[[1]])
    is.dag(net_list[[1]])








    node_size = 10
    df_graph_random$graph[[4]] %>%
        ggraph(layout = "circle") +
        geom_node_point(size = node_size, shape = 21, fill = "gray", colour = "black", stroke = node_size/5) +
        geom_edge_link(aes(color = InteractionType), width = node_size/10,
                       arrow = arrow(length = unit(node_size/2, "mm"), type = "closed", angle = 30, ends = "last"),
                       start_cap = circle(node_size/2+1, "mm"),
                       end_cap = circle(node_size/2+1, "mm")) +
        scale_edge_color_manual(values = interaction_color) +
        theme_void() +
        theme(legend.position = "none", plot.margin = margin(0,0,0,0)) +
        labs(x = "")


}





















