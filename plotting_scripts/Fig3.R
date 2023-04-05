library(tidyverse)
library(cowplot)
library(tidygraph)
library(ggraph)
library(igraph)
source(here::here("analysis/00-metadata.R"))

communities <- read_csv(paste0(folder_data, "output/communities_remained.csv"), show_col_types = F)
isolates <- read_csv(paste0(folder_data, "output/isolates_remained.csv"), show_col_types = F)
pairs <- read_csv(paste0(folder_data, "output/pairs_remained.csv"), show_col_types = F)
pairs_freq <- read_csv(paste0(folder_data, "temp/93a-pairs_freq.csv"), show_col_types = F)
load(paste0(folder_data, "temp/38-communities_network.Rdata"))


# Plot individual networks ----
set.seed(3)
node_size = 3
edge_width = .8
plot_network_hierarchy <- function(net, n_rank = 10, n_break = 10) {
    #net <- communities_network$Network[[13]]
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
#        filter(outcome != "5-inconclusive", !is.na(outcome)) %>%
#        filter(outcome == "1-exclusion" | outcome == "2-exclusion" | outcome == "3-coexistence" | outcome == "4-coexistence") %>%
        filter(outcome == "1-exclusion" | outcome == "2-exclusion") %>%
        mutate(Temp = sample(c(-1, 1), size = n(), replace = T))
    strength_angle <- as_tibble(temp)$Temp * 0.06

    temp %>%
        ggraph(layout = "nicely") +
        geom_hline(yintercept = c(-n_rank:-1), color = "grey90") +
        geom_node_point(size = node_size, shape = 21, fill = "grey", stroke = node_size/5, color = "black") +
        geom_node_text(aes(label = Isolate), size = 3) +
        geom_edge_arc(strength = strength_angle, alpha = 0.5,
                      aes(color = outcome), width = edge_width,
                      arrow = arrow(length = unit(edge_width, "mm"), type = "closed", angle = 30, ends = "last"),
                      start_cap = circle(node_size/2, "mm"),
                      end_cap = circle(node_size/2, "mm")) +
        scale_edge_color_manual(values = outcome_colors, labels = outcome_labels) +
        scale_x_continuous(limits = c(0.1, 0.9), expand = c(0,0)) +
        scale_y_continuous(limits = c(-n_break-1, 0), breaks = -n_break:-1, labels = n_break:1) +
        theme_void() +
        theme(
            legend.position = "none",
            legend.title = element_blank(),
            axis.title = element_blank(),
            strip.text = element_blank(),
            plot.margin = unit(c(0,0,0,0),"mm")
        ) +
        labs(y = "Rank")

}
communities_network_hierarchy <- communities_network %>%
    #mutate(Community = factor(Community, Community)) %>%
    arrange(CommunityLabel) %>%
    rowwise() %>%
    mutate(NetworkHierarchyPlot = plot_network_hierarchy(Network, n_rank = CommunitySize) %>% list())
#plot_network_hierarchy(communities_network_hierarchy$Network[[13]])

# Check the isolate identity in hierarchy ----
p_net_hierarchy_list <- communities_network_hierarchy$NetworkHierarchyPlot
p_net_hierarchy_list[[12]] <- p_net_hierarchy_list[[12]] +
    scale_y_continuous(limits = c(-10-1, 0), breaks = -10:-1, labels = 10:1, position = "right") +
    theme(axis.title.y = element_text(color = 1, size = 10, angle = 270, margin = margin(l = 2, unit = "mm")),
          axis.text.y = element_text(color = 1, size = 10, margin = margin(l = 1, unit = "mm")))
p_axistitle <- ggdraw() + draw_label("Community", fontface = 'plain', x = .5, hjust = .5) + theme(plot.margin = margin(5, 0, 5, 7))
p_temp <- plot_grid(plotlist = p_net_hierarchy_list,
                    rel_widths = c(communities_network_hierarchy$CommunitySize / max(communities_network_hierarchy$CommunitySize))^1.5,
                    labels = 1:12, label_fontface = "plain", label_x = c(rep(0.5, 11), 0.45), hjust = c(rep(.5, 11), 1),
                    nrow = 1, axis = "tb", align = "h") + paint_white_background()
p1 <- plot_grid(p_axistitle, p_temp, ncol = 1, rel_heights = c(.1, 1)) + paint_white_background()


# line legend
edge_width = 0.5
p_legend <- get_legend({
    p_net_hierarchy_list[[12]] +
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

p1 <- ggdraw(p1) +
    draw_plot(p_legend, x = 0.4, y = 0.15, width = 0.1, height = 0.1,  hjust = 0, vjust = 0)


# Bionomial distribution of nontranistivity ----
subset_exclusion_network <- function (net) {
    net %>%
        activate(edges) %>%
        filter(outcome %in% c("1-exclusion", "2-exclusion"))
}

communities_network <- communities_network %>%
    rowwise() %>%
    mutate(ExclusionNetwork = list(subset_exclusion_network(Network))) %>%
    mutate(TriadCensus = sum(triad_census(Network))) %>%
    mutate(TriadCensusFull = sum(triad_census(Network) %>% `[`(c(9, 10, 12:16)))) %>%
    mutate(TriadCensusFullExclusion = sum(triad_census(ExclusionNetwork) %>% `[`(c(9, 10, 12:16)))) %>%
    mutate(NTCensusFullExclusion = triad_census(ExclusionNetwork) %>% `[`(c(10))) %>%
    select(-CommunityPairSize, -Network, -ExclusionNetwork)


sum(communities_network$TriadCensus) # 284 possible triads with 0, 1, 2, or 3 links in exclusion-coexistence network
sum(communities_network$TriadCensusFull) # 200 fully connected nodes (3 links) in exclusion-coexistence network
sum(communities_network$TriadCensusFullExclusion) # 77 fully connected nodes in exclusion-only network
n_triads <- sum(communities_network$TriadCensusFullExclusion)
n_NT <- sum(communities_network$NTCensusFullExclusion) # Number of nontransitive triad

set.seed(2)
N = 10000
tb <- tibble(
    BootstrapID = 1:N,
    NumberExpectedRPS = rbinom(N, n_triads, prob = 1/4)
)

range(tb$NumberExpectedRPS)

p2 <- tb %>%
    ggplot() +
    geom_histogram(aes(x = NumberExpectedRPS), color = "black", fill = "white", binwidth = 1) +
    annotate("point", x = n_NT, y = 0, shape = 21, size = 2, stroke = 1, color = "maroon") +
    #scale_x_continuous(limits = c(0, max(tb$NumberExpectedRPS)+2)) +
    annotate("text", x = 0, y = 1020, label = paste0("Binomial(n= ", n_triads, ", p = 1/4)"), size = 2.5, hjust = 0) +
    theme_classic() +
    theme(
        plot.background = element_blank()
    ) +
    guides() +
    labs(x = "number of expected nontransitivity", y = "count")

#
p <- ggdraw(plot_grid(p1, labels = "A")) +
    draw_plot(plot_grid(p2, labels = "B") , x = 0, y = 0, width = 0.3, height = 0.5, hjust = 0, vjust = 0)


ggsave(here::here("plots/Fig3.png"), p, width = 10, height = 5)
#ggsave(here::here("plots/Fig3_2.svg"), p, width = 10, height = 5)



