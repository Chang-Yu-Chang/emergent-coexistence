# Figure 1: species coexistence is emergent

library(tidyverse)
library(tidymodels)
library(cowplot)
library(ggraph)
library(tidygraph)
#library(ggridges)
#library(ggsci)
library(flextable)
library(officer)

#library(ggprism)
source(here::here("plotting_scripts/network_functions.R"))

sequences_abundance <- read_csv(here::here("data/temp/sequences_abundance.csv"))
communities <- read_csv(here::here("data/output/communities.csv"))
pairs <- read_csv(here::here("data/output/pairs.csv")) %>% mutate(InteractionType = ifelse(InteractionType == "neutrality", "coexistence", InteractionType))
pairs_example_outcomes <- read_csv(here::here("data/output/pairs_example_outcomes.csv"))
pairs_example_outcomes_finer <- read_csv(here::here("data/output/pairs_example_outcomes_finer.csv"))
pairs_freq <- read_csv(here::here("data/output/pairs_freq.csv"))
networks_motif <- read_csv(here::here("data/output/networks_motif.csv"))
networks_motif_randomized <- read_csv(here::here("data/output/networks_motif_randomized.csv"))
networks_motif_randomized_percentile <- read_csv(here::here("data/output/networks_motif_randomized_percentile.csv")) %>% mutate(Community = factor(Community, communities$Community)) %>% arrange(Community)
load(here::here("data/output/network_community.Rdata"))
load(here::here("data/output/motif_list.Rdata"))
interaction_color <- assign_interaction_color()

# Figure 1A: two hypotheses; two network structures
load(here::here("data/output/motif_list.Rdata"))
node_size = 5; text_size = 3
p1 <- motif_list[[7]] %>%
    ggraph(layout = "circle") +
    geom_node_point(shape = 21, size = node_size, fill = "gray", stroke = node_size/5) +
    geom_edge_link(aes(color = InteractionType), width = node_size/5,
                   arrow = arrow(length = unit(node_size/2, "mm"), type = "closed", angle = 30, ends = "last"),
                   start_cap = circle(node_size/2+1, "mm"),
                   end_cap = circle(node_size/2+1, "mm")) +
    geom_text(x = 0.25, y = 0.433, label = "coexistence", angle = -30, vjust = -1.5, fontface = "bold", size = text_size, color = interaction_color["coexistence"]) +
    geom_text(x = -1.1, y = 0, label = "All-coexistence", angle = 90, size = text_size*1.2, color = 1) +
    scale_edge_color_manual(values = interaction_color) +
    scale_color_manual(values = c("black" = "black", "white" = NA)) +
    scale_shape_manual(values = c("real" = 21, "fake" = NA)) +
    scale_x_continuous(limits = c(-1.2, 1.2)) +
    scale_y_continuous(limits = c(-1.2, 1.2)) +
    theme_minimal() +
    theme(legend.position = "none", panel.grid = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
          plot.background = element_rect(fill = NA, color = NA), plot.margin = margin(0,0,0,0, unit = "pt")) +
    labs(x = "", y = "")

p2 <- motif_list[[1]] %>%
    ggraph(layout = "circle") +
    geom_node_point(shape = 21, size = node_size, fill = "gray", stroke = node_size/5) +
    geom_edge_link(aes(color = InteractionType), width = node_size/5,
                   arrow = arrow(length = unit(node_size/2, "mm"), type = "closed", angle = 30, ends = "last"),
                   start_cap = circle(node_size/2+1, "mm"),
                   end_cap = circle(node_size/2+1, "mm")) +
    geom_text(x = 0.25, y = 0.433, label = "exclusion", angle = -30, vjust = -1.5, fontface = "bold", size = text_size, color = interaction_color["exclusion"]) +
    geom_text(x = -1.1, y = 0, label = "Nontransitivity", angle = 90, size = text_size*1.2, color = 1) +
    scale_edge_color_manual(values = interaction_color) +
    scale_color_manual(values = c("black" = "black", "white" = NA)) +
    scale_shape_manual(values = c("real" = 21, "fake" = NA)) +
    scale_x_continuous(limits = c(-1.2, 1.2)) +
    scale_y_continuous(limits = c(-1.2, 1.2)) +
    theme_minimal() +
    theme(legend.position = "none", panel.grid = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
          plot.background = element_rect(fill = NA, color = NA), plot.margin = margin(0,0,0,0, unit = "pt")) +
    labs(x = "", y = "")

p_A <- plot_grid(p1, p2, ncol = 1) + theme(panel.background = element_rect(color = NA, fill = "white"))
ggsave(here::here("plots/Fig1A-two_hypothese.png"), p_A, width = 1.5, height = 3)


# Figure 1B: pairwise outcomes
## Plot pairs example dynamics
p_pairs_example_outcomes <- pairs_example_outcomes %>%
    filter(InteractionType != "neutrality") %>%
    left_join(pairs_freq, by = c("Community", "Isolate1", "Isolate2")) %>%
    mutate(Isolate1InitialODFreq = factor(Isolate1InitialODFreq)) %>%
    mutate(InteractionType = factor(InteractionType, c("exclusion", "coexistence"))) %>%
    ggplot(aes(x = Time, y = Isolate1MeasuredFreq, color = Isolate1InitialODFreq, group = Isolate1InitialODFreq)) +
    geom_point(size = 1) +
    geom_line(size = .5) +
    scale_y_continuous(breaks = c(0, .5, 1), limits = c(0,1)) +
    facet_grid(.~InteractionType) +
    theme_bw() +
    theme(panel.spacing = unit(2, "mm"), strip.text.x = element_blank(),
          panel.border = element_rect(color = 1, fill = NA, size = 1),
          axis.title = element_text(size = 10), axis.text = element_text(color = 1, size = 8)) +
    guides(color = F) +
    labs(x = "Time", y = "Frequency")
## The frequencies of coexistence vs. exclusion
temp <- pairs %>% filter(Assembly == "self_assembly") %>%
    mutate(InteractionType = ifelse(InteractionType == "neutrality", "coexistence", InteractionType)) %>%
    mutate(InteractionType = factor(InteractionType, c("exclusion", "coexistence"))) %>%
    group_by(InteractionType) %>% summarize(Count = n()) %>% ungroup() %>% mutate(Fraction = Count / sum(Count))
p_pairs_interaction <- temp %>%
    ggplot() +
    geom_col(aes(x = InteractionType, y = Count, fill = InteractionType), color = 1) +
    geom_text(x = Inf, y = Inf, label = paste0("n = ", sum(temp$Count)), vjust = 1, hjust = 1.5) +
    geom_text(aes(x = InteractionType, y = Count, label = paste0(round(Fraction, 3) * 100,"%")), nudge_y = 10) +
    scale_fill_manual(values = assign_interaction_color(level = "simple")) +
    scale_y_continuous(limits = c(0, 150), expand = c(0,0)) +
    theme_classic() +
    theme(axis.title.x = element_blank(), axis.title.y = element_text(size = 10),
          axis.text = element_text(size = 8, color = "black"),
          legend.position = "none") +
    labs(x = "", y = "Number of pairs", fill = "")
p_B <- plot_grid(p_pairs_interaction, p_pairs_example_outcomes, ncol = 1, axis = "lf", align = "h", rel_heights = c(2,1))
ggsave(here::here("plots/Fig1B-pairwise_competition.png"), p_B, width = 3, height = 4)


# Figure 1C: nontransitive motif count
networks_motif <- read_csv(here::here("data/output/networks_motif.csv")) %>% filter(str_detect(Community, "C\\d"))
networks_motif_randomized <- read_csv(here::here("data/output/networks_motif_randomized.csv")) %>% filter(str_detect(Community, "C\\d")) %>% mutate(Community = ordered(Community, communities$Community))

## Diagram of nontransitive motif
load(here::here("data/output/motif_list.Rdata"))
p_motif1 <- plot_competitive_network(motif_list[[1]], node_size = 3)

##
networks_motif_total <- networks_motif %>%
        group_by(Motif) %>%
        summarize(Count = sum(Count)) %>%
        filter(Motif == 1)

p_motif_count <- networks_motif_randomized %>%
    filter(Motif == 1) %>%
    group_by(Motif, Replicate) %>%
    summarize(Count = sum(Count)) %>%
    group_by(Motif) %>%
    # 5% and 95% percentiles in randomized networks
    mutate(p5 = quantile(Count, 0.05), p95 = quantile(Count, 0.95), ColoredTails = ifelse(Count <= p5, "tail", ifelse(Count >= p95, "head", "body"))) %>%
    ggplot() +
    geom_vline(xintercept = 0, color = 1) +
    geom_hline(yintercept = 0, color = 1) +
    geom_histogram(aes(y = Count, x = after_stat(count / max(count)), fill = ColoredTails), alpha = .3, color = 1) +
    geom_point(data = filter(networks_motif, Motif == 1), aes(x = 0, y = Count, color = "observed network"), pch = 1, size = 2, stroke = 2, inherit.aes = F) +
    scale_fill_manual(values = c("head" = "#FF0000A0", "body" = "#A0A0A0A0", "tail" = "#FF0000A0"),
                      labels = c("head" = "top 5%", "body" = "middle", "tail" = "bottom 5%")) +
    scale_color_manual(values = c("observed network" = "red")) +
    coord_flip() +
    facet_grid(.~Motif) +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 3)) +
    theme_cowplot() +
    theme(panel.background = element_rect(color = 1, size = 1.5, fill = NA),
          panel.spacing = unit(0, "mm"), strip.text = element_blank(),
          legend.position = "top",
          axis.title = element_text(size = 10), axis.text = element_text(size = 8),
          plot.background = element_rect(fill = "white", color = NA)) +
    #guides(fill = guide_legend(title = "randomized network"), color = guide_legend(title = "")) +
    guides(fill = F, color = F) +
    labs(x = "Probability density", y = "Count of nontransitive motif")

p_C <- p_motif_count
ggsave(here::here("plots/Fig1C-nontransitive_motif_counts.png"), p_C, width = 2, height = 3)


# Figure 1D: hierarchy metrics
communities_hierarchy <- read_csv(here::here("data/output/communities_hierarchy.csv")) %>% mutate(Community = factor(Community, communities$Community))
p_D <- communities_hierarchy %>%
    pivot_longer(cols = matches("\\d"), names_to = c("Variable", "Metric"), names_pattern = "(.*)(\\d)") %>%
    pivot_wider(names_from = Variable) %>%
    mutate(Significance = as.logical(Significance)) %>%
    #filter(Metric == 1) %>%
    ggplot(aes(x = Metric, y = HierarchyScore)) +
    geom_boxplot(width = .5, lwd = .8) +
    geom_jitter(shape = 1, size = 2, width = .2, stroke = .8) +
    scale_x_discrete (position = "bottom", labels = c("1" = "Higgins2017", "2" = "Rank-following\npairs", "3" = "Fraction of transitive motifs")) +
    scale_y_continuous(limits = c(0.5,1), breaks = c(0.5, 0.75, 1)) +
    #facet_wrap(.~Metric, scales = "free", nrow = 1) +
    theme_classic() +
    theme(panel.grid.major.x = element_line(color = "grey", linetype = 2),
          axis.title.x = element_blank(), strip.text = element_blank(),
          axis.text = element_text(size = 8, color = 1),
          axis.ticks.x = element_blank(),
          panel.border = element_rect(fill = NA, color = 1, size = 1.5)) +
    labs(x = "Metrics", y = "Hierarchy score")
p_D
ggsave(here::here("plots/Fig1D-hierarchy_metrics.png"), p_D, width = 5, height = 5)

#
p_aligned <- plot_grid(p_C, p_D, nrow = 1, scale = .9,  labels = LETTERS[3:4], axis = "lrtb", align = "hv")
p <- plot_grid(p_A, p_B, p_aligned, nrow = 1, scale = c(.7, .9, 1), rel_widths = c(1,1.5,3), labels = c(LETTERS[1:2],"")) + theme(plot.background = element_rect(fill = "white", color = NA))
ggsave(here::here("plots/Fig1.png"), p, width = 9, height = 3)





# Figure S1: cartoon of self-assembly experiment and isolate abundance
p1 <- ggdraw() + draw_image(here::here("plots/cartoons/Fig1A.png")) + theme(plot.background = element_rect(fill = "white", color = NA))
color_sets <- tibble(Color = c("yellow", "deepskyblue3", "blue", "darkorchid2", "firebrick", "orange2", "grey"),
                     Family = c("Aeromonadaceae", "Enterobacteriaceae", "Moraxellaceae", "Pseudomonadaceae","Comamonadaceae","Alcaligenaceae", "Sphingobacteriaceae"))
temp <- sequences_abundance %>%
    filter(AlignmentType == "local") %>%
    filter(AllowMismatch == Inf) %>%
    filter(BasePairMismatch <= 4) %>%
    mutate(Community = ordered(Community,  communities$Community))

p2 <- temp %>%
    ggplot() +
    geom_bar(aes(x = Community, y = RelativeAbundance, fill = Family), size = .5, color = "grey30", position = "stack", stat = "identity", col = 1) +
    theme_bw() +
    scale_fill_manual(values = setNames(color_sets$Color, color_sets$Family)) +
    scale_y_continuous(breaks = c(0, .5, 1), expand = c(0,0), limits = c(0, 1)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5), panel.grid = element_blank(), axis.ticks.x = element_blank(),
          axis.title.x = element_blank(), panel.border = element_rect(color = 1, size = 1)) +
    labs(y = "Relative abundance")

## Stats
temp %>%
    group_by(Community) %>%
    summarize(Total = sum(RelativeAbundance)) %>%
    summarize(Mean = mean(Total))

p_S1 <- plot_grid(p1, p2, nrow = 1, rel_widths = c(1,2), labels = LETTERS[1:2], scale = .9) + theme(plot.background = element_rect(fill = "white", color = NA))
ggsave(here::here("plots/FigS1-glucose_community_bar.png"), p_S1, width = 8, height = 3)

# Figure S2: pairwise competition cartoon
p_S2 <- ggdraw() + draw_image(here::here("plots/cartoons/Fig1C.png")) + theme(plot.background = element_rect(fill = "white", color = NA))
ggsave(here::here("plots/FigS2-pairwise_experiment_cartoon.png"), p_S2, width = 4, height = 3)

# Figure S3: finer grain pairwise coexistence
## Plot pairs example dynamics
p_pairs_example_outcomes_finer <- pairs_example_outcomes_finer %>%
    left_join(pairs_freq, by = c("Community", "Isolate1", "Isolate2")) %>%
    mutate(Isolate1InitialODFreq = factor(Isolate1InitialODFreq)) %>%
    mutate(InteractionTypeFiner = factor(InteractionTypeFiner, names(assign_interaction_color(level = "finer")))) %>%
    ggplot(aes(x = Time, y = Isolate1MeasuredFreq, color = Isolate1InitialODFreq, group = Isolate1InitialODFreq)) +
    geom_point(size = 2) +
    geom_line(size = 1) +
    scale_y_continuous(breaks = c(0, .5, 1), limits = c(0,1)) +
    facet_grid(.~InteractionTypeFiner) +
    theme_bw() +
    theme(panel.spacing = unit(2, "mm"), strip.text.x = element_blank(),
          panel.border = element_rect(color = 1, fill = NA, size = 1)) +
    guides(color = "none") +
    labs(x = "Time", y = "Frequency")
## The frequencies of coexistence vs. exclusion
temp <- pairs %>% filter(Assembly == "self_assembly") %>%
    mutate(InteractionTypeFiner = factor(InteractionTypeFiner, names(assign_interaction_color(level = "finer")))) %>%
    group_by(InteractionTypeFiner) %>% summarize(Count = n()) %>% ungroup() %>% mutate(Fraction = Count / sum(Count))
p_pairs_interaction_finer <- temp %>%
    ggplot() +
    geom_col(aes(x = InteractionTypeFiner, y = Count, fill = InteractionTypeFiner), color = 1) +
    geom_text(x = -Inf, y = Inf, label = paste0("n = ", sum(temp$Count)), vjust = 1, hjust = -0.1) +
    geom_text(aes(x = InteractionTypeFiner, y = Count, label = paste0(round(Fraction, 3) * 100,"%")), nudge_y = 6) +
    scale_fill_manual(values = assign_interaction_color(level = "finer")) +
    scale_x_discrete(breaks = c("competitive exclusion", "stable coexistence", "neutrality", "mutual exclusion", "frequency-dependent coexistence"),
                     labels = c("competitive\nexclusion", "stable\ncoexistence", "neutrality", "mutual\nexclusion", "frequency-dependent\ncoexistence")) +
    scale_y_continuous(limits = c(0, 150), expand = c(0,0)) +
    theme_classic() +
    theme(axis.title.x = element_blank(), axis.title.y = element_text(size = 10),
          axis.text.x = element_text(size = 10, color = "black"), axis.text.y = element_text(color = "black"),
          legend.position = "none") +
    labs(x = "", y = "Number of pairs", fill = "")
p_S3 <- plot_grid(p_pairs_interaction_finer, p_pairs_example_outcomes_finer, ncol = 1, axis = "lf", align = "h", rel_heights = c(2,1))
ggsave(here::here("plots/FigS3-pairwise_outcomes_finer.png"), p_S3, width = 7, height = 4)




# Figure S4: total motif count and example motifs
if (FALSE) {
## 13 networks
p_net_list <- lapply(net_list, function(x) plot_competitive_network(x, node_size = 3) + theme(plot.background = element_rect(fill = NA)))
temp <- get_legend(p_net_list[[1]] + theme_bw() + theme(legend.title = element_blank(), legend.text = element_text(size = 10)))
p1 <- plot_grid(plotlist = c(p_net_list[1:13], list(temp)), nrow = 3, labels = c(communities$Community[1:13], ""), label_size = 10)

}
## motif diagram
load(here::here("data/output/motif_list.Rdata"))
p_motif_list <- lapply(motif_list, function(x) plot_competitive_network(x, node_size = 3))
p2 <- plot_grid(plotlist = p_motif_list, nrow = 1)

## motif counts
networks_motif <- read_csv(here::here("data/output/networks_motif.csv")) %>% filter(str_detect(Community, "C\\d"))
networks_motif_randomized <- read_csv(here::here("data/output/networks_motif_randomized.csv")) %>% filter(str_detect(Community, "C\\d")) %>% mutate(Community = ordered(Community, communities$Community))
networks_motif_randomized_total_percentile <- networks_motif_randomized %>%
    group_by(Replicate, Motif) %>%
    summarize(Count = sum(Count)) %>%
    group_by(Motif) %>%
    arrange(Motif, Count) %>%
    slice(c(1000 * 0.05, 1000 * 0.95)) %>%
    mutate(Percentile = c("p5", "p95")) %>%
    select(Motif, Count, Percentile)
networks_motif_total <- networks_motif %>%
    group_by(Motif) %>%
    summarize(Count = sum(Count))
p3 <- networks_motif_randomized %>%
    group_by(Motif, Replicate) %>%
    summarize(Count = sum(Count)) %>%
    group_by(Motif) %>%
    # 5% and 95% percentiles in randomized networks
    mutate(p5 = quantile(Count, 0.05), p95 = quantile(Count, 0.95), ColoredTails = ifelse(Count <= p5, "tail", ifelse(Count >= p95, "head", "body"))) %>%
    ggplot() +
    geom_vline(xintercept = 0, color = 1) +
    geom_histogram(aes(y = Count, fill = ColoredTails), binwidth = 1, color = NA) +
    geom_point(data = networks_motif_total, x = 0, aes(y = Count, color = "observed network"), pch = 1, size = 2, stroke = 2, inherit.aes = F) +
    scale_fill_manual(values = c("head" = "#FF0000A0", "body" = "#A0A0A0A0", "tail" = "#0000FFA0"),
                      labels = c("head" = "top 5%", "body" = "middle", "tail" = "bottom 5%")) +
    scale_color_manual(values = c("observed network" = "red")) +
    facet_grid(.~Motif) +
    #scale_y_discrete(expand = c(0, 0)) +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 2)) +
    theme_cowplot() +
    theme(panel.background = element_rect(color = 1, size = 1), panel.spacing = unit(0, "mm"),
          strip.background = element_rect(color = 1, fill = NA, size = 1)) +
    guides(fill = guide_legend(title = "randomized network"), color = guide_legend(title = "")) +
    labs(x = "Probability density", y = "Motif count")

#p_S4 <- plot_grid(p1, p2, p3, ncol = 1, rel_heights = c(4,1,2), axis = "lr", align = "v", labels = LETTERS[1:3]) + theme(plot.background = element_rect(fill = "white", color = NA))
p_S4 <- plot_grid(p2, p3, ncol = 1, rel_heights = c(1,3), axis = "lr", align = "v", labels = LETTERS[1:2]) + theme(plot.background = element_rect(fill = "white", color = NA))
ggsave(here::here("plots/FigS4-total_motif_counts.png"), p_S4, width = 10, height = 3)




# Figure S5: network matrices and motifs
## Matrix and graph
networks_motif <- read_csv(here::here("data/output/networks_motif.csv")) %>% filter(str_detect(Community, "C\\d"))
networks_motif_randomized <- read_csv(here::here("data/output/networks_motif_randomized.csv")) %>% filter(str_detect(Community, "C\\d")) %>% mutate(Community = ordered(Community, communities$Community))
p_net_matrix_list <- lapply(net_list, function(x) plot_adjacent_matrix(x) + theme(plot.margin = grid::unit(c(5,0,3,0), "mm")))
p_net_list <- lapply(net_list, function(x) plot_competitive_network(x, node_size = 2) + theme(plot.background = element_rect(fill = NA)))
p_list <- rep(list(NA), length(p_net_list))
for (i in 1:length(net_list)) p_list[[i]] <- ggdraw(p_net_matrix_list[[i]]) + draw_plot(plot = p_net_list[[i]], x = -.1, y = -.1, width = 0.7, height = 0.7)
## Motif count
plot_motif_count <- function (x = 1) {
    motif_randomized_subset <- networks_motif_randomized_percentile %>%
        filter(Community %in% communities$Community[x]) %>%
        mutate(Community = factor(Community, communities$Community[x]))
    motif_community_subset <- networks_motif %>%
        filter(Community %in% communities$Community[x]) %>%
        mutate(Community = factor(Community, communities$Community[x]))

    ggplot() +
        # 5% and 95% percentiles in randomized networks
        geom_point(data = motif_randomized_subset, aes(x = Motif, y = Count, group = Motif, color = "randomized network")) +
        geom_segment(data = motif_randomized_subset %>% pivot_wider(id_cols = c(Community, Motif), names_from = Percentile, values_from = Count),
                     aes(x = Motif, xend = Motif, y = p5, yend = p95, color = "randomized network")) +
        # Observations
        geom_point(data = motif_community_subset, aes(x = Motif, y = Count, color = "observed network")) +
        scale_x_continuous(breaks = 1:7) +
        scale_color_manual(values = c("observed network" = "red", "randomized network" = "black"))+
        #facet_wrap(Community ~., scale = "free_y", nrow = 1)  +
        theme_classic() +
        theme(panel.background = element_rect(color = 1, size = 1), legend.position = "none")
}
p_motif_count_list <- rep(list(NA), length(p_net_list))
for (i in 1:length(p_net_list)) {
    if (i %in% c(1, 6)) p_motif_count_list[[i]] <- plot_motif_count(i) + theme(axis.title.x = element_blank())
    if (i %in% c(2:5, 7:10)) p_motif_count_list[[i]] <- plot_motif_count(i) + theme(axis.title = element_blank())
    if (i == 11) p_motif_count_list[[i]] <- plot_motif_count(i)
    if (i %in% 12:13) p_motif_count_list[[i]] <- plot_motif_count(i) + theme(axis.title.y = element_blank())
}

## Example matrix
plot_example_matrix <- function(graph) {
    graph_ranked <- graph %>%
        activate(nodes) %>%
        select(Isolate, PlotRank) %>%
        activate(edges) %>%
        mutate(
            fromRank = .N()$PlotRank[match(from, .N()$Isolate)],
            toRank = .N()$PlotRank[match(to, .N()$Isolate)])
    n_nodes <- igraph::vcount(graph_ranked)
    interaction_color <- assign_interaction_color(level = "matrix")

    graph_ranked %>%
        filter(fromRank <= toRank) %>%
        bind_edges(tibble(from = 1:n_nodes, to = 1:n_nodes, fromRank = 1:n_nodes, toRank = 1:n_nodes, InteractionType = "self")) %>%
        as_tibble() %>%
        ggplot() +
        geom_tile(aes(x = factor(toRank), y = ordered(fromRank, n_nodes:1), fill = InteractionType), width = 0.9, height = 0.9) +
        scale_x_discrete(position = "top", labels = c("top\nrank", rep("", n_nodes-2), "bottom\nrank")) +
        scale_y_discrete(position = "right", labels = c("bottom\nrank", rep("", n_nodes-2), "top\nrank")) +
        scale_fill_manual(breaks = c("exclusion", "coexistence", "self"), values = c(assign_interaction_color(), "self" = "black")) +
        theme_bw() +
        theme(axis.ticks = element_blank(), axis.title = element_blank(), legend.title = element_blank(),
              axis.text = element_text(size = 10, color = 1),
              panel.border = element_blank(), panel.grid = element_blank())
}
p_example_matrix <- plot_example_matrix(net_list[[1]]) + guides(fill = F) + theme(plot.background = element_rect(fill = "grey90"), panel.background = element_rect(fill = "grey90"))
## Get legend for matrix
shared_legend_matrix <- get_legend(plot_example_matrix(net_list[[1]]) + theme(legend.text = element_text(size = 15), legend.justification = "right", plot.background = element_rect(fill = "grey90"), panel.background = element_rect(fill = "grey90")))
## Get legend for line
p_temp <- plot_motif_count(1) + theme(legend.position = "right", legend.title = element_blank(), legend.text = element_text(size = 15))
shared_legend_line <- cowplot::get_legend(p_temp)
p1 <- plot_grid(plotlist = c(list(p_example_matrix), list(shared_legend_matrix), list(NULL), list(shared_legend_line), rep(list(NULL), 1)), nrow = 1) + theme(plot.background = element_rect(fill = "white", color = NA))
p2 <- list(p_list[1:5], p_motif_count_list[1:5],
             p_list[6:10], p_motif_count_list[6:10],
             p_list[11:13], rep(list(NULL), 2),
             p_motif_count_list[11:13]) %>%
    unlist(recursive = F) %>%
    plot_grid(plotlist = ., labels = c(communities$Community[1:5], rep("", 5), communities$Community[6:10], rep("", 5), communities$Community[11:13], rep("", 7)),
              ncol = 5, axis = "tbrl", align = "v")
p_S5 <- plot_grid(p1, p2, ncol = 1, rel_heights = c(1,6)) + theme(plot.background = element_rect(fill = "white", color = NA))
ggsave(here::here("plots/FigS5-networks_matrix.png"), p_S5, width = 10, height = 12)


# Figure S6: hierarchy score for communities, compared to 1000 randomized ones
## 1. Higgins2017
## 2. Violation of ranks
## 3. Fraction of transitive motifs
communities_hierarchy <- read_csv(here::here("data/output/communities_hierarchy.csv")) %>% mutate(Community = factor(Community, communities$Community))
communities_hierarchy_randomized <- read_csv(here::here("data/output/communities_hierarchy_randomized.csv")) %>% mutate(Community = factor(Community, communities$Community))
networks_motif <- read_csv(here::here("data/output/networks_motif.csv")) %>% filter(Motif == 2) %>% filter(str_detect(Community, "C\\d"))
networks_motif_randomized <- read_csv(here::here("data/output/networks_motif_randomized.csv")) %>% filter(Motif == 2) %>% filter(str_detect(Community, "C\\d")) %>% mutate(Community = ordered(Community, communities$Community))

b = 1000
networks_motif_pvalue <- networks_motif_randomized %>%
    group_by(Community) %>%
    # Add a bottom row to each community
    bind_rows(communities %>% select(Community) %>% mutate(Replicate = b+1, Motif = 2, Count = -1, Fraction = -.001) %>% filter(str_detect(Community, "C\\d"))) %>%
    arrange(Community, desc(Fraction)) %>%
    left_join(networks_motif %>% select(Community, Motif, FractionObserved = Fraction)) %>%
    mutate(Percentile = 0:b / b) %>%
    filter(FractionObserved > Fraction) %>%
    slice(1) %>%
    mutate(Significance = Percentile < 0.05) %>%
    select(Community, HierarchyScore3 = FractionObserved, Percentile3 = Percentile, Significance3 = Significance) %>%
    mutate(Community = factor(Community, communities$Community)) %>%
    arrange(Community)

communities_hierarchy <- communities_hierarchy %>% left_join(networks_motif_pvalue) %>%
    pivot_longer(cols = matches("\\d"), names_to = c("Variable", "Metric"), names_pattern = "(.*)(\\d)") %>%
    pivot_wider(names_from = Variable) %>%
    mutate(Significance = as.logical(Significance))

 plot_top <- function(x) {
    communities_hierarchy %>%
        filter(Metric == x) %>%
        ggplot(aes(x = Metric, y = HierarchyScore)) +
        geom_boxplot(width = .3) +
        geom_jitter(shape = 1, size = 2, width = .1) +
        scale_x_discrete (position = "top", labels = c("1" = "Higgins et al 2017", "2" = "Rank-following pairs", "3" = "Fraction of transitive motifs")) +
        scale_y_continuous(limits = c(0,1), breaks = c(0, 0.5, 1)) +
        facet_wrap(.~Metric, scales = "free", nrow = 1) +
        theme_classic() +
        theme(axis.title.x = element_blank(), strip.text = element_blank(),
              axis.text.x = element_text(size = 13), axis.ticks.x = element_blank(),
              panel.border = element_rect(fill = NA, color = 1, size = 1.5)) +
        labs(x = "Metrics", y = "Score")
}
plot_bottom <- function(x) {
    communities_hierarchy %>%
        filter(Metric == x) %>%
        ggplot() +
        geom_rect(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, aes(fill = Significance), alpha = .5) +
        geom_histogram(data = communities_hierarchy_randomized, aes(y = HierarchyScore1), color = 1, fill = "white") +
        geom_hline(aes(yintercept = HierarchyScore), color = "red") +
        geom_text(aes(x = Inf, y = -Inf, label = paste0("p=", Percentile)), vjust = -1, hjust = 1.5) +
        scale_x_continuous(breaks = scales::pretty_breaks(n = 3)) +
        scale_y_continuous(breaks = c(0.5, 0.75, 1)) +
        scale_fill_manual(values = c("TRUE" = "#DB7469", "FALSE" = "white")) +
        facet_wrap(Community ~., ncol = 2, scales = "free_x", dir = "v") +
        theme_classic() +
        theme(panel.grid.major = element_blank(), legend.position = "none", panel.border = element_rect(color = 1, fill = NA, size = 1.5)) +
        labs(x = "Count", y = "Hierarchy score")
}
p_top <- plot_grid(plot_top(1), plot_top(2), plot_top(3), nrow = 1, labels = LETTERS[1:3], align = "hv")
p_bottom <- plot_grid(plot_bottom(1), plot_bottom(2), plot_bottom(3), nrow = 1)

if (FALSE) {
    p1 <- communities_hierarchy %>%
        filter(Metric == 1) %>%
        ggplot() +
        geom_rect(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, aes(fill = Significance), alpha = .5) +
        geom_histogram(data = communities_hierarchy_randomized, aes(y = HierarchyScore1), color = 1, fill = "white") +
        geom_hline(aes(yintercept = HierarchyScore), color = "red") +
        geom_text(aes(x = Inf, y = -Inf, label = paste0("p=", Percentile)), vjust = -1, hjust = 1.5) +
        scale_x_continuous(breaks = scales::pretty_breaks(n = 3)) +
        scale_y_continuous(breaks = c(0.5, 0.75, 1)) +
        scale_fill_manual(values = c("TRUE" = "#DB7469", "FALSE" = "white")) +
        facet_wrap(Community ~., ncol = 2, scales = "free_x", dir = "v") +
        theme_classic() +
        theme(panel.grid.major = element_blank(), legend.position = "none", panel.border = element_rect(color = 1, fill = NA, size = 1.5)) +
        labs(x = "Count", y = "Hierarchy score") +
        ggtitle("Higgins et al 2017")

    p2 <- communities_hierarchy %>%
        filter(Metric == 2) %>%
        ggplot() +
        geom_rect(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, aes(fill = Significance), alpha = .5) +
        geom_histogram(data = communities_hierarchy_randomized, aes(y = HierarchyScore2), color = 1, fill = "white") +
        geom_hline(aes(yintercept = HierarchyScore), color = "red") +
        geom_text(aes(x = Inf, y = -Inf, label = paste0("p=", Percentile)), vjust = -1, hjust = 1.5) +
        scale_x_continuous(breaks = scales::pretty_breaks(n = 3)) +
        scale_y_continuous(breaks = c(0, 0.5, 1)) +
        scale_fill_manual(values = c("TRUE" = "#DB7469", "FALSE" = "white")) +
        facet_wrap(Community ~., ncol = 2, scales = "free_x", dir = "v") +
        theme_classic() +
        theme(panel.grid.major = element_blank(), legend.position = "none", panel.border = element_rect(color = 1, fill = NA, size = 1.5)) +
        labs(x = "Count", y = "Hierarchy score") +
        ggtitle("Violation of ranks")

    p3 <- communities_hierarchy %>%
        filter(Metric == 3) %>%
        ggplot() +
        geom_rect(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, aes(fill = Significance), alpha = .5) +
        geom_histogram(data = networks_motif_randomized, aes(y = Fraction), color = 1, fill = "white") +
        geom_hline(aes(yintercept = HierarchyScore), color = "red") +
        geom_text(aes(x = Inf, y = Inf, label = paste0("p=", Percentile)), vjust = 1.5, hjust = 1.5) +
        scale_x_continuous(breaks = scales::pretty_breaks(n = 3)) +
        scale_y_continuous(breaks = c(0,0.5,1)) +
        scale_fill_manual(values = c("TRUE" = "#DB7469", "FALSE" = "white")) +
        facet_wrap(Community ~., ncol = 2, scales = "free_x", dir = "v") +
        theme_classic() +
        theme(panel.grid.major = element_blank(), legend.position = "none", panel.border = element_rect(color = 1, fill = NA, size = 1.5)) +
        labs(x = "Count", y = "Hierarchy score") +
        ggtitle("Fraction of transitive motifs")
}
#p_lower <- plot_grid(p1, p2, p3, nrow = 1, labels = LETTERS[2:4])
p_S6 <- plot_grid(p_top, p_bottom, ncol = 1, rel_heights = c(1, 4))

ggsave(here::here("plots/FigS6-hierarchy.png"), p_S6, width = 9, height = 13)



# Figure S7: diagonal analysis
communities <- read_csv(here::here("data/output/communities.csv")) %>% filter(str_detect(Community, "C\\d"))
networks_diag <- read_csv(here::here("data/output/networks_diag.csv"))
networks_diag_randomized <- read_csv(here::here("data/output/networks_diag_randomized.csv"))
networks_diag_sum <- networks_diag %>%
    group_by(DistanceToDiagonal) %>%
    summarize(ObservedCountCoexistenceSum = sum(CountCoexistence)) %>%
    bind_rows(tibble(DistanceToDiagonal = 7:11, ObservedCountCoexistenceSum = 0))
networks_diag_randomized_sum <- networks_diag_randomized %>%
    group_by(Replicate, DistanceToDiagonal) %>%
    summarize(CountCoexistenceSum = sum(CountCoexistence))
## Statistics
stat_diag <- networks_diag_randomized_sum %>%
    group_by(DistanceToDiagonal) %>%
    # find percentile
    bind_rows(tibble(Replicate = 0, DistanceToDiagonal = 1:11, CountCoexistenceSum = 0)) %>%
    left_join(networks_diag_sum) %>%
    arrange(DistanceToDiagonal, desc(CountCoexistenceSum)) %>%
    mutate(Percentile = (1:n())/n()) %>%
    filter(CountCoexistenceSum <= ObservedCountCoexistenceSum) %>%
    group_by(DistanceToDiagonal) %>%
    slice(1) %>%
    select(DistanceToDiagonal, Percentile) %>%
    # Asterisk
    mutate(Significance = ifelse(Percentile < 0.001 | Percentile > 0.999, "***",
                                 ifelse(Percentile < 0.01 | Percentile > 0.99, "**",
                                        ifelse(Percentile < 0.05 | Percentile > 0.95, "**",
                                               ifelse(Percentile > 0.05 & Percentile < 0.95, "n.s.", NA)))))

p_S7 <- networks_diag_sum %>%
    ggplot() +
    # Random networks
    geom_boxplot(data = networks_diag_randomized_sum, aes(x = DistanceToDiagonal, y = CountCoexistenceSum, group = DistanceToDiagonal, color = "randomized network"),
                 outlier.size = 1) +
    geom_jitter(data = networks_diag_randomized_sum, aes(x = DistanceToDiagonal, y = CountCoexistenceSum, group = DistanceToDiagonal, color = "randomized network"),
                size = .1, alpha = 0.5, width = .3) +
    # Observed networks
    geom_point(aes(x = DistanceToDiagonal, y = ObservedCountCoexistenceSum, group = DistanceToDiagonal, color = "observed network"), size = 2) +
    geom_line(aes(x = DistanceToDiagonal, y = ObservedCountCoexistenceSum, color = "observed network")) +
    # Asterisk
    geom_text(data = stat_diag, aes(x = DistanceToDiagonal, y = Inf, label = Significance), vjust = 2) +
    scale_x_continuous(breaks = 1:11) +
    scale_y_continuous(limits = c(0, 40)) +
    scale_color_manual(values = c("observed network" = "red", "randomized network" = "black"))+
    theme_classic() +
    theme(legend.position = "top", legend.title = element_blank(),
          panel.border = element_rect(fill = NA, color = 1, size = 1.5)) +
    #labs(x = "Distance to diagonal (|i-j|)", y = "Count of pairwise coexistence")
    labs(x = "Difference in rankings", y = "Count of pairwise coexistence")

ggsave(here::here("plots/FigS7-diagonal_analysis.png"), p_S7, width = 4, height = 4)


# Figure SXX. Raw pair frequencies

## Number of pairs using CASEU
pairs_freq %>%
    filter(str_detect(Community, "C\\d+R\\d+")) %>%
    filter(Time == "T8") %>%
    #group_by(Community, Isolate1, Isolate2) %>%
    unite("Pair", Community, Isolate1, Isolate2, sep = "_") %>%
    select(Pair, Isolate1InitialODFreq, RawDataType) %>%
    pivot_wider(names_from = Isolate1InitialODFreq, values_from = RawDataType, names_prefix = "f") %>%
    mutate(Method = ifelse(f5 == "Sanger" | f50 == "Sanger" | f95 == "Sanger", "sanger",
                           ifelse(f5 == "CFU" & f50 == "CFU" & f95 == "CFU", "cfu", NA))) %>%
    group_by(Method) %>%
    summarize(Count = n())

# Figure SXX. raw data for pair



# Table S1. Interaction tables
ft1 <- read_csv(here::here("data/output/pairs_interaction_table.csv")) %>%
    mutate(InteractionType = ifelse(InteractionType == "neutrality", "coexistence", InteractionType)) %>%
    arrange(InteractionType, InteractionTypeFiner, FromRare, FromMedium, FromAbundant) %>%
    setNames(c("From 5%", "From 50%", "From 95%", "Outcome", "Finer outcome", "Count")) %>%
    flextable() %>%
    width(j = 1:3, width = 1) %>%
    width(j = 5, width = 2.5)

save_as_image(ft1, here::here("plots/TableS1.png"))

read_csv(here::here("data/output/pairs_interaction_table.csv")) %>%
    pull(Count) %>% sum


if (FALSE) {
    # Figure1 1E: coexistence pair count for both random assembly and self-assembly
    temp <- pairs %>%
        mutate(Assembly = factor(Assembly, c("random_assembly", "across_community", "self_assembly"))) %>%
        filter(Assembly %in% c("random_assembly", "self_assembly")) %>%
        group_by(Assembly, InteractionType) %>% summarize(Count = n()) %>% mutate(Fraction = Count / sum(Count))
    n_size <- temp %>% summarize(Count = sum(Count))
    p_E <- temp %>%
        ggplot() +
        geom_col(aes(x = Assembly, y = Count, fill = InteractionType), position = "fill", color = 1, width = 0.8) +
        scale_fill_manual(values = assign_interaction_color(level = "simple")) +
        scale_x_discrete(labels = c("random_assembly" = paste0("random assembly\nn=", pull(filter(n_size, Assembly == "random_assembly"), Count)),
                                    "self_assembly" = paste0("self assembly\nn=", pull(filter(n_size, Assembly == "self_assembly"), Count)))) +
        scale_y_continuous(expand = c(0,0), breaks = c(0, .5, 1)) +
        theme_classic() +
        theme(axis.title.x = element_blank(), legend.position = "top", axis.text.x = element_text(size = 10)) +
        labs(x = "", y  = "Fraction", fill = "")
    ggsave(here::here("plots/Fig1E-pairwsie_competition_assembly.png"), p_E, width = 3, height = 3)

    ## Stat: does assembly explain for coexistence ratio?
    observed_stat <- pairs %>%
        filter(Assembly %in% c("random_assembly", "self_assembly")) %>%
        chisq_test(InteractionType ~ Assembly) %>% pull(statistic)
    null_stat <- pairs %>%
        filter(Assembly %in% c("random_assembly", "self_assembly")) %>%
        specify(InteractionType ~ Assembly, success = "coexistence") %>%
        hypothesize(null = "independence") %>%
        generate(reps = 1000, type = "permute") %>%
        calculate(stat = "Chisq", order = c("random_assembly", "self_assembly"))
    null_stat %>%
        get_p_value(obs_stat = observed_stat, direction = "two-sided")




    # Figure I: fraction of pairwise coexistence across communities
    pairs_count <- pairs %>%
        filter(Assembly == "self_assembly") %>%
        group_by(Community) %>%
        summarize(Count = n())

    p_I <- pairs %>%
        mutate(Community = factor(Community, communities$Community)) %>%
        filter(Assembly == "self_assembly") %>%
        ggplot() +
        geom_bar(aes(x = Community, fill = InteractionType), color = 1, position = position_fill(), size = .5) +
        geom_text(data = pairs_count, aes(x = Community, y = 1, label = Count), vjust = -.5) +
        geom_rect(xmin = -Inf, xmax = Inf, ymin = 0, ymax = 1, color = grey(0.1), fill = NA, size = .5) +
        scale_fill_manual(values = assign_interaction_color()) +
        scale_y_continuous(breaks = c(0,.5,1), limit = c(0, 1.15), expand = c(0,0)) +
        theme_classic() +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
              axis.line.y = element_blank(),
              legend.title = element_blank(), legend.position = "top") +
        labs(y = "Fraction")
    ggsave(here::here("plots/Fig1I-pairs_counts.png"), p_I, width = 4, height = 3)


    # Figure 1K: relative abundance within pairs
    pairs_freq_FN <- pairs_freq %>%
        filter(Isolate1InitialODFreq == 50) %>%
        filter(Time == "T8") %>%
        left_join(pairs) %>%
        filter(PairFermenter == "FN", InteractionType == "coexistence") %>%
        mutate(FractionFermenter = ifelse(Fermenter1, Isolate1MeasuredFreq, 1 - Isolate1MeasuredFreq),
               FractionRespirator = ifelse(Fermenter1, 1 - Isolate1MeasuredFreq, Isolate1MeasuredFreq)) %>%
        mutate(RF_Ratio = FractionRespirator / FractionFermenter) %>%
        select(Community, Isolate1, Isolate2, Isolate1InitialODFreq, Isolate1MeasuredFreq,
               PairFermenter, Fermenter1, Fermenter2,
               FractionFermenter, FractionRespirator, RF_Ratio)

    p1 <- pairs_freq_FN %>%
        ggplot() +
        geom_boxplot(aes(x = PairFermenter, y = RF_Ratio), color = 1) +
        geom_jitter(aes(x = PairFermenter, y = RF_Ratio), color = 1, shape = 1, size = 2, width = 0.2) +
        scale_y_log10(limits = c(0.01, 1), minor_breaks = rep(1:9, 4)*(10^rep(-2:1, each = 9)), guide = "prism_minor") +
        theme_classic() +
        theme(axis.title.x = element_blank(), axis.text.x = element_blank(),
              axis.ticks.y = element_line(size = .8), axis.ticks.length.y = unit(2, "mm")) +
        labs(y = "R/F")

    p2 <- pairs_freq_FN %>%
        arrange(FractionFermenter) %>%
        mutate(Pair = 1:n()) %>%
        select(Pair, FractionFermenter, FractionRespirator) %>%
        pivot_longer(cols = starts_with("Fraction"), names_prefix = "Fraction", names_to = "Fermenter", values_to = "Fraction") %>%
        ggplot() +
        geom_col(aes(x = Pair, y = Fraction, fill = Fermenter), color = 1) +
        scale_x_continuous(breaks = 1:22, expand = c(0,0)) +
        scale_y_continuous(breaks = c(0,.5, 1), expand = c(0,0)) +
        scale_fill_npg() +
        theme_classic() +
        theme(legend.position = "right", legend.title = element_blank())
    p_K <- plot_grid(p1, p2, nrow = 1, rel_widths = c(1,3), labels = LETTERS[1:2], axis = "tb", align = "h")
    ggsave(here::here("plots/Fig1K-pair_relative_abundance.png"), p_K, width = 7, height = 3)


}










