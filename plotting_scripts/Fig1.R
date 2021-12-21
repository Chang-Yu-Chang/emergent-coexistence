# Figure 1: species coexistence is emergent

library(tidyverse)
library(tidymodels)
library(cowplot)
library(ggraph)
library(tidygraph)
library(ggsci)
library(ggprism)
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


# Figure 1A: cartoon of experiment
p_A <- ggdraw() + draw_image(here::here("plots/cartoons/Fig1A.png")) + theme(plot.background = element_rect(fill = "white"))
ggsave(here::here("plots/Fig1A-assembly_cartoon.png"), p_A, width = 3, height = 2)

# Figure 1B: isolate abundance
## color set from Sylvie
color_sets <- tibble(Color = c("yellow", "deepskyblue3", "blue", "darkorchid2", "firebrick", "orange2", "grey"),
                     Family = c("Aeromonadaceae", "Enterobacteriaceae", "Moraxellaceae", "Pseudomonadaceae","Comamonadaceae","Alcaligenaceae", "Sphingobacteriaceae"))
p_B <- sequences_abundance %>%
    filter(AlignmentType == "local") %>%
    filter(AllowMismatch == Inf) %>%
    filter(BasePairMismatch <= 4) %>%
    mutate(Community = ordered(Community,  communities$Community)) %>%
    ggplot() +
    geom_bar(aes(x = Community, y = RelativeAbundance, fill = Family), position = "stack", stat = "identity", col = 1) +
    theme_bw() +
    scale_fill_manual(values = setNames(color_sets$Color, color_sets$Family)) +
    scale_y_continuous(breaks = c(0, .5, 1), expand = c(0,0), limits = c(0, 1)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5), panel.grid = element_blank(), axis.ticks.x = element_blank(),
          axis.title.x = element_blank(), panel.border = element_rect(color = 1, size = 1)) +
    labs(y = "Relative abundance")
ggsave(here::here("plots/Fig1B-relative_abundance.png"), p_B, width = 5, height = 3)

# Figure 1C: pairwise comeptiton cartoon
p_C <- ggdraw() + draw_image(here::here("plots/cartoons/Fig1C.png")) + theme(plot.background = element_rect(fill = "white"))
ggsave(here::here("plots/Fig1C-pairwise_competition_cartoon.png"), p_C, width = 4, height = 3)

# Figure 1D
## Plot pairs example dynamics
p_pairs_example_outcomes <- pairs_example_outcomes %>%
    filter(InteractionType != "neutrality") %>%
    left_join(pairs_freq, by = c("Community", "Isolate1", "Isolate2")) %>%
    mutate(Isolate1InitialODFreq = factor(Isolate1InitialODFreq)) %>%
    mutate(InteractionType = factor(InteractionType, c("coexistence", "exclusion"))) %>%
    ggplot(aes(x = Time, y = Isolate1MeasuredFreq, color = Isolate1InitialODFreq, group = Isolate1InitialODFreq)) +
    geom_point(size = 2) +
    geom_line(size = 1) +
    scale_y_continuous(breaks = c(0, .5, 1), limits = c(0,1)) +
    facet_grid(.~InteractionType) +
    theme_bw() +
    theme(axis.title.x = element_blank(), panel.spacing = unit(2, "mm"), strip.text.x = element_blank(),
          panel.border = element_rect(color = 1, fill = NA, size = 1)) +
    guides(color = "none") +
    labs(y = "Frequency")
## The frequencies of coexistence vs. exclusion
temp <- pairs %>% filter(Assembly == "self_assembly") %>%
    mutate(InteractionType = ifelse(InteractionType == "neutrality", "coexistence", InteractionType)) %>%
    mutate(InteractionType = factor(InteractionType, c("coexistence", "exclusion"))) %>%
    group_by(InteractionType) %>% summarize(Count = n()) %>% ungroup() %>% mutate(Fraction = Count / sum(Count))
p_pairs_interaction <- temp %>%
    ggplot() +
    geom_col(aes(x = InteractionType, y = Count, fill = InteractionType), color = 1) +
    geom_text(x = -Inf, y = Inf, label = paste0("n = ", sum(temp$Count)), vjust = 1, hjust = -0.1) +
    geom_text(aes(x = InteractionType, y = Count, label = paste0(round(Fraction, 3) * 100,"%")), nudge_y = 6) +
    scale_fill_manual(values = assign_interaction_color(level = "simple")) +
    scale_y_continuous(limits = c(0, 150), expand = c(0,0)) +
    theme_classic() +
    theme(axis.title.x = element_blank(), axis.title.y = element_text(size = 10),
          axis.text.x = element_text(size = 10, color = "black"), axis.text.y = element_text(color = "black"),
          legend.position = "none") +
    labs(x = "", y = "Number of pairs", fill = "")
p_D <- plot_grid(p_pairs_interaction, p_pairs_example_outcomes, ncol = 1, axis = "lf", align = "h", rel_heights = c(2,1))
ggsave(here::here("plots/Fig1D-pairwise_competition.png"), p_D, width = 3, height = 4)

# Figur1 1E; Similar coexistence count for both random assembly and self-assembly
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


# Figure 1F: hierarchy score for communities, compared to 1000 randomized ones
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

p <- communities_hierarchy %>%
    ggplot(aes(x = Metric, y = HierarchyScore)) +
    geom_boxplot() +
    geom_jitter(shape = 1, size = 2) +
    scale_x_discrete (labels = c("1" = "Higgins et al 2017", "2" = "Violation\nof ranks", "3" = "Fraction of\ntransitive motifs")) +
    facet_wrap(.~Metric, scales = "free", nrow = 1) +
    theme_classic() +
    theme(axis.title.x = element_blank(), strip.text = element_blank(),
          panel.border = element_rect(fill = NA, color = 1, size = 1.5)) +
    labs(x = "Metrics", y = "Hierarchy score")

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

p_lower <- plot_grid(p1, p2, p3, nrow = 1, labels = LETTERS[2:4])
p_F <- plot_grid(p, p_lower, ncol = 1, labels = c("A", ""), rel_heights = c(1, 4))

ggsave(here::here("plots/Fig1F-hierarchy.png"), p_F, width = 9, height = 13)



# Figure 1G: network and motifs
networks_motif <- read_csv(here::here("data/output/networks_motif.csv")) %>% filter(str_detect(Community, "C\\d"))
networks_motif_randomized <- read_csv(here::here("data/output/networks_motif_randomized.csv")) %>% filter(str_detect(Community, "C\\d")) %>% mutate(Community = ordered(Community, communities$Community))
p_net_matrix_list <- lapply(net_list, function(x) plot_adjacent_matrix(x))
p_net_list <- lapply(net_list, function(x) plot_competitive_network(x, node_size = 2) + theme(plot.background = element_rect(fill = NA)))
p_list <- rep(list(NA), length(p_net_list))
for (i in 1:length(net_list)) p_list[[i]] <- ggdraw(p_net_matrix_list[[i]]) + draw_plot(plot = p_net_list[[i]], x = -.1, y = -.1, width = 0.7, height = 0.7)
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
        facet_wrap(Community ~., scale = "free_y", nrow = 1)  +
        theme_classic() +
        theme(panel.background = element_rect(color = 1, size = 1), legend.position = "none")
}
p_motif_count_list <- rep(list(NA), length(p_net_list))
for (i in 1:length(p_net_list)) {
    if (i == 1) p_motif_count_list[[i]] <- plot_motif_count(i) + theme(axis.title.x = element_blank())
    if (i %in% 2:7) p_motif_count_list[[i]] <- plot_motif_count(i) + theme(axis.title = element_blank())
    if (i == 8) p_motif_count_list[[i]] <- plot_motif_count(i)
    if (i %in% 9:13) p_motif_count_list[[i]] <- plot_motif_count(i) + theme(axis.title.y = element_blank())
}
shared_legend <- {plot_motif_count(1) +
        theme(legend.position = "right", legend.title = element_blank(), legend.text = element_text(size = 10))} %>%
    cowplot::get_legend()
p_G <- list(p_list[1:7], p_motif_count_list[1], p_motif_count_list[2:7],
            p_list[8:13], list(NULL), p_motif_count_list[8:13], list(shared_legend)) %>%
    unlist(recursive = F) %>%
    plot_grid(plotlist = ., ncol = 7, align = "v", axis = "lrtb") + theme(plot.background = element_rect(fill = "white"))
ggsave(here::here("plots/Fig1G-networks.png"), p_G, width = 16, height = 7)


# Figure 1H: total motif count and example motifs
networks_motif <- read_csv(here::here("data/output/networks_motif.csv")) %>% filter(str_detect(Community, "C\\d"))
networks_motif_randomized <- read_csv(here::here("data/output/networks_motif_randomized.csv")) %>% filter(str_detect(Community, "C\\d")) %>% mutate(Community = ordered(Community, communities$Community))
load(here::here("data/output/motif_list.Rdata"))
p_motif_list <- lapply(motif_list, function(x) plot_competitive_network(x, node_size = 3))
p_motif_example <- plot_grid(plotlist = p_motif_list, nrow = 1)
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
p_motif_count <- ggplot() +
    # 5% and 95% percentiles in randomized networks
    geom_point(data = networks_motif_randomized_total_percentile, aes(x = Motif, y = Count, group = Motif, color = "randomized network")) +
    geom_segment(data = networks_motif_randomized_total_percentile %>% pivot_wider(id_cols = Motif, names_from = Percentile, values_from = Count),
                 aes(x = Motif, xend = Motif, y = p5, yend = p95, color = "randomized network")) +
    # Observations
    geom_point(data = networks_motif_total, aes(x = Motif, y = Count, color = "observed network")) +
    scale_x_continuous(breaks = 1:7, position = "top") +
    scale_color_manual(values = c("observed network" = "red", "randomized network" = "black"))+
    facet_grid(.~Motif, scales = "free_x") +
    theme_classic() +
    theme(panel.background = element_rect(color = 1, size = 1),
          panel.spacing = unit(0, "mm"), strip.text = element_blank(),
          axis.title.x = element_blank(),
          legend.position = "bottom", legend.title = element_blank(),
          legend.text = element_text(size = 12), axis.text = element_text(size = 12))
p_H <- plot_grid(p_motif_example, p_motif_count, nrow = 2, rel_heights = c(1,3), axis = "lr", align = "v") + theme(plot.background = element_rect(fill = "white"))
ggsave(here::here("plots/Fig1H-motif_counts.png"), p_H, width = 8, height = 4)

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

# Figure J: finer grain pairwise coexistence
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
    theme(axis.title.x = element_blank(), panel.spacing = unit(2, "mm"), strip.text.x = element_blank(),
          panel.border = element_rect(color = 1, fill = NA, size = 1)) +
    guides(color = "none") +
    labs(y = "Frequency")
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
p_J <- plot_grid(p_pairs_interaction_finer, p_pairs_example_outcomes_finer, ncol = 1, axis = "lf", align = "h", rel_heights = c(2,1))
ggsave(here::here("plots/Fig1J-pairwise_competition_finer.png"), p_J, width = 7, height = 4)


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

# Figure 1L: diagonal analysis
communities <- read_csv(here::here("data/output/communities.csv")) %>% filter(str_detect(Community, "C\\d"))
networks_diag <- read_csv(here::here("data/output/networks_diag.csv"))
networks_diag_randomized <- read_csv(here::here("data/output/networks_diag_randomized.csv"))
networks_diag_sum <- networks_diag %>%
    group_by(DistanceToDiagonal) %>%
    summarize(CountCoexistenceSum = sum(CountCoexistence)) %>%
    bind_rows(tibble(DistanceToDiagonal = 7:11, CountCoexistenceSum = 0))
networks_diag_randomized_sum <- networks_diag_randomized %>%
    group_by(Replicate, DistanceToDiagonal) %>%
    summarize(CountCoexistenceSum = sum(CountCoexistence))

p_L <- networks_diag_sum %>%
    ggplot() +
    # Random networks
    geom_boxplot(data = networks_diag_randomized_sum, aes(x = DistanceToDiagonal, y = CountCoexistenceSum, group = DistanceToDiagonal, color = "randomized network"),
                 outlier.size = 1) +
    geom_jitter(data = networks_diag_randomized_sum, aes(x = DistanceToDiagonal, y = CountCoexistenceSum, group = DistanceToDiagonal, color = "randomized network"),
                size = .1, alpha = 0.5, width = .3) +
    # Observed networks
    geom_point(aes(x = DistanceToDiagonal, y = CountCoexistenceSum, group = DistanceToDiagonal, color = "observed network"), size = 2) +
    geom_line(aes(x = DistanceToDiagonal, y = CountCoexistenceSum, color = "observed network")) +
    scale_x_continuous(breaks = 1:11) +
    scale_color_manual(values = c("observed network" = "red", "randomized network" = "black"))+
    theme_classic() +
    theme(legend.position = "top", legend.title = element_blank(),
          panel.border = element_rect(fill = NA, color = 1, size = 1.5)) +
    #labs(x = "Distance to diagonal (|i-j|)", y = "Count of pairwise coexistence")
    labs(x = "Difference in rankings", y = "Count of pairwise coexistence")

ggsave(here::here("plots/Fig1L-diagonal_analysis.png"), p_L, width = 4, height = 4)


# Figure 1M






