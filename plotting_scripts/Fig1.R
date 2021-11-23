# Figure 1 for pairwise network

library(tidyverse)
library(tidymodels)
library(cowplot)
library(ggraph)
library(tidygraph)
library(ggsci)
source(here::here("plotting_scripts/network_functions.R"))

sequences_abundance <- read_csv(here::here("data/temp/sequences_abundance.csv"))
communities <- read_csv(here::here("data/output/communities.csv"))
pairs <- read_csv(here::here("data/output/pairs.csv")) %>% mutate(InteractionType = ifelse(InteractionType == "neutrality", "coexistence", InteractionType))
pairs_example_outcomes <- read_csv(here::here("data/output/pairs_example_outcomes.csv"))
pairs_freq <- read_csv(here::here("data/output/pairs_freq.csv"))
interaction_type <- c("exclusion", "coexistence", "neutrality", "mutual exclusion", "frequency-dependent\ncoexistence")
interaction_color = c("#DB7469", "#557BAA", "#8650C4", "red", "blue")
names(interaction_color) <- interaction_type
networks_motif <- read_csv(here::here("data/output/networks_motif.csv"))


# Figure 1A
p_A <- ggdraw() + draw_image(here::here("plots/cartoons/Fig1A.png"), )

# Figure 1B: isolate abundance
# color set from Sylvie
color_sets <- tibble(Color = c("yellow", "deepskyblue3", "blue", "darkorchid2", "firebrick", "orange2", "grey"),
                     Family = c("Aeromonadaceae", "Enterobacteriaceae", "Moraxellaceae", "Pseudomonadaceae","Comamonadaceae","Alcaligenaceae", "Sphingobacteriaceae"))
plot_abundance <- function(df, label_x = "Community", label_y = "RelativeAbundance", fill = "CommunityESVID") {
    ggplot(df) +
        geom_bar(aes_string(x = label_x, y = label_y, fill = fill),
                 position = "stack", stat = "identity", col = 1) +
        theme_bw()
}

p_B <- sequences_abundance %>%
    filter(AlignmentType == "local") %>%
    filter(AllowMismatch == Inf) %>%
    filter(BasePairMismatch <= 4) %>%
    mutate(Community = ordered(Community,  communities$Community)) %>%
    plot_abundance(label_x = "Community", label_y = " RelativeAbundance", fill = "Family") +
    scale_fill_manual(values = setNames(color_sets$Color, color_sets$Family)) +
    scale_x_discrete(expand=c(0,0)) +
    scale_y_continuous(expand=c(0,0), limits = c(0,1), breaks = seq(0,1, .25)) +
    theme(axis.text.x = element_text(angle = 90), panel.grid.major = element_blank(), axis.ticks.x = element_blank()) +
    labs(x = "", y = "Relative abundance")

# Figure 1C
p_C <- ggdraw() + draw_image(here::here("plots/cartoons/Fig1C.png"))


# Figure 1D
## Plot pairs example dynamics
p_pairs_example_outcomes <- pairs_example_outcomes %>%
    filter(InteractionType != "neutrality") %>%
    left_join(pairs_freq, by = c("Community", "Isolate1", "Isolate2")) %>%
    mutate(Isolate1InitialODFreq = factor(Isolate1InitialODFreq),
           InteractionType = factor(InteractionType, interaction_type)) %>%
    ggplot(aes(x = Time, y = Isolate1MeasuredFreq, color = Isolate1InitialODFreq, group = Isolate1InitialODFreq)) +
    geom_point(size = 2) +
    geom_line(size = 1) +
    scale_y_continuous(breaks = seq(0,1,0.5)) +
    facet_grid(.~InteractionType) +
    theme_bw() +
    theme(axis.title.x = element_blank(), panel.spacing = unit(0, "mm"), strip.text.x = element_blank(),
          panel.border = element_rect(color = 1, fill = NA)) +
    guides(color = "none") +
    labs(x = "transfer", y = "frequency")
## The frequencies of coexistence vs. exclusion
temp <- pairs %>% filter(Assembly == "self_assembly") %>%
    mutate(InteractionType = ifelse(InteractionType == "neutrality", "coexistence", InteractionType)) %>%
    mutate(InteractionType = factor(InteractionType, c("coexistence", "exclusion"))) %>%
    group_by(InteractionType) %>% summarize(Count = n()) %>% ungroup() %>% mutate(Fraction = Count / sum(Count))
p_pairs_interaction <- temp %>%
    ggplot() +
    geom_col(aes(x = InteractionType, y = Count, fill = InteractionType), color = 1) +
    geom_text(x = -Inf, y = Inf, label = paste0("n = ", sum(temp$Count)), vjust = 1, hjust = -0.1) +
    geom_text(aes(x = InteractionType, y = Count, label = paste0(round(Fraction, 3) * 100,"%")), nudge_y = 5) +
    scale_fill_manual(values = interaction_color, breaks = c("exclusion", "coexistence")) +
    scale_y_continuous(limits = c(0, 150), expand = c(0,0)) +
    theme_classic() +
    theme(axis.title.x = element_blank(), legend.position = "top",
          axis.text.x = element_text(size = 10, color = "black"), axis.text.y = element_text(color = "black"),
          axis.title.y = element_text(size = 10)) +
    labs(x = "", y = "Number of pairs", fill = "")
ggsave(here::here("plots/Fig1-self_assembly.png"), p_pairs_interaction, width = 3, height = 4)
p_D <- plot_grid(p_pairs_interaction, p_pairs_example_outcomes, ncol = 1, axis = "lf", align = "h", rel_heights = c(1, 0.5))

# Similar coexistence count for both random assembly and self-assembly
temp <- pairs %>%
    mutate(Assembly = factor(Assembly, c("random_assembly", "across_community", "self_assembly"))) %>%
    filter(Assembly %in% c("random_assembly", "self_assembly")) %>%
    group_by(Assembly, InteractionType) %>% summarize(Count = n()) %>% mutate(Fraction = Count / sum(Count))
n_size <- temp %>% summarize(Count = sum(Count))
p_pairs_interaction_random <- temp %>%
    ggplot() +
    geom_col(aes(x = Assembly, y = Count, fill = InteractionType), position = "fill", color = 1, width = 0.8) +
    #geom_text(aes(x = Assembly, y = Fraction, label = paste0(round(Fraction, 3) * 100,"%")), nudge_y = 5) +
    scale_fill_manual(values = interaction_color, breaks = c("exclusion", "coexistence")) +
    scale_x_discrete(labels = c("random_assembly" = paste0("random assembly\nn=", pull(filter(n_size, Assembly == "random_assembly"), Count)),
                                "self_assembly" = paste0("self assembly\nn=", pull(filter(n_size, Assembly == "self_assembly"), Count)))) +
    scale_y_continuous(expand = c(0,0), breaks = c(0, .5, 1)) +
    theme_classic() +
    theme(axis.title.x = element_blank(), legend.position = "right", axis.text.x = element_text(size = 10)) +
    labs(x = "", y  = "Fraction", fill = "")
ggsave(here::here("plots/Fig1-random_assembly.png"), p_pairs_interaction_random, width = 4.5, height = 4.5)
## Stat
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



# Nontransitivity
p1 <- networks_motif %>%
    mutate(Community = factor(Community)) %>% mutate(Motif = as.character(Motif)) %>%
    group_by(Community) %>%
    filter(Motif %in% 2) %>%
    ggplot(aes(x = Motif, y = Fraction)) +
    geom_boxplot() +
    geom_jitter(shape = 21, size = 2, width = 0.2, height = 0, color = "red") +
    geom_text(x = -Inf, y = Inf, label = paste0("n = ", length(unique(networks_motif$Community))), vjust = 1, hjust = -0.1) +
    scale_x_discrete(labels = c("Motif1" = "Nontransitive", "Motif2" = "Transitive")) +
    scale_y_continuous(breaks = c(0, 0.5, 1), limits = c(-0.05,1), expand = c(0,0)) +
    #scale_color_npg(labels = c("communityPairs" = "Community", "randomNetworks" = "Species pool")) +
    theme_classic() +
    theme(legend.position = "right", legend.title = element_blank(), strip.text = element_blank(),
          panel.spacing = unit(0, "pt"),
          axis.text.x = element_blank(),
          axis.title.y = element_text(size = 11), axis.text.y = element_text(size = 10),
          legend.text = element_text(size = 10)) +
    labs(x = "", y = "Hierarchy")
p1
ggsave("../plots/Fig1-motif.png", plot = p1, width = 2, height = 3)



  #
p_top <- plot_grid(p_A, p_B, nrow = 1, rel_widths = c(1, 2), scale = .9)
p_middle <- plot_grid(p_C, p_D, nrow = 1, rel_widths = c(2, 1), axis = "tb", align = "vh", scale = .9)
p <- plot_grid(p_top, p_middle, ncol = 1, rel_heights = c(1, 2))
ggsave(here::here("plots/Fig1.pdf"), p, width = 12, height = 10)









#ggsave(here::here("plots/Fig1.png"), p, width = 10, height = 3)


if (FALSE) {
    # Figure 1D
    load(here::here("data/output/network_community.Rdata")) # Load observed networks net_list
    load(here::here("data/output/example_motif_list.Rdata")) # Load example motif graphs example_motifs
    network_pairs <- net_list %>%
        lapply(function(x){activate(x, edges) %>% as_tibble}) %>%
        bind_rows(.id = "Community") %>%
        filter(!(InteractionType == "coexistence" & from > to)) %>%
        mutate(InteractionType = ifelse(InteractionType == "neutrality", "coexistence", InteractionType)) %>%
        group_by(InteractionType) %>%
        summarize(Count = n())

    network_pairs %>%
        ggplot(aes(x = InteractionType, y = Count, fill = InteractionType)) +
        geom_col(color = 1) +
        geom_text(aes(y = Count + 5, label = paste0(round(Count/sum(Count), 2) * 100, "%")), position = position_dodge(width = 0.9)) +
        annotate("text", x = 0.5, y = 150, hjust = "inward", vjust = "inward", label = paste0("n=", sum(network_pairs$Count))) +
        scale_y_continuous(expand = c(0,0), limits = c(0, 150)) +
        scale_fill_manual(values = c(exclusion="#DB7469", coexistence="#557BAA")) +
        theme_cowplot() +
        guides(fill = F) +
        labs(x = "", y = "Number of pairwise competition")



















    isolates <- fread("../data/output/isolates.csv")
    pairs <- fread("../data/output/pairs.csv")
    pairs_melted <- fread("../data/output/pairs_melted.csv")
    isolates_random <- fread("../data/output/isolates_random.csv")
    pairs_random <- fread("../data/output/pairs_random.csv")
    pairs_random_melted <- fread("../data/output/pairs_random_melted.csv")
    communities <- fread("../data/output/communities.csv")
    community_names_ordered_by_size <- communities %>% arrange(CommunitySize) %>% pull(Community)
    community_abundance <- fread("../data/temp/communities_abundance.csv")
    #simulated_motif_counts <- fread("../data/temp/simulated_motif_counts.txt")
    #observed_motif_counts <- fread("../data/temp/observed_motif_counts.txt")
    #random_motif_counts <- fread("../data/temp/random_motif_counts.txt")
    #random_motif_counts_percentile <- fread("../data/temp/random_motif_counts_percentile.txt")
    network_motif <- fread(here::here("data/output/networks_motif.csv"))
    networks_motif_randomized <- fread(here::here("data/output/networks_motif_randomized.csv"))
    load("../data/output/network_community.Rdata") # Load observed networks net_list
    load("../data/output/example_motif_list.Rdata") # Load example motif graphs example_motifs
    interaction_color <- assign_interaction_color()

    # Panel A cartoon for experiment
    # p_A <- ggdraw() + draw_image("../data/experimental_scheme/Fig1A.png")


    # Panel B: pairwise coexistence vs pairwise competition
    network_pairs <- net_list %>%
        lapply(function(x){activate(x, edges) %>% as_tibble}) %>%
        bind_rows(.id = "Community") %>%
        filter(!(InteractionType == "coexistence" & from > to)) %>%
        mutate(InteractionType = ifelse(InteractionType == "neutrality", "coexistence", InteractionType)) %>%
        group_by(InteractionType) %>%
        summarize(Count = n())

    p_B <- network_pairs %>%
        ggplot(aes(x = InteractionType, y = Count, fill = InteractionType)) +
        geom_col(color = 1) +
        geom_text(aes(y = Count + 5, label = paste0(round(Count/sum(Count), 2) * 100, "%")), position = position_dodge(width = 0.9)) +
        annotate("text", x = 0.5, y = 150, hjust = "inward", vjust = "inward", label = paste0("n=", sum(network_pairs$Count))) +
        scale_y_continuous(expand = c(0,0), limits = c(0, 150)) +
        scale_fill_manual(values = c(exclusion="#DB7469", coexistence="#557BAA")) +
        theme_cowplot() +
        guides(fill = F) +
        labs(x = "", y = "Number of pairwise competition")

    ggsave("../plots/Fig1B.png", plot = p_B, width = 3, height = 4)

    # Panel C: motif counts pooled
    p_motifs <- plot_grid(plotlist = lapply(example_motif_list, function(x)plot_competitive_network(x, node_size=3)), nrow = 1)

    randomized_motif_counts_aggregated <- networks_motif_randomized %>%
        mutate(Motif = as.numeric(sub("Motif", "", Motif))) %>%
        group_by(Motif, Randomization) %>% summarize(Count = sum(CountMotif)) %>%
        summarize(p5 = quantile(Count, probs = 0.05), p95 = quantile(Count, probs = 0.95))
    observed_motif_counts_aggregated <- networks_motif %>%
        mutate(Motif = as.numeric(sub("Motif", "", Motif))) %>%
        group_by(Motif) %>%
        summarize(ObservedTotalCount = sum(CountMotif)) %>%
        left_join(randomized_motif_counts_aggregated) %>%
        mutate(Enrichment = ifelse(ObservedTotalCount < p5 , "underrepresented",
                                   ifelse(ObservedTotalCount > p95, "overrepresented", "Nonsignificant")))
    p_motif_counts <- observed_motif_counts_aggregated %>%
        ggplot() +
        geom_rect(aes(xmin = Motif -.5, xmax = Motif + .5, ymin = -Inf, ymax = Inf, fill = Enrichment), alpha = 0.5) +
        geom_segment(data = randomized_motif_counts_aggregated, aes(x = Motif, xend = Motif, y = p5, yend = p95), size = 1) +
        geom_point(data = randomized_motif_counts_aggregated, aes(x = Motif, y = p5, color = "random [5th nd 95th percentile]"), size = 2) +
        geom_point(data = randomized_motif_counts_aggregated, aes(x = Motif, y = p95, color = "random [5th nd 95th percentile]"), size = 2) +
        geom_point(aes(x = Motif, y = ObservedTotalCount, color = "observation"), size = 2) +
        scale_x_continuous(breaks = 1:7, expand = c(0,0)) +
        scale_fill_manual(values = c(overrepresented="#A6FFA1", underrepresented="#DB7469")) +
        scale_color_manual(values = c(`observation`="red", "random [5th nd 95th percentile]"="black")) +
        facet_grid(.~Motif, scale = "free_x") +
        theme_cowplot() +
        theme(legend.position = "bottom", legend.title = element_blank(),
              strip.background = element_blank(), strip.text = element_blank(), panel.spacing = unit(0, "cm")) +
        guides(fill = F) +
        labs(x = "", y = "Motif Count")
    p_C <- plot_grid(p_motifs, p_motif_counts, ncol = 1, axis = "tblr", align = "v", rel_heights = c(2,6))

    ggsave("../plots/Fig1C.png", p_C, width = 8, height = 4)

    # Panel D: number of intransitive triad loops
    p_net_list <- net_list %>% lapply(function(x) plot_competitive_network(x, node_size = 2))
    p_motif1 <- plot_competitive_network(example_motif_list[[1]], node_size = 2)
    for (i in 13:1) p_net_list[[i+1]] <- p_net_list[[i]]
    p_net_list[[1]] <- p_motif1
    p_community_graph <- plot_grid(plotlist = p_net_list, nrow = 1, greedy = F, scale = 1.2)

    p_intransitivity <- observed_motif_counts %>%
        filter(Motif == 1) %>%
        select(Community, Motif, RelativeMotifCount) %>%
        bind_rows(tibble(Community = "Motif1", Motif = 1, RelativeMotifCount = 1)) %>%
        mutate(Community = factor(Community, levels = c("Motif1", communities$Community))) %>%
        ggplot() +
        geom_point(aes(x = Community, y = RelativeMotifCount), size = 2, shape = 21) +
        scale_y_continuous(limits = c(0,1)) +
        facet_grid(.~Community, scales = "free_x") +
        theme_cowplot() +
        theme(axis.text.x = element_blank(), panel.spacing.x = unit(0, "cm"), strip.background = element_blank(), strip.text = element_blank()) +
        labs(x = "", y = "Intransitivity")

    p_D <- plot_grid(p_community_graph, p_intransitivity, ncol = 1, axis = "tblr", align = "v", rel_heights = c(2,5))
    ggsave("../plots/Fig1D.png", p_D, width = 12, height = 3)

    # p <- plot_grid(p_B, p_C, ncol = 2, axis = "tblr", align = "hv")
    # ggsave("../plots/Fig1.png", p, width = 10, height = 10)



    # Panel XX: pairwise coexistence
    df_pairs <- bind_rows(select(pairs, Community, Isolate1, Isolate2, InteractionType, From, To, Family1, Family2, Fermenter1, Fermenter2),
                          select(pairs_random, Community, Isolate1, Isolate2, InteractionType, From, To, Family1, Family2, Fermenter1, Fermenter2)) %>%
        mutate(Treatment = ifelse(grepl("C\\d+R\\d+", Community), "within community",
                                  ifelse(grepl("AcrAss", Community,), "across community", "natural isolates"))) %>%
        mutate(Treatment = factor(Treatment, c("natural isolates", "across community", "within community"))) %>%
        mutate(PairFermenter = ifelse(Fermenter1==T & Fermenter2==T, "FF", ifelse(Fermenter1==T & Fermenter2==F | Fermenter1==F & Fermenter2==T, "FR", ifelse(Fermenter1==F & Fermenter2==F, "RR", NA))))

    p <- df_pairs %>%
        group_by(Treatment, InteractionType) %>%
        summarize(Count = n()) %>%
        group_by(Treatment) %>%
        mutate(Fraction = Count / sum(Count), SampleSize = sum(Count)) %>%
        ggplot() +
        geom_col(aes(x = Treatment, y = Fraction, fill = InteractionType), color = 1) +
        geom_text(aes(x = Treatment, label = paste0("n = ", SampleSize)), y = 1, vjust = 1.5) +
        scale_y_continuous(limits = c(0,1), breaks = c(0, 0.5, 1), expand = c(0,0)) +
        scale_fill_manual(values = interaction_color) +
        theme_cowplot() +
        theme(legend.position = "top", legend.title = element_blank(), axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
        labs(x = "", y = "Fraction")
    ggsave("../plots/Fig1S1.png", plot = p, width = 5, height = 5)

    # Panel SXX:
    p <- df_pairs %>%
        filter(!is.na(PairFermenter)) %>%
        group_by(Treatment, PairFermenter, InteractionType) %>%
        summarize(Count = n()) %>%
        group_by(Treatment, PairFermenter) %>%
        mutate(Fraction = Count / sum(Count), SampleSize = sum(Count)) %>%
        ggplot() +
        geom_col(aes(x = Treatment, y = Fraction, fill = InteractionType), color = 1) +
        geom_text(aes(x = Treatment, label = paste0("n = ", SampleSize)), y = 1, vjust = 1.5) +
        scale_y_continuous(limits = c(0,1), breaks = c(0, 0.5, 1), expand = c(0,0)) +
        scale_fill_manual(values = interaction_color) +
        facet_grid(PairFermenter~.) +
        theme_cowplot() +
        theme(legend.position = "top", legend.title = element_blank(), axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
        labs(x = "", y = "Fraction")
    ggsave("../plots/Fig1S2.png", plot = p, width = 5, height = 10)


    # Panel XX: OD CASEU
    p <- pairs_random_melted %>%
        filter(Time == "T0", RawDataType == "Sanger") %>%
        #filter(Community == "AcrAss1") %>%
        #filter(Community == "RanAss1") %>%
        ggplot(aes(x = Isolate1InitialODFreq, y = Isolate1MeasuredFreq)) +
        geom_jitter(width = .02, shape = 21, aes(color = Community), size = 2) +
        geom_abline(slope = 1, intercept = 0) +
        scale_color_discrete(labels = c("AcrAss1"="across community 1", "RanAss1"="random community 1", "RanAss2"="random community 2")) +
        theme_cowplot() +
        theme(legend.position = "right") +
        panel_border(color = 1) +
        labs(x = "Frequency by OD", y = "Frequency by CASEU")
    p
    ggsave("../plots/Fig1S3.png", plot = p, width = 7, height = 5)


}


if (FALSE) {
    ## Motif demo
    colors_grey <- grey(seq(1,0, length.out = length(example_motif_list)))
    names(colors_grey) = 1:7
    p_motif_list <- rep(list(NA), 7)
    for (i in 1:length(example_motif_list)) {
        p_motif_list[[i]] <- example_motif_list[[i]] %>%
            plot_competitive_network(node_size = 5) +
            theme(panel.background = element_rect(fill = colors_grey[i], color = NA),
                  plot.margin = unit(c(0, 0, 0, 0), "cm"))
    }
    p_motifs <- plot_grid(plotlist = p_motif_list, nrow = 1)

    ## Motif count
    motif_type <- c("others", "0-scored")
    motif_color = c("#DB7469", "#557BAA")
    names(motif_color) <- motif_type

    p1 <- simulated_motif_counts %>%
        mutate(Motif = factor(Motif)) %>%
        group_by(CommunitySize, ProbPairCoexistence, Motif) %>%
        summarize(MeanCount = mean(Count)) %>%
        group_by(CommunitySize, ProbPairCoexistence) %>%
        mutate(SumMeanCount = sum(MeanCount), RelativeMeanCount = MeanCount/SumMeanCount) %>%
        filter(CommunitySize == 12) %>%
        ggplot() +
        geom_area(aes(x = ProbPairCoexistence, y = RelativeMeanCount, fill = Motif), color = 1) +
        scale_x_continuous(expand = c(0,0), breaks = c(0, 0.5, 1)) +
        scale_y_continuous(expand = c(0,0), breaks = c(0, 0.5, 1)) +
        scale_fill_manual(values = colors_grey) +
        theme_cowplot() +
        panel_border(color = 1) +
        theme(legend.title = element_blank(), legend.position = "none", legend.direction = "horizontal") +
        labs(x = "Probability of pairwise coexistence", y = "Relative motif count")

    p_B <- plot_grid(p_motifs, p1, ncol = 1, rel_heights = c(1,5))
}


if (FALSE){
    simulated_motif_counts_mean <- simulated_motif_counts %>%
        group_by(CommunitySize, ProbPairCoexistence, Motif) %>%
        summarize(MeanCount = mean(Count))

    simulated_motif_counts_mean %>%
        ggplot() +
        geom_area(aes(x = ProbPairCoexistence, y = MeanCount, fill = Motif), color = 1) +
        scale_x_continuous(expand = c(0,0)) +
        scale_y_continuous(expand = c(0,0)) +
        facet_grid(CommunitySize ~ ., scales = "free_y") +
        theme_cowplot()


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


