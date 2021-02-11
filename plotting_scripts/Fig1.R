# Experimental data

library(tidygraph)
library(ggraph)
library(tidyverse)
library(data.table)
library(cowplot)
source("network_functions.R")
whether_plot_supp <- T

# Read data
communities <- fread("../data/output/communities.csv")
community_names_ordered_by_size <- communities %>% arrange(CommunitySize) %>% pull(Community)
community_abundance <- fread("../data/temp/communities_abundance.csv")
#simulated_motif_counts <- fread("../data/temp/simulated_motif_counts.txt")
#observed_motif_counts <- fread("../data/temp/observed_motif_counts.txt")
random_motif_counts <- fread("../data/temp/random_motif_counts.txt")
random_motif_counts_percentile <- fread("../data/temp/random_motif_counts_percentile.txt")
load("../data/temp/graph_list.Rdata") # Load observed networks graph_list
load("../data/temp/example_motif_list.Rdata") # Load example motif graphs example_motifs

# Panel A cartoon for experiment
# p_A <- ggdraw() + draw_image("../data/experimental_scheme/Fig1A.png")


# Panel B: pairwise coexistence vs pairwise competition
network_pairs <- graph_list %>%
    lapply(function(x){activate(x, edges) %>% as_tibble}) %>%
    bind_rows(.id = "Community") %>%
    filter(!(InteractionType == "coexistence" & from > to)) %>%
    mutate(InteractionType = ifelse(InteractionType == "neutrality", "coexistence", InteractionType)) %>%
    group_by(InteractionType) %>%
    summarize(Count = n())

p_B <- network_pairs %>%
    ggplot(aes(x = InteractionType, y = Count, fill = InteractionType)) +
    geom_col() +
    geom_text(aes(y = Count + 5, label = paste0(round(Count/sum(Count), 2) * 100, "%")), position = position_dodge(width = 0.9)) +
    annotate("text", x = 0.5, y = 150, hjust = "inward", vjust = "inward", label = paste0("n=", sum(network_pairs$Count))) +
    scale_y_continuous(expand = c(0,0), limits = c(0, 150)) +
    scale_fill_manual(values = c(exclusion="#DB7469", coexistence="#557BAA")) +
    theme_cowplot() +
    guides(fill = F) +
    labs(x = "", y = "Number of pairwise competition")

ggsave("../plots/Fig1B.png", plot = p_B, width = 3, height = 4)

# Panel SS: pairwise coexistence in



# Panel C: motif counts pooled
p_motifs <- plot_grid(plotlist = lapply(example_motif_list, function(x)plot_competitive_network(x, node_size=3)), nrow = 1)
p_motifs

random_motif_counts_aggregated <- random_motif_counts %>%
    group_by(Motif, Seed) %>% summarize(Count = sum(Count)) %>%
    summarize(p5 = quantile(Count, probs = 0.05), p95 = quantile(Count, probs = 0.95))

observed_motif_counts_aggregated <- observed_motif_counts %>%
    group_by(Motif) %>%
    summarize(ObservedTotalCount = sum(Count)) %>%
    left_join(random_motif_counts_aggregated) %>%
    mutate(Enrichment = ifelse(ObservedTotalCount < p5 , "underrepresented",
                               ifelse(ObservedTotalCount > p95, "overrepresented", "Nonsignificant")))

p_motif_counts <- observed_motif_counts_aggregated %>%
    ggplot() +
    geom_rect(aes(xmin = Motif -.5, xmax = Motif + .5, ymin = -Inf, ymax = Inf, fill = Enrichment), alpha = 0.5) +
    geom_segment(aes(x = Motif, xend = Motif, y = p5, yend = p95), size = 1) +
    geom_point(aes(x = Motif, y = p5, color = "random [5th nd 95th percentile]"), size = 2) +
    geom_point(aes(x = Motif, y = p95, color = "random [5th nd 95th percentile]"), size = 2) +
    geom_point(data = observed_motif_counts_aggregated, aes(x = Motif, y = ObservedTotalCount, color = "observation"), size = 2) +
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
p_graph_list <- graph_list %>% lapply(function(x) plot_competitive_network(x, node_size = 2))
p_motif1 <- plot_competitive_network(example_motif_list[[1]], node_size = 2)
for (i in 13:1) p_graph_list[[i+1]] <- p_graph_list[[i]]
p_graph_list[[1]] <- p_motif1
p_community_graph <- plot_grid(plotlist = p_graph_list, nrow = 1, greedy = F, scale = 1.2)

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



if (whether_plot_supp) {
    ## 0-scored motif vs others
    p2 <- simulated_motif_counts %>%
        mutate(MotifType = ifelse(Motif %in% c(1, 7), "0-scored", "others")) %>%
        group_by(CommunitySize, Seed, ProbPairCoexistence, MotifType) %>%
        summarize(Count = sum(Count)) %>%
        group_by(CommunitySize, ProbPairCoexistence, MotifType) %>%
        summarize(MeanCount = mean(Count)) %>%
        filter(CommunitySize == 12) %>%
        ggplot() +
        geom_area(aes(x = ProbPairCoexistence, y = MeanCount, fill = MotifType), color = 1) +
        scale_x_continuous(expand = c(0,0), breaks = c(0, 0.5, 1)) +
        scale_y_continuous(expand = c(0,0), breaks = c(0, 100, 200)) +
        scale_fill_manual(values = motif_color) +
        theme_cowplot() +
        theme(legend.title = element_blank(),
              legend.position = "top") +
        labs(x = "Fraction of pairwise coexistence", y = "Mean motif count")

    ggsave("../plots/Fig1B_example.png", plot = p2, width = 4, height = 4)


    # Possible supp
    # Panel C
    ps1 <- simulated_motif_counts %>%
        mutate(Motif = factor(Motif)) %>%
        group_by(CommunitySize, ProbPairCoexistence, Motif) %>%
        summarize(MeanCount = mean(Count)) %>%
        group_by(CommunitySize, ProbPairCoexistence) %>%
        mutate(SumMeanCount = sum(MeanCount), RelativeMeanCount = MeanCount/SumMeanCount) %>%
        ggplot() +
        geom_area(aes(x = ProbPairCoexistence, y = RelativeMeanCount, fill = Motif), color = 1) +
        scale_x_continuous(expand = c(0,0), breaks = c(0, 0.5, 1)) +
        scale_y_continuous(expand = c(0,0), breaks = c(0, 0.5, 1)) +
        scale_fill_manual(values = colors_grey) +
        facet_grid(.~CommunitySize) +
        guides(fill = guide_legend(nrow = 1)) +
        theme_cowplot() +
        theme(legend.position = "top",
              legend.direction = "horizontal") +
        panel_border(color = 1) +
        labs(x = "Probability of pairwise coexistence", y = "Relative motif count")

    ggsave("../plots/FigS1.png", plot = ps1, width = 10, height = 4)

    # Panel XX: motif count as a function of community size
    p_motifs <- example_motif_list %>%
        lapply(function(x) plot_competitive_network(x, node_size = 3)) %>%
        plot_grid(plotlist = ., nrow = 1)

    summary_network_motifs <- graph_list %>%
        lapply(summarize_network_motif) %>%
        bind_rows(.id = "Community")

    ps2 <- summary_network_motifs %>%
        ggplot(aes(x = CommunitySize, y = RelativeMotifCount)) +
        geom_jitter(size = 3, shape = 21, width = 0.1) +
        geom_smooth(method = "lm", formula = y ~ x) +
        scale_x_continuous(limits = c(2, 13), breaks = c(2, 7, 12)) +
        scale_y_continuous(limits = c(-0.001,1.001), breaks = c(0, 0.5, 1)) +
        facet_wrap(.~Motif, nrow = 1) +
        theme_cowplot() +
        theme(strip.background = element_blank(), strip.text = element_blank()) +
        panel_border(color = "black") +
        labs(x = "Community size", y = "Relative motif count")

    p <- plot_grid(p_motifs, ps2, ncol = 1, axis = "rl", align = "hv", rel_heights = c(2,5))
    ggsave("../plots/FigS2.png", plot = p, width = 10, height = 4)

    # Panel: relative abundance of communities
    community_abundance_subset <- community_abundance %>%
        mutate(CommunityESVID = factor(CommunityESVID)) %>%
        filter(Community %in% community_names_ordered_by_size, Transfer == 12) %>%
        filter(!grepl("Glu-T12-X11-Rep1", SampleID), !grepl("Glu-T12-X1-Rep2-AA", SampleID),
               !grepl("Glu-T12-X1-Rep4-AA", SampleID), !grepl("Glu-T12-X1-Rep6-AA", SampleID),
               !grepl("Glu-T12-X4-Rep1", SampleID), !grepl("Glu-T12-X7-Rep1", SampleID)) %>%
        as_tibble()

    library(rRDPData)
    predict_RDP_taxonomy <- function(sequence_set) {
        pred <- predict(rdp(), sequence_set, confidence = 0) #  Predict 16s
        conf_score <- attr(pred, "confidence") %>% as.data.frame()  # Confidence score
        colnames(pred) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus")
        colnames(conf_score) <- paste0(c("Domain", "Phylum", "Class", "Order", "Family", "Genus"), "Score")
        pred$ID <- rownames(pred) %>% as.numeric()
        conf_score$ID <- rownames(pred) %>% as.numeric()
        return(pred)
    }
    esv_set <- DNAStringSet(community_abundance_subset$ESV)
    names(esv_set) <- community_abundance_subset$CommunityESVID # Rename the sequence
    esv_RDP <- predict_RDP_taxonomy(esv_set) %>%
        mutate(CommunityESVID = factor(ID)) %>%
        select(CommunityESVID, Family, Genus)

    family_name <- c("Pseudomonadaceae", "Enterobacteriaceae", "Aeromonadaceae",  "Sphingobacteriaceae", "Xanthomonadaceae", "Moraxellaceae", "Other")
    #        "Enterococcaceae", "Alcaligenaceae", "Oxalobacteraceae", "Comamonadaceae", "Porphyromonadaceae", "Flavobacteriaceae", "Nocardiaceae", "Sphingomonadaceae", "Brucellaceae")
    family_color <- c("#E21F27", "#397EB8", "#4fb049", "#984f9f", "#a94624", "#8fd1c6", "#989898")
    names(family_color) <- family_name
    genus_name <- c("Pseudomonas", "Klebsiella", "Citrobacter", "Enterobacter", "Raoultella", "Stenotrophomonas", "Aeromonas")
    #, "Yersinia", "Buttiauxella", "Enterococcus", "Erwinia", "Sphingobacterium", "Bordetella", "Acinetobacter", "Salmonella", "Xanthomonas", "Pedobacter", "Pantoea", "Herbaspirillum", "Comamonas", "Dysgonomonas", "Delftia", "Flavobacterium", "Rhodococcus", "Novosphingobium", "Ochrobactrum", "Achromobacter", "Providencia"
    genus_color <- c("#7d1416", "#3eb7c4", "#fdf8ce", "#225fa9", "#a2d6b3", "#f57e2d", "#8fd1c6")
    names(genus_color) <- genus_name

    p_family <- left_join(community_abundance_subset, esv_RDP) %>%
        filter(RelativeAbundance >= 0.01) %>%
        mutate(Family = as.character(Family)) %>%
        select(Family, Genus, Community, Transfer, CommunityESVID, RelativeAbundance) %>%
        mutate(Family = ifelse(Family %in% family_name, Family, "Other") %>% ordered(family_name)) %>%
        mutate(Community = ordered(Community, community_names_ordered_by_size)) %>%
        ggplot() +
        geom_bar(aes(x = Community, y = RelativeAbundance, fill = Family), stat = "identity", color = "grey20") +
        scale_x_discrete(expand = c(0,0)) +
        scale_y_continuous(expand = c(0,0), breaks = seq(0, 1, by = 0.2)) +
        scale_fill_manual(values = family_color) +
        theme_cowplot() +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
        labs(x = "community", y = "relative abundance")

    p_genus <- left_join(community_abundance_subset, esv_RDP) %>%
        filter(RelativeAbundance >= 0.01) %>%
        mutate(Genus = as.character(Genus)) %>%
        select(Family, Genus, Community, Transfer, CommunityESVID, RelativeAbundance) %>%
        mutate(Genus = ifelse(Genus == "Salmonella", "Klebsiella", Genus)) %>%
        mutate(Genus = ifelse(Genus %in% genus_name, Genus, "Other") %>% ordered(genus_name)) %>%
        mutate(Community = ordered(Community, community_names_ordered_by_size)) %>%
        ggplot() +
        geom_bar(aes(x = Community, y = RelativeAbundance, fill = Genus), stat = "identity", color = "grey20") +
        scale_x_discrete(expand = c(0,0)) +
        scale_y_continuous(expand = c(0,0), breaks = seq(0, 1, by = 0.2)) +
        scale_fill_manual(values = genus_color) +
        theme_cowplot() +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
        labs(x = "community", y = "relative abundance")

    ps0 <- plot_grid(p_family, p_genus, ncol = 1, axis = "lr", align = "hv")
    ggsave("../plots/FigS0.png", plot = ps0, width = 6, height = 8)

    # Panel: adjacent matrix
    graph_list_ordered_by_size <- rep(list(NA), length(graph_list))
    for (i in 1:length(graph_list)) graph_list_ordered_by_size[[i]] <- graph_list[[community_names_ordered_by_size[i]]]
    names(graph_list_ordered_by_size) <- community_names_ordered_by_size
    p_mats <- graph_list_ordered_by_size %>% lapply(function(x) {plot_adjacent_matrix(x) + ggtitle("")})
    p_mats[[14]] <- (plot_adjacent_matrix(graph_list_ordered_by_size[[13]]) +
                         theme(legend.position = "right", legend.title = element_blank(), legend.text = element_text(size = 25)) +
                         guides(fill = guide_legend())) %>%
        ggpubr::get_legend() %>%
        ggpubr::as_ggplot()
    ps3 <- plot_grid(plotlist = p_mats, nrow = 2)#, labels = community_names_ordered_by_size)

    ggsave("../plots/FigS3.png", plot = ps3, width = 18, height = 6)

    random_motif_counts_percentile <- random_motif_counts_percentile %>%
        mutate(Motif = factor(Motif)) %>%
        mutate(Community = ordered(Community, levels = community_names_ordered_by_size))
    colors <- c("observed" = "red", "random [5th and 95th percentiles]" = "black")

    ps4 <- summary_network_motifs %>%
        mutate(Community = ordered(Community, levels = community_names_ordered_by_size)) %>%
        ggplot() +
        geom_point(data = random_motif_counts_percentile, aes(x = Motif, y = Count, color = "random [5th and 95th percentiles]")) +
        geom_segment(data = pivot_wider(random_motif_counts_percentile, names_from = Percentile, values_from = Count),
                     aes(x = Motif, xend = Motif, y = p05, yend = p95, color = "random [5th and 95th percentiles]")) +
        geom_point(aes(color = "observed", x = Motif, y = Count)) +
        scale_color_manual(values = colors) +
        facet_wrap(Community~., scales = "free_y", nrow = 2) +
        theme_cowplot() +
        theme(legend.position = "bottom") +
        panel_border(color = "black") +
        labs(color = "")

    ggsave("../plots/FigS4.png", ps4, width = 10, height = 4)

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


