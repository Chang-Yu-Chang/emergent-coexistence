#' Script for supplement figures

library(tidyverse)
library(cowplot)
library(broom)
library(infer) # For tidyverse statistics
library(ggsignif) # For adding asterisks to boxplots
library(RColorBrewer) # For color
library(ggsci) # For color
library(ggraph)
library(tidygraph)
library(officer)
library(flextable)
source(here::here("analysis/00-metadata.R"))

isolates <- read_csv(paste0(folder_data, "output/isolates.csv"), show_col_types = F)
pairs <- read_csv(paste0(folder_data, "output/pairs.csv"), show_col_types = F)
communities <- read_csv(paste0(folder_data, "temp/00c-communities.csv"), show_col_types = F)
accuracy <- read_csv(paste0(folder_data, "temp/91-accuracy.csv"), show_col_types = F)
pairs_freq <- read_csv(paste0(folder_data, "temp/93-pairs_freq.csv"), show_col_types = F)
load(paste0(folder_data, "temp/95-communities_network.Rdata"))
communities_hierarchy <- read_csv(paste0(folder_data, "temp/95-communities_hierarchy.csv"), show_col_types = F)

# 0. Clean up column factors ----
# Arrange communities by size
communities <- communities %>%
    arrange(CommunitySize) %>%
    mutate(Community = factor(Community, Community))

# Clean up the pairs data ----
pairs <- pairs %>%
    # Remove no-colony pairs
    unite(col = "Pair", Community, Isolate1, Isolate2, sep = "_", remove = F) %>%
    filter(!(Pair %in% pairs_no_colony)) %>%
    # Remove low-accuracy model pairs
    filter(AccuracyMean > 0.9)

# Figure S4 model accuracy----
accuracy_to_plot <- accuracy %>%
    # Remove pairs that have cocultures with no colony
    unite(col = "Pair", Community, Isolate1, Isolate2, sep = "_", remove = F) %>%
    filter(!(Pair %in% pairs_no_colony))

p1 <- accuracy_to_plot %>%
    ggplot() +
    geom_histogram(aes(x = Accuracy), color = 1, binwidth = 0.01, breaks = seq(0.6,1,.01), fill = NA) +
    geom_text(x = -Inf, y = Inf, label = paste0("N=", nrow(accuracy_to_plot)), vjust = 2, hjust = -1) +
    geom_vline(xintercept = 0.9, color = "red", linetype = 2) +
    theme_classic() +
    labs(x = "Accuracy", y = "Count")


accuracy_to_plot_count <- accuracy_to_plot %>%
    group_by(AccuracyPassThreshold) %>%
    count(name = "Count")
p2 <- accuracy_to_plot_count %>%
    ggplot() +
    geom_col(aes(x = AccuracyPassThreshold, y = Count), color = 1, fill = NA, width = .8) +
    geom_text(aes(x = AccuracyPassThreshold, y = Count, label = paste0("n=", Count)), vjust = -1) +
    geom_text(x = -Inf, y = Inf, label = paste0("N=", nrow(accuracy_to_plot)), vjust = 2, hjust = -1) +
    scale_x_discrete(labels = c("FALSE" = "Accuracy<0.9", "TRUE" = "Accuracy>0.9")) +
    scale_y_continuous(limits = c(0, 560)) +
    theme_classic() +
    labs(x = "")
p <- plot_grid(p1, p2, nrow = 1, labels = LETTERS[1:2], scale = 0.9) + paint_white_background()
ggsave(here::here("plots/FigS4-random_forest_accuracy.png"), p, width = 8, height = 4)



# Figure S5 machine vs. human ----
pairs_to_include <- pairs %>%
    select(Community, Isolate1, Isolate2) %>%
    mutate(Include = T)
pairs_T8_combined <- read_csv(paste0(folder_data, "temp/92-pairs_T8_combined.csv"), show_col_types = F) %>%
    left_join(pairs_to_include) %>%
    filter(!is.na(Include))

p1 <- pairs_T8_combined %>%
    ggplot() +
    geom_abline(slope = 1, intercept = 0, color = "red", linetype = 2) +
    geom_point(aes(x = TotalCount_human, y = TotalCount_machine), shape = 21, size = 2) +
    geom_text(x = -Inf, y = Inf, label = paste0("N=", nrow(pairs_T8_combined)), vjust = 2, hjust = -1) +
    scale_x_log10() +
    scale_y_log10() +
    theme_classic() +
    labs(x = "Segmentation CFU count", y = "Manual CFU count")


p2 <- pairs_T8_combined %>%
    filter(!is.na(Isolate1Count_human)) %>%
    ggplot() +
    geom_abline(slope = 1, intercept = 0, color = "red", linetype = 2) +
    geom_hline(yintercept = c(0,1), color = gray(.8), linetype = 2) +
    geom_vline(xintercept = c(0,1), color = gray(.8), linetype = 2) +
    geom_point(aes(x = Isolate1CFUFreq_human, y = Isolate1CFUFreq_machine), shape = 21, size = 2, stroke = .4) +
    geom_text(x = 0.5, y = 0.9, label = paste0("N=", pairs_T8_combined %>% filter(!is.na(Isolate1Count_human)) %>% nrow())) +
    scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 21)) +
    theme_classic() +
    labs(x = "Random Forest CFU frequency", y = "Manual CFU frequency")

p <- plot_grid(p1, p2, nrow = 1, axis = "tblr", align = "h", scale = .9, labels = c("A", "B")) +
    paint_white_background()

ggsave(here::here("plots/FigS5-human_machine_comparison.png"), p, width = 8, height = 4)

## cocultures without human results
pairs_T8_combined %>%
    filter(is.na(Isolate1Count_human)) %>%
    nrow()
## R-squared
pairs_T8_combined %>%
    filter(!is.na(Isolate1Count_human)) %>%
    lm(TotalCount_human ~ TotalCount_machine, data = .) %>%
    summary()
pairs_T8_combined %>%
    filter(!is.na(Isolate1Count_human)) %>%
    lm(Isolate1CFUFreq_human ~ Isolate1CFUFreq_machine, data = .) %>%
    summary()


# Figure S6 pairwise competition between highly abundant species ----
# This csv keeps all Sangers that match to one ESV with up to 0-2 mismatch
isolates_abundance_loose <- read_csv(paste0(folder_data, "temp/13-isolates_abundance_loose.csv"), show_col_types = F)

# Highly abundant isolate species have >0.05 relative ESV abundance
isolates_abundant <- isolates_abundance_loose %>%
    select(Community, Isolate, RelativeAbundance) %>%
    filter(RelativeAbundance > 0.05)

isolates_abundant_richness <- isolates_abundant %>%
    group_by(Community) %>%
    summarize(Richness = n()) %>%
    left_join(communities) %>%
    arrange(CommunityLabel)

#
pairs_abundant <- pairs %>%
    select(PairID, Community, Isolate1, Isolate2, InteractionType) %>%
    left_join(rename(isolates_abundant, Isolate1 = Isolate, RelativeAbundance1 = RelativeAbundance)) %>%
    left_join(rename(isolates_abundant, Isolate2 = Isolate, RelativeAbundance2 = RelativeAbundance)) %>%
    filter(!is.na(RelativeAbundance1), !is.na(RelativeAbundance2)) %>%
    left_join(pairs)

pairs_abundant_richness <- pairs_abundant %>%
    group_by(Community) %>%
    summarize(Richness = n()) %>%
    left_join(communities) %>%
    arrange(CommunityLabel)

# Stats
pairs_abundant %>%
    group_by(InteractionType) %>%
    count(name = "Count") %>%
    ungroup() %>%
    mutate(Fraction = Count / sum(Count))

#
p <- pairs_abundant %>%
    filter(!is.na(FitnessFunction)) %>%
    group_by(Community, InteractionType) %>%
    count(name = "Count") %>%
    group_by(Community) %>% mutate(Fraction = Count / sum(Count), TotalCount = sum(Count)) %>%
    ungroup() %>%
    mutate(Community = factor(Community, communities$Community)) %>%
    arrange(Community) %>%
    left_join(communities, by = "Community") %>%
    mutate(CommunityLabel = factor(CommunityLabel)) %>%
    replace_na(list(InteractionType = "unknown")) %>%
    ggplot() +
    geom_col(aes(x = CommunityLabel, fill = InteractionType, y = Fraction), color = 1, width = .8, linewidth = .5) +
    annotate("text", x = 1:13, y = 1.15, label = isolates_abundant_richness$Richness, size = 4) +
    annotate("text", x = 14, y = 1.15, label = "n. of species", size = 4, hjust = 0) +
    annotate("segment", x = .5, xend = 18, y = 1.1, yend = 1.1, color = "black") +
    annotate("text", x= 1:13, y = 1.05, label = pairs_abundant_richness$Richness, size = 4) +
    annotate("text", x = 14, y = 1.05, label = "n. of tested pairs", size = 4, hjust = 0) +
    scale_fill_manual(values = assign_interaction_color(), breaks = c("coexistence", "exclusion", "unknown")) +
    scale_x_discrete(breaks = 1:13, expand = c(0.01, 0)) +
    scale_y_continuous(breaks = c(0,.5,1), limit = c(0, 1.3), expand = c(0,0)) +
    coord_cartesian(xlim = c(0.5, 13.5), ylim = c(0, 1), clip = "off") +
    theme_classic() +
    theme(legend.text = element_text(size = 12),
          axis.text = element_text(color = 1, size = 12),
          #axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          axis.title = element_text(color = 1, size = 12),
          legend.title = element_blank(),
          legend.position = "right",
          plot.margin = unit(c(2,.5,.5,.5), "cm")
    ) +
    labs(x = "Community", y = "Fraction")

ggsave(here::here("plots/FigS6-pairwise_competition_abundant.png"), p, width = 10, height = 4)



# Figure S7 16S mismatch versus probability of coexistence ----
# Data with mismatch
pairs_mismatch <- pairs %>%
    filter(!is.na(Mismatch)) %>%
    mutate(MismatchGroup = case_when(
        Mismatch < 90 ~ "<90",
        Mismatch >= 90 ~ ">=90"
    )) %>%
    filter(InteractionType != "unknown") %>%
    mutate(InteractionType = factor(InteractionType, c("exclusion", "coexistence"))) %>%
    mutate(InteractionTypeBinary = ifelse(InteractionType == "coexistence", 1, 0))

pairs_mismatch_group <- pairs_mismatch %>%
    group_by(MismatchGroup, InteractionType) %>%
    summarize(Count = n()) %>%
    group_by(MismatchGroup) %>%
    mutate(Fraction = Count / sum(Count), TotalCount = sum(Count))

pairs_function_group <- pairs_mismatch %>%
    group_by(PairFermenter, InteractionType) %>%
    summarize(Count = n()) %>%
    group_by(PairFermenter)


# Statistics
extract_statistics_equation <- function (model) {
    p_value <- summary(model)$coefficients[2,4]
    # For glm model
    r_square <- with(summary(model), 1 - deviance/null.deviance)
    paste0(
        "\n", case_when(
            p_value < 0.001 ~ "p<0.001",
            p_value < 0.0001 ~ "p<0.0001",
            TRUE ~ paste0("p=", as.character(round(p_value, 3)))
        ),
        "\nR-squared=", round(r_square,4))
}
model1 <- glm(InteractionType ~ Mismatch, family = "binomial", data = pairs_mismatch)

pairs_mismatch %>%
    group_by(Mismatch, InteractionType) %>%
    chisq_test(response = MismatchGroup)

# scatterplot, categorical lm
p1 <- pairs_mismatch %>%
    ggplot(aes(x = Mismatch, y = InteractionTypeBinary)) +
    geom_point(shape = 21, size = 3) +
    geom_smooth(method = "glm", method.args = list(family = "binomial")) +
    annotate("text", x = 100, y = 0.8, label = extract_statistics_equation(model1), hjust = 0) +
    scale_y_continuous(labels = c("0" = "exclusion", "1" = "coexistence"), breaks = c(0,1), expand = c(0,.2)) +
    theme_classic() +
    labs(x = "# of nucleotide difference in 16S", y = "")


# boxplot by mismatch
p2 <- pairs_mismatch %>%
    group_by(MismatchGroup, InteractionType) %>%
    summarize(Count = n()) %>%
    mutate(Fraction = Count / sum(Count), TotalCount = sum(Count)) %>%
    ggplot() +
    geom_col(aes(x = MismatchGroup, y = Fraction, fill = InteractionType), color = 1) +
    geom_text(aes(x = MismatchGroup, label = paste0("n=", TotalCount)), y = 0.9) +
    scale_fill_manual(values = interaction_color) +
    scale_y_continuous(expand = c(0,0)) +
    theme_classic() +
    theme(legend.position = "right") +
    guides(fill = guide_legend(title = "")) +
    labs(x = "# of nucleotide difference in 16S")

p <- plot_grid(p1, p2, nrow = 1, labels = LETTERS[1:2], scale = 0.9, axis = "tb", align = "h") + paint_white_background()
ggsave(here::here("plots/FigS7-mismatch_vs_coexistence.png"), p, width = 9, height = 4)


if (FALSE) {

# Figure S7a 16S mismatch versus rank difference
## Data with mismatch
pairs_mismatch <- pairs %>%
    left_join(rename_with(select(isolates, ExpID, Rank), ~paste0(.x, "1")), by = "ExpID1") %>%
    left_join(rename_with(select(isolates, ExpID, Rank), ~paste0(.x, "2")), by = "ExpID2") %>%
    left_join(select(communities, Community, CommunitySize), by = "Community") %>%
    mutate(RankDifference = abs(Rank1 - Rank2)) %>%
    mutate(RankDifferenceStandardized = RankDifference / CommunitySize) %>%
    filter(!is.na(Mismatch)) %>%
    mutate(MismatchGroup = case_when(
        Mismatch < 90 ~ "<90",
        Mismatch >= 90 ~ ">=90"
    ))

## Linear regression
extract_statistics_equation <- function (model) {
    p_value <- summary(model)$coefficients[2,4]
    paste0(#"y=", round(tidy(model)$estimate[1], 2), "+", round(tidy(model)$estimate[2], 2), "x",
           "\n", case_when(
               p_value < 0.001 ~ "p<0.001",
               p_value < 0.0001 ~ "p<0.0001",
               TRUE ~ paste0("p=", as.character(round(p_value, 3)))
           ),
           "\nR-squared=", round(summary(model1)$r.squared,2))
}
model1 <- lm(RankDifference ~ Mismatch, data = pairs_mismatch)
tidy(model1)
model2 <- lm(RankDifferenceStandardized ~ Mismatch, data = pairs_mismatch)
summary(model2)
tidy(model2)

## t test
model3 <- pairs_mismatch %>%
    t_test(formula = RankDifference ~ MismatchGroup, order = c("<90", ">=90"), alternative = "two-sided")
model4 <- pairs_mismatch %>%
    t_test(formula = RankDifferenceStandardized ~ MismatchGroup, order = c("<90", ">=90"), alternative = "two-sided")


## Rank difference raw
p1 <- pairs_mismatch %>%
    ggplot() +
    geom_point(aes(x = Mismatch, y = RankDifference), shape = 21, size = 2, stroke = 1) +
    geom_smooth(formula = y~x, aes(x = Mismatch, y = RankDifference), method = "lm") +
    annotate("text", x = 10, y = 11, label = extract_statistics_equation(model1), hjust = 0) +
    scale_y_continuous(limits = c(0,12), breaks = 1:12) +
    theme_classic() +
    labs(x = "# of nucleotide difference in 16S", y = expression(paste("|", R[x] - R[y], "|")))

## Rank difference, standardized
p2 <- pairs_mismatch %>%
    ggplot() +
    geom_point(aes(x = Mismatch, y = RankDifferenceStandardized), shape = 21, size = 2, stroke = 1) +
    geom_smooth(formula = y~x, aes(x = Mismatch, y = RankDifferenceStandardized), method = "lm") +
    annotate("text", x = 1, y = 0.9, label = extract_statistics_equation(model2), hjust = 0) +
    scale_y_continuous(limits = c(0,1), breaks = c(0, 0.5, 1)) +
    theme_classic() +
    labs(x = "# of nucleotide difference in 16S", y = expression(paste("|", R[x] - R[y], "|", "/", N[c])))

## Discrete data
p3 <- pairs_mismatch %>%
    ggplot(aes(x = MismatchGroup, y = RankDifference)) +
    geom_boxplot(width = .5) +
    geom_signif(comparisons = list(c("<90", ">=90")), map_signif_level = TRUE, ) +
    geom_jitter(shape = 21, height = 0, width = .2) +
    scale_y_continuous(limits = c(0,12), breaks = 1:12) +
    theme_classic() +
    labs(x = "# of nucleotide difference in 16S", y = expression(paste("|", R[x] - R[y], "|")))

p4 <- pairs_mismatch %>%
    ggplot(aes(x = MismatchGroup, y = RankDifferenceStandardized)) +
    geom_boxplot(width = .5) +
    geom_signif(comparisons = list(c("<90", ">=90")), map_signif_level = TRUE) +
    geom_jitter(shape = 21, height = 0, width = .2) +
    scale_y_continuous(limits = c(0,1), breaks = c(0, 0.5, 1)) +
    theme_classic() +
    labs(x = "# of nucleotide difference in 16S", y = expression(paste("|", R[x] - R[y], "|", "/", N[c])))


p <- plot_grid(p1, p2, p3, p4, nrow = 2, labels = LETTERS[1:4], scale = 0.9, align = "hv") + paint_white_background()
ggsave(here::here("plots/FigS7-mismatch_vs_rank.png"), p, width = 9, height = 8)
}


# Figure S8 ESV abundance vs. strain ranks ----
## Isolate rank data
isolates_rank <- isolates %>%
    left_join(select(communities, Community, CommunitySize), by = "Community") %>%
    mutate(RankStandardized = Rank / CommunitySize) %>%
    mutate(RankRelativeAbundanceStandardized = RankRelativeAbundance / CommunitySize) %>%
    filter(!is.na(RelativeAbundance))

## Analysis
cor.test(isolates_rank$RankRelativeAbundance, isolates_rank$Rank, method = c("kendall")) %>%
    broom::tidy()

cor.test(isolates_rank$RelativeAbundance, isolates_rank$Rank, method = c("kendall")) %>%
    broom::tidy()


##
p1 <- isolates_rank %>%
    ggplot(aes(x = RelativeAbundance, y = Rank)) +
    geom_point(shape = 21, size = 2, stroke = 1, position = position_jitter(width = 0.15, height = 0.15)) +
    geom_smooth(formula = y~x, method = "lm") +
    #annotate("text", x = 1, y = 12, label = extract_statistics_equation(model1), hjust = 0) +
    #scale_x_continuous(breaks = 1:12) +
    scale_y_continuous(breaks = 1:12) +
    theme_classic() +
    labs(x = "Ranked ESV abundance", y = "Competitive rank")

p <- plot_grid(p1, p2, nrow = 1, labels = LETTERS[1:2], scale = 0.9, align = "hv") + paint_white_background()
ggsave(here::here("plots/FigS8-abundance_vs_rank.png"), p, width = 8, height = 4)



# Figure S9 ESV abundance vs. relative frequencies in pairwise coexistence ----
pairs_freq_ESV <- pairs_freq %>%
    left_join(pairs, by = join_by(Community, Isolate1, Isolate2)) %>%
    filter(InteractionType == "coexistence") %>%
    group_by(PairID, ExpID1, ExpID2) %>%
    summarize(meanIsolate1CFUFreqMean = mean(Isolate1CFUFreqMean)) %>%
    left_join(rename_with(select(isolates, ExpID, RelativeAbundance), ~paste0(.x, "1")), by = "ExpID1") %>%
    left_join(rename_with(select(isolates, ExpID, RelativeAbundance), ~paste0(.x, "2")), by = "ExpID2") %>%
    mutate(RelativeRelativeAbundance1 = RelativeAbundance1 / (RelativeAbundance1 + RelativeAbundance2)) %>%
    mutate(RelativeRelativeAbundance2 = RelativeAbundance2 / (RelativeAbundance1 + RelativeAbundance2)) %>%
    filter(!is.na(RelativeAbundance1))

p <- pairs_freq_ESV %>%
    ggplot(aes(x = RelativeAbundance1, y = meanIsolate1CFUFreqMean)) +
    geom_point(shape = 21, size = 2, stroke = .5) +
    geom_smooth(formula = y~x, method = "lm") +
    scale_x_continuous(limits = c(0,1)) +
    scale_y_continuous(limits = c(0,1)) +
    theme_classic() +
    labs(x = "ESV abundance", y = "relative CFU frequency")

#p <- plot_grid(p1, p2, nrow = 1, labels = LETTERS[1:2], scale = 0.9, align = "hv") + paint_white_background()
ggsave(here::here("plots/FigS9-abundance_vs_pair_frequency.png"), p, width = 4, height = 3)

## Analysis
lm(meanIsolate1CFUFreqMean ~ RelativeAbundance1, data = pairs_freq_ESV) %>% summary()




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
        geom_edge_link(aes(color = InteractionType), width = edge_width/2) +
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

ggsave(here::here("plots/FigS10-community_graph.png"), p, width = 13, height = 4)

# Figure SXX h score robustness ----
# library(igraph)
# g <- make_ring(10)
# transitivity(g)
# plot(g)
# g2 <- sample_gnp(100, 10/100)
# plot(g2)
# transitivity(g2)   # this is about 10/1000



# Table S1 image object features ----
features_example <- read_csv(paste0(folder_pipeline, "images/D-07-feature/merged/D_T8_C1R2_1.csv"), show_col_types = F)

features <- tibble(Feature = names(features_example)) %>%
    filter(!Feature == "ObjectID") %>%
    mutate(`Feature type` = case_when(
        str_sub(Feature, 1, 1) == "s" ~ "shape",
        str_sub(Feature, 1, 1) == "m" ~ "moment",
        str_sub(Feature, 1, 1) == "b" & str_sub(Feature, 1, 6) != "b.tran" ~ "intensity",
        str_sub(Feature, 1, 6) == "b.tran" ~ "transect"
    )) %>%
    mutate(Description =
               c("area size (in pixels)",
                 "perimeter (in pixels)",
                 "mean radius (in pixels)",
                 "standard deviation of the mean radius (in pixels)",
                 "min radius (in pixels)",
                 "max radius (in pixels)",
                 "center of mass x (in pixels)",
                 "center of mass y (in pixels)",
                 "elliptical fit major axis (in pixels)",
                 "elliptical eccentricity defined by sqrt(1-minoraxis^2/majoraxis^2). Circle eccentricity is 0 and straight line eccentricity is 1.",
                 "object angle (in radians)",
                 "green channel mean intensity",
                 "green channel standard deviation intensity",
                 "green channel mad intensity",
                 "green channel 1st quantile intensity",
                 "green channel 5th quantile intensity",
                 "green channel 10th quantile intensity",
                 "green channel 20th quantile intensity",
                 "green channel 50th quantile intensity",
                 "green channel 80th quantile intensity",
                 "green channel 90th quantile intensity",
                 "green channel 95th quantile intensity",
                 "green channel 99th quantile intensity",
                 "green channel transect mean intensity",
                 "green channel transect standard deviation intensity",
                 "green channel transect mad intensity",
                 "green channel the inmost transect pixel intensity",
                 "green channel the outmost transect pixel intensity",
                 "green channel difference between b.center and b.periphery",
                 "green channel pixel intensity at the 5% on the scaled transect",
                 "green channel pixel intensity at the 10% on the scaled transect",
                 "green channel pixel intensity at the 50% on the scaled transect",
                 "green channel pixel intensity at the 90% on the scaled transect",
                 "green channel pixel intensity at the 95% on the scaled transect",
                 "red channel mean intensity",
                 "red channel standard deviation intensity",
                 "red channel mad intensity",
                 "blue channel mean intensity",
                 "blue channel standard deviation intensity",
                 "blue channel mad intensity"
               )) %>%
    mutate(` ` = 1:n()) %>%
    select(` `, everything())

ft <- features %>%
    flextable() %>%
    width(j = 1, width = 1) %>%
    width(j = 2:3, width = 2) %>%
    width(j = 4, width = 5) %>%
    vline(j = 1, border = NULL, part = "all")

save_as_image(ft, here::here("plots/TableS1-features.png"), webshot = "webshot2")


# Table S2 community overview and labels ----

temp1 <- isolates %>%
    distinct(Community, Batch) %>%
    filter(!(Community == "C11R1" & Batch == "B2")) %>%
    rows_update(tibble(Community = "C11R1", Batch = "B2 and C"), by = "Community")
temp2 <- pairs %>%
    group_by(Community) %>%
    filter(!is.na(InteractionType)) %>%
    filter(AccuracyMean > 0.9) %>%
    count(name = "ActualPairs") %>%
    ungroup()
ft <- communities %>%
    left_join(temp1, by = "Community") %>%
    left_join(temp2, by = "Community") %>%
    select(CommunityLabel, Community, Batch, CommunitySize, CommunityPairSize, ActualPairs) %>%
    mutate(CommunityLabel = as.character(CommunityLabel)) %>%
    rename(`Community` = CommunityLabel, `Internal label` = Community, `Number of isolates` = CommunitySize,
           `Number of tested pairs` = CommunityPairSize, `Number of applicable pairs` = ActualPairs) %>%
    janitor::adorn_totals() %>%
    flextable() %>%
    width(j = 1:2, width = 1.1) %>%
    width(j = 3, width = 0.8) %>%
    width(j = 4:6, width = 1.3) %>%
    hline(i = 13, border = fp_border(color = "black", style = "solid", width = 2))

save_as_image(ft, here::here("plots/TableS2-communities.png"), webshot = "webshot2")

# Table S3 list of isolates, images used for monocultures, CFU, OD ----
isolates$OD620
isolates_epsilon <- read_csv(paste0(folder_data, "temp/06-isolates_epsilon.csv"), show_col_types = F) %>%
    select(Batch, Community, Isolate, Time, image_name, ColonyCount, OD620, Epsilon) %>%
    mutate(Community = factor(Community, communities$Community)) %>%
    mutate(OD620 = round(OD620, 3)) %>%
    mutate(Epsilon = format(Epsilon, scientific = T, digits = 2)) %>%
    arrange(Batch, Community) %>%
    rename(`Image name` = image_name, `Colony count` = ColonyCount)


t1 <- isolates_epsilon %>% slice(1:35)
t2 <- isolates_epsilon %>% slice(36:n())

ft1 <- t1 %>%
    flextable() %>%
    width(j = c(1:4, 6), width = .8) %>%
    width(j = 5, width = 1) %>%
    width(j = 6, width = 1.5) %>%
    highlight(i = which(t1$Time %in% c("T0", "T1")), j = 4, color = "yellow")
ft2 <- t2 %>%
    slice(1:35) %>%
    flextable() %>%
    width(j = c(1:4, 6), width = .8) %>%
    width(j = 5, width = 1) %>%
    width(j = 6, width = 1.5) %>%
    highlight(i = which(t2$Time %in% c("T0", "T1")), j = 4, color = "yellow")

save_as_image(ft1, here::here("plots/TableS3-monoculture_1.png"), webshot = "webshot2")
save_as_image(ft2, here::here("plots/TableS3-monoculture_2.png"), webshot = "webshot2")



# Table S4 fitness function ----
make_interaction_type <- function () {
    #' This function generates the fitness function table.
    #' There are a total of 27 possibilities
    interaction_type <- tibble(
        FromRare = rep(c(1, -1, 0), each = 9),
        FromMedium = rep(rep(c(1, -1, 0), each = 3), 3),
        FromAbundant = rep(c(1, -1, 0), 9),
        InteractionType = NA,
        InteractionTypeFiner = NA
    )
    # Assign interaction types to combinations of frequency changes signs
    interaction_type$InteractionType[c(1,10,13,14)] <- "exclusion"
    interaction_type$InteractionType[c(2:6,8,11,20,23, 9,18,21,24,25,26,27)] <- "coexistence"

    # Assign finer interaction types to combinations of frequency changes signs
    interaction_type$InteractionTypeFiner[c(1,14)] <- "competitive exclusion"
    interaction_type$InteractionTypeFiner[c(10,13)] <- "mutual exclusion"
    interaction_type$InteractionTypeFiner[c(2,5,8)] <- "stable coexistence"
    interaction_type$InteractionTypeFiner[c(4,6,11,20)] <- "frequency-dependent coexistence"
    interaction_type$InteractionTypeFiner[c(3)] <- "coexistence at 95%"
    interaction_type$InteractionTypeFiner[c(23)] <- "coexistence at 5%"
    interaction_type$InteractionTypeFiner[c(9,18,21,24:26)] <- "neutrality"
    interaction_type$InteractionTypeFiner[c(9,18,21,24:26)] <- "2-freq neutrality"
    interaction_type$InteractionTypeFiner[c(27)] <- "3-freq neutrality"
    interaction_type <- interaction_type %>%  mutate(FitnessFunction = paste(FromRare, FromMedium, FromAbundant, sep = "_"))

    # Fill in unknown
    interaction_type %>%
        replace_na(list(InteractionType = "unknown", InteractionTypeFiner = "unknown"))
}
reformat <- function(x) {
    x <- ifelse(x == "1", "+", x)
    x <- ifelse(x == "-1", "-", x)
    x <- factor(x, c("+", "-", "0"))
    return(x)
}
interaction_type <- make_interaction_type()

pairs_interaction <- pairs %>%
    filter(!is.na(InteractionType)) %>%
    filter(AccuracyMean > 0.9) %>%
    group_by(InteractionType, InteractionTypeFiner, FitnessFunction) %>%
    count(name = "Count")

ft <- interaction_type %>%
    left_join(pairs_interaction, by = c("InteractionType", "InteractionTypeFiner", "FitnessFunction")) %>%
    replace_na(list(Count = 0)) %>%
    mutate(across(starts_with("From"), reformat)) %>%
    select(-FitnessFunction) %>%
    rename(`Outcome` = InteractionType, `Finer outcome` = InteractionTypeFiner, `From rare` = FromRare, `From medium` = FromMedium, `From abundant` = FromAbundant) %>%
    mutate(` ` = 1:n()) %>%
    select(` `, everything()) %>%
    janitor::adorn_totals() %>%
    flextable() %>%
    width(j = 1, width = .5) %>%
    width(j = 2:5, width = 1.3) %>%
    width(j = 6, width = 3) %>%
    vline(j = 1, border = NULL, part = "all") %>%
    hline(i = 27, border = fp_border(color = "black", style = "solid", width = 2))

save_as_image(ft, here::here("plots/TableS4-fitness_function.png"), webshot = "webshot2")




























