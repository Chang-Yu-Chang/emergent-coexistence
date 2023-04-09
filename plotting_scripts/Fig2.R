library(tidyverse)
library(cowplot)
source(here::here("processing_scripts/00-metadata.R"))

communities <- read_csv(paste0(folder_data, "output/communities_remained.csv"), show_col_types = F)
pairs <- read_csv(paste0(folder_data, "output/pairs_remained.csv"), show_col_types = F)
pairs_freq <- read_csv(paste0(folder_data, "temp/25-pairs_freq.csv"), show_col_types = F)
pairs_boots <- read_csv(paste0(folder_data, "temp/07-pairs_boots.csv"), show_col_types = F)


# Figure 2A: cartoon----
pA <- ggdraw() + draw_image(here::here("plots/cartoons/Fig2A.png")) + paint_white_background()

# Figure 2B: pairwise outcomes per community ----
pB <- pairs %>%
    group_by(Community, outcome) %>%
    count(name = "Count") %>%
    # Total count
    group_by(Community) %>% mutate(Fraction = Count / sum(Count), TotalCount = sum(Count)) %>%
    left_join(communities, by = "Community") %>%
    #mutate(outcome = factor(outcome, c("5-inconclusive", "1-exclusion", "2-exclusion", "3-coexistence", "4-coexistence"))) %>%
    ungroup() %>%
    ggplot() +
    geom_col(aes(x = CommunityLabel, fill = outcome, y = Fraction), width = .8, linewidth = .5, position = position_stack(reverse = T)) +
    # Number of distinct ESVs
    annotate("text", x = 13, y = 1.25, label = "n. of ESVs", size = 4, hjust = 0) +
    annotate("text", x = 1:12, y = 1.25, label = communities$ESVRichness, size = 4) +
    annotate("segment", x = .5, xend = 16, y = 1.2, yend = 1.2, color = "black") +
    # Number of isolates
    annotate("text", x = 13, y = 1.15, label = "n. of isolates", size = 4, hjust = 0) +
    annotate("text", x = 1:12, y = 1.15, label = communities$CommunitySize, size = 4) +
    annotate("segment", x = .5, xend = 16, y = 1.1, yend = 1.1, color = "black") +
    # Number of tested pairs
    annotate("text", x = 13, y = 1.05, label = "n. of tested pairs", size = 4, hjust = 0) +
    geom_text(aes(x = CommunityLabel, y = 1.05, label = TotalCount), size = 4) +
    scale_fill_manual(values = outcome_colors, breaks = names(outcome_colors), labels = outcome_labels) +
    scale_x_continuous(breaks = 1:12, expand = c(0.01, 0)) +
    scale_y_continuous(breaks = seq(0,1,0.2), limit = c(0, 1.45), expand = c(0,0)) +
    coord_cartesian(xlim = c(0.5, 12.5), ylim = c(0, 1), clip = "off") +
    theme_classic() +
    theme(
        legend.text = element_text(size = 12),
        legend.title = element_blank(),
        legend.key.width = unit(8, "mm"),
        legend.key.height = unit(8, "mm"),
        legend.spacing.x = unit(11, "mm"),
        legend.position = "bottom",
        panel.border = element_rect(color = 1, fill = NA),
        axis.text = element_text(color = 1, size = 12),
        axis.title = element_text(color = 1, size = 12),
        plot.margin = unit(c(20, 30, 5, 5), "mm")
    ) +
    guides(fill = guide_legend(byrow = T, nrow = 1)) +
    labs(x = "community", y = "fraction")

pB_legend <- get_legend(pB)

# Figure 2C: the frequency plot ----
flip_winner_species_freq <- function (pairs_freq, pairs_mean_eq_mean) {
    #' Flip the frequency if the mean equilibrium frequency is > 0.5
    pairs_freq <- pairs_freq %>%
        left_join(pairs_mean_eq_mean)
    temp_index <- which(pairs_freq$FlipFreq)
    if (length(temp_index) !=0) {
        pairs_freq_flipped <- pairs_freq[temp_index, ] %>%
            rename(temp = Isolate1, Isolate1 = Isolate2) %>%
            rename(Isolate2 = temp) %>%
            mutate(
                Isolate1CFUFreqMean = 1-Isolate1CFUFreqMean,
                Isolate1InitialODFreq = 100-Isolate1InitialODFreq,
                Isolate1CFUFreqMedian = 1-Isolate1CFUFreqMedian,
                Isolate1CFUFreqPercentile5 = 1-Isolate1CFUFreqPercentile5,
                Isolate1CFUFreqPercentile95 = 1-Isolate1CFUFreqPercentile95
            )

        bind_rows(pairs_freq[-temp_index, ], pairs_freq_flipped) %>%
            arrange(Time, PairID) %>%
            return()
    } else if (length(temp_index) == 0) {
        return(pairs_freq)
    }
}
flip_winner_species_freq_shade <- function (pairs_mean_eq_measures, pairs_mean_eq_mean) {
    #' Flip the frequency if the mean equilibrium frequency is > 0.5
    pairs_mean_eq_measures <- pairs_mean_eq_measures %>%
        left_join(select(pairs_mean_eq_mean, Community, Isolate1, Isolate2, FlipFreq))
    temp_index <- which(pairs_mean_eq_measures$FlipFreq)
    if (length(temp_index) !=0) {
        pairs_freq_flipped <- pairs_mean_eq_measures[temp_index, ] %>%
            rename(temp = Isolate1, Isolate1 = Isolate2) %>%
            rename(Isolate2 = temp) %>%
            mutate(
                MeanMeanIsolate1CFUFreq = 1-MeanMeanIsolate1CFUFreq,
                MeanIsolate1CFUFreqPercentile5 = 1-MeanIsolate1CFUFreqPercentile5,
                MeanIsolate1CFUFreqPercentile95 = 1-MeanIsolate1CFUFreqPercentile95
            )

        bind_rows(pairs_mean_eq_measures[-temp_index, ], pairs_freq_flipped) %>%
            arrange(PairID) %>%
            return()
    } else if (length(temp_index) == 0) {
        return(pairs_mean_eq_measures)
    }
}

pairs_mean_eq_mean <- pairs_boots %>%
    filter(Time == "T8") %>%
    group_by(Community, Isolate1, Isolate2, BootstrapID) %>%
    # Mean of the three equilibrium frequencies
    summarize(MeanIsolate1CFUFreq = mean(Isolate1CFUFreq)) %>%
    arrange(Community, Isolate1, Isolate2, MeanIsolate1CFUFreq) %>%
    summarize(MeanMeanIsolate1CFUFreq = mean(MeanIsolate1CFUFreq)) %>%
    # Flip the frequency if the mean equilibirum freq > 0.5
    mutate(FlipFreq = ifelse(MeanMeanIsolate1CFUFreq > 0.5, T, F))


pairs_mean_eq_measures <- read_csv(paste0(folder_data, "temp/25-pairs_mean_eq_measures.csv"), show_col_types = F) %>%
    right_join(select(pairs, PairID, outcome)) %>%
    filter(outcome %in% c("3-coexistence", "4-coexistence")) %>%
    flip_winner_species_freq_shade(pairs_mean_eq_mean)
pairs_freq <- mutate(pairs_freq, -PairID) %>%
    left_join(select(pairs, PairID, outcome, AccuracyMean)) %>%
    remove_ineligible_pairs()

line_size = 1
plot_category_freq <- function (pairs_freq, pairs_mean_eq_measures, outcome_category = "1-exclusion") {
    pairs_mean_eq_measures <- pairs_mean_eq_measures %>%
        filter(outcome %in% outcome_category)
    pairs_freq %>%
        filter(outcome %in% outcome_category) %>%
        # Order by outcome
        flip_winner_species_freq(pairs_mean_eq_mean) %>%
        mutate(Isolate1InitialODFreq = factor(Isolate1InitialODFreq, c(95,50,5))) %>%
        mutate(Time = as.character(Time) %>% str_replace("T", "") %>% as.numeric()) %>%
        mutate(Nudge = case_when(Isolate1InitialODFreq == 5 ~ .3, Isolate1InitialODFreq == 50 ~.6, Isolate1InitialODFreq == 95 ~ .9)) %>%
        ggplot() +
        geom_rect(aes(fill = outcome), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, alpha = .03, linewidth = .8) +
        geom_hline(linewidth = line_size/3, yintercept = c(0,1), linetype = 2, color = "black") +
        # CIs of eq freq
        geom_rect(data = pairs_mean_eq_measures, aes(ymin = MeanIsolate1CFUFreqPercentile5, ymax = MeanIsolate1CFUFreqPercentile95),
                  xmin = -Inf, xmax = Inf, fill = grey(0.85), color = NA) +
        # Mean of eq freq
        geom_hline(data = pairs_mean_eq_measures, aes(yintercept = MeanMeanIsolate1CFUFreq),
                   linewidth = line_size/2, color = grey(0.95), linetype = 2) +
        # Frequency change
        geom_line(aes(x = Time+Nudge, y = Isolate1CFUFreqMean, color = Isolate1InitialODFreq, group = Isolate1InitialODFreq), linewidth = line_size) +
        #geom_point(aes(x = Time+Nudge, y = Isolate1CFUFreqMean, color = Isolate1InitialODFreq, group = Isolate1InitialODFreq), size = line_size*2, shape = 21) +
        geom_segment(aes(x = Time+Nudge, xend = Time+Nudge, y = Isolate1CFUFreqPercentile5, yend = Isolate1CFUFreqPercentile95, color = factor(Isolate1InitialODFreq)), linewidth = line_size) +
        scale_x_continuous(breaks = c(0,8), limits = c(-3,11)) +
        scale_y_continuous(breaks = c(0, 1), expand = c(0, 0.1)) +
        scale_color_manual(values = frequency_color, label = c("95%", "50%", "5%")) +
        scale_fill_manual(values = outcome_colors, labels = outcome_labels) +
        facet_wrap(.~PairID, nrow = 8, dir = "v") +
        theme_classic() +
        theme(
            panel.spacing = unit(0, "mm"),
            panel.border = element_blank(),
            panel.grid.minor.y = element_blank(),
            axis.line = element_blank(),
            legend.position = "none",
            legend.background = element_blank(),
            legend.title = element_text(size = 20),
            legend.text = element_text(size = 15),
            legend.key.size = unit(10, "mm"),
            legend.spacing.y = unit(5, "mm"),
            strip.text = element_blank(),
            plot.title = element_text(color = outcome_colors[outcome_category], size = 13, face = "bold", margin = margin(0,0,0,0, "mm")),
            plot.background = element_blank()
        ) +
        guides(
            fill = guide_legend(title = "Pairwise competiton outcome", override.aes = list(alpha = 0.8), order = 1),
            color = guide_legend(title = "Inital frequencies", order = 2)
        ) +
        labs(x = "transfer", y = "frequency") +
        ggtitle(outcome_labels[which(names(outcome_colors) == outcome_category)])
}
outcome_labels <- c("competitive exclusion",
                    "on the path to\ncompetitive exclusion",
                    "stable coexistence\n(mutual invasibility)",
                    "coexistence \nwithout\nevidence of\nmutual\ninvasibility",
                    "inconclusive")

p1 <- plot_category_freq(pairs_freq, pairs_mean_eq_measures, outcome_category = "1-exclusion")
p2 <- plot_category_freq(pairs_freq, pairs_mean_eq_measures, outcome_category = "2-exclusion") + theme(axis.title.y = element_blank())
p3 <- plot_category_freq(pairs_freq, pairs_mean_eq_measures, outcome_category = "3-coexistence") + theme(axis.title.y = element_blank())
p4 <- plot_category_freq(pairs_freq, pairs_mean_eq_measures, outcome_category = "4-coexistence") + theme(axis.title.y = element_blank())
p5 <- plot_category_freq(pairs_freq, pairs_mean_eq_measures, outcome_category = "5-inconclusive") + theme(axis.title.y = element_blank())
pC_legend <- get_legend(plot_category_freq(pairs_freq, pairs_mean_eq_measures, outcome_category = names(outcome_colors)) + theme(legend.position = "bottom", legend.text = element_text(size = 13), legend.title = element_text(size = 13)) + guides(fill = "none", color = guide_legend(title = "initial frequencies", override.aes = list(linewidth = 2))))


pC <- plot_grid(
    plot_grid(
        p1, p2, p3, p4, p5, nrow = 1, axis = "tb", align = "h",
        rel_widths = c(5, 9, 3, 1.5, 2), scale = .99
        #labels = outcome_labels, label_x = 0.05, hjust = 0, label_y = .97, vjust = 0, label_size = 13
    ),
    pC_legend, ncol = 1, rel_heights = c(10, .5)
) + paint_white_background()

# Assemble panels ----
p_top <- plot_grid(pA, pB + guides(fill = "none"), nrow = 1, labels = c("A", "B"), scale = c(1, 0.95), rel_widths = c(1.3, 1), axis = "b")
p <- plot_grid(p_top, NULL, pC, ncol = 1, labels = c("", "", "C"),
               scale = c(1, 1, .95), rel_heights = c(1, 0, 2)) + paint_white_background()
ggsave(here::here("plots/Fig2.png"), p, width = 15, height = 10)

# Save vector based
p_top <- plot_grid(pA, pB + guides(fill = "none"), nrow = 1, labels = c("A", "B"), scale = c(1, 0.95), rel_widths = c(1.3, 1), axis = "b")
p <- plot_grid(p_top, NULL, pC, ncol = 1, labels = c("", "", "C"),
               scale = c(1, 1, .95), rel_heights = c(1, 0, 2)) + paint_white_background()
ggsave(here::here("plots/Fig2.pdf"), p, width = 15, height = 10)


# Stat
pairs %>%
    group_by(outcome) %>%
    count(name = "Count") %>%
    ungroup() %>%
    mutate(Fraction = Count / sum(Count), TotalCount = sum(Count))





