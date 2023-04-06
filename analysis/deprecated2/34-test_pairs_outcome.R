library(tidyverse)
library(cowplot)
library(ggsci)
source(here::here("analysis/00-metadata.R"))

communities <- read_csv(paste0(folder_data, "output/communities_remained.csv"), show_col_types = F)
pairs_ID <- read_csv(paste0(folder_data, "temp/00c-pairs_ID.csv"), show_col_types = F)
pairs <- read_csv(paste0(folder_data, "output/pairs_remained.csv"), show_col_types = F) # 145 pairs
pairs_freq <- read_csv(paste0(folder_data, "temp/93a-pairs_freq.csv"), show_col_types = F)
pairs_freq <- pairs_freq %>% left_join(pairs) %>% remove_ineligible_pairs()

pairs_T0_boots <- read_csv(paste0(folder_data, "temp/06-pairs_T0_boots.csv"), show_col_types = F) # bootstraps using T0 mean and sd
pairs_T8_boots <- read_csv(paste0(folder_data, "temp/06-pairs_T8_boots.csv"), show_col_types = F) # bootstraps using random forest object probabilities=

pairs_boots <- bind_rows(pairs_T0_boots, pairs_T8_boots) %>%
    select(Community, Isolate1, Isolate2, Isolate1InitialODFreq, Time, BootstrapID, Isolate1CFUFreq)

# Check the pairwise outcome ----
p <- pairs %>%
    group_by(Community, outcome) %>%
    count(name = "Count") %>%
    # Total count
    group_by(Community) %>% mutate(Fraction = Count / sum(Count), TotalCount = sum(Count)) %>%
    left_join(communities, by = "Community") %>%
    replace_na(list(outcome = "inconclusive")) %>%
    ungroup() %>%
    ggplot() +
    geom_col(aes(x = CommunityLabel, fill = outcome, y = Fraction), color = 1, width = .8, linewidth = .5) +
    annotate("text", x = 1:13, y = 1.15, label = communities$CommunitySize, size = 4) +
    annotate("text", x = 14, y = 1.15, label = "n. of species", size = 4, hjust = 0) +
    annotate("segment", x = .5, xend = 18, y = 1.1, yend = 1.1, color = "black") +
    geom_text(aes(x = CommunityLabel, y = 1.05, label = TotalCount), size = 4) +
    annotate("text", x = 14, y = 1.05, label = "n. of tested pairs", size = 4, hjust = 0) +
    scale_fill_manual(values = outcome_colors) +
    scale_x_continuous(breaks = 1:13, expand = c(0.01, 0)) +
    scale_y_continuous(breaks = seq(0,1,0.2), limit = c(0, 1.3), expand = c(0,0)) +
    coord_cartesian(xlim = c(0.5, 13.5), ylim = c(0, 1), clip = "off") +
    theme_classic() +
    theme(legend.text = element_text(size = 10),
          legend.title = element_blank(),
          legend.key.size = unit(.5, "cm"),
          legend.spacing.y = unit(.3, "cm"),
          legend.position = "right",
          panel.border = element_rect(color = 1, fill = NA),
          axis.text = element_text(color = 1, size = 10),
          axis.title = element_text(color = 1, size = 10),
          plot.margin = unit(c(1,.5,.5,.5), "cm")
    ) +
    guides(fill = guide_legend(byrow = TRUE)) +
    labs(x = "community", y = "fraction")

ggsave(paste0(folder_data, "temp/34-01-pairs_outcome.png"), p, width = 6, height = 3)

pairs %>%
    group_by(outcome) %>%
    count(name = "Count") %>%
    ungroup() %>%
    mutate(Fraction = Count / sum(Count), TotalCount = sum(Count))


# Test plot Figure 3----

# Find 5% and 95% CIs of for each T0 and T8 freqeunecy
pairs_boots_percentile <- pairs_boots %>%
    group_by(Community, Isolate1, Isolate2, Isolate1InitialODFreq, Time) %>%
    arrange(Isolate1CFUFreq) %>%
    mutate(Percentile = paste0("P", 0.1 * (1:1000))) %>%
    filter(Percentile %in% c("P5", "P95")) %>%
    ungroup() %>%
    select(Community, Isolate1, Isolate2, Isolate1InitialODFreq, Time, Isolate1CFUFreq, Measure = Percentile)

pairs_boots_mean <- pairs_boots %>%
    group_by(Community, Isolate1, Isolate2, Isolate1InitialODFreq, Time) %>%
    summarize(Isolate1CFUFreqMean = mean(Isolate1CFUFreq),
              Isolate1CFUFreqMedian = median(Isolate1CFUFreq)) %>%
    ungroup() %>%
    pivot_longer(cols = starts_with("Isolate1CFUFreq"), names_prefix = "Isolate1CFUFreq",
                 names_to = "Measure", values_to = "Isolate1CFUFreq")


# Check if the percentile makes sense for each pairs ----
pairs_boots_comm <- pairs_boots %>%
    filter(Community == "C1R2", Isolate1 == 1, Isolate2 == 2)
pairs_boots_percentile_comm <- pairs_boots_percentile %>%
    filter(Community == "C1R2", Isolate1 == 1, Isolate2 == 2)
pairs_boots_mean_comm <- pairs_boots_mean %>%
    filter(Community == "C1R2", Isolate1 == 1, Isolate2 == 2)

p <- pairs_boots_comm %>%
    ggplot() +
    geom_histogram(aes(x = Isolate1CFUFreq), color = "black", fill = "white", binwidth = 0.01) +
    geom_vline(data = pairs_boots_percentile_comm, aes(xintercept = Isolate1CFUFreq, color = Measure), linewidth = 1) +
    geom_vline(data = pairs_boots_mean_comm, aes(xintercept = Isolate1CFUFreq, color = Measure), linewidth = 1) +
    scale_color_manual(values = RColorBrewer::brewer.pal(6, "Set2")) +
    facet_grid(Isolate1InitialODFreq ~ Time, scales = "free_y") +
    coord_flip() +
    #theme_bw() +
    theme_classic() +
    theme() +
    guides() +
    labs()

ggsave(paste0(folder_data, "temp/34-02-check_percentiles.png"), p, width = 6, height = 6)


# Check the means ----
line_size = 1
p <- pairs_boots_mean %>%
    filter(Measure == "Mean") %>%
    # Order by outcome
    left_join(pairs) %>% arrange(outcome, PairID) %>% mutate(PairID = factor(PairID, unique(PairID))) %>%
    drop_na(outcome) %>%
    mutate(Isolate1InitialODFreq = factor(Isolate1InitialODFreq, c(95,50,5))) %>%
    ggplot() +
    geom_rect(aes(fill = outcome), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, alpha = .03, linewidth = .8) +
    geom_hline(linewidth = line_size/2, yintercept = c(0,1), linetype = 1, color = "black") +
    geom_line(aes(x = Time, y = Isolate1CFUFreq, color = Isolate1InitialODFreq, group = Isolate1InitialODFreq), linewidth = line_size) +
    geom_point(aes(x = Time, y = Isolate1CFUFreq, color = Isolate1InitialODFreq, group = Isolate1InitialODFreq), size = line_size*2, shape = 21) +
    scale_y_continuous(breaks = c(0, .5, 1)) +
    scale_color_manual(values = frequency_color, label = c("95%", "50%", "5%")) +
    scale_fill_manual(values = outcome_colors) +
    facet_wrap(.~PairID, nrow = 10, dir = "v") +
    theme_classic() +
    theme(panel.spacing = unit(2, "mm"),
          panel.border = element_blank(),
          panel.grid.minor.y = element_blank(),
          strip.background = element_rect(color = NA, fill = "white"),
          plot.background = element_rect(color = NA, fill = "white"),
          plot.margin = margin(0,0,0,0, "mm")) +
    guides(alpha = "none") +
    labs(x = "Time", y = "Frequency")

ggsave(paste0(folder_data, "temp/34-03-check_means.png"), p, width = 20, height = 10)


# Check the percentiles on plot ----
pairs_boots_percentile_wider <- pairs_boots_percentile %>%
    pivot_wider(names_from = Measure, values_from = Isolate1CFUFreq) %>%
    left_join(pairs) %>% arrange(outcome, PairID) %>% mutate(PairID = factor(PairID, unique(PairID))) %>%
    drop_na(outcome)
line_size = 1
p <- pairs_boots_mean %>%
    filter(Measure == "Mean") %>%
    # Order by outcome
    left_join(pairs) %>% arrange(outcome, PairID) %>% mutate(PairID = factor(PairID, unique(PairID))) %>%
    drop_na(outcome) %>%
    #filter(PairID %in% 1:10) %>%
    mutate(Isolate1InitialODFreq = factor(Isolate1InitialODFreq, c(95,50,5))) %>%
    ggplot() +
    geom_rect(aes(fill = outcome), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, alpha = .03, linewidth = .8) +
    geom_hline(linewidth = line_size/2, yintercept = c(0,1), linetype = 1, color = "black") +
    geom_line(aes(x = Time, y = Isolate1CFUFreq, color = Isolate1InitialODFreq, group = Isolate1InitialODFreq), linewidth = line_size) +
    geom_point(aes(x = Time, y = Isolate1CFUFreq, color = Isolate1InitialODFreq, group = Isolate1InitialODFreq), size = line_size*2, shape = 21) +
    geom_segment(data = pairs_boots_percentile_wider, linewidth = line_size,
                 aes(x = Time, xend = Time, y = P5, yend = P95, color = factor(Isolate1InitialODFreq))) +
    scale_y_continuous(breaks = c(0, .5, 1)) +
    scale_color_manual(values = frequency_color, label = c("95%", "50%", "5%")) +
    scale_fill_manual(values = outcome_colors) +
    facet_wrap(.~PairID, nrow = 10, dir = "v") +
    theme_classic() +
    theme(panel.spacing = unit(2, "mm"),
          panel.border = element_blank(),
          panel.grid.minor.y = element_blank(),
          strip.background = element_rect(color = NA, fill = "white"),
          plot.background = element_rect(color = NA, fill = "white"),
          plot.margin = margin(0,0,0,0, "mm")) +
    guides(alpha = "none") +
    labs(x = "Time", y = "Frequency")

ggsave(paste0(folder_data, "temp/34-04-check_percentiles.png"), p, width = 20, height = 10)


# Add the mean of the three frequnecies ----
pairs_boots_mean_three <- pairs_boots_mean %>%
    filter(Measure == "Mean", Time == "T8") %>%
    group_by(Community, Isolate1, Isolate2) %>%
    summarize(Isolate1CFUFreq_mean_three = mean(Isolate1CFUFreq)) %>%
    left_join(pairs) %>% arrange(outcome, PairID) %>% mutate(PairID = factor(PairID, unique(PairID))) %>%
    drop_na(outcome)

line_size = 1
p <- pairs_boots_mean %>%
    filter(Measure == "Mean") %>%
    # Order by outcome
    left_join(pairs) %>% arrange(outcome, PairID) %>% mutate(PairID = factor(PairID, unique(PairID))) %>%
    drop_na(outcome) %>%
    mutate(Isolate1InitialODFreq = factor(Isolate1InitialODFreq, c(95,50,5))) %>%
    ggplot() +
    geom_rect(aes(fill = outcome), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, alpha = .03, linewidth = .8) +
    geom_hline(data = pairs_boots_mean_three, aes(yintercept = Isolate1CFUFreq_mean_three), linewidth = line_size/2, color = "black", linetype = 2) +
    geom_hline(linewidth = line_size/2, yintercept = c(0,1), linetype = 1, color = "black") +
    geom_line(aes(x = Time, y = Isolate1CFUFreq, color = Isolate1InitialODFreq, group = Isolate1InitialODFreq), linewidth = line_size) +
    geom_point(aes(x = Time, y = Isolate1CFUFreq, color = Isolate1InitialODFreq, group = Isolate1InitialODFreq), size = line_size*2, shape = 21) +
    geom_segment(data = pairs_boots_percentile_wider, linewidth = line_size,
                 aes(x = Time, xend = Time, y = P5, yend = P95, color = factor(Isolate1InitialODFreq))) +
    scale_y_continuous(breaks = c(0, .5, 1)) +
    scale_color_manual(values = frequency_color, label = c("95%", "50%", "5%")) +
    scale_fill_manual(values = outcome_colors) +
    facet_wrap(.~PairID, nrow = 10, dir = "v") +
    theme_classic() +
    theme(panel.spacing = unit(2, "mm"),
          panel.border = element_blank(),
          panel.grid.minor.y = element_blank(),
          strip.background = element_rect(color = NA, fill = "white"),
          plot.background = element_rect(color = NA, fill = "white"),
          plot.margin = margin(0,0,0,0, "mm")) +
    guides(alpha = "none") +
    labs(x = "Time", y = "Frequency")

ggsave(paste0(folder_data, "temp/34-05-check_mean_three.png"), p, width = 20, height = 10)


# Find CIs and menas for means of equilibrium frequencies ----
pairs_mean_eq <- pairs_boots %>% # 185000 rows
    filter(Time == "T8") %>%
    group_by(Community, Isolate1, Isolate2, BootstrapID) %>%
    # Mean of the three equilibrium frequencies
    summarize(MeanIsolate1CFUFreq = mean(Isolate1CFUFreq))
## 5th and 95th percentile
pairs_mean_eq_percentile <- pairs_mean_eq %>% # 185 rows
    arrange(Community, Isolate1, Isolate2, MeanIsolate1CFUFreq) %>%
    mutate(Percentile = paste0("Percentile", 0.1 * (1:1000))) %>%
    filter(Percentile %in% c("Percentile5", "Percentile95")) %>%
    ungroup() %>%
    select(Community, Isolate1, Isolate2, MeanIsolate1CFUFreq, Measure = Percentile) %>%
    pivot_wider(names_from = Measure, names_prefix = "MeanIsolate1CFUFreq", values_from = MeanIsolate1CFUFreq)
## mean
pairs_mean_eq_mean <- pairs_mean_eq %>% # 185 rows
    arrange(Community, Isolate1, Isolate2, MeanIsolate1CFUFreq) %>%
    #group_by(Community, Isolate1, Isolate2) %>%
    summarize(MeanMeanIsolate1CFUFreq = mean(MeanIsolate1CFUFreq))

pairs_mean_eq_measures <- pairs_mean_eq_mean %>% left_join(pairs_mean_eq_percentile) %>%
    left_join(pairs_ID)

# Check if the mean and percentiles make sense
pairs_mean_eq_measures_comm <- pairs_mean_eq_measures %>%
    filter(Community == "C1R2") %>%
    left_join(pairs_ID)
p <- pairs_mean_eq %>%
    filter(Community == "C1R2") %>%
    left_join(pairs_ID) %>%
    ggplot() +
    geom_histogram(aes(x = MeanIsolate1CFUFreq), color = "black", fill = "white") +
    # Mean
    geom_vline(data = pairs_mean_eq_measures_comm, aes(xintercept = MeanMeanIsolate1CFUFreq, color = "mean")) +
    # Percentile 5th
    geom_vline(data = pairs_mean_eq_measures_comm, aes(xintercept = MeanIsolate1CFUFreqPercentile5, color = "5th")) +
    # Percentile 95th
    geom_vline(data = pairs_mean_eq_measures_comm, aes(xintercept = MeanIsolate1CFUFreqPercentile95, color = "95th")) +
    theme_classic() +
    facet_wrap(.~PairID, nrow = 3, dir = "v", scales = "free") +
    theme(
        panel.spacing = unit(2, "mm"),
        panel.border = element_blank(),
        panel.grid.minor.y = element_blank(),
        legend.position = "right",
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 15),
        legend.key.size = unit(10, "mm"),
        legend.spacing.y = unit(5, "mm"),
        #strip.text = element_blank(),
        plot.background = element_rect(color = NA, fill = "white")
    ) +
    guides(
        fill = guide_legend(title = "Pairwise competiton outcome", override.aes = list(alpha = 0.8), order = 1),
        color = guide_legend(title = "Inital frequencies", order = 2)
    ) +
    labs()

ggsave(paste0(folder_data, "temp/34-06-check_mean_eq.png"), p, width = 6, height = 6)

#
write_csv(pairs_mean_eq_measures, paste0(folder_data, "temp/34-pairs_mean_eq_measures.csv"))


# Chekc the new figure design. No flipping freqnuecy---
pairs_mean_eq_mean <- pairs_boots %>%
    filter(Time == "T8") %>%
    group_by(Community, Isolate1, Isolate2, BootstrapID) %>%
    # Mean of the three equilibrium frequencies
    summarize(MeanIsolate1CFUFreq = mean(Isolate1CFUFreq)) %>%
    arrange(Community, Isolate1, Isolate2, MeanIsolate1CFUFreq) %>%
    summarize(MeanMeanIsolate1CFUFreq = mean(MeanIsolate1CFUFreq)) %>%
    # Flip the frequency if the mean equilibirum freq > 0.5
    mutate(FlipFreq = ifelse(MeanMeanIsolate1CFUFreq > 0.5, T, F))

pairs_mean_eq_measures <- read_csv(paste0(folder_data, "temp/34-pairs_mean_eq_measures.csv"), show_col_types = F) %>%
    right_join(select(pairs, PairID, outcome)) %>%
    filter(outcome %in% c("3-coexistence", "4-coexistence"))

line_size = 1
p <- pairs_freq %>%
    # Order by outcome
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
    scale_y_continuous(breaks = c(0, .5, 1), expand = c(0, 0.1), labels = c("0", "0.5", "1")) +
    scale_color_manual(values = frequency_color, label = c("95%", "50%", "5%")) +
    scale_fill_manual(values = outcome_colors, labels = outcome_labels) +
    facet_wrap(.~PairID, nrow = 10, dir = "v") +
    theme_classic() +
    theme(
        panel.spacing = unit(0, "mm"),
        panel.border = element_blank(),
        panel.grid.minor.y = element_blank(),
        legend.position = "right",
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 15),
        legend.key.size = unit(10, "mm"),
        legend.spacing.y = unit(5, "mm"),
        strip.text = element_blank(),
        plot.background = element_rect(color = NA, fill = "white")
    ) +
    guides(
        fill = guide_legend(title = "Pairwise competiton outcome", override.aes = list(alpha = 0.8), order = 1),
        color = guide_legend(title = "Inital frequencies", order = 2)
    ) +
    labs(x = "transfer", y = "frequency")

ggsave(paste0(folder_data, "temp/34-07-test_new_fig3.png"), p, width = 20, height = 10)

# Check the new figure design. Flip frequnecy if the mean eq freq > 0.5 ----
pairs_mean_eq_mean <- pairs_boots %>%
    filter(Time == "T8") %>%
    group_by(Community, Isolate1, Isolate2, BootstrapID) %>%
    # Mean of the three equilibrium frequencies
    summarize(MeanIsolate1CFUFreq = mean(Isolate1CFUFreq)) %>%
    arrange(Community, Isolate1, Isolate2, MeanIsolate1CFUFreq) %>%
    summarize(MeanMeanIsolate1CFUFreq = mean(MeanIsolate1CFUFreq)) %>%
    # Flip the frequency if the mean equilibirum freq > 0.5
    mutate(FlipFreq = ifelse(MeanMeanIsolate1CFUFreq > 0.5, T, F))

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

pairs_mean_eq_measures <- read_csv(paste0(folder_data, "temp/34-pairs_mean_eq_measures.csv"), show_col_types = F) %>%
    right_join(select(pairs, PairID, outcome)) %>%
    filter(outcome %in% c("3-coexistence", "4-coexistence")) %>%
    flip_winner_species_freq_shade(pairs_mean_eq_mean)

line_size = 1
p <- pairs_freq %>%
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
    scale_y_continuous(breaks = c(0, .5, 1), expand = c(0, 0.1), labels = c("0", "0.5", "1")) +
    scale_color_manual(values = frequency_color, label = c("95%", "50%", "5%")) +
    scale_fill_manual(values = outcome_colors, labels = outcome_labels) +
    facet_wrap(.~PairID, nrow = 10, dir = "v") +
    theme_classic() +
    theme(
        panel.spacing = unit(0, "mm"),
        panel.border = element_blank(),
        panel.grid.minor.y = element_blank(),
        legend.position = "right",
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 15),
        legend.key.size = unit(10, "mm"),
        legend.spacing.y = unit(5, "mm"),
        strip.text = element_blank(),
        plot.background = element_rect(color = NA, fill = "white")
    ) +
    guides(
        fill = guide_legend(title = "Pairwise competiton outcome", override.aes = list(alpha = 0.8), order = 1),
        color = guide_legend(title = "Inital frequencies", order = 2)
    ) +
    labs(x = "transfer", y = "frequency")

ggsave(paste0(folder_data, "temp/34-08-test_new_fig3_flipped.png"), p, width = 20, height = 10)



# Check the new figure design, divide the panels ----
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
        facet_wrap(.~PairID, nrow = 10, dir = "v") +
        theme_classic() +
        theme(
            panel.spacing = unit(0, "mm"),
            panel.border = element_blank(),
            panel.grid.minor.y = element_blank(),
            axis.line = element_blank(),
            legend.position = "none",
            legend.title = element_text(size = 20),
            legend.text = element_text(size = 15),
            legend.key.size = unit(10, "mm"),
            legend.spacing.y = unit(5, "mm"),
            strip.text = element_blank(),
            plot.background = element_rect(color = outcome_colors[outcome_category], fill = "white", linewidth = 3)
        ) +
        guides(
            fill = guide_legend(title = "Pairwise competiton outcome", override.aes = list(alpha = 0.8), order = 1),
            color = guide_legend(title = "Inital frequencies", order = 2)
        ) +
        labs(x = "transfer", y = "frequency")
}

p1 <- plot_category_freq(pairs_freq, pairs_mean_eq_measures, outcome_category = "1-exclusion")
p2 <- plot_category_freq(pairs_freq, pairs_mean_eq_measures, outcome_category = "2-exclusion")
p3 <- plot_category_freq(pairs_freq, pairs_mean_eq_measures, outcome_category = "3-coexistence")
p4 <- plot_category_freq(pairs_freq, pairs_mean_eq_measures, outcome_category = "4-coexistence")
p5 <- plot_category_freq(pairs_freq, pairs_mean_eq_measures, outcome_category = "5-inconclusive")
p_legend <- get_legend(plot_category_freq(pairs_freq, pairs_mean_eq_measures, outcome_category = names(outcome_colors)) + theme(legend.position = "bottom") + guides(fill = "none", color = guide_legend(title = "Initial frequencies", override.aes = list(linewidth = 2))))

p <- plot_grid(
    plot_grid(
        p1, p2, p3, p4, p5, nrow = 1, axis = "tb", align = "h",
        rel_widths = c(4, 3, 3, 3, 4), scale = 0.9,
        labels = outcome_labels, label_x = 0.05, hjust = 0
    ),
    p_legend, ncol = 1, rel_heights = c(10, 1)
) + paint_white_background()
ggsave(paste0(folder_data, "temp/34-09-test_new_fig3_panels.png"), p, width = 20, height = 10)

# Check the new figure design, divide the panels, ver 2 ----
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
        #geom_rect(aes(fill = outcome), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, alpha = .03, linewidth = .8) +
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
        facet_wrap(.~PairID, nrow = 10, dir = "v") +
        theme_classic() +
        theme(
            panel.spacing = unit(0, "mm"),
            panel.border = element_blank(),
            panel.grid.minor.y = element_blank(),
            axis.line = element_blank(),
            legend.position = "none",
            legend.title = element_text(size = 20),
            legend.text = element_text(size = 15),
            legend.key.size = unit(10, "mm"),
            legend.spacing.y = unit(5, "mm"),
            strip.text = element_blank(),
            plot.background = element_rect(color = outcome_colors[outcome_category], fill = "white", linewidth = 3)
        ) +
        guides(
            fill = guide_legend(title = "Pairwise competiton outcome", override.aes = list(alpha = 0.8), order = 1),
            color = guide_legend(title = "Inital frequencies", order = 2)
        ) +
        labs(x = "transfer", y = "frequency")
}

p1 <- plot_category_freq(pairs_freq, pairs_mean_eq_measures, outcome_category = "1-exclusion")
p2 <- plot_category_freq(pairs_freq, pairs_mean_eq_measures, outcome_category = "2-exclusion")
p3 <- plot_category_freq(pairs_freq, pairs_mean_eq_measures, outcome_category = "3-coexistence")
p4 <- plot_category_freq(pairs_freq, pairs_mean_eq_measures, outcome_category = "4-coexistence")
p5 <- plot_category_freq(pairs_freq, pairs_mean_eq_measures, outcome_category = "5-inconclusive")
p_legend <- get_legend(plot_category_freq(pairs_freq, pairs_mean_eq_measures, outcome_category = names(outcome_colors)) + theme(legend.position = "bottom") + guides(fill = "none", color = guide_legend(title = "Initial frequencies", override.aes = list(linewidth = 2))))

p <- plot_grid(
    plot_grid(
        p1, p2, p3, p4, p5, nrow = 1, axis = "tb", align = "h",
        rel_widths = c(4, 3, 3, 3, 4), scale = 0.9,
        labels = outcome_labels, label_x = 0.05, hjust = 0
    ),
    p_legend, ncol = 1, rel_heights = c(10, 1)
) + paint_white_background()
ggsave(paste0(folder_data, "temp/34-10-test_new_fig3_panels_v2.png"), p, width = 20, height = 10)

