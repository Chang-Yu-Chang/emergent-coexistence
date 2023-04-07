library(tidyverse)
library(cowplot)
source(here::here("analysis/00-metadata.R"))

pairs_ID <- read_csv(paste0(folder_data, "temp/00c-pairs_ID.csv"), show_col_types = F)
pairs_boots <- read_csv(paste0(folder_data, "temp/07-pairs_boots.csv"), show_col_types = F)
pairs_freq <- read_csv(paste0(folder_data, "temp/25-pairs_freq.csv"), show_col_types = F)
#pairs_outcome <- read_csv(paste0(folder_data, "temp/26-pairs_outcome.csv"), show_col_types = F)


# 1. Check if the percentile makes sense for each pairs ----
pairs_boots_comm <- pairs_boots %>%
    filter(Community == "C1R2", Isolate1 == 1, Isolate2 == 2)
pairs_freq_comm <- pairs_freq %>%
    filter(Community == "C1R2", Isolate1 == 1, Isolate2 == 2)

p <- pairs_boots_comm %>%
    ggplot() +
    geom_histogram(aes(x = Isolate1CFUFreq), color = "black", fill = "white", binwidth = 0.01) +
    geom_vline(data = pairs_freq_comm, aes(xintercept = Isolate1CFUFreqMean, color = "mean"), linewidth = 1) +
    geom_vline(data = pairs_freq_comm, aes(xintercept = Isolate1CFUFreqPercentile5, color = "5th"), linewidth = 1) +
    geom_vline(data = pairs_freq_comm, aes(xintercept = Isolate1CFUFreqPercentile95, color = "95th"), linewidth = 1) +
    scale_color_manual(values = RColorBrewer::brewer.pal(6, "Set2")) +
    facet_grid(Isolate1InitialODFreq ~ Time, scales = "free_y") +
    coord_flip() +
    theme_classic() +
    theme() +
    guides() +
    labs()

ggsave(paste0(folder_data, "temp/26a-01-check_percentiles.png"), p, width = 6, height = 6)

if (FALSE) {

# 2. Check the means ----
line_size = 1
p <- pairs_freq %>%
    # Order by outcome
    left_join(pairs_outcome) %>%
    arrange(outcome, PairID) %>%
    mutate(PairID = factor(PairID, unique(PairID))) %>%
    drop_na(outcome) %>%
    mutate(Isolate1InitialODFreq = factor(Isolate1InitialODFreq, c(95,50,5))) %>%
    ggplot() +
    geom_rect(aes(fill = outcome), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, alpha = .03, linewidth = .8) +
    geom_hline(linewidth = line_size/2, yintercept = c(0,1), linetype = 1, color = "black") +
    geom_line(aes(x = Time, y = Isolate1CFUFreqMean, color = Isolate1InitialODFreq, group = Isolate1InitialODFreq), linewidth = line_size) +
    geom_point(aes(x = Time, y = Isolate1CFUFreqMean, color = Isolate1InitialODFreq, group = Isolate1InitialODFreq), size = line_size*2, shape = 21) +
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

ggsave(paste0(folder_data, "temp/26a-02-check_means.png"), p, width = 20, height = 10)


# 3. Check the percentiles on plot ----
line_size = 1
p <- pairs_freq %>%
    # Order by outcome
    left_join(pairs_outcome) %>% arrange(outcome, PairID) %>% mutate(PairID = factor(PairID, unique(PairID))) %>%
    drop_na(outcome) %>%
    #filter(PairID %in% 1:10) %>%
    mutate(Isolate1InitialODFreq = factor(Isolate1InitialODFreq, c(95,50,5))) %>%
    ggplot() +
    geom_rect(aes(fill = outcome), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, alpha = .03, linewidth = .8) +
    geom_hline(linewidth = line_size/2, yintercept = c(0,1), linetype = 1, color = "black") +
    geom_line(aes(x = Time, y = Isolate1CFUFreqMean, color = Isolate1InitialODFreq, group = Isolate1InitialODFreq), linewidth = line_size) +
    geom_point(aes(x = Time, y = Isolate1CFUFreqMean, color = Isolate1InitialODFreq, group = Isolate1InitialODFreq), size = line_size*2, shape = 21) +
    geom_segment(aes(x = Time, xend = Time, y = Isolate1CFUFreqPercentile5, yend = Isolate1CFUFreqPercentile95, color = factor(Isolate1InitialODFreq)),  linewidth = line_size) +
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

ggsave(paste0(folder_data, "temp/26a-03-check_percentiles.png"), p, width = 20, height = 10)


# 4. Add the mean of the three frequenceis ----
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

}
