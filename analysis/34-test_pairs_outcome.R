library(tidyverse)
library(cowplot)
library(ggsci)
source(here::here("analysis/00-metadata.R"))

#pairs_ID <- read_csv(paste0(folder_data, "temp/00c-pairs_ID.csv"), show_col_types = F)
pairs <- read_csv(paste0(folder_data, "output/pairs.csv"), show_col_types = F)
communities <- read_csv(paste0(folder_data, "temp/00c-communities.csv"), show_col_types = F)

pairs <- pairs %>%
    unite(col = "Pair", Community, Isolate1, Isolate2, sep = "_", remove = F) %>%
    filter(!(Pair %in% pairs_no_colony)) %>%
    filter(AccuracyMean > 0.9)

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


# Figure 3----
pairs_T0_boots <- read_csv(paste0(folder_data, "temp/06-pairs_T0_boots.csv"), show_col_types = F) # bootstraps using T0 mean and sd
pairs_T8_boots <- read_csv(paste0(folder_data, "temp/06-pairs_T8_boots.csv"), show_col_types = F) # bootstraps using random forest object probabilities=

pairs_boots <- bind_rows(pairs_T0_boots, pairs_T8_boots) %>%
    select(Community, Isolate1, Isolate2, Isolate1InitialODFreq, Time, BootstrapID, Isolate1CFUFreq)

# Find 5% and 95%
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


# Combine the 5th and 95th percentile
# pairs_freq <- bind_rows(pairs_boots_mean, pairs_boots_percentile) %>%
#     left_join(pairs_ID)

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





