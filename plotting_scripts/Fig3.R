library(tidyverse)
library(cowplot)
source(here::here("analysis/00-metadata.R"))

communities <- read_csv(paste0(folder_data, "output/communities_remained.csv"), show_col_types = F)
pairs <- read_csv(paste0(folder_data, "output/pairs_remained.csv"), show_col_types = F) # 145 pairs
pairs_freq <- read_csv(paste0(folder_data, "temp/93a-pairs_freq.csv"), show_col_types = F)
pairs_freq <- pairs_freq %>% left_join(pairs) %>% remove_ineligible_pairs()

#
pairs_T0_boots <- read_csv(paste0(folder_data, "temp/06-pairs_T0_boots.csv"), show_col_types = F) # bootstraps using T0 mean and sd
pairs_T8_boots <- read_csv(paste0(folder_data, "temp/06-pairs_T8_boots.csv"), show_col_types = F) # bootstraps using random forest object probabilities=

pairs_boots <- bind_rows(pairs_T0_boots, pairs_T8_boots) %>%
    select(Community, Isolate1, Isolate2, Isolate1InitialODFreq, Time, BootstrapID, Isolate1CFUFreq)

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
ggsave(here::here("plots/Fig3.png"), p, width = 16, height = 8)


# Stat
pairs %>%
    group_by(outcome) %>%
    count(name = "Count") %>%
    ungroup() %>%
    mutate(Fraction = Count / sum(Count), TotalCount = sum(Count))





