library(tidyverse)
library(cowplot)
source(here::here("analysis/00-metadata.R"))

communities <- read_csv(paste0(folder_data, "temp/00c-communities.csv"), show_col_types = F)
pairs <- read_csv(paste0(folder_data, "output/pairs.csv"), show_col_types = F)
pairs_freq <- read_csv(paste0(folder_data, "temp/93a-pairs_freq.csv"), show_col_types = F)

pairs_freq <- pairs_freq %>% left_join(pairs) %>% remove_ineligible_pairs()
pairs <- remove_ineligible_pairs(pairs)

pairs_freq_mean_three <- pairs_freq %>%
    filter(Time == "T8") %>%
    flip_winner_species_freq() %>%
    group_by(PairID, Community, Isolate1, Isolate2) %>%
    summarize(Isolate1CFUFreq_mean_three = mean(Isolate1CFUFreqMean))



line_size = 1
p <- pairs_freq %>%
    # Order by outcome
    flip_winner_species_freq %>%
    mutate(Isolate1InitialODFreq = factor(Isolate1InitialODFreq, c(95,50,5))) %>%
    mutate(Time = as.character(Time) %>% str_replace("T", "") %>% as.numeric()) %>%
    mutate(Nudge = case_when(Isolate1InitialODFreq == 5 ~ .3, Isolate1InitialODFreq == 50 ~.6, Isolate1InitialODFreq == 95 ~ .9)) %>%
    ggplot() +
    geom_rect(aes(fill = outcome), color = 1, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, alpha = .03, linewidth = .8) +
    geom_hline(linewidth = line_size/2, yintercept = c(0,1), linetype = 1, color = "black") +
    geom_hline(data = pairs_freq_mean_three, aes(yintercept = Isolate1CFUFreq_mean_three), linewidth = line_size/2, color = "black", linetype = 2) +
    geom_line(aes(x = Time+Nudge, y = Isolate1CFUFreqMean, color = Isolate1InitialODFreq, group = Isolate1InitialODFreq), linewidth = line_size) +
    #geom_point(aes(x = Time+Nudge, y = Isolate1CFUFreqMean, color = Isolate1InitialODFreq, group = Isolate1InitialODFreq), size = line_size*2, shape = 21) +
    geom_segment(aes(x = Time+Nudge, xend = Time+Nudge, y = Isolate1CFUFreqPercentile5, yend = Isolate1CFUFreqPercentile95, color = factor(Isolate1InitialODFreq)), linewidth = line_size) +
    scale_x_continuous(breaks = c(0,8), limits = c(-3,11)) +
    scale_y_continuous(breaks = c(0, .5, 1)) +
    scale_color_manual(values = frequency_color, label = c("95%", "50%", "5%")) +
    scale_fill_manual(values = outcome_colors, labels = outcome_labels) +
    facet_wrap(.~PairID, nrow = 10, dir = "v") +
    theme_classic() +
    theme(
        panel.spacing = unit(2, "mm"),
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



