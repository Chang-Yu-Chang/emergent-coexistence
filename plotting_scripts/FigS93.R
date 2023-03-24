library(tidyverse)
library(cowplot)
library(gridExtra)
source(here::here("analysis/00-metadata.R"))
#source(here::here("plotting_scripts/Fig2.R"))

pairs <- read_csv(paste0(folder_data, "output/pairs.csv"), show_col_types = F)
pairs_freq <- read_csv(paste0(folder_data, "temp/93-pairs_freq.csv"), show_col_types = F)
load(paste0(folder_data, "temp/95-communities_network.Rdata"))

# Clean up pairs data
pairs <- pairs %>%
    # Remove no-colony pairs
    unite(col = "Pair", Community, Isolate1, Isolate2, sep = "_", remove = F) %>%
    filter(!(Pair %in% pairs_no_colony)) %>%
    # Remove low-accuracy model pairs
    filter(AccuracyMean > 0.9)

pairs_freq <- pairs_freq %>%
    left_join(distinct(pairs, PairID, Community, Isolate1, Isolate2)) %>%
    select(PairID, everything()) %>%
    filter(!is.na(PairID))

# Negative frequency dependent of isolate
p <- pairs_freq %>%
    #filter(Community == "C8R4") %>%
    select(PairID, Isolate1, Isolate2, Isolate1InitialODFreq, Time, Isolate1CFUFreqMean) %>%
    #group_by(PairID, Isolate1, Isolate2) %>%
    #pivot_wider(names_from = Isolate1InitialODFreq, names_prefix = "F", values_from = Isolate1CFUFreqMean) %>%
    pivot_wider(names_from = Time, values_from = Isolate1CFUFreqMean) %>%
    mutate(Fitness = log(T8/T0)) %>%
    ggplot() +
    geom_rect(data = pairs, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, aes(fill = InteractionTypeFiner), alpha = .4) +
    geom_smooth(aes(x = T0, y = Fitness), method = "lm", formula = y~x, se = F, linewidth = .3) +
    geom_point(aes(x = T0, y = Fitness), shape = 21, stroke = .5, size = 1) +
    geom_hline(yintercept = 0, linetype = 2) +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 3), limits = c(0,1), expand = c(0,0)) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 3), expand = c(0, 0.5)) +
    #scale_y_log10() +
    scale_fill_manual(values = assign_interaction_color(level = "finer")) +
    facet_wrap(.~PairID, ncol = 10, scales = "free_y") +
    theme_classic() +
    theme(
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 5),
        strip.text = element_text(size = 8),
        strip.background = element_blank(),
        legend.position = "bottom",
        panel.spacing = unit(0, "mm"),
        panel.border = element_rect(color = 1, linewidth = 0.5, fill = NA)
    ) +
    guides() +
    labs()

ggsave(here::here("plots/FigS93-pairs_negative_freq_dep.png"), p, width = 8, height = 12)



interaction_color

get_interaction_legend <- function (pairs) {
    panel_colors <- c(interaction_color[1], rep(interaction_color[2], 6), interaction_color[3]) %>% setNames(pairs_interaction_finer$InteractionTypeFiner)
    panel_fills <- assign_interaction_color(level = "finer")[-3]
    temp <- pairs %>%
        ggplot() +
        geom_tile(aes(x = Isolate1, y = Isolate2, color = InteractionType, fill = InteractionTypeFiner), height = .8, width = .8, alpha = .9, size = 1) +
        scale_color_manual(values = panel_colors,
                           breaks = names(panel_colors),
                           labels = paste0(pairs_interaction_finer$InteractionTypeFiner, " (", round(pairs_interaction_finer$Fraction, 3) * 100,"%)")) +
        scale_fill_manual(values = panel_fills,
                          breaks = names(panel_fills),
                          labels = paste0(pairs_interaction_finer$InteractionTypeFiner, " (", round(pairs_interaction_finer$Fraction, 3) * 100,"%)")) +
        theme(legend.position = "right",
              legend.spacing.y = unit(0.2, "cm"),
              legend.text = element_text(size = 12),
              legend.background = element_blank(),
              legend.margin = margin(0,0,0,0)) +
        guides(fill = guide_legend(byrow = TRUE), color = "none") +
        labs(color = "", fill = "") +
        paint_white_background()

    return(get_legend(temp))
}
plot_example_freq <- function(pairs_freq, line_size = 0.5) {
    axis_lower <- -0.5
    axis_upper <- 1.5
    pairs_freq %>%
        mutate(Time = case_when(Time == "T0" ~ 0,Time == "T8" ~ 1)) %>%
        mutate(Isolate1InitialODFreq = factor(Isolate1InitialODFreq, c(95,50,5))) %>%
        ggplot() +
        geom_rect(aes(color = InteractionType, fill = InteractionTypeFiner),
                  xmin = axis_lower, xmax = axis_upper, ymin = axis_lower+0.2, ymax = axis_upper-0.2,
                  alpha = .08, linewidth = .8) +
        geom_hline(linewidth = line_size, yintercept = c(0,1), linetype = 1, color = "grey90") +
        geom_line(aes(x = Time, y = Isolate1CFUFreqMean, color = Isolate1InitialODFreq, group = Isolate1InitialODFreq),
                  linewidth = line_size*2) +
        geom_point(aes(x = Time, y = Isolate1CFUFreqMean, color = Isolate1InitialODFreq, group = Isolate1InitialODFreq),
                   size = line_size*2, ) +
        geom_segment(aes(x = Time, xend = Time,
                         y = Isolate1CFUFreqMean + 2*Isolate1CFUFreqSd,
                         yend = Isolate1CFUFreqMean - 2*Isolate1CFUFreqSd,
                         color = Isolate1InitialODFreq),
                     linewidth = line_size*2, ) +
        scale_x_continuous(expand = c(0, 0), limits = c(axis_lower, axis_upper)) +
        scale_y_continuous(expand = c(0, 0), breaks = c(0, .5, 1), limits = c(axis_lower+0.2, axis_upper-0.2)) +
        #scale_color_manual(values = c(frequency_color, interaction_color), label = c("95%", "50%", "5%", "exclusion", "coexistence")) +
        scale_color_manual(values = frequency_color, label = c("95%", "50%", "5%")) +
        scale_fill_manual(values = assign_interaction_color(level = "finer")) +
        theme_classic() +
        theme(panel.spacing = unit(2, "mm"),
              panel.border = element_blank(),
              panel.grid.minor.y = element_blank(),
              axis.line = element_blank(),
              axis.title = element_blank(), axis.text = element_blank(),
              axis.ticks = element_blank(),
              plot.background = element_blank(),
              plot.title = element_text(size = 5, margin = margin(0,0,0,0)),
              plot.margin = margin(0,0,0,0, "mm")) +
        guides(color = "none", fill = "none") +
        labs(x = "Time", y = "Frequency") +
        #ggtitle(unique(pairs_freq$PairID)) +
        NULL
}

pairs_freq %>%
    left_join(pairs) %>%
    filter(PairID == 3) %>%
    plot_example_freq()












