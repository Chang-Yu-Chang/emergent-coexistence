library(tidyverse)
library(cowplot)
library(gridExtra)
source(here::here("analysis/00-metadata.R"))
source(here::here("plotting_scripts/Fig2.R"))

pairs <- read_csv(paste0(folder_data, "output/pairs.csv"), show_col_types = F)
pairs_freq <- read_csv(paste0(folder_data, "temp/93-pairs_freq.csv"), show_col_types = F)
load(paste0(folder_data, "temp/95-communities_network.Rdata"))

# Clean up pairs data ----
pairs <- pairs %>%
    # Remove no-colony pairs
    unite(col = "Pair", Community, Isolate1, Isolate2, sep = "_", remove = F) %>%
    filter(!(Pair %in% pairs_no_colony)) %>%
    # Remove low-accuracy model pairs
    filter(AccuracyMean > 0.9)

#
get_interaction_legend <- function (pairs) {
    panel_colors <- c(interaction_color[1], rep(interaction_color[2], 6), interaction_color[3]) %>% setNames(pairs_interaction_finer$InteractionTypeFiner)
    panel_fills <- assign_interaction_color(level = "finer")[-3]
    temp <- pairs %>%
        #mutate(InteractionType = factor(InteractionType)) %>%
        #mutate(InteractionTypeFiner = factor(InteractionTypeFiner)) %>%
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

# Legend for fill
pairs_interaction_finer <- pairs %>%
    group_by(InteractionType, InteractionTypeFiner) %>%
    count(name = "Count") %>%
    ungroup() %>%
    mutate(Fraction = Count / sum(Count)) %>%
    mutate(InteractionTypeFiner = factor(InteractionTypeFiner, interaction_type_finer)) %>%
    arrange(InteractionTypeFiner)
p_legend_fill <- pairs %>%
    get_interaction_legend()

# Append competition outcome to frequencies
pairs_example_freq <- pairs %>%
    filter(!is.na(FitnessFunction)) %>%
    mutate(InteractionType = factor(InteractionType, c("exclusion", "coexistence"))) %>%
    mutate(InteractionTypeFiner = factor(InteractionTypeFiner, interaction_type_finer)) %>%
    arrange(InteractionType, InteractionTypeFiner) %>%
    left_join(pairs_freq, multiple = "all") %>%
    select(PairID, InteractionType, InteractionTypeFiner, Isolate1InitialODFreq, Time, Isolate1CFUFreqMean, Isolate1CFUFreqSd) %>%
    mutate(Isolate1InitialODFreq = factor(Isolate1InitialODFreq), Time = factor(Time, c("T0", "T8")))

temp_list <- pairs_example_freq %>%
    arrange(InteractionTypeFiner, PairID) %>%
    group_split(InteractionTypeFiner, PairID) %>%
    lapply(plot_example_freq)

p_example <- pairs_example_freq %>%
    filter(PairID == 3) %>%
    plot_example_freq(line_size = 0.5) +
    theme(axis.text = element_text(size = 10), axis.title = element_text(size = 10)) +
    scale_x_continuous(breaks = c(0,1), labels = c("start", "end"), limits = c(-0.5, 1.5)) +
    labs(y = "CFU frequency", color = "Inital OD frequency")

## Grid layout
# m <- matrix(c(1:nrow(pairs), rep(NA, 10-(nrow(pairs)%%10))), nrow = 10)
m <- matrix(c(1:nrow(pairs), rep(NA, 15-(nrow(pairs)%%15))), nrow = 15)
p_waffle <- arrangeGrob(grobs = temp_list, layout_matrix = m)
p_waffle <- plot_grid(NULL, p_waffle, rel_widths = c(1, 2.3), scale = c(1, 0.95)) + paint_white_background()

ss = .3
p <- ggdraw(p_waffle) +
    draw_plot(p_example, x = 0.01, y = .7, width = ss*0.6, height = ss*0.8, hjust = 0, vjust = 0) +
    draw_plot(p_legend_color, x = .18, y = .75, width = ss/2, height = ss/2, hjust = 0, vjust = 0) +
    draw_plot(p_legend_fill, x = .01, y = .2, width = ss, height = ss, hjust = 0, vjust = 0) +
    theme(panel.background = element_blank(), plot.background = element_rect(color = NA, fill = "white"),
          plot.margin = unit(c(0,0,0,0), "mm"))

ggsave(here::here("plots/Fig3.png"), p, width = 12, height = 8)
