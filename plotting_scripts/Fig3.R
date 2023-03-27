library(tidyverse)
library(cowplot)
library(gridExtra)
source(here::here("analysis/00-metadata.R"))
source(here::here("plotting_scripts/Fig2.R"))

pairs <- read_csv(paste0(folder_data, "output/pairs.csv"), show_col_types = F)
#pairs_freq <- read_csv(paste0(folder_data, "temp/93a-pairs_freq.csv"), show_col_types = F)
pairs_interaction <- read_csv(paste0(folder_data, "temp/93a-pairs_interaction.csv"), show_col_types = F)
#pairs_accuracy <- read_csv(paste0(folder_data, "temp/91-pairs_accuracy.csv"), show_col_types = F)

flip_winner_species_pairs <- function (pairs) {
    temp_index <- which(pairs$FlipLoser)

    pairs_flipped <- pairs[temp_index, ] %>%
        rename(temp = Isolate1, Isolate1 = Isolate2) %>%
        rename(Isolate2 = temp) %>%
        mutate(FitnessFunction = "-1_-1_-1") %>%
        mutate(F5 = 1-F5, F50 = 1-F50, F95 = 1-F95)

    bind_rows(pairs[-temp_index, ], pairs_flipped) %>%
        arrange(PairID)
}
flip_winner_species_freq <- function (pairs_freq) {
    temp_index <- which(pairs_freq$FlipLoser)

    pairs_freq_flipped <- pairs_freq[temp_index, ] %>%
        rename(temp = Isolate1, Isolate1 = Isolate2) %>%
        rename(Isolate2 = temp) %>%
        #mutate(FitnessFunction = "-1_-1_-1") %>%
        mutate(Isolate1CFUFreqMean = 1-Isolate1CFUFreqMean) %>%
        mutate(Isolate1InitialODFreq = 100-Isolate1InitialODFreq)

    bind_rows(pairs_freq[-temp_index, ], pairs_freq_flipped) %>%
        arrange(Time, PairID, PairFreqID)
}

# Clean up pairs data
pairs <- pairs_interaction %>%
    left_join(pairs_accuracy) %>%
    # Remove no-colony pairs
    unite(col = "Pair", Community, Isolate1, Isolate2, sep = "_", remove = F) %>%
    filter(!(Pair %in% pairs_no_colony)) %>%
    select(-Pair) %>%
    # Remove low-accuracy model pairs
    filter(AccuracyMean > 0.9) %>%
    # Flip isolate 1 if it's a winner
    flip_winner_species_pairs

# pairs <- pairs %>%
#     # Remove no-colony pairs
#     unite(col = "Pair", Community, Isolate1, Isolate2, sep = "_", remove = F) %>%
#     filter(!(Pair %in% pairs_no_colony)) %>%
#     select(-Pair) %>%
#     # Remove low-accuracy model pairs
#     filter(AccuracyMean > 0.9) %>%
#     # Flip isolate 1 if it's a winner
#     flip_winner_species_pairs

pairs_freq <- pairs_freq %>%
    left_join(select(pairs, PairID, FlipLoser)) %>%
    select(PairID, everything()) %>%
    filter(!is.na(PairID)) %>%
    flip_winner_species_freq()


interaction_colors <- c("1-exclusion" = "red",
                        "2-exclusion" = "pink",
                        "3-coexistence" = "blue",
                        "4-coexistence" = "lightblue",
                        "5-inconclusive" = "#808080")


get_interaction_legend <- function (pairs) {
    #panel_colors <- c(interaction_color[1], rep(interaction_color[2], 6), interaction_color[3]) %>% setNames(pairs_interaction_frac$InteractionType)
    #panel_fills <- assign_interaction_color(level = "finer")[-3]
    panel_fills <- interaction_colors
    temp <- pairs %>%
        ggplot() +
        geom_tile(aes(x = Isolate1, y = Isolate2, fill = InteractionType), height = .8, width = .8, alpha = .9, size = 1) +
        scale_fill_manual(values = panel_fills,
                          breaks = names(panel_fills),
                          labels = paste0(pairs_interaction_frac$InteractionType, " (", round(pairs_interaction_frac$Fraction, 3) * 100,"%)")) +
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
        geom_rect(aes(fill = InteractionType),
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
        scale_fill_manual(values = interaction_colors) +
        #scale_fill_manual(values = assign_interaction_color(level = "finer")) +
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

# Append competition outcome to frequencies
pairs_example_freq <- pairs %>%
    filter(!is.na(FitnessFunction)) %>%
    #mutate(InteractionType = factor(InteractionType, c("exclusion", "coexistence"))) %>%
    #mutate(InteractionTypeFiner = factor(InteractionTypeFiner, interaction_type_finer)) %>%
    arrange(InteractionType) %>%
    left_join(pairs_freq, multiple = "all") %>%
    select(PairID, InteractionType, Isolate1InitialODFreq, Time, Isolate1CFUFreqMean, Isolate1CFUFreqSd) %>%
    mutate(Isolate1InitialODFreq = factor(Isolate1InitialODFreq), Time = factor(Time, c("T0", "T8")))

temp_list <- pairs_example_freq %>%
    arrange(InteractionType, PairID) %>%
    group_split(InteractionType, PairID) %>%
    lapply(plot_example_freq)


# pairs_example_freq %>%
#     arrange(InteractionType, PairID) %>%
#     filter(PairID == 3) %>%
#     plot_example_freq()
#



# Legend for fill
pairs_interaction_frac <- pairs %>%
    group_by(InteractionType) %>%
    count(name = "Count") %>%
    ungroup() %>%
    mutate(Fraction = Count / sum(Count), TotalCount = sum(Count)) %>%
    #mutate(InteractionType = factor(InteractionType, interaction_type_finer)) %>%
    arrange(InteractionType)
p_legend_fill <- pairs %>%
    get_interaction_legend()


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





if (FALSE) {

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

























}
