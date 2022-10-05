#' This scripts plots the competition outcomes

library(tidyverse)
library(cowplot)
library(gridExtra)
folder_script <- "~/Desktop/lab/emergent-coexistence/output/check/"
folder_main <- "~/Dropbox/lab/emergent-coexistence/data/raw/plate_scan/emergent_coexistence_plate_scan_check/"

communities <- read_csv("~/Dropbox/lab/emergent-coexistence/data/output/communities.csv", col_types = cols())
pairs_freq_ID <- read_csv(paste0(folder_main, "meta/00-pairs_freq_ID.csv"), show_col_types = F)
pairs_ID <- distinct(pairs_freq_ID, Batch, Community, Isolate1, Isolate2) %>% mutate(PairID = 1:n())
pairs_outcome_classifiedT8 <- read_csv(paste0(folder_main, "meta/95-pairs_outcome_classifiedT8.csv"), show_col_types = F) # pairwise outcome data
pairs_outcome_bootstrappedT8 <- read_csv(paste0(folder_main, "meta/95-pairs_outcome_bootstrappedT8.csv"), show_col_types = F)
#pairs_freq <- read_csv("~/Dropbox/lab/emergent-coexistence/data/raw/plate_scan/emergent_coexistence_plate_scan_check/result_pairwise_competition_arranged.csv", show_col_types = F) # human-eye results
pairs_freq <- read_csv(paste0(folder_main, "meta/95-pairs_freq.csv"), show_col_types = F) # frequency data
pairs_freq_boots <- read_csv(paste0(folder_main, "meta/95-pairs_freq_boots.csv"), show_col_types = F)


paint_white_background <- function(x) theme(plot.background = element_rect(color = NA, fill = "white"))

# 1. Configure the column types ----
# 1.1 Sort the communities by size
pairs_freq_ID <- pairs_freq_ID %>% mutate(Community = factor(Community, communities$Community))
pairs_outcome_classifiedT8 <- pairs_outcome_classifiedT8 %>% mutate(Community = factor(Community, communities$Community))
pairs_outcome_bootstrappedT8 <- pairs_outcome_bootstrappedT8 %>% mutate(Community = factor(Community, communities$Community))
pairs_freq <- pairs_freq %>% mutate(Community = factor(Community, communities$Community))
pairs_freq_boots <- pairs_freq_boots %>% mutate(Community = factor(Community, communities$Community))

# 1.2 Add Pair ID
pairs_outcome_classifiedT8 <- pairs_outcome_classifiedT8 %>% left_join(pairs_ID)
pairs_outcome_bootstrappedT8 <- pairs_outcome_bootstrappedT8 %>% left_join(pairs_ID)
pairs_freq <- pairs_freq %>% left_join(pairs_ID)
pairs_freq_boots <- pairs_freq_boots %>% left_join(pairs_ID)


# 2. Plot the competition outcomes ----
assign_interaction_color <- function (level = "simple") {
    if (level == "simple") {
        interaction_type <- c("exclusion", "coexistence")
        interaction_color <- c("#DB7469", "#557BAA")
        names(interaction_color) <- interaction_type
        return(interaction_color)
    }
    if (level == "matrix") {
        interaction_type <- c("exclusion", "coexistence", "exclusion violating rank", "bistability", "neutrality", "self", "undefined")
        interaction_color <- c("#DB7469", "#557BAA", "#8CB369", "#EECF6D", "#8650C4", "black", "grey80")
        names(interaction_color) <- interaction_type
        return(interaction_color)
    }
    if (level == "finer") {
        interaction_type <- c("competitive exclusion", "stable coexistence", "mutual exclusion", "frequency-dependent coexistence", "neutrality", "exclusion violating rank")
        interaction_color <- c("#DB7469", "#557BAA", "#FFBC42", "#B9FAF8", "#8650C4", "#8CB369")
        names(interaction_color) <- interaction_type
        return(interaction_color)
    }
}
fill_names <- assign_interaction_color()

# 2.1 Bootstrapped T8 ----
pairs_outcome_bootstrappedT8 %>%
    ungroup() %>%
    #unite(col = "FitnessFunction", FromRare, FromMedium, FromAbundant, sep = "_") %>%
    group_by(FitnessFunction, InteractionType, InteractionTypeFiner) %>%
    count(name = "Count")

## Exclusion vs. coexistence
p1 <- pairs_outcome_bootstrappedT8 %>%
    mutate(InteractionType = factor(InteractionType, factor(c("coexistence", "exclusion")))) %>%
    group_by(Community, InteractionType, .drop = F) %>%
    count(name = "Count") %>%
    group_by(Community) %>%
    mutate(Fraction = Count / sum(Count), TotalCount = sum(Count)) %>%
    ggplot() +
    geom_col(aes(x = Community, y = Fraction, fill = InteractionType), color = 1) +
    geom_text(aes(x = Community, y = .9, label = paste0("n=",TotalCount))) +
    scale_fill_manual(values = fill_names) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, vjust = .5),
          legend.title = element_blank(),
          legend.position = "right") +
    ggtitle("1000 bootstraps based on random-forest-predicted object probabilities")

# 2.2 Classified T8 ----
pairs_outcome_classifiedT8 %>%
    ungroup() %>%
    #unite(col = "FitnessFunction", FromRare, FromMedium, FromAbundant, sep = "_") %>%
    group_by(FitnessFunction, InteractionType, InteractionTypeFiner) %>%
    count(name = "Count")

## Exclusion vs. coexistence
p2 <- pairs_outcome_classifiedT8 %>%
    mutate(InteractionType = factor(InteractionType, factor(c("coexistence", "exclusion")))) %>%
    group_by(Community, InteractionType, .drop = F) %>%
    count(name = "Count") %>%
    group_by(Community) %>%
    mutate(Fraction = Count / sum(Count), TotalCount = sum(Count)) %>%
    ggplot() +
    geom_col(aes(x = Community, y = Fraction, fill = InteractionType), color = 1) +
    geom_text(aes(x = Community, y = .1, label = paste0("n=",TotalCount))) +
    scale_fill_manual(values = fill_names) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, vjust = .5),
          legend.title = element_blank(),
          legend.position = "right") +
    ggtitle("Random forest classified objects")

p <- plot_grid(p1, p2, nrow = 2, align = "v", axis = "lr")

ggsave(paste0(folder_main, "meta/95a-competition_outcome.png"), p, width = 8, height = 8)


# 3. Waffle plot ----
count_interaction_finer <- function(pairs) {
pairs_interaction_finer <- pairs %>%
    mutate(InteractionType = factor(InteractionType, c("exclusion", "coexistence"))) %>%
    mutate(InteractionTypeFiner = factor(InteractionTypeFiner, c("competitive exclusion", "mutual exclusion", "stable coexistence", "frequency-dependent coexistence", "neutrality"))) %>%
    group_by(InteractionType, InteractionTypeFiner, .drop = 0) %>%
    count(name = "Count") %>% ungroup() %>%
    mutate(Fraction = Count / sum(Count)) %>%
    mutate(Label = str_replace(InteractionTypeFiner, " ", "\n")) %>%
    filter(
        (InteractionType == "coexistence" & InteractionTypeFiner == "stable coexistence") |
            (InteractionType == "coexistence" & InteractionTypeFiner == "frequency-dependent coexistence") |
            (InteractionType == "coexistence" & InteractionTypeFiner == "neutrality") |
            (InteractionType == "exclusion" & InteractionTypeFiner == "competitive exclusion") |
            (InteractionType == "exclusion" & InteractionTypeFiner == "mutual exclusion")
    )


}
get_interaction_legend <- function (pairs) {
    temp <- pairs %>%
        ggplot() +
        geom_tile(aes(x = Isolate1, y = Isolate2, fill = InteractionTypeFiner), height = .8, width = .8, alpha = .9) +
        scale_fill_manual(values = assign_interaction_color(level = "finer"),
                          breaks = c("competitive exclusion", "mutual exclusion", "stable coexistence", "frequency-dependent coexistence", "neutrality"),
                          labels = paste0(pairs_interaction_finer$InteractionTypeFiner, " (", round(pairs_interaction_finer$Fraction, 3) * 100,"%)")) +
        theme(legend.title = element_blank(),
              legend.position = "right",
              legend.spacing.y = unit("2", "mm"),
              legend.text = element_text(size = 12)) +
        guides(fill = guide_legend(byrow = T)) +
        paint_white_background()

    return(get_legend(temp))

}
plot_example_freq <- function(pairs_freq) {
    pairs_freq %>%
        mutate(Isolate1InitialODFreq = factor(Isolate1InitialODFreq, c(95,50,5))) %>%
        ggplot() +
        geom_rect(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, aes(fill = InteractionTypeFiner), alpha = .1) +
        geom_hline(size = .2, yintercept = c(0,1), linetype = 1, color = "grey90") +
        geom_line(size = .4, aes(x = Time, y = Isolate1CFUFreqMean, color = Isolate1InitialODFreq, group = Isolate1InitialODFreq)) +
        geom_point(size = .2, aes(x = Time, y = Isolate1CFUFreqMean, color = Isolate1InitialODFreq, group = Isolate1InitialODFreq)) +
        geom_segment(size = .4, aes(x = Time, xend = Time,
                                    y = Isolate1CFUFreqMean + 2*Isolate1CFUFreqSd,
                                    yend = Isolate1CFUFreqMean - 2*Isolate1CFUFreqSd,
                                    color = Isolate1InitialODFreq)) +
        scale_y_continuous(breaks = c(0, .5, 1), limits = c(-.1, 1.1)) +
        scale_color_manual(values = frequency_color, label = c("95%", "50%", "5%")) +
        scale_fill_manual(values = assign_interaction_color(level = "finer")) +
        theme_bw() +
        theme(panel.spacing = unit(2, "mm"),
              panel.border = element_rect(color = 1, fill = NA, size = .5),
              panel.grid = element_blank(),
              panel.grid.minor.y = element_blank(),
              axis.title = element_blank(), axis.text = element_blank(),
              axis.ticks = element_blank(),
              plot.background = element_blank(),
              plot.title = element_text(size = 5, margin = margin(0,0,0,0)),
              plot.margin = margin(0,0,0,0, "mm")) +
        guides(color = "none", fill = "none") +
        labs(x = "Time", y = "Frequency") +
        ggtitle(unique(pairs_freq$PairID))
}
# 3.1 boostrapped T8 ----
pairs_interaction_finer <- count_interaction_finer(pairs_outcome_bootstrappedT8)
p_legend_fill <- get_interaction_legend(pairs_outcome_bootstrappedT8)

# Append competiton outcome to frequencies
pairs_example_freq <- pairs_outcome_bootstrappedT8 %>%
    mutate(InteractionType = factor(InteractionType, c("exclusion", "coexistence"))) %>%
    mutate(InteractionTypeFiner = ordered(InteractionTypeFiner, c("competitive exclusion", "mutual exclusion", "stable coexistence", "frequency-dependent coexistence", "neutrality"))) %>%
    arrange(InteractionType, InteractionTypeFiner) %>%
    left_join(pairs_freq_boots) %>%
    select(PairID, InteractionType, InteractionTypeFiner, Isolate1InitialODFreq, Time, Isolate1CFUFreqMean, Isolate1CFUFreqSd) %>%
    mutate(Isolate1InitialODFreq = factor(Isolate1InitialODFreq), Time = factor(Time, c("T0", "T8")))

frequency_color <- c( "95"="#292F36", "50"="#9F87AF", "5"="#7D7C7C")
pairs_example_freq %>%
    filter(PairID == 1) %>%
    plot_example_freq()

temp_list <- pairs_example_freq %>%
    arrange(InteractionTypeFiner, PairID) %>%
    group_split(InteractionTypeFiner, PairID) %>%
    lapply(plot_example_freq)

## Grid layout
m <- matrix(c(1:186, rep(NA, 4)), nrow = 10)
p_waffle <- arrangeGrob(grobs = temp_list, layout_matrix = m)
p_waffle <- plot_grid(p_waffle, NULL, rel_widths = c(3, 1), scale = c(.9, 1)) + paint_white_background()

#p_right <- plot_grid(p_legend_fill, p_legend_color, nrow = 2)
ss = .3
p <- ggdraw(p_waffle) +
    draw_plot(p_legend_fill, x = .86, y = .7, width = ss/2, height = ss/2, hjust = 0.5, vjust = .5) +
    #draw_plot(p_legend_color, x = .77, y = .2, width = ss/2, height = ss/2, hjust = 0.5, vjust = .5) +
    theme(panel.background = element_blank(), plot.background = element_rect(color = NA, fill = "white"),
          plot.margin = unit(c(0,0,0,0), "mm"))

ggsave(paste0(folder_main, "meta/95a-waffle_boots.png"), p, width = 12, height = 6)
#ggsave(here::here("plots/Fig3.png"), p, width = 12, height = 4)



# 3.2 no boostrapped T8 ----
pairs_interaction_finer <- count_interaction_finer(pairs_outcome_classifiedT8)
p_legend_fill <- get_interaction_legend(pairs_outcome_classifiedT8)

# Append competiton outcome to frequencies
pairs_example_freq <- pairs_outcome_classifiedT8 %>%
    mutate(InteractionType = factor(InteractionType, c("exclusion", "coexistence"))) %>%
    mutate(InteractionTypeFiner = ordered(InteractionTypeFiner, c("competitive exclusion", "mutual exclusion", "stable coexistence", "frequency-dependent coexistence", "neutrality"))) %>%
    arrange(InteractionType, InteractionTypeFiner) %>%
    left_join(pairs_freq) %>%
    select(PairID, InteractionType, InteractionTypeFiner, Isolate1InitialODFreq, Time, Isolate1CFUFreqMean, Isolate1CFUFreqSd) %>%
    mutate(Isolate1InitialODFreq = factor(Isolate1InitialODFreq), Time = factor(Time, c("T0", "T8")))

pairs_example_freq %>%
    filter(PairID == 1) %>%
    plot_example_freq()

temp_list <- pairs_example_freq %>%
    arrange(InteractionTypeFiner, PairID) %>%
    group_split(InteractionTypeFiner, PairID) %>%
    lapply(plot_example_freq)

## Grid layout
m <- matrix(c(1:186, rep(NA, 4)), nrow = 10)
p_waffle <- arrangeGrob(grobs = temp_list, layout_matrix = m)
p_waffle <- plot_grid(p_waffle, NULL, rel_widths = c(3, 1), scale = c(.9, 1)) + paint_white_background()

#p_right <- plot_grid(p_legend_fill, p_legend_color, nrow = 2)
ss = .3
p <- ggdraw(p_waffle) +
    draw_plot(p_legend_fill, x = .86, y = .7, width = ss/2, height = ss/2, hjust = 0.5, vjust = .5) +
    #draw_plot(p_legend_color, x = .77, y = .2, width = ss/2, height = ss/2, hjust = 0.5, vjust = .5) +
    theme(panel.background = element_blank(), plot.background = element_rect(color = NA, fill = "white"),
          plot.margin = unit(c(0,0,0,0), "mm"))

ggsave(paste0(folder_main, "meta/95a-waffle_classification.png"), p, width = 12, height = 6)




#
#
#
#
#
#
# pairs_outcome_classifiedT8 %>%
#     group_by(Community, InteractionType) %>%
#     count(name = "Count") %>%
#     group_by(Community) %>% mutate(Fraction = Count / sum(Count)) %>% ungroup() %>%
#     mutate(Community = factor(Community, communities$Community)) %>%
#     arrange(Community) %>%
#     left_join(communities, by = "Community") %>%
#     mutate(CommunityLabel = factor(CommunityLabel)) %>%
#     ggplot() +
#     geom_col(aes(x = CommunityLabel, fill = InteractionType, y = Fraction), color = 1, width = .8, size = .5) +
#     geom_text(data = communities, aes(x = CommunityLabel, y = .1, label = paste0("n=", CommunityPairSize)), vjust = -.5, size = 3) +
#     scale_fill_manual(values = assign_interaction_color()) +
#     scale_y_continuous(breaks = c(0,.5,1), limit = c(0, 1), expand = c(0,0)) +
#     theme_classic() +
#     theme(legend.text = element_text(size = 12),
#           axis.text = element_text(color = 1, size = 12),
#           axis.title = element_text(color = 1, size = 12),
#           legend.title = element_blank(),
#           legend.position = "top") +
#     guides(fill = guide_legend(reverse = T)) +
#     labs(x = "Community", y = "Fraction")











