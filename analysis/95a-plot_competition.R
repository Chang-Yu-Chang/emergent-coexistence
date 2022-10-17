#' This scripts plots the competition outcomes

library(tidyverse)
library(cowplot)
library(gridExtra)
folder_script <- "~/Desktop/lab/emergent-coexistence/analysis/"
folder_main <- "~/Dropbox/lab/emergent-coexistence/plate_scan_pipeline/"

communities <- read_csv("~/Dropbox/lab/emergent-coexistence/data/output/communities.csv", col_types = cols())
isolates <- read_csv(paste0(folder_main, "meta/97-isolates.csv"), show_col_types = F)
pairs_freq_ID <- read_csv(paste0(folder_main, "meta/00-pairs_freq_ID.csv"), show_col_types = F)
pairs_ID <- distinct(pairs_freq_ID, Batch, Community, Isolate1, Isolate2) %>% mutate(PairID = 1:n())
pairs_interaction <- read_csv(paste0(folder_main, "meta/95-pairs_interaction.csv"), show_col_types = F)
pairs_freq <- read_csv(paste0(folder_main, "meta/95-pairs_freq.csv"), show_col_types = F) # frequency data
accuracy <- read_csv(paste0(folder_main, "meta/92-accuracy.csv"), show_col_types = F)

#
paint_white_background <- function(x) theme(plot.background = element_rect(color = NA, fill = "white"))

# 1. Configure the column types ----
# 1.1 Add Pair ID
pairs_interaction <- pairs_interaction %>% left_join(pairs_ID)
pairs_freq <- pairs_freq %>% left_join(pairs_ID)

# 1.2 Sort the communities by size
pairs_freq_ID <- pairs_freq_ID %>% mutate(Community = factor(Community, communities$Community))
pairs_interaction <- pairs_interaction %>% mutate(Community = factor(Community, communities$Community))
pairs_freq <- pairs_freq %>% mutate(Community = factor(Community, communities$Community))

# 1.3 append random forest model accuracy values
pairs_accuracy <- accuracy %>%
    select(image_name_pair, Accuracy) %>%
    left_join(pairs_freq_ID) %>%
    group_by(Community, Isolate1, Isolate2) %>%
    arrange(Community, Isolate1, Isolate2) %>%
    summarize(Accuracy = mean(Accuracy), Count = n())
pairs_interaction %>%
    left_join(pairs_accuracy) %>%
    filter(Accuracy < 0.9)

## Remove low accuracy pairs
pairs_interaction <- pairs_interaction %>%
    left_join(pairs_accuracy) %>%
    filter(Accuracy > 0.9)

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
        interaction_type <- c("competitive exclusion", "stable coexistence",
                              "mutual exclusion", "frequency-dependent coexistence",
                              "coexistence at 5%", "coexistence at 95%",
                              "2-freq neutrality", "3-freq neutrality")
        #interaction_type <- c("competitive exclusion", "stable coexistence", "mutual exclusion", "frequency-dependent coexistence", "neutrality", "exclusion violating rank")
        interaction_color <- c("#DB7469", "#557BAA",
                               "#FFBC42", "#B9FAF8",
                               "lightslateblue", "cyan",
                               "#8650C4", "maroon")
        names(interaction_color) <- interaction_type
        return(interaction_color)
    }
}
fill_names <- assign_interaction_color()

# Bootstrapped T8
## Exclusion vs. coexistence
p1 <- pairs_interaction %>%
    filter(!is.na(FitnessFunction)) %>%
    mutate(InteractionType = factor(InteractionType, factor(c("coexistence", "exclusion")))) %>%
    mutate(Community = factor(Community, communities$Community)) %>%
    group_by(Community, InteractionType, .drop = F) %>%
    count(name = "Count") %>%
    group_by(Community) %>%
    mutate(Fraction = Count / sum(Count), TotalCount = sum(Count)) %>%
    ggplot() +
    geom_col(aes(x = Community, y = Fraction, fill = InteractionType), color = 1) +
    geom_text(aes(x = Community, y = .9, label = paste0("n=", TotalCount))) +
    scale_fill_manual(values = c(fill_names, "NA" = "grey")) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, vjust = .5),
          legend.title = element_blank(),
          legend.position = "right") +
    ggtitle("")

p2 <- pairs_interaction %>%
    filter(!is.na(FitnessFunction)) %>%
    #mutate(InteractionType = factor(InteractionType, factor(c("coexistence", "exclusion")))) %>%
    mutate(Community = factor(Community, communities$Community)) %>%
    group_by(Community, InteractionTypeFiner, .drop = F) %>%
    count(name = "Count") %>%
    group_by(Community) %>%
    mutate(Fraction = Count / sum(Count), TotalCount = sum(Count)) %>%
    ggplot() +
    geom_col(aes(x = Community, y = Fraction, fill = InteractionTypeFiner), color = 1) +
    #geom_text(aes(x = Community, y = .9, label = paste0("n=", TotalCount))) +
    scale_fill_manual(values = c(assign_interaction_color("finer"), "NA" = "grey")) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, vjust = .5),
          legend.title = element_blank(),
          legend.position = "right") +
    ggtitle("")

p <- plot_grid(p1, p2, nrow = 2, axis = "rl", align = "v")
ggsave(paste0(folder_main, "meta/95a-competition_outcome.png"), p, width = 8, height = 8)


# 3. Waffle plot ----
interaction_type_finer <- c("competitive exclusion", "stable coexistence",
                            "mutual exclusion", "frequency-dependent coexistence",
                            "coexistence at 5%", "coexistence at 95%",
                            "2-freq neutrality", "3-freq neutrality")

make_interaction_type <- function () {
    #' This function generates the fitness function table.
    #' There are a total of 27 possibilities
    interaction_type <- tibble(
        FromRare = rep(c(1, -1, 0), each = 9),
        FromMedium = rep(rep(c(1, -1, 0), each = 3), 3),
        FromAbundant = rep(c(1, -1, 0), 9),
        InteractionType = NA,
        InteractionTypeFiner = NA
    )
    ## Assign interaction types to combinations of frequency changes signs
    interaction_type$InteractionType[c(1,10,13,14)] <- "exclusion"
    interaction_type$InteractionType[c(2:6,8,11,20,23, 9,18,21,24,25,26,27)] <- "coexistence"

    ## Assign finer interaction types to combinations of frequency changes signs
    interaction_type$InteractionTypeFiner[c(1,14)] <- "competitive exclusion"
    interaction_type$InteractionTypeFiner[c(10,13)] <- "mutual exclusion"
    interaction_type$InteractionTypeFiner[c(2,5,8)] <- "stable coexistence"
    interaction_type$InteractionTypeFiner[c(4,6,11,20)] <- "frequency-dependent coexistence"
    interaction_type$InteractionTypeFiner[c(3)] <- "coexistence at 95%"
    interaction_type$InteractionTypeFiner[c(23)] <- "coexistence at 5%"
    interaction_type$InteractionTypeFiner[c(9,18,21,24:26)] <- "neutrality"
    interaction_type$InteractionTypeFiner[c(9,18,21,24:26)] <- "2-freq neutrality"
    interaction_type$InteractionTypeFiner[c(27)] <- "3-freq neutrality"
    interaction_type <- interaction_type %>%  mutate(FitnessFunction = paste(FromRare, FromMedium, FromAbundant, sep = "_"))
}
interaction_type_table <- make_interaction_type() %>%
    filter(!is.na(InteractionType)) %>%
    distinct(InteractionType, InteractionTypeFiner) %>%
    mutate(InteractionTypeFiner = factor(InteractionTypeFiner, interaction_type_finer)) %>%
    arrange(InteractionTypeFiner)
count_interaction_finer <- function(pairs_interaction) {
    pairs_interaction %>%
        mutate(InteractionType = factor(InteractionType, c("exclusion", "coexistence"))) %>%
        mutate(InteractionTypeFiner = factor(InteractionTypeFiner, interaction_type_finer)) %>%
        group_by(InteractionType, InteractionTypeFiner, .drop = 0) %>%
        count(name = "Count") %>% ungroup() %>%
        mutate(Fraction = Count / sum(Count)) %>%
        mutate(Label = str_replace(InteractionTypeFiner, " ", "\n")) %>%
        right_join(interaction_type_table) %>%
        arrange(InteractionTypeFiner)
}
get_interaction_legend <- function (pairs) {
    temp <- pairs %>%
        ggplot() +
        geom_tile(aes(x = Isolate1, y = Isolate2, fill = InteractionTypeFiner), height = .8, width = .8, alpha = .9) +
        scale_fill_manual(values = assign_interaction_color(level = "finer"),
                          breaks = names(assign_interaction_color(level = "finer")),
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

#
pairs_interaction_finer <- count_interaction_finer(pairs_interaction)
p_legend_fill <- get_interaction_legend(pairs_interaction)

# Append competition outcome to frequencies
pairs_example_freq <- pairs_interaction %>%
    mutate(InteractionType = factor(InteractionType, c("exclusion", "coexistence"))) %>%
    mutate(InteractionTypeFiner = factor(InteractionTypeFiner, interaction_type_finer)) %>%
    arrange(InteractionType, InteractionTypeFiner) %>%
    left_join(pairs_freq) %>%
    select(PairID, InteractionType, InteractionTypeFiner, Isolate1InitialODFreq, Time, Isolate1CFUFreqMean, Isolate1CFUFreqSd) %>%
    mutate(Isolate1InitialODFreq = factor(Isolate1InitialODFreq), Time = factor(Time, c("T0", "T8")))

frequency_color <- c( "95"="#292F36", "50"="#9F87AF", "5"="#7D7C7C")

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



