#' This scripts generate the figures for pool pairs and community pairs


library(tidyverse)
library(cowplot)
source(here::here("analysis/00-metadata.R"))

# 0. parameters ----
input_parameters <- read_csv(here::here("simulation/01-input_parameters.csv"), col_types = cols())
input_poolPairs <- read_csv(here::here("simulation/03a-input_poolPairs.csv"), col_types = cols())
input_withinCommunityPairs <- read_csv(here::here("simulation/03b-input_withinCommunityPairs.csv"), col_types = cols())

# Generate family-species and class-resource table for matching
sa <- input_parameters$sa[1]
ma <- input_parameters$ma[1]
sal <- tibble(Family = paste0("F", c(rep(0, sa), rep(1, sa))), Species = paste0("S", 0:(sa * 2 - 1)))
mal <- tibble(Class = paste0("T", c(rep(0, ma), rep(1, ma))), Resource = paste0("R", 0:(ma * 2 - 1)))

# 1. Pool pairs ----
poolPairs_N_freq <- read_csv(paste0(folder_simulation, "12-aggregated_pairs/poolPairs_N_freq.csv"), col_types = cols()) %>%
    mutate(Time = factor(Time, c("init", "end"))) %>%
    mutate(Pair = factor(Pair, paste0("P", 1:1000))) %>%
    mutate(InitialFrequency = factor(InitialFrequency, c(5, 50, 95)))

# Barplot
monocultureSets_richness <- read_csv(paste0(folder_simulation, "11-aggregated/monocultureSets_richness.csv"), col_types = cols()) %>%
    mutate(Community = factor(Community, paste0("W", 0:19))) %>%
    mutate(PairSize = choose(Richness, 2))
poolPairs_N_outcome <- read_csv(paste0(folder_simulation, "12-aggregated_pairs/poolPairs_N_outcome.csv"), col_types = cols()) %>%
    mutate(Community = factor(Community, paste0("W", 0:19)))

p1 <- poolPairs_N_outcome %>%
    group_by(Community, InteractionType, .drop = F) %>%
    summarize(Count = n()) %>%
    arrange(Community, InteractionType) %>%
    mutate(Fraction = Count / sum(Count)) %>%
    left_join(monocultureSets_richness) %>%
    ggplot() +
    geom_col(aes(x = Community, y = Fraction, fill = InteractionType), color = 1) +
    annotate("text", x = 1:20, y = 1.15, label = monocultureSets_richness$Richness, size = 4) +
    annotate("text", x = 21, y = 1.15, label = c("n. of species"), size = 4, hjust = 0) +
    annotate("segment", x = .5, xend = 25, y = 1.1, yend = 1.1, color = "black") +
    annotate("text", x = 21, y = 1.05, label = c("n. of tested pairs"), size = 4, hjust = 0) +
    annotate("text", x = 1:20, y = 1.05, label = monocultureSets_richness$PairSize, size = 4) +
    scale_fill_manual(values = interaction_color) +
    scale_x_discrete(breaks = paste0("W", 0:19)) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 1.2)) +
    coord_cartesian(xlim = c(0.5, 20.5), ylim = c(0, 1), clip = "off") +
    theme_classic() +
    theme(plot.margin = unit(c(2,.5,.5,.5), "cm"))  +
    guides(color = "none", fill = guide_legend(title = "")) +
    labs(x = "Set")

#ggsave(here::here("simulation/plots/23-poolPairs.png"), p, width = 8, height = 4)

# 2. Community pairs ----
withinCommunityPairs_N_freq <- read_csv(paste0(folder_simulation, "12-aggregated_pairs/withinCommunityPairs_N_freq.csv"), col_types = cols()) %>%
    mutate(Time = factor(Time, c("init", "end"))) %>%
    mutate(Pair = factor(Pair, paste0("P", 1:1000))) %>%
    mutate(InitialFrequency = factor(InitialFrequency, c(5, 50, 95)))

# # Line plot
# withinCommunityPairs_N_freq %>%
#     filter(Community == "W6") %>%
#     ggplot(aes(x = Time, y = Frequency1, color = InitialFrequency, group = InitialFrequency)) +
#     geom_line(linewidth = 1) +
#     geom_point(size = 1, shape = 21) +
#     scale_color_manual(values = c("F0" = "#8A89C0", "F1" = "#FFCB77")) +
#     scale_linetype_manual(values = c("F0" = 1, "F1" = 2)) +
#     #scale_y_log10() +
#     facet_wrap(~Pair) +
#     theme_classic() +
#     theme(panel.border = element_rect(color = 1, fill = NA)) +
#     guides(alpha = "none", color = "none") +
#     labs()


# Barplot ----
communities_richness <- read_csv(paste0(folder_simulation, "11-aggregated/communities_richness.csv"), col_types = cols()) %>%
    mutate(Community = factor(Community, paste0("W", 0:19))) %>%
    mutate(PairSize = choose(Richness, 2))
withinCommunityPairs_N_outcome <- read_csv(paste0(folder_simulation, "12-aggregated_pairs/withinCommunityPairs_N_outcome.csv"), col_types = cols()) %>%
    mutate(Community = factor(Community, paste0("W", 0:19)))

p2 <- withinCommunityPairs_N_outcome %>%
    group_by(Community, InteractionType, .drop = F) %>%
    summarize(Count = n()) %>%
    arrange(Community, InteractionType) %>%
    mutate(Fraction = Count / sum(Count)) %>%
    left_join(communities_richness) %>%
    ggplot() +
    geom_col(aes(x = Community, y = Fraction, fill = InteractionType), color = 1) +
    #geom_text(aes(x = Community, label = Richness), y = 0.9) +
    annotate("text", x = 1:20, y = 1.15, label = communities_richness$Richness, size = 4) +
    annotate("text", x = 21, y = 1.15, label = c("n. of species"), size = 4, hjust = 0) +
    annotate("segment", x = .5, xend = 25, y = 1.1, yend = 1.1, color = "black") +
    annotate("text", x = 21, y = 1.05, label = c("n. of tested pairs"), size = 4, hjust = 0) +
    annotate("text", x = 1:20, y = 1.05, label = communities_richness$PairSize, size = 4) +
    scale_fill_manual(values = interaction_color) +
    scale_x_discrete(breaks = paste0("W", 0:19)) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 1.2)) +
    coord_cartesian(xlim = c(0.5, 20.5), ylim = c(0, 1), clip = "off") +
    theme_classic() +
    theme(plot.margin = unit(c(2,.5,.5,.5), "cm"))  +
    guides(color = "none", fill = guide_legend(title = "")) +
    labs(x = "Community")


#ggsave(here::here("simulation/plots/23-withinCommunityPairs.png"), p, width = 8, height = 4)

p <- plot_grid(p1, p2, nrow = 2, align = "v", axis = "rl", labels = c("Pool pairs", "Within-community pairs"), hjust = 0, label_x = 0.01) + paint_white_background()
ggsave(here::here("simulation/plots/23-pairs.png"), p, width = 8, height = 8)





















