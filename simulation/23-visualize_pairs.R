#' This scripts generate the figures for pool pairs and community pairs
#'
#' N: consumer abundance
#' R: resource abundance

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
df_poolPairs_N_freq <- read_csv(paste0(folder_simulation, "12-aggregated_pairs/df_poolPairs_N_freq.csv"), col_types = cols()) %>%
    mutate(Time = factor(Time, c("init", "end"))) %>%
    mutate(Pair = factor(Pair, paste0("P", 1:1000))) %>%
    mutate(InitialFrequency = factor(InitialFrequency, c(5, 50, 95)))

# Line plot
df_poolPairs_N_freq %>%
    filter(Community == "W0") %>%
    ggplot(aes(x = Time, y = Frequency1, color = InitialFrequency, group = InitialFrequency)) +
    geom_line(linewidth = 1) +
    geom_point(size = 1, shape = 21) +
    scale_color_manual(values = c("F0" = "#8A89C0", "F1" = "#FFCB77")) +
    scale_linetype_manual(values = c("F0" = 1, "F1" = 2)) +
    #scale_y_log10() +
    facet_wrap(~Pair) +
    theme_classic() +
    theme(panel.border = element_rect(color = 1, fill = NA)) +
    guides(alpha = "none", color = "none") +
    labs()

#ggsave(here::here("simulation/plots/23-monoculture_line.png"), p1, width = 5, height = 4)


# Barplot ----
df_network_richness <- read_csv(paste0(folder_simulation, "11-aggregated/df_network_richness.csv"), col_types = cols()) %>%
    mutate(Community = factor(Community, paste0("W", 0:19)))
df_poolPairs_N_outcome <- read_csv(paste0(folder_simulation, "12-aggregated_pairs/df_poolPairs_N_outcome.csv"), col_types = cols()) %>%
    mutate(Community = factor(Community, paste0("W", 0:19)))

p <- df_poolPairs_N_outcome %>%
    group_by(Community, InteractionType, .drop = F) %>%
    summarize(Count = n()) %>%
    arrange(Community, InteractionType) %>%
    mutate(Fraction = Count / sum(Count)) %>%
    left_join(df_network_richness) %>%
    ggplot() +
    geom_col(aes(x = Community, y = Fraction, fill = InteractionType)) +
    geom_text(aes(x = Community, label = Richness), y = 0.9) +
    scale_fill_manual(values = interaction_color) +
    scale_x_discrete(breaks = paste0("W", 0:19)) +
    theme_classic()

ggsave(here::here("simulation/plots/23-poolPairs.png"), p, width = 8, height = 4)

# 2. Community pairs ----
df_withinCommunityPairs_N_freq <- read_csv(paste0(folder_simulation, "12-aggregated_pairs/df_withinCommunityPairs_N_freq.csv"), col_types = cols()) %>%
    mutate(Time = factor(Time, c("init", "end"))) %>%
    mutate(Pair = factor(Pair, paste0("P", 1:1000))) %>%
    mutate(InitialFrequency = factor(InitialFrequency, c(5, 50, 95)))

# Line plot
df_withinCommunityPairs_N_freq %>%
    filter(Community == "W6") %>%
    ggplot(aes(x = Time, y = Frequency1, color = InitialFrequency, group = InitialFrequency)) +
    geom_line(linewidth = 1) +
    geom_point(size = 1, shape = 21) +
    scale_color_manual(values = c("F0" = "#8A89C0", "F1" = "#FFCB77")) +
    scale_linetype_manual(values = c("F0" = 1, "F1" = 2)) +
    #scale_y_log10() +
    facet_wrap(~Pair) +
    theme_classic() +
    theme(panel.border = element_rect(color = 1, fill = NA)) +
    guides(alpha = "none", color = "none") +
    labs()

#ggsave(here::here("simulation/plots/23-lines.png"), p1, width = 5, height = 4)


# Barplot ----
df_communities_richness <- read_csv(paste0(folder_simulation, "11-aggregated/df_communities_richness.csv"), col_types = cols()) %>%
    mutate(Community = factor(Community, paste0("W", 0:19)))
df_withinCommunityPairs_N_outcome <- read_csv(paste0(folder_simulation, "12-aggregated_pairs/df_withinCommunityPairs_N_outcome.csv"), col_types = cols()) %>%
    mutate(Community = factor(Community, paste0("W", 0:19)))

p <- df_withinCommunityPairs_N_outcome %>%
    group_by(Community, InteractionType, .drop = F) %>%
    summarize(Count = n()) %>%
    arrange(Community, InteractionType) %>%
    mutate(Fraction = Count / sum(Count)) %>%
    left_join(df_communities_richness) %>%
    ggplot() +
    geom_col(aes(x = Community, y = Fraction, fill = InteractionType)) +
    geom_text(aes(x = Community, label = Richness), y = 0.9) +
    scale_fill_manual(values = interaction_color) +
    scale_x_discrete(breaks = paste0("W", 0:19)) +
    theme_classic()

ggsave(here::here("simulation/plots/23-withinCommunityPairs.png"), p, width = 8, height = 4)























