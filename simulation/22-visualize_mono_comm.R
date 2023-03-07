#' This scripts generate the figures for monocultures and communities
#'
#' N: consumer abundance
#' R: resource abundance

library(tidyverse)
library(cowplot)
source(here::here("analysis/00-metadata.R"))

# 0. parameters ----
input_parameters <- read_csv(here::here("simulation/01-input_parameters.csv"), col_types = cols())
input_monocultures <- read_csv(here::here("simulation/02a-input_monocultures.csv"), col_types = cols())
input_communities <- read_csv(here::here("simulation/02b-input_communities.csv"), col_types = cols())

# Generate family-species and class-resource table for matching
sa <- input_parameters$sa[1]
ma <- input_parameters$ma[1]
sal <- tibble(Family = paste0("F", c(rep(0, sa), rep(1, sa))), Species = paste0("S", 0:(sa * 2 - 1)))
mal <- tibble(Class = paste0("T", c(rep(0, ma), rep(1, ma))), Resource = paste0("R", 0:(ma * 2 - 1)))


# 1. Communities ----
communities_abundance <- read_csv(paste0(folder_simulation, "11-aggregated/communities_abundance.csv"), col_types = cols()) %>%
    mutate(Community = ordered(Community, paste0("W", 0:(input_communities$n_wells[1]-1)))) %>%
    mutate(Time = ordered(Time, c("init", paste0("T", 1:20), "end"))) %>%
    arrange(Community, Time)

# Line plot
p1 <- communities_abundance %>%
    filter(Abundance != 0) %>%
    ggplot(aes(x = Time, y = Abundance, color = Family, group = Species)) +
    geom_line(linewidth = .2) +
    geom_point(size = 1, shape = 21) +
    scale_color_manual(values = c("F0" = "#8A89C0", "F1" = "#FFCB77")) +
    scale_linetype_manual(values = c("F0" = 1, "F1" = 2)) +
    facet_wrap(.~Community, ncol = 5) +
    theme_classic() +
    theme(panel.border = element_rect(color = 1, fill = NA)) +
    guides(alpha = "none") +
    labs()
ggsave(here::here("simulation/plots/22-community_line.png"), p1, width = 12, height = 10)

# Barplot over time
p2 <- communities_abundance %>%
    filter(Abundance != 0) %>%
    ggplot(aes(x = Time, y = Abundance, fill = Family, color = Species)) +
    geom_col() +
    scale_fill_manual(values = c("F0" = "#8A89C0", "F1" = "#FFCB77")) +
    scale_color_manual(values = rep("black", length(sal$Species))) +
    facet_wrap(.~Community, ncol = 5) +
    theme_classic() +
    guides(color = "none") +
    labs()

ggsave(here::here("simulation/plots/22-community_bar.png"), p2, width = 12, height = 10)

# Barplot over time, standard
p3 <- communities_abundance %>%
    filter(Abundance != 0) %>%
    ggplot(aes(x = Time, y = Abundance, fill = Family, color = Species)) +
    geom_col(position = "fill") +
    scale_fill_manual(values = c("F0" = "#8A89C0", "F1" = "#FFCB77")) +
    scale_color_manual(values = rep("black", length(sal$Species))) +
    facet_wrap(.~Community, ncol = 5) +
    theme_classic() +
    guides(color = "none") +
    labs()
ggsave(here::here("simulation/plots/22-community_bar_standard.png"), p3, width = 12, height = 10)



# Barplot final time point
communities_abundance_abundant <- communities_abundance %>%
    filter(Time == max(Time)) %>%
    filter(Abundance != 0) %>%
    group_by(Community) %>%
    mutate(TotalAbundance = sum(Abundance)) %>%
    filter(Abundance > 0.01 * sum(Abundance))

communities_abundance_summmary <- communities_abundance_abundant %>%
    summarize(Richness = n())
p4 <- communities_abundance_abundant %>%
    filter(Time == max(Time)) %>%
    filter(Abundance != 0) %>%
    group_by(Community) %>%
    filter(Abundance > 0.01 * sum(Abundance)) %>%
    mutate(RelativeAbundance = Abundance / sum(Abundance)) %>%
    ggplot() +
    geom_col(aes(x = Community, y = RelativeAbundance, fill = Family), color = 1) +
    geom_text(data = communities_abundance_summmary, aes(x = Community, label = Richness), y = 1.1) +
    scale_fill_manual(values = c("F0" = "#8A89C0", "F1" = "#FFCB77")) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 1.1)) +
    coord_cartesian(ylim = c(0, 1), clip = "off") +
    theme_classic() +
    theme(plot.margin = unit(c(1,.5,.5,.5), "cm")) +
    guides(color = "none") +
    labs()
ggsave(here::here("simulation/plots/22-community_bar_final.png"), p4, width = 9, height = 3)




# 2. Monocultures ----
monocultures_abundance <- read_csv(paste0(folder_simulation, "11-aggregated/monocultures_abundance.csv"), col_types = cols()) %>%
    mutate(Well = ordered(Well, paste0("W", 0:(input_monocultures$n_wells[1]-1)))) %>%
    mutate(Time = ordered(Time, c("init", paste0("T", 1:20)))) %>%
    arrange(Well, Time)

# Line plot
p1 <- monocultures_abundance %>%
    ggplot(aes(x = Time, y = Abundance, color = Family, group = Species)) +
    geom_line(linewidth = .2) +
    geom_point(size = 1, shape = 21) +
    scale_color_manual(values = c("F0" = "#8A89C0", "F1" = "#FFCB77")) +
    scale_linetype_manual(values = c("F0" = 1, "F1" = 2)) +
    theme_classic() +
    theme(panel.border = element_rect(color = 1, fill = NA)) +
    guides(alpha = "none") +
    labs()
ggsave(here::here("simulation/plots/22-monoculture_line.png"), p1, width = 5, height = 4)






