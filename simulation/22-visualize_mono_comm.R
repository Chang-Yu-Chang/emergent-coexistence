#' This scripts generate the figures for monocultures and communities
#'
#' N: consumer abundance
#' R: resource abundance

library(tidyverse)
library(cowplot)
source(here::here("analysis/00-metadata.R"))
source(here::here("simulation/01-generate_input.R"))

# 0. parameters ----
input_parameters <- read_csv(here::here("simulation/01-input_parameters.csv"), col_types = cols())
input_monocultures <- read_csv(here::here("simulation/02a-input_monocultures.csv"), col_types = cols())
input_communities <- read_csv(here::here("simulation/02b-input_communities.csv"), col_types = cols())
input_communitiesWithoutCrossfeeding <- read_csv(here::here("simulation/02c-input_communitiesWithoutCrossfeeding.csv"), col_types = cols())


# 1. Monocultures ----
monocultures_abundance <- read_csv(paste0(folder_simulation, "aggregated/12-monocultures_abundance.csv"), col_types = cols()) %>%
    mutate(Well = ordered(Well, paste0("W", 0:(input_monocultures$n_wells[1]-1)))) %>%
    mutate(Time = ordered(Time, c("init", paste0("T", 1:20)))) %>%
    arrange(Well, Time)

# Line plot
p1 <- monocultures_abundance %>%
    filter(Abundance > 0) %>%
    #mutate(UniqueWell = paste0(Well, Species)) %>%
    ggplot(aes(x = Time, y = Abundance, color = Family, group = Species)) +
    geom_line(linewidth = .2) +
    geom_point(size = 1, shape = 21) +
    #scale_color_manual(values = c("F0" = "#8A89C0", "F1" = "#FFCB77"), labels = c("F0" = "fermenter", "F1" = "repirator")) +
    #scale_linetype_manual(values = c("F0" = 1, "F1" = 2)) +
    theme_classic() +
    theme(panel.border = element_rect(color = 1, fill = NA)) +
    guides(alpha = "none", color = guide_legend(title = "")) +
    labs()


monocultures_abundance_Tinit <- monocultures_abundance %>%
    filter(Time == min(Time)) %>%
    filter(Abundance > 0) %>%
    select(Well, Family, Species)

monocultures_abundance_Tend <- monocultures_abundance %>%
    filter(Time == max(Time)) %>%
    filter(Abundance > 0) %>%
    select(Well, Family, Species, Abundance)

p2 <- monocultures_abundance_Tinit %>%
    left_join(monocultures_abundance_Tend) %>%
    replace_na(list(Abundance = 0)) %>%
    mutate(Alive = ifelse(Abundance > 0, "alive", "nah")) %>%
    group_by(Family, Alive) %>%
    summarize(Count = n()) %>%
    mutate(TotalCount = sum(Count),
           Fraction = Count / TotalCount) %>%
    ggplot(aes(x = Family, y = Fraction, fill = Alive)) +
    geom_col(color = 1) +
    geom_text(aes(x = Family, label = paste0("n=", TotalCount)), y = Inf, vjust = 1) +
    #scale_color_manual(values = c("F0" = "#8A89C0", "F1" = "#FFCB77"), labels = c("F0" = "fermenter", "F1" = "repirator")) +
    #scale_linetype_manual(values = c("F0" = 1, "F1" = 2)) +
    theme_classic() +
    theme(panel.border = element_rect(color = 1, fill = NA)) +
    guides(alpha = "none", color = guide_legend(title = "")) +
    labs()
p <- plot_grid(p1, p2, scale = 0.9) + paint_white_background()
ggsave(here::here("simulation/plots/22-monoculture-01-line.png"), p, width = 10, height = 4)


# 2. Communities ----
communities_abundance <- read_csv(paste0(folder_simulation, "aggregated/12-communities_abundance.csv"), col_types = cols()) %>%
    mutate(Community = ordered(Community, paste0("W", 0:(input_communities$n_wells[1]-1)))) %>%
    mutate(Time = ordered(Time, c("init", paste0("T", 1:20), "end"))) %>%
    arrange(Community, Time)

# Barplot over time, absolute
p <- communities_abundance %>%
    filter(Abundance != 0) %>%
    filter(Time != "end") %>%
    ggplot(aes(x = Time, y = Abundance, fill = Family, color = Species)) +
    geom_col() +
    #scale_fill_manual(values = c("F0" = "#8A89C0", "F1" = "#FFCB77"), labels = c("F0" = "fermenter", "F1" = "repirator")) +
    scale_color_manual(values = rep("black", length(sal$Species))) +
    facet_wrap(.~Community, ncol = 5) +
    theme_classic() +
    guides(color = "none", fill = guide_legend(title = "")) +
    labs()

ggsave(here::here("simulation/plots/22-communities-01-bar_abs.png"), p, width = 12, height = 10)

# Barplot over time, standard
temp <- tibble(Community = paste0("W", 0:19),
       Community_new = factor(paste0("community ", 1:20), paste0("community ", 1:20)))
p <- communities_abundance %>%
    filter(Abundance != 0) %>%
    filter(Time != "end") %>%
    left_join(temp) %>%
    ggplot(aes(x = Time, y = Abundance, fill = Family, color = Species)) +
    geom_col(position = "fill", linewidth = .2) +
    #scale_fill_manual(values = c("F0" = "#8A89C0", "F1" = "#FFCB77"), labels = c("F0" = "fermenter", "F1" = "repirator")) +
    scale_color_manual(values = rep("black", length(sal$Species))) +
    scale_x_discrete(breaks = c("init", paste0("T", seq(5,20, 5))), labels = c(0, seq(5, 20, 5))) +
    #scale_x_discrete(breaks = c("init", paste0("T", seq(5,20, 5))), labels = c("T0", paste0("T", seq(5,20, 5)))) +
    scale_y_continuous(breaks = seq(0,1,0.2), expand = c(0,0)) +
    facet_wrap(.~Community_new, ncol = 4) +
    theme_classic() +
    theme(
       strip.background = element_rect(color = NA, fill = NA)
    ) +
    guides(color = "none", fill = guide_legend(title = "Family")) +
    labs(x = "transfer", y = "relative abundance")

ggsave(here::here("simulation/plots/22-communities-02-bar_fraction.png"), p, width = 10, height = 10)


# Barplot final time point
communities_abundance_abundant <- communities_abundance %>%
    filter(Time == max(Time)) %>%
    filter(Abundance != 0) %>%
    group_by(Community) %>%
    mutate(TotalAbundance = sum(Abundance)) %>%
    filter(Abundance > 0.01 * sum(Abundance))
temp <- tibble(Community = paste0("W", 0:19),
               Community_new = factor(1:20, 1:20))

communities_richness <- read_csv(paste0(folder_simulation, "aggregated/12-communities_richness.csv"), col_types = cols())
communities_abundance_richness <- communities_abundance_abundant %>%
    filter(Time == max(Time)) %>%
    group_by(Community, .drop = F) %>%
    summarize(Richness = n()) %>%
    left_join(temp)


p <- communities_abundance_abundant %>%
    filter(Time == max(Time)) %>%
    group_by(Community) %>%
    mutate(RelativeAbundance = Abundance / sum(Abundance)) %>%
    left_join(temp) %>%
    ggplot() +
    geom_col(aes(x = Community_new, y = RelativeAbundance, fill = Family), color = 1) +
    geom_text(data = communities_abundance_richness, aes(x = Community_new, label = Richness), y = 1.1) +
    annotate("text", x = 1:20, y = 1.15, label = communities_richness$Richness, size = 4) +
    annotate("text", x = 21, y = 1.1, label = c("n. of species"), size = 4, hjust = 0) +
    #scale_fill_manual(values = c("F0" = "#8A89C0", "F1" = "#FFCB77"), labels = c("F0" = "fermenter", "F1" = "repirator")) +
    scale_y_continuous(breaks = seq(0,1, 0.2), expand = c(0,0), limits = c(0,1.1)) +
    coord_cartesian(xlim = c(0.5, 20.5), ylim = c(0, 1), clip = "off") +
    theme_classic() +
    theme(plot.margin = unit(c(1,.5,.5,.5), "cm")) +
    guides(color = "none", fill = guide_legend(title = "Family")) +
    labs(x = "community", y = "relative abundance")

ggsave(here::here("simulation/plots/22-communities-03-bar_final.png"), p, width = 9, height = 3)

# 3. Communities without crossfeeding -----
communitiesWithoutCrossfeeding_abundance <- read_csv(paste0(folder_simulation, "aggregated/12-communitiesWithoutCrossfeeding_abundance.csv"), col_types = cols()) %>%
    mutate(Community = ordered(Community, paste0("W", 0:(input_communitiesWithoutCrossfeeding$n_wells[1]-1)))) %>%
    mutate(Time = ordered(Time, c("init", paste0("T", 1:20), "end"))) %>%
    arrange(Community, Time)

# Barplot over time, absolute
p <- communitiesWithoutCrossfeeding_abundance %>%
    filter(Abundance != 0) %>%
    filter(Time != "end") %>%
    ggplot(aes(x = Time, y = Abundance, fill = Family, color = Species)) +
    geom_col() +
    #scale_fill_manual(values = c("F0" = "#8A89C0", "F1" = "#FFCB77"), labels = c("F0" = "fermenter", "F1" = "repirator")) +
    scale_color_manual(values = rep("black", length(sal$Species))) +
    facet_wrap(.~Community, ncol = 5) +
    theme_classic() +
    guides(color = "none", fill = guide_legend(title = "")) +
    labs()

ggsave(here::here("simulation/plots/22-communitiesWithoutCrossfeeding-01-bar_abs.png"), p, width = 12, height = 10)

# Barplot over time, fraction
p <- communitiesWithoutCrossfeeding_abundance %>%
    filter(Abundance != 0) %>%
    filter(Time != "end") %>%
    ggplot(aes(x = Time, y = Abundance, fill = Family, color = Species)) +
    geom_col(position = "fill") +
    #scale_fill_manual(values = c("F0" = "#8A89C0", "F1" = "#FFCB77"), labels = c("F0" = "fermenter", "F1" = "repirator")) +
    scale_color_manual(values = rep("black", length(sal$Species))) +
    facet_wrap(.~Community, ncol = 5) +
    theme_classic() +
    guides(color = "none", fill = guide_legend(title = "")) +
    labs()
ggsave(here::here("simulation/plots/22-communitiesWithoutCrossfeeding-02-bar_fraction.png"), p, width = 12, height = 10)


# Barplot final time point
communitiesWithoutCrossfeeding_abundance_abundant <- communitiesWithoutCrossfeeding_abundance %>%
    filter(Time == max(Time)) %>%
    filter(Abundance != 0) %>%
    group_by(Community) %>%
    mutate(TotalAbundance = sum(Abundance)) %>%
    filter(Abundance > 0.01 * sum(Abundance))

#communitiesWithoutCrossfeeding_richness <- read_csv(paste0(folder_simulation, "11-aggregated/communitiesWithoutCrossfeeding_richness.csv"), col_types = cols())
communitiesWithoutCrossfeeding_abundance_richness <- communitiesWithoutCrossfeeding_abundance_abundant %>%
    filter(Time == max(Time)) %>%
    group_by(Community, .drop = F) %>%
    summarize(Richness = n())

p <- communitiesWithoutCrossfeeding_abundance_abundant %>%
    filter(Time == max(Time)) %>%
    group_by(Community) %>%
    mutate(RelativeAbundance = Abundance / sum(Abundance)) %>%
    ggplot() +
    geom_col(aes(x = Community, y = RelativeAbundance, fill = Family), color = 1) +
    geom_text(data = communitiesWithoutCrossfeeding_abundance_richness, aes(x = Community, label = Richness), y = 1.1) +
    #annotate("text", x = 1:20, y = 1.15, label = communitiesWithoutCrossfeeding_richness$Richness, size = 4) +
    annotate("text", x = 21, y = 1.1, label = c("n. of species"), size = 4, hjust = 0) +
    #scale_fill_manual(values = c("F0" = "#8A89C0", "F1" = "#FFCB77"), labels = c("F0" = "fermenter", "F1" = "repirator")) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 1.1)) +
    coord_cartesian(xlim = c(0.5, 20.5), ylim = c(0, 1), clip = "off") +
    theme_classic() +
    theme(plot.margin = unit(c(1,.5,.5,.5), "cm")) +
    guides(color = "none", fill = guide_legend(title = "")) +
    labs(y = "Relative abundance")

ggsave(here::here("simulation/plots/22-communitiesWithoutCrossfeeding-03-bar_final.png"), p, width = 9, height = 3)

