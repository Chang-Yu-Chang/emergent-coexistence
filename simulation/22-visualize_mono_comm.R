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

mcrm_family_colors <- RColorBrewer::brewer.pal(10, "Set3") %>% setNames(paste0("F", 0:9))
#n_timepoints <- input_communities$n_timesteps[1]
n_timesteps <- input_communities$n_timesteps[1]
n_timepoints <- input_communities$n_timepoints[1]
time_ids <- tibble(Time = c("init", paste0("T", 1:10000)), TimeID = 0:10000)

# 1. Monocultures ----
monocultures_abundance <- read_csv(paste0(folder_simulation, "aggregated/12-monocultures_abundance.csv"), col_types = cols()) %>%
    mutate(Well = ordered(Well, paste0("W", 0:(input_monocultures$n_wells[1]-1)))) %>%
    mutate(Time = ordered(Time, c("init", paste0("T", 1:n_timepoints)))) %>%
    arrange(Well, Time)

# Line plot
p1 <- monocultures_abundance %>%
    filter(Abundance > 0) %>%
    left_join(time_ids) %>%
    #mutate(UniqueWell = paste0(Well, Species)) %>%
    ggplot(aes(x = TimeID, y = Abundance, color = Family, group = Species)) +
    geom_line(linewidth = .5) +
    scale_color_manual(values = mcrm_family_colors) +
    #geom_point(size = 1, shape = 21) +
    #scale_color_manual(values = c("F0" = "#8A89C0", "F1" = "#FFCB77"), labels = c("F0" = "fermenter", "F1" = "repirator")) +
    #scale_linetype_manual(values = c("F0" = 1, "F1" = 2)) +
    theme_classic() +
    theme(panel.border = element_rect(color = 1, fill = NA)) +
    guides(alpha = "none", color = guide_legend(title = "")) +
    labs()

if (FALSE) {
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
    labs(x = "time steps")

}
#p <- plot_grid(p1, p2, scale = 0.9) + paint_white_background()
p <- p1
ggsave(here::here("simulation/plots/22-monoculture-01-line.png"), p, width = 5, height = 4)


# 2. Communities ----
communities_abundance <- read_csv(paste0(folder_simulation, "aggregated/12-communities_abundance.csv"), col_types = cols()) %>%
    mutate(Community = ordered(Community, paste0("W", 0:(input_communities$n_wells[1]-1)))) %>%
    mutate(Time = ordered(Time, c("init", paste0("T", 1:n_timepoints), "end"))) %>%
    arrange(Community, Time)

# Barplot over time, standard
temp <- tibble(Community = paste0("W", 0:19),
       Community_new = factor(paste0("community ", 1:20), paste0("community ", 1:20)))
p <- communities_abundance %>%
    filter(Abundance != 0) %>%
    filter(Time != "end") %>%
    left_join(temp) %>%
    ggplot(aes(x = Time, y = Abundance, fill = Family, color = Species)) +
    geom_col(position = "fill", linewidth = .2) +
    scale_fill_manual(values = mcrm_family_colors) +
    scale_color_manual(values = rep("black", length(sal$Species))) +
    scale_x_discrete(labels = n_timesteps/1e6 * seq(0, n_timepoints, 1)) +
    scale_y_continuous(breaks = seq(0,1,0.2), expand = c(0,0)) +
    facet_wrap(.~Community_new, ncol = 4) +
    theme_classic() +
    theme(strip.background = element_rect(color = NA, fill = NA)) +
    guides(color = "none", fill = guide_legend(title = "Family")) +
    labs(x = "time unit [10^6]", y = "relative abundance")

ggsave(here::here("simulation/plots/22-communities-02-bar_fraction.png"), p, width = 10, height = 10)


# Barplot final time point
communities_richness <- read_csv(paste0(folder_simulation, "aggregated/12-communities_richness.csv"), col_types = cols()) %>%
    #mutate(Community = factor(Community, paste0("W", 0:(nrow(input_withinCommunityPairs)-1)))) %>%
    mutate(PairSize = choose(Richness, 2)) %>%
    arrange(desc(Richness)) %>%
    slice(1:20) %>%
    mutate(CommunityLabel = 1:20)

p <- communities_abundance %>%
    filter(Time == max(Time), Abundance > 0) %>%
    group_by(Community) %>%
    mutate(RelativeAbundance = Abundance / sum(Abundance)) %>%
    left_join(communities_richness) %>%
    ggplot() +
    geom_col(aes(x = CommunityLabel, y = RelativeAbundance, fill = Family), color = 1) +
    geom_text(data = communities_richness, aes(x = CommunityLabel, label = Richness), y = 1.1) +
    annotate("text", x = 1:20, y = 1.15, label = communities_richness$Richness, size = 4) +
    annotate("text", x = 21, y = 1.1, label = c("n. of species"), size = 4, hjust = 0) +
    scale_fill_manual(values = mcrm_family_colors) +
    scale_x_continuous(breaks = 1:20, expand = c(0,0.1)) +
    scale_y_continuous(breaks = seq(0,1, 0.2), expand = c(0,0), limits = c(0,1.1)) +
    coord_cartesian(xlim = c(0.5, 20.5), ylim = c(0, 1), clip = "off") +
    theme_classic() +
    theme(plot.margin = unit(c(10,10,5,5), "mm"),
          panel.border = element_rect(color = 1, fill = NA)) +
    guides(color = "none", fill = guide_legend(title = "Family", ncol = 1)) +
    labs(x = "community", y = "relative abundance")

ggsave(here::here("simulation/plots/22-communities-04-bar_final.png"), p, width = 9, height = 3)


# 3. Communities without crossfeeding -----
communitiesWithoutCrossfeeding_abundance <- read_csv(paste0(folder_simulation, "aggregated/12-communitiesWithoutCrossfeeding_abundance.csv"), col_types = cols()) %>%
    mutate(Community = ordered(Community, paste0("W", 0:(input_communitiesWithoutCrossfeeding$n_wells[1]-1)))) %>%
    mutate(Time = ordered(Time, c("init", paste0("T", 1:n_timepoints), "end"))) %>%
    arrange(Community, Time)

# Barplot final time point
communitiesWithoutCrossfeeding_richness <- read_csv(paste0(folder_simulation, "aggregated/12-communitiesWithoutCrossfeeding_richness.csv"), col_types = cols()) %>%
    mutate(PairSize = choose(Richness, 2)) %>%
    arrange(desc(Richness)) %>%
    slice(1:20) %>%
    mutate(CommunityLabel = 1:20)

p <- communitiesWithoutCrossfeeding_abundance %>%
    filter(Time == max(Time), Abundance > 0) %>%
    group_by(Community) %>%
    mutate(RelativeAbundance = Abundance / sum(Abundance)) %>%
    left_join(communitiesWithoutCrossfeeding_richness) %>%
    ggplot() +
    geom_col(aes(x = CommunityLabel, y = RelativeAbundance, fill = Family), color = 1) +
    geom_text(data = communitiesWithoutCrossfeeding_richness, aes(x = CommunityLabel, label = Richness), y = 1.1) +
    annotate("text", x = 1:20, y = 1.15, label = communitiesWithoutCrossfeeding_richness$Richness, size = 4) +
    annotate("text", x = 21, y = 1.1, label = c("n. of species"), size = 4, hjust = 0) +
    scale_fill_manual(values = mcrm_family_colors) +
    scale_x_continuous(breaks = 1:20, expand = c(0,0.1)) +
    scale_y_continuous(breaks = seq(0,1, 0.2), expand = c(0,0), limits = c(0,1.1)) +
    coord_cartesian(xlim = c(0.5, 20.5), ylim = c(0, 1), clip = "off") +
    theme_classic() +
    theme(plot.margin = unit(c(10,10,5,5), "mm"),
          panel.border = element_rect(color = 1, fill = NA)) +
    guides(color = "none", fill = guide_legend(title = "Family", ncol = 2)) +
    labs(x = "community", y = "relative abundance")

ggsave(here::here("simulation/plots/22-communitiesWithoutCrossfeeding-03-bar_final.png"), p, width = 9, height = 3)

