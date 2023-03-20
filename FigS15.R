library(tidyverse)
library(cowplot)
source(here::here("analysis/00-metadata.R"))
source(here::here("simulation/01-generate_input.R"))

# 0. parameters ----
input_parameters <- read_csv(here::here("simulation/01-input_parameters.csv"), col_types = cols())
input_monocultures <- read_csv(here::here("simulation/02a-input_monocultures.csv"), col_types = cols())
input_communities <- read_csv(here::here("simulation/02b-input_communities.csv"), col_types = cols())
input_communitiesWithoutCrossfeeding <- read_csv(here::here("simulation/02c-input_communitiesWithoutCrossfeeding.csv"), col_types = cols())


# 2. Communities ----
communities_abundance <- read_csv(paste0(folder_simulation, "aggregated/12-communities_abundance.csv"), col_types = cols()) %>%
    mutate(Community = ordered(Community, paste0("W", 0:(input_communities$n_wells[1]-1)))) %>%
    mutate(Time = ordered(Time, c("init", paste0("T", 1:20), "end"))) %>%
    arrange(Community, Time)

# Barplot over time, standard
temp <- tibble(Community = paste0("W", 0:19),
               Community_new = factor(paste0("community ", 1:20), paste0("community ", 1:20)))
p1 <- communities_abundance %>%
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
    facet_wrap(.~Community_new, ncol = 5) +
    theme_classic() +
    theme(strip.background = element_rect(color = NA, fill = NA)) +
    guides(color = "none", fill = guide_legend(title = "Family")) +
    labs(x = "transfer", y = "relative abundance")

#ggsave(here::here("simulation/plots/22-communities-02-bar_fraction.png"), p, width = 10, height = 10)


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


p2 <- communities_abundance_abundant %>%
    filter(Time == max(Time)) %>%
    group_by(Community) %>%
    mutate(RelativeAbundance = Abundance / sum(Abundance)) %>%
    left_join(temp) %>%
    ggplot() +
    geom_col(aes(x = Community_new, y = RelativeAbundance, fill = Family), color = 1) +
    geom_text(data = communities_abundance_richness, aes(x = Community_new, label = Richness), y = 1.1) +
    annotate("text", x = 1:20, y = 1.15, label = communities_richness$Richness, size = 4) +
    annotate("text", x = 21, y = 1.1, label = c("n. of species"), size = 4, hjust = 0) +
    scale_x_discrete(expand = c(0,.1)) +
    scale_y_continuous(breaks = seq(0,1, 0.2), expand = c(0,0), limits = c(0,1.1)) +
    coord_cartesian(xlim = c(0.5, 20.5), ylim = c(0, 1), clip = "off") +
    theme_classic() +
    theme(plot.margin = unit(c(1,.5,.5,.5), "cm"),
          panel.border = element_rect(color = 1, fill = NA)) +
    guides(color = "none", fill = guide_legend(title = "Family")) +
    labs(x = "community", y = "relative abundance")

#ggsave(here::here("simulation/plots/22-communities-03-bar_final.png"), p, width = 9, height = 3)

if (FALSE) {
# 3. Communities without crossfeeding -----
communitiesWithoutCrossfeeding_abundance <- read_csv(paste0(folder_simulation, "aggregated/12-communitiesWithoutCrossfeeding_abundance.csv"), col_types = cols()) %>%
    mutate(Community = ordered(Community, paste0("W", 0:(input_communitiesWithoutCrossfeeding$n_wells[1]-1)))) %>%
    mutate(Time = ordered(Time, c("init", paste0("T", 1:20), "end"))) %>%
    arrange(Community, Time)

# Barplot over time, fraction
p3 <- communitiesWithoutCrossfeeding_abundance %>%
    filter(Abundance != 0) %>%
    filter(Time != "end") %>%
    filter(Community == "W0") %>%
    left_join(temp) %>%
    ggplot(aes(x = Time, y = Abundance, fill = Family, color = Species)) +
    geom_col(position = "fill", linewidth = .2) +
    scale_color_manual(values = rep("black", length(sal$Species))) +
    scale_x_discrete(breaks = c("init", paste0("T", seq(5,20, 5))), labels = c(0, seq(5, 20, 5)), expand = c(0,0)) +
    scale_y_continuous(breaks = seq(0,1,0.2), expand = c(0,0)) +
    theme_classic() +
    theme(strip.background = element_rect(color = NA, fill = NA),
          panel.border = element_rect(color = 1, fill = NA)) +
    guides(color = "none", fill = guide_legend(title = "Family")) +
    labs(x = "transfer", y = "relative abundance") +
    ggtitle("community 1 without cross-feeding")


}

# Assemble panels ----
p <- plot_grid(
    p1 + guides(fill = "none"),
    #plot_grid(p3 + guides(fill = "none"), p2 + guides(fill = "none"), nrow = 1, rel_widths = c(1, 3), scale = c(.9, .95), labels = c("B", "C")),
    p2,
    ncol = 1, rel_heights = c(3,1), labels = c("A", "B"), scale = c(.95, .95)) +
    paint_white_background()
ggsave(here::here("plots/FigS15-in_silico_assembly.png"), p, width = 10, height = 12)










