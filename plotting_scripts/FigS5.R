library(tidyverse)
library(cowplot)
library(broom)
source(here::here("processing_scripts/00-metadata.R"))

factorize_communities <- function (x) x %>% mutate(Community = factor(Community, paste0("C", rep(1:12, each = 8), "R", rep(1:8, 12))))
fitness_stable <- read_csv(paste0(folder_data, "temp/15-fitness_stable.csv"), show_col_types = F) %>% factorize_communities
eq_freq_stable <- read_csv(paste0(folder_data, "temp/15-eq_freq_stable.csv"), show_col_types = F) %>% factorize_communities
fitness_transient2 <- read_csv(paste0(folder_data, "temp/15-fitness_transient2.csv"), show_col_types = F) %>% factorize_communities
eq_freq_transient2 <- read_csv(paste0(folder_data, "temp/15-eq_freq_transient2.csv"), show_col_types = F) %>% factorize_communities

fitness_stable <- fitness_stable %>% left_join(eq_freq_stable) %>% replace_na(list(Significance = "p>=0.05"))
fitness_transient2 <- fitness_transient2 %>% left_join(eq_freq_transient2) %>% replace_na(list(Significance = "p>=0.05"))


fitness_mean <- bind_rows(
    fitness_stable %>% mutate(ESVType = "stable"),
    fitness_transient2 %>% mutate(ESVType = "transient")
) %>%
    group_by(ESVType, Community, ESV_ID) %>%
    summarize(MeanFitness = mean(Fitness), NumberPoint = n(), SdFitness = sd(Fitness))
table(fitness_mean$ESVType) # 99 stable ESVs and 110 transient ESVs


eq_freq2 <- bind_rows(eq_freq_stable, eq_freq_transient2)
nrow(eq_freq2) # 99+110 = 209 rows

fitness_eq_freq2 <- left_join(fitness_mean, eq_freq2)
fitness_eq_freq2_filtered <- fitness_eq_freq2 %>% filter(Slope < 0)
nrow(fitness_eq_freq2_filtered) # 175 ESVs
table(fitness_eq_freq2_filtered$ESVType) # 95 stable ESVs and 80 transient ESVs


#
p1 <- fitness_eq_freq2_filtered %>%
    arrange(ESVType) %>%
    ggplot() +
    geom_vline(xintercept = 0, linewidth = 0.1, linetype = 2) +
    geom_hline(yintercept = 0, linewidth = 0.1, linetype = 2) +
    geom_point(aes(x = MeanFitness, y = PredictedEqAbundance, color = ESVType), shape = 21, size = 3, stroke = .8) +
    #geom_segment(aes(x = MeanFitness - 2*SdFitness, xend = MeanFitness + 2*SdFitness, y = Xintercept, yend = Xintercept, color = ESVType)) +
    scale_color_manual(values = c(stable = "firebrick1", transient = grey(0.7)), label = c(stable = "stable ESV", transient = "transient ESV")) +
    theme_classic() +
    theme(
        panel.border = element_rect(color = 1, fill = NA, linewidth = 1),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 15),
        legend.position = "bottom",
        legend.title = element_blank(),
        legend.margin = margin(0,3,1,1, unit = "mm"),
        legend.box.margin = margin(0,0,0,0, unit = "mm"),
        legend.text = element_text(size = 15),
        legend.key.size = unit(15, "mm"),
        legend.spacing.x = unit(1, "mm"),
        legend.spacing.y = unit(2, "mm"),
        legend.background = element_rect(color = NA, fill = NA),

    ) +
    guides(color = guide_legend(title = "", override.aes = list(stroke = 1))) +
    labs(x = "average fitness value", y = "equilibrium frequency\npredicted from assembly dynamics")

# Distribution of fitness value with error bar
temp_order <- fitness_eq_freq2_filtered %>%
    arrange(desc(ESVType), desc(Community)) %>%
    pull(CommunityESV)

p2 <- fitness_eq_freq2_filtered %>%
    mutate(CommunityESV = factor(CommunityESV, temp_order)) %>%
    ggplot() +
    geom_point(aes(x = CommunityESV, y = MeanFitness, color = ESVType), shape = 21, size = 2, stroke = 1) +
    geom_segment(aes(x = CommunityESV, xend = CommunityESV, y = MeanFitness - SdFitness, yend = MeanFitness + SdFitness, color = ESVType)) +
    geom_hline(yintercept = 0, linewidth = 0.5, linetype = 2) +
    scale_color_manual(values = c(stable = "firebrick1", transient = grey(0.7)), label = c(stable = "stable ESV", transient = "transient ESV")) +
    coord_flip() +
    scale_x_discrete(position = "top") +
    theme_classic() +
    theme(
        panel.border = element_rect(fill = NA, color = 1, linewidth = 1),
        panel.grid.major.y = element_line(color = grey(0.9)),
        legend.position = "none",
        axis.title = element_text(size = 15),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 6, hjust = 0)
        #axis.text.x = element_text(size = 6, angle = 45, hjust = 0)
    ) +
    guides() +
    labs(x = "", y = "invasion fitness")

# Assemble panels
p <- plot_grid(
    plot_grid(p1, NULL, ncol = 1, rel_heights = c(1, 1.5)), p2,
    labels = LETTERS[1:2], scale = c(.95, .95),
    nrow = 1, align = "h", axis = "tb", rel_widths = c(1,1.2)) +
    paint_white_background()
ggsave(here::here("plots/FigS5.png"), p, width = 10, height = 13)

