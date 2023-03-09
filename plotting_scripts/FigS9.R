#' CURRENT THIS SCRIPT IS almost IDENTICAL TO 17-ESV_fitness_change.R

library(tidyverse)
library(cowplot)
source(here::here("analysis/00-metadata.R"))

communities <- read_csv(paste0(folder_data, "temp/00c-communities.csv"), show_col_types = F)
communities <- communities %>% mutate(Community = factor(Community, Community))
communities_abundance <- read_csv(paste0(folder_data, "raw/community_ESV/CSS_Emergent_Simplicity_Timeseries.csv"), show_col_types = F)

# 0. Clean up communities_abundance data ----
communities_abundance <- communities_abundance %>%
    mutate(Family = ifelse(Relative_Abundance < 0.01 | !(Family %in% names(family_colors)), "Others", Family)) %>%
    mutate(Community = paste0("C", Inoculum, "R", Replicate)) %>%
    filter(Community %in% communities$Community) %>%
    mutate(Family = factor(Family, names(family_colors))) %>%
    arrange(Community, Family, ESV) %>%
    mutate(CommunityESV = paste0(Community, "-", ESV))

# Community relative abundance over time ----
p <- communities_abundance %>%
    arrange(Transfer, Family) %>%
    ggplot() +
    geom_col(aes(x = Transfer, y = Relative_Abundance, fill = Family, alpha = ESV), linewidth = .3) +
    scale_fill_manual(values = family_colors) +
    scale_alpha_discrete(range = c(0.5, 1)) +
    scale_x_continuous(breaks = 0:12, expand = c(0,.3)) +
    scale_y_continuous(expand = c(0, 0), breaks = seq(0,1, 0.2)) +
    facet_wrap(~Community, scale = "free") +
    theme_classic() +
    theme(panel.grid = element_blank(),
          axis.title = element_text(size = 12),
          axis.text = element_text(color = 1, size = 12),
          strip.background = element_rect(color = NA, fill = NA),
          strip.text = element_text(size = 12),
          panel.spacing = unit(1, units = "cm"),
          panel.border = element_rect(color = 1, linewidth = 1, fill = NA)) +
    guides(alpha = "none") +
    labs(x = "Transfer", y = "Relative abundance")
ggsave(here::here("plots/FigS9-communities_abundance.png"), p, width = 10, height = 6)











