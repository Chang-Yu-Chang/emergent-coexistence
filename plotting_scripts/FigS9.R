#' CURRENT THIS SCRIPT IS almost IDENTICAL TO 17-ESV_fitness_change.R

library(tidyverse)
library(cowplot)
source(here::here("analysis/00-metadata.R"))

communities <- read_csv(paste0(folder_data, "temp/00c-communities.csv"), show_col_types = F)
communities <- communities %>% mutate(Community = factor(Community, Community))
communities_abundance <- read_csv(paste0(folder_data, "raw/community_ESV/CSS_Emergent_Simplicity_Timeseries.csv"), show_col_types = F)

# 0. Clean up communities_abundance data ----
communities_abundance <- communities_abundance %>%
    #mutate(Family = ifelse(Relative_Abundance < 0.01 | !(Family %in% names(family_colors)), "Others", Family)) %>%
    mutate(Community = paste0("C", Inoculum, "R", Replicate)) %>%
    filter(Community %in% communities$Community) %>%
    mutate(Family = ifelse(!(Family %in% names(family_colors[-1])), "Others", Family)) %>%
    mutate(Family = factor(Family, names(family_colors))) %>%
    arrange(Community, Family, ESV) %>%
    mutate(CommunityESV = paste0(Community, "-", ESV)) %>%
    #filter(Family != "Others") %>%
    {.}

# Communities abundance at T12
communities_abundance_T1 <- communities_abundance %>%
    filter(Transfer == 1)
communities_abundance_T12 <- communities_abundance %>%
    mutate(Family = ifelse(Relative_Abundance < 0.01 | !(Family %in% names(family_colors)), "Others", Family)) %>%
    filter(Relative_Abundance > 0.01) %>%
    filter(Transfer == 12)

# Subset species that are present both at T1 and at T12
communities_abundance_sp <- communities_abundance %>%
    filter(CommunityESV %in% communities_abundance_T12$CommunityESV) %>%
    filter(CommunityESV %in% communities_abundance_T1$CommunityESV) %>%
    filter(Family %in% names(family_colors)) %>%
    filter(Transfer != 0)

# Richness
communities_abundance_sp_richness <- communities_abundance_sp %>%
    filter(Transfer == 12) %>%
    group_by(Community) %>%
    summarize(Count = n())

# Community relative abundance over time ----
p1 <- communities_abundance %>%
    arrange(Transfer, Family) %>%
    ggplot() +
    geom_col(aes(x = Transfer, y = Relative_Abundance, fill = Family, alpha = ESV), linewidth = .3) +
    scale_fill_manual(values = family_colors) +
    scale_alpha_discrete(range = c(0.5, 1)) +
    scale_x_continuous(expand = c(0,.3), breaks = seq(0, 12, 2)) +
    scale_y_continuous(expand = c(0, 0), breaks = seq(0,1, 0.2)) +
    facet_wrap(~Community, ncol = 1, scale = "free") +
    theme_classic() +
    theme(panel.grid = element_blank(),
          axis.title = element_text(size = 12),
          axis.text = element_text(color = 1, size = 12),
          strip.background = element_rect(color = NA, fill = NA),
          strip.text = element_text(size = 12),
          panel.spacing = unit(1, units = "cm"),
          panel.border = element_rect(color = 1, linewidth = 1, fill = NA)) +
    guides(alpha = "none", fill = "none") +
    labs(x = "Transfer", y = "Relative abundance")


# Species relative abundance over time ----
plot_species_abundance <- function (species_abundance) {
    species_abundance %>%
        filter(CommunityESV %in% communities_abundance_sp$CommunityESV) %>%
        ggplot() +
        geom_col(aes(x = Transfer, y = Relative_Abundance, fill = Family, alpha = ESV), linewidth = .3) +
        geom_hline(yintercept = 0.01, linetype = 2) +
        scale_fill_manual(values = family_colors) +
        scale_alpha_discrete(range = c(0.5, 1)) +
        scale_x_continuous(breaks = seq(0,12,2)) +
        facet_wrap(~ESV, nrow = 1, scales = "free_y") +
        theme_classic() +
        theme(panel.border = element_rect(color = 1, fill = NA),
              strip.background = element_rect(color = NA, fill = NA),
              plot.margin = unit(rep(0.5, 4), "cm"),
              #plot.background = element_rect(color = 1, linewidth = 1, fill = NA)
        ) +
        guides(alpha = "none", fill = "none") +
        labs(x = "Transfer", y = "Relative abundance")
}

temp <- communities_abundance_sp_richness %>%
    rename(comm = Community) %>%
    rowwise() %>%
    mutate(CommunityAbundance = communities_abundance %>% filter(Community == comm) %>% list) %>%
    mutate(CommunityAbundnacePlot = list(plot_species_abundance(CommunityAbundance)))

p2 <- plot_grid(
    plot_grid(temp$CommunityAbundnacePlot[[1]], NULL,
              nrow = 1, rel_widths = c(2, 2), scale = 1, label_x = 0
              #labels = c(communities_abundance_sp_richness$Community[1])
              ),
    plot_grid(temp$CommunityAbundnacePlot[[2]], NULL,
              nrow = 1, rel_widths = c(3, 1), scale = 1, label_x = 0
              #labels = c(communities_abundance_sp_richness$Community[2])
              ),
    plot_grid(temp$CommunityAbundnacePlot[[3]], NULL,
              nrow = 1, rel_widths = c(1, 0), scale = 1, label_x = 0
              #labels = c(communities_abundance_sp_richness$Community[3])
              ),
    plot_grid(temp$CommunityAbundnacePlot[[4]], NULL,
              nrow = 1, rel_widths = c(3, 1), scale = 1, label_x = 0
              #labels = c(communities_abundance_sp_richness$Community[4])
              ),
    ncol = 1,
    scale = 1
) + paint_white_background()


p <- plot_grid(p1, p2, nrow = 1, rel_widths = c(1, 3), labels = c("A", "B")) + paint_white_background()
ggsave(here::here("plots/FigS9-communities_abundance.png"), p, width = 10, height = 10)


