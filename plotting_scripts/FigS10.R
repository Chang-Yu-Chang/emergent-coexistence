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

# Communities abundance at T12
communities_abundance_T1 <- communities_abundance %>%
    filter(Transfer == 1)
communities_abundance_T12 <- communities_abundance %>%
    mutate(Family = ifelse(Relative_Abundance < 0.01 | !(Family %in% names(family_colors)), "Others", Family)) %>%
    filter(Relative_Abundance > 0.01) %>%
    filter(Transfer == 12)

# Subset species that are present both at T1 and at T12
communities_abundance_sp <- communities_abundance %>%
    filter(CommunityESV %in% communities_abundance_T12$CommunityESV,
           CommunityESV %in% communities_abundance_T1$CommunityESV) %>%
    filter(Family %in% names(family_colors)) %>%
    filter(Transfer != 0)

# Richness
communities_abundance_sp_richness <- communities_abundance_sp %>%
    filter(Transfer == 12) %>%
    group_by(Community) %>%
    summarize(Count = n())

# Calculate Malthusian fitness
communities_abundance_fitness <- communities_abundance_sp %>%
    select(CommunityESV, Transfer, Relative_Abundance) %>%
    arrange(CommunityESV, Transfer) %>%
    pivot_wider(names_from = Transfer, names_prefix = "T", values_from = Relative_Abundance) %>%
    mutate(Fitness = log(T12/T1)) %>%
    select(CommunityESV, Fitness)


# Relative abundacne of each species over time
if (FALSE) {

p <- communities_abundance %>%
    filter(CommunityESV %in% communities_abundance_sp$CommunityESV) %>%
    group_by(Community) %>%
    ggplot() +
    geom_col(aes(x = Transfer, y = Relative_Abundance, fill = Family, alpha = ESV), linewidth = .3) +
    geom_hline(yintercept = 0.01, linetype = 2) +
    scale_fill_manual(values = family_colors) +
    scale_alpha_discrete(range = c(0.5, 1)) +
    #facet_grid(CommunityESV ~ Community, scales = "free_y") +
    facet_wrap(CommunityESV ~ Community, scales = "free_y") +
    theme_classic() +
    theme(panel.border = element_rect(color = 1, fill = NA, linewidth = 1)) +
    guides(alpha = "none") +
    labs(x = "Transfer", y = "Relative abundance")
}
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
              plot.background = element_rect(color = 1, linewidth = 1, fill = NA)
              ) +
        guides(alpha = "none", fill = "none") +
        labs(x = "Transfer", y = "Relative abundance")
}

temp <- communities_abundance_sp_richness %>%
    rename(comm = Community) %>%
    rowwise() %>%
    mutate(CommunityAbundance = communities_abundance %>% filter(Community == comm) %>% list) %>%
    mutate(CommunityAbundnacePlot = list(plot_species_abundance(CommunityAbundance)))

# p1 <- plot_grid(
#     plot_grid(temp$CommunityAbundnacePlot[[1]], temp$CommunityAbundnacePlot[[3]],
#               ncol = 1, rel_heights = c(2,4), scale = 1, label_x = 0,
#               labels = c(communities_abundance_sp_richness$Community[c(1,3)])),
#     plot_grid(temp$CommunityAbundnacePlot[[2]], temp$CommunityAbundnacePlot[[4]],
#               ncol = 1, rel_heights = c(3,3), scale = 1, label_x = 0,
#               labels = c(communities_abundance_sp_richness$Community[c(2,4)])),
#     nrow = 1,
#     rel_heights = c(1,1),
#     scale = 0.9
# ) +
#     paint_white_background()

p1 <- plot_grid(
    plot_grid(temp$CommunityAbundnacePlot[[1]], NULL, temp$CommunityAbundnacePlot[[3]],
              nrow = 1, rel_widths = c(2, 0.05, 4), scale = 1, label_x = 0,
              labels = c(communities_abundance_sp_richness$Community[c(1,3)])),
    plot_grid(temp$CommunityAbundnacePlot[[2]], NULL, temp$CommunityAbundnacePlot[[4]],
              nrow = 1, rel_widths = c(3, 0.05, 3), scale = 1, label_x = 0,
              labels = c(communities_abundance_sp_richness$Community[c(2,4)])),
    ncol = 1,
    scale = 0.95
) + paint_white_background()

# x_T1 vs. log(x_T12/x_T1)
p2 <- communities_abundance_fitness %>%
    left_join(communities_abundance_T1) %>%
    ggplot(aes(x = Relative_Abundance, y = Fitness)) +
    geom_point(shape = 21, size = 3, stroke = 1) +
    geom_hline(yintercept = 0, linetype = 2) +
    theme_classic() +
    labs(x = expression(x[T1]), y = expression(log(x[T12]/x[T1])))


# log(x[Ti+1]/x[i]) over transfers

communities_abundance_time <- communities_abundance_sp %>%
    select(CommunityESV, Transfer, Relative_Abundance) %>%
    arrange(CommunityESV, Transfer) %>%
    group_by(CommunityESV) %>%
    # time step change
    mutate(Fitness = log(lead(Relative_Abundance) / Relative_Abundance))

communities_abundance_time %>%
    ggplot(aes(x = Relative_Abundance, y = Fitness)) +
    geom_point(shape = 21, size = 3, stroke = 1) +
    geom_hline(yintercept = 0, linetype = 2) +
    #facet_wrap(~CommunityESV) +
    #facet_wrap(~Transfer, scales = "free") +
    theme_classic() +
    theme(panel.border = element_rect(color = 1, fill = NA)) +
    labs(x = "x_T", y = "log(x_{Ti+1}/x_{Ti})")

p3 <- communities_abundance_time %>%
    filter(Transfer != 12) %>%
    ggplot() +
    geom_boxplot(aes(x = Transfer, y = Fitness, group = Transfer), outlier.color = NA) +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_point(aes(x = Transfer, y = Fitness), position = position_jitter(width = 0.1, height = 0),
               size = 2, shape = 21) +
    scale_x_continuous(breaks = 1:12) +
    theme_classic() +
    theme() +
    labs(x = "Transfer", y = expression(log(x[Ti+1]/x[i])))


#

p <- plot_grid(
    p1,
    plot_grid(p2, p3, NULL, rel_widths= c(1, 1.5, 0.5), nrow = 1, scale = 1),
    rel_heights = c(1.5, 1),
    ncol = 1,
    scale = c(1),
    labels = c("A", "B")
) + paint_white_background()

ggsave(here::here("plots/FigS10-species_abundance.png"), p, width = 15, height = 10)

if (FALSE) {
    # x_Ti vs. log(x_{Ti+1}/x_{Ti})
    communities_abundance_sp %>%
        select(CommunityESV, Transfer, Relative_Abundance) %>%
        arrange(CommunityESV, Transfer) %>%
        group_by(CommunityESV) %>%
        # time step change
        mutate(Fitness = log(lead(Relative_Abundance) / Relative_Abundance)) %>%
        ggplot(aes(x = Relative_Abundance, y = Fitness)) +
        geom_point(shape = 21, size = 3, stroke = 1) +
        geom_hline(yintercept = 0, linetype = 2) +
        #facet_wrap(~CommunityESV) +
        #facet_wrap(~Transfer, scales = "free") +
        theme_classic() +
        theme(panel.border = element_rect(color = 1, fill = NA)) +
        labs(x = "x_T", y = "log(x_{Ti+1}/x_{Ti})")
    #ggsave(paste0(folder_data, "temp/17-fitness_changes.png"), p, width = 15, height = 10)

    ## Statistics
    library(broom)
    communities_abundance_sp %>%
        select(CommunityESV, Transfer, Relative_Abundance) %>%
        arrange(CommunityESV, Transfer) %>%
        group_by(CommunityESV) %>%
        # time step change
        mutate(Fitness = log(lead(Relative_Abundance) / Relative_Abundance))  %>%
        filter(!is.na(Fitness)) %>%
        #
        nest(data = c(-CommunityESV)) %>%
        mutate(
            fit = map(data, ~ lm(Fitness ~ Relative_Abundance, data = .x)),
            tidied = map(fit, tidy)
        ) %>%
        unnest(tidied) %>%
        select(-data, -fit) %>%
        filter(term == "Relative_Abundance") %>%
        filter(p.value < 0.05)


}

