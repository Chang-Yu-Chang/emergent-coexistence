library(tidyverse)
library(cowplot)
library(broom)
source(here::here("analysis/00-metadata.R"))

factorize_communities <- function (x) x %>% mutate(Community = factor(Community, paste0("C", rep(1:12, each = 8), "R", rep(1:8, 12))))
fitness2 <- read_csv(paste0(folder_data, "temp/15-fitness2.csv"), show_col_types = F) %>% factorize_communities
ESV_eq_freq2 <- read_csv(paste0(folder_data, "temp/15-ESV_eq_freq2.csv"), show_col_types = F) %>% factorize_communities
fitness_stable <- fitness2 %>% filter(ESVType == "stable") %>% left_join(ESV_eq_freq2) %>% replace_na(list(Significance = "p>=0.05"))
ESV_eq_freq_stable <- ESV_eq_freq2 %>% filter(ESVType == "stable")
fitness_transient2 <- fitness2 %>% filter(ESVType == "transient") %>% left_join(ESV_eq_freq2) %>% replace_na(list(Significance = "p>=0.05"))
ESV_eq_freq_transient2 <- ESV_eq_freq2 %>% filter(ESVType == "transient")


fitness_mean <- bind_rows(
    fitness_stable %>% mutate(ESVType = "stable"),
    fitness_transient2 %>% mutate(ESVType = "transient")
) %>%
    group_by(ESVType, Community, ESV_ID) %>%
    summarize(MeanFitness = mean(Fitness), NumberPoint = n(), SdFitness = sd(Fitness))
table(fitness_mean$ESVType) # 99 stable ESVs and 110 transient ESVs

nrow(ESV_eq_freq2) # 99+110 = 209 rows

fitness_eq_freq2 <- left_join(fitness_mean, ESV_eq_freq2)
fitness_eq_freq2_filtered <- fitness_eq_freq2 %>% filter(Slope < 0)
nrow(fitness_eq_freq2_filtered) # 175 ESVs
table(fitness_eq_freq2_filtered$ESVType) # 95 stable ESVs and 80 transient ESVs


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

# Distribution of fitness value with error bar----
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
    plot_grid(p1, NULL, ncol = 1, rel_heights = c(1, 2)), p2,
    labels = LETTERS[1:2], scale = c(.95, .95),
    nrow = 1, align = "h", axis = "tb", rel_widths = c(1,1.2)) +
    paint_white_background()
ggsave(here::here("plots/FigS5-transient_negative_fitness.png"), p, width = 10, height = 13)



if (FALSE) {

communities <- read_csv(paste0(folder_data, "output/communities_remained.csv"), show_col_types = F) %>%
    mutate(Community = factor(Community, Community))
communities_abundance <- read_csv(paste0(folder_data, "raw/community_ESV/Emergent_Comunity_Data.csv"), show_col_types = F) %>%
    filter(Carbon_Source == "Glucose" | Carbon_Source == "Original") %>%
    mutate(Community = factor(paste0("C", Inoculum, "R", Replicate), paste0("C", rep(1:12, each = 8), "R", rep(1:8, 12))))%>%
    arrange(Community, Family, Transfer, ESV)
communities_abundance_temporal <- communities_abundance %>%
    filter(Transfer != 0) %>%
    # Filter for those that has temporal data
    filter(Inoculum %in% c(2,6) | Replicate == 4) %>%
    select(Community, Transfer, ESV_ID, Relative_Abundance) %>%
    arrange(Community, Transfer, ESV_ID)

# A complete list of ESV-Transfer
communities_abundance_temporal_complete <- communities_abundance_temporal %>%
    distinct(Community, ESV_ID) %>%
    slice(rep(1:n(), each = 12)) %>%
    mutate(Transfer = rep(1:12, n()/12))
nrow(distinct(communities_abundance_temporal_complete, Community, ESV_ID)) # 755 unique ESVs
nrow(communities_abundance_temporal_complete) # 755 ESVs * 12 transfers = 9060 rows

# Calculate fitness
communities_abundance_fitness <- communities_abundance_temporal %>%
    right_join(communities_abundance_temporal_complete) %>%
    group_by(Community, ESV_ID) %>%
    arrange(Community, ESV_ID, Transfer) %>%
    mutate(Fitness = log(lead(Relative_Abundance) / Relative_Abundance))

# Stable ESVs ----
ESV_stable <- communities_abundance_temporal %>%
    filter(Transfer %in% c(9:12)) %>%
    pivot_wider(id_cols = c(Community, ESV_ID), names_from = Transfer, names_prefix = "T", values_from = Relative_Abundance) %>%
    filter(
        # Species with complete presence at T8-12
        (!is.na(T9) & !is.na(T10) & !is.na(T11) & !is.na(T12) ) |
            # Community C4R4, C9R4 have missing time points at T11. For these communities, add back the species
            (Community %in% c("C4R4", "C9R4") & (!is.na(T9) & !is.na(T10) & !is.na(T12))) |
            # Community C5R4 have missing time points at T9. For these communities, add back the species
            (Community %in% c("C5R4") & (!is.na(T10) & !is.na(T11) & !is.na(T12)))
    ) %>%
    distinct(Community, ESV_ID) %>%
    mutate(CommunityESV = paste0(Community, ESV_ID))

communities_abundance_fitness_stable <- communities_abundance_fitness %>%
    drop_na() %>%
    mutate(CommunityESV = paste0(Community, ESV_ID)) %>%
    filter(CommunityESV %in% ESV_stable$CommunityESV)

# Correlation
cor_test_long <- function (long_data) cor.test(long_data$Relative_Abundance, long_data$Fitness, method = "spearman", exact = F, alternative = c("less"))
tb_cor_stable <- communities_abundance_fitness_stable %>%
    nest(data = c(-Community, -ESV_ID)) %>%
    mutate(fit = map(data, ~ cor_test_long(.x)),
           tidied = map(fit, tidy)) %>%
    unnest(tidied)
ESV_sig_stable <- tb_cor_stable %>%
    select(Community, ESV_ID, estimate, p.value) %>%
    filter(p.value < 0.05, estimate < 0) %>%
    mutate(CommunityESV = paste0(Community, ESV_ID)) %>%
    mutate(Significance = "p<0.05")

nrow(tb_cor_stable) # 99 stable ESVs
nrow(ESV_sig_stable) # 52 shows significant negative correlation
range(ESV_sig_stable$estimate) # range of correlation [-0.95, -0.53]
tb_cor_stable %>% filter(Community %in% c("C4R4", "C5R4", "C9R4")) %>% nrow() # 9 ESVs from C4R4, C5R4, C9R4 have only three time points

# Empirical equilibrium frequency
eq_freq_stable <- communities_abundance_temporal %>%
    mutate(CommunityESV = paste0(Community, ESV_ID)) %>%
    filter(Transfer %in% 9:12) %>%
    filter(CommunityESV %in% ESV_stable$CommunityESV) %>%
    arrange(Community, ESV_ID, CommunityESV, Transfer) %>%  # 387 rows. 99 * 4 - 9 = 387. 9 ESVs from C4R4, C5R4, C9R4 have only three time points
    group_by(Community, ESV_ID, CommunityESV) %>%
    summarize(EquilibriumAbundance = mean(Relative_Abundance)) # 99 ESVs

# Linear model predicted equilibrium frequency
xintercept_stable <- communities_abundance_fitness_stable %>%
    nest(data = c(-Community, -ESV_ID, -CommunityESV)) %>%
    mutate(fit = map(data, ~ lm(Fitness ~ Relative_Abundance, data = .x)),
           tidied = map(fit, tidy)) %>%
    unnest(tidied) %>%
    select(Community, ESV_ID, CommunityESV, term, estimate) %>%
    pivot_wider(names_from = term, values_from = estimate) %>%
    # X intercept
    mutate(Xintercept = -`(Intercept)` / Relative_Abundance) %>%
    select(Community, ESV_ID, CommunityESV, Xintercept, Slope = Relative_Abundance, Yintercept = `(Intercept)`) %>%
    mutate(ESVType = "stable") %>%
    # Mark significant correlation
    left_join(select(tb_cor_stable, Community, ESV_ID, estimate, p.value))


# Ephemeral ESVs ----
ESV_ephemeral2 <- communities_abundance_fitness %>%
    drop_na() %>%
    mutate(CommunityESV = paste0(Community, ESV_ID)) %>%
    filter(!(CommunityESV %in% ESV_stable$CommunityESV)) %>%
    group_by(CommunityESV, Community, ESV_ID) %>%
    count() %>%
    # Include the ESVs with >=3 data points
    filter(n>=3) %>%
    # Remove the artifact C10R4 Stenotrophomonas
    filter(CommunityESV != "C10R4Stenotrophomonas")

communities_abundance_fitness_ephemeral2 <- communities_abundance_fitness %>%
    drop_na() %>%
    mutate(CommunityESV = paste0(Community, ESV_ID)) %>%
    filter(CommunityESV %in% ESV_ephemeral2$CommunityESV)

# Correlation
tb_cor_ephemeral2 <- communities_abundance_fitness_ephemeral2 %>%
    nest(data = c(-Community, -ESV_ID)) %>%
    mutate(fit = map(data, ~ cor_test_long(.x)),
           tidied = map(fit, tidy)) %>%
    unnest(tidied)

ESV_sig_ephemeral2 <- tb_cor_ephemeral2 %>%
    select(Community, ESV_ID, estimate, p.value) %>%
    filter(p.value < 0.05, estimate < 0) %>%
    mutate(CommunityESV = paste0(Community, ESV_ID)) %>%
    mutate(Significance = "p<0.05")

nrow(tb_cor_ephemeral2) # 110 ephemeral ESVs
nrow(ESV_sig_ephemeral2) # 21 ephemeral ESVs shows significant negative correlation
range(ESV_sig_ephemeral2$estimate) # range of correlation [-1, -0.714]

# Linear model predicted equilibrium frequency
xintercept_ephemeral2 <- communities_abundance_fitness_ephemeral2 %>%
    #filter(CommunityESV %in% ESV_sig_ephemeral$CommunityESV) %>%
    nest(data = c(-Community, -ESV_ID, -CommunityESV)) %>%
    mutate(fit = map(data, ~ lm(Fitness ~ Relative_Abundance, data = .x)),
           tidied = map(fit, tidy)) %>%
    unnest(tidied) %>%
    select(Community, ESV_ID, CommunityESV, term, estimate) %>%
    pivot_wider(names_from = term, values_from = estimate) %>%
    # X intercept
    mutate(Xintercept = -`(Intercept)` / Relative_Abundance) %>%
    select(Community, ESV_ID, CommunityESV, Xintercept, Slope = Relative_Abundance, Yintercept = `(Intercept)`) %>%
    mutate(ESVType = "transient") %>%
    # Mark significant correlation
    left_join(select(tb_cor_ephemeral2, Community, ESV_ID, estimate, p.value))


# Fitness ----
communities_abundance_fitness_mean <- bind_rows(
    communities_abundance_fitness_stable %>% mutate(ESVType = "stable"),
    communities_abundance_fitness_ephemeral2 %>% mutate(ESVType = "transient")
) %>%
    group_by(ESVType, Community, ESV_ID) %>%
    summarize(MeanFitness = mean(Fitness), NumberPoint = n(), SdFitness = sd(Fitness))
table(communities_abundance_fitness_mean$ESVType) # 99 stable ESVs and 110 transient ESVs

# Linear model predicted equilibrium frequency ----
xintercept <- bind_rows(xintercept_stable, xintercept_ephemeral2)
table(xintercept$ESVType) # 99 stable ESVs and 110 transient ESVs

fitness_eq_freq <- left_join(communities_abundance_fitness_mean, xintercept)
nrow(fitness_eq_freq) # 99+110 = 209 rows

fitness_eq_freq_filtered <- fitness_eq_freq %>%
    filter(Slope < 0)
nrow(fitness_eq_freq_filtered) # 175 ESVs
table(fitness_eq_freq_filtered$ESVType) # 95 stable ESVs and 80 transient ESVs


# Plots ----
p1 <- fitness_eq_freq_filtered %>%
    arrange(ESVType) %>%
    ggplot() +
    geom_vline(xintercept = 0, linewidth = 0.1, linetype = 2) +
    geom_hline(yintercept = 0, linewidth = 0.1, linetype = 2) +
    geom_point(aes(x = MeanFitness, y = Xintercept, color = ESVType), shape = 21, size = 3, stroke = .8) +
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

# Distribution of fitness value with error bar----
temp_order <- fitness_eq_freq_filtered %>%
    arrange(desc(ESVType), desc(Community)) %>%
    pull(CommunityESV)

p2 <- fitness_eq_freq_filtered %>%
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

# predicted equilibirum fitness
# p3 <- fitness_eq_freq_filtered %>%
#     mutate(CommunityESV = factor(CommunityESV, temp_order)) %>%
#     ggplot() +
#     geom_point(aes(x = CommunityESV, y = Xintercept, color = ESVType), shape = 21, size = 2, stroke = 1) +
#     #geom_segment(aes(x = CommunityESV, xend = CommunityESV, y = MeanFitness - SdFitness, yend = MeanFitness + SdFitness, color = ESVType)) +
#     geom_hline(yintercept = 0, linewidth = 0.5, linetype = 2) +
#     scale_color_manual(values = c(stable = "firebrick1", transient = grey(0.7)), label = c(stable = "stable ESV", transient = "transient ESV")) +
#     coord_flip() +
#     scale_x_discrete(position = "top") +
#     theme_classic() +
#     theme(
#         panel.border = element_rect(fill = NA, color = 1, linewidth = 1),
#         panel.grid.major.y = element_line(color = grey(0.9)),
#         legend.position = "none",
#         axis.title = element_text(size = 15),
#         axis.text.x = element_text(size = 10),
#         axis.text.y = element_text(size = 8, hjust = 0)
#         #axis.text.x = element_text(size = 6, angle = 45, hjust = 0)
#     ) +
#     guides() +
#     labs(x = "ESV", y = "predicted equilibrium frequency")

p <- plot_grid(
    plot_grid(p1, NULL, ncol = 1, rel_heights = c(1, 2)), p2,
    labels = LETTERS[1:2], scale = c(.95, .95),
    nrow = 1, align = "h", axis = "tb", rel_widths = c(1,1.2)) +
    paint_white_background()
ggsave(here::here("plots/FigS5-transient_negative_fitness.png"), p, width = 10, height = 13)



}
