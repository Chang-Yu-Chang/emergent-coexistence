library(tidyverse)
library(cowplot)
library(broom)
source(here::here("analysis/00-metadata.R"))

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


# FigS3: Abundance vs. fitness for stable ESVs ----
# ESVs that have data from T9-12. Three excpetions: C4R4 and C9R4 do not have T11, and C5R4 do not have T9
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

# Empirical equilibrium frequency
communities_eq_freq_stable <- communities_abundance_temporal %>%
    mutate(CommunityESV = paste0(Community, ESV_ID)) %>%
    filter(Transfer %in% 9:12) %>%
    filter(CommunityESV %in% ESV_stable$CommunityESV) %>%
    arrange(Community, ESV_ID, CommunityESV, Transfer) %>%  # 387 rows. 99 * 4 - 9 = 387. 9 ESVs from C4R4, C5R4, C9R4 have only three time points
    group_by(Community, ESV_ID, CommunityESV) %>%
    summarize(EquilibriumAbundance = mean(Relative_Abundance)) # 99 ESVs

tb_cor_stable %>% filter(Community %in% c("C4R4", "C5R4", "C9R4")) %>% nrow() # 9 ESVs from C4R4, C5R4, C9R4 have only three time points

# Linear model predicted equilibrium frequency
xintercept_stable <- communities_abundance_fitness_stable %>%
    #filter(CommunityESV %in% ESV_sig_stable$CommunityESV) %>%
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

p <- communities_abundance_fitness_stable %>%
    ggplot() +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_point(aes(x = Relative_Abundance, y = Fitness), shape = 21) +
    geom_smooth(data = communities_abundance_fitness_stable %>% left_join(ESV_sig_stable) %>% replace_na(list(Significance = "p>=0.05")),
                aes(x = Relative_Abundance, y = Fitness, color = Significance), method = stats::lm, formula = y ~ x, se = F) +
    # Linear model predicted
    geom_vline(data = xintercept_stable, aes(xintercept = Xintercept), color = "green", linetype = 2) +
    # Mean of T9-12
    geom_vline(data = communities_eq_freq_stable, aes(xintercept = EquilibriumAbundance), color = "navyblue", linetype = 2) +
    scale_color_manual(values = c("p<0.05" = "pink", "p>=0.05" = grey(0.8))) +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 3)) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 3)) +
    facet_wrap(Community~ESV_ID, scales = "free", ncol = 9) +
    theme_classic() +
    theme(
        axis.text = element_text(size = 8, angle = 30, hjust = 1),
        axis.title = element_text(size = 15),
        strip.text = element_text(size = 8),
        panel.border = element_rect(color = 1, fill = NA),
        legend.position = "top"
    ) +
    guides(color = "none") +
    labs(x = expression(x[i]), y = expression(log(x[i+1]/x[i])))

ggsave(here::here("plots/FigS3-species_fitness_all_transfers_linear.png"), p, width = 12, height = 15)


# FigS4: ESVs that are ephemeral ----
ESV_ephemeral <- communities_abundance_fitness %>%
    drop_na() %>%
    mutate(CommunityESV = paste0(Community, ESV_ID)) %>%
    filter(!(CommunityESV %in% ESV_stable$CommunityESV)) %>%
    group_by(CommunityESV, Community, ESV_ID) %>%
    count() %>%
    filter(n>=5) %>%
    # Remove the artifact C10R4 Stenotrophomonas
    filter(CommunityESV != "C10R4Stenotrophomonas")

communities_abundance_fitness_ephemeral <- communities_abundance_fitness %>%
    drop_na() %>%
    mutate(CommunityESV = paste0(Community, ESV_ID)) %>%
    filter(CommunityESV %in% ESV_ephemeral$CommunityESV)

# Correlation
tb_cor_ephemeral <- communities_abundance_fitness_ephemeral %>%
    nest(data = c(-Community, -ESV_ID)) %>%
    mutate(fit = map(data, ~ cor_test_long(.x)),
           tidied = map(fit, tidy)) %>%
    unnest(tidied)

ESV_sig_ephemeral <- tb_cor_ephemeral %>%
    select(Community, ESV_ID, estimate, p.value) %>%
    filter(p.value < 0.05, estimate < 0) %>%
    mutate(CommunityESV = paste0(Community, ESV_ID)) %>%
    mutate(Significance = "p<0.05")

nrow(tb_cor_ephemeral) # 46 ephemeral ESVs
nrow(ESV_sig_ephemeral) # 10 ephemeral ESVs shows significant negative correlation
range(ESV_sig_ephemeral$estimate) # range of correlation [-1, -0.714]

# Linear model predicted equilibrium frequency
xintercept_ephemeral <- communities_abundance_fitness_ephemeral %>%
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
    left_join(select(tb_cor_ephemeral, Community, ESV_ID, estimate, p.value))


p <- communities_abundance_fitness_ephemeral %>%
    ggplot() +
    geom_smooth(data = communities_abundance_fitness_ephemeral %>% left_join(ESV_sig_ephemeral) %>% replace_na(list(Significance = "p>=0.05")),
                aes(x = Relative_Abundance, y = Fitness, color = Significance),
                method = stats::lm, formula = y ~ x, se = F) +
    geom_point(aes(x = Relative_Abundance, y = Fitness), shape = 21) +
    geom_hline(yintercept = 0, linetype = 2) +
    scale_color_manual(values = c("p<0.05" = "pink", "p>=0.05" = grey(0.8))) +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 3)) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 3)) +
    facet_wrap(Community~ESV_ID, scales = "free", ncol = 9) +
    theme_classic() +
    theme(
        axis.text = element_text(size = 8, angle = 30, hjust = 1),
        axis.title = element_text(size = 15),
        strip.text = element_text(size = 8),
        panel.border = element_rect(color = 1, fill = NA),
        legend.position = "top"
    ) +
    guides(color = "none") +
    labs(x = expression(x[i]), y = expression(log(x[i+1]/x[i])))

ggsave(here::here("plots/FigS4-species_fitness_all_transfers_linear_ephemeral.png"), p, width = 12, height = 9)


# Comparison of linear model predicted equilibrium frequency between stable and transient ESVs ----
xintercept <- bind_rows(xintercept_stable, xintercept_ephemeral)
table(xintercept$ESVType) # 99 stable ESVs and 46 transient ESVs

xintercept %>%
    filter(p.value < 0.05, estimate < 0) %>%
    pull(ESVType) %>%
    table() # For siginficant negaative correlation, 52 stable ESVs, 10 transient ESVs

p <- xintercept %>%
    filter(p.value < 0.05, estimate < 0) %>%
    ggplot() +
    geom_histogram(aes(x = Xintercept, fill = ESVType), color = 1, binwidth = 0.05) +
    scale_fill_manual(values = RColorBrewer::brewer.pal(3, "Set3"), label = c(stable = "stable ESV", transient = "transient ESV")) +
    scale_x_continuous(limits = c(-0.1,1.1), breaks = seq(0,1,.1), expand = c(0,0)) +
    scale_y_continuous(expand = c(0,0)) +
    theme_classic() +
    theme(
        legend.position = c(0.8, 0.8)
    ) +
    guides(fill = guide_legend(title = "")) +
    labs(x = "equilibirum frequency")

ggsave(here::here("plots/FigS3a-ESV_eq_freq.png"), p, width = 4, height = 3)



# Correlation plot: for all stable ESVs, color by siginficance ----
# x: the equilibrium frequency calculated as mean of T9-T12 vs.
# y: the equilibrium frequency predicted by the negative linear model

eq_freq_stable <- communities_eq_freq_stable %>%
    left_join(xintercept_stable) %>%
    drop_na(Xintercept)

eq_freq_ephemeral <- communities_eq_freq_ephemeral %>%
    left_join(xintercept_ephemeral) %>%
    drop_na(Xintercept)

#

p <- eq_freq_stable %>%
    mutate(ShowNegFreqDep = case_when(
        estimate < 0 & p.value < 0.05 ~ "neg freq dep",
        T ~ "no evidence of neg freq dep"
    )) %>%
    ggplot() +
    geom_abline(intercept = 0, slope = 1, linetype = 2, color = "red") +
    geom_point(aes(x = EquilibriumAbundance, y = Xintercept, color = ShowNegFreqDep), shape = 21, size = 2, stroke = 1) +
    scale_color_manual(values = c("neg freq dep" = "black", "no evidence of neg freq dep" = grey(0.8))) +
    scale_x_continuous(limits = c(0,1)) +
    scale_y_continuous(limits = c(0,1)) +
    theme_classic() +
    theme(
        panel.border = element_rect(color = 1, fill = NA, linewidth = 0.5),
        #legend.position = c(0.7, 0.2),
        legend.background = element_rect(color = NA, fill = NA)
    ) +
    guides() +
    labs(x = "mean of T9-T12", y = "predicetd by linear model")

ggsave(here::here("plots/FigS3b-ESV_eq_freq_predicted.png"), p, width = 6, height = 4)

cor.test(eq_freq_stable$EquilibriumAbundance, eq_freq_stable$Xintercept, method = "pearson") %>%
    tidy()


# Correlation plot: for all stable ESVs, color by communities used in this study ----
# x: the equilibrium frequency calculated as mean of T9-T12 vs.
# y: the equilibrium frequency predicted by the negative linear model

eq_freq_stable_comm <- eq_freq_stable %>%
    mutate(InThisStudy = case_when(
        Community %in% communities$Community ~ "four communities in current study",
        T ~ "other 22 communities"
    )) %>%
    arrange(desc(InThisStudy))

eq_freq_stable_comm %>%
    group_by(InThisStudy, Community) %>%
    count() %>%
    pull(InThisStudy) %>%
    table() # 4 communities in current study, 22 other communities

eq_freq_stable_comm %>%
    group_by(InThisStudy, Community) %>%
    count() # 18 red points. Or 18 ESVs in the current study


p <- eq_freq_stable_comm %>%
    ggplot() +
    geom_abline(intercept = 0, slope = 1, linetype = 2, color = "black") +
    geom_point(aes(x = EquilibriumAbundance, y = Xintercept, color = InThisStudy), shape = 21, size = 2, stroke = 1) +
    scale_color_manual(values = c("four communities in current study" = "red", "other 22 communities" = grey(0.8))) +
    scale_x_continuous(limits = c(0,1)) +
    scale_y_continuous(limits = c(0,1)) +
    theme_classic() +
    theme(
        panel.border = element_rect(color = 1, fill = NA, linewidth = 0.5),
        #legend.position = c(0.7, 0.2),
        legend.background = element_rect(color = NA, fill = NA)
    ) +
    guides() +
    labs(x = "mean of T9-T12", y = "predicetd by linear model")
# There are 16 red points because two points have negative predicted ESV eq freq from the linear model
ggsave(here::here("plots/FigS3c-ESV_eq_freq_predicted_comm.png"), p, width = 6, height = 4)

cor.test(eq_freq_stable$EquilibriumAbundance, eq_freq_stable$Xintercept, method = "pearson") %>%
    tidy()

table(eq_freq_stable_comm$Slope < 0) # 95 ESVs have slope <0, 4 ESVs have slope > 0
table(eq_freq_stable_comm$Xintercept > 0) # 94 ESVs have the predicted x intercept >0 , 5 ESVs have x intercept < 0

eq_freq_stable_comm %>%
    mutate(NegativeSlope = Slope < 0) %>%
    mutate(PositiveXintercept = Xintercept > 0) %>%
    group_by(NegativeSlope, PositiveXintercept) %>%
    count()

# For both stable and ephemeral ESVs, compare average fitness value vs. equilibirum freq ----
eq_freq_stable
xintercept





