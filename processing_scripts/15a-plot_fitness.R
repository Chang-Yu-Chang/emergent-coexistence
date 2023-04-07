#' This script checks and plot the fitness data
library(tidyverse)
library(cowplot)
library(broom)
source(here::here("analysis/00-metadata.R"))

factorize_communities <- function (x) x %>% mutate(Community = factor(Community, paste0("C", rep(1:12, each = 8), "R", rep(1:8, 12))))
communities <- read_csv(paste0(folder_data, "temp/00c-communities.csv"), show_col_types = F) %>% mutate(Community = factor(Community, Community))
fitness_stable <- read_csv(paste0(folder_data, "temp/15-fitness_stable.csv"), show_col_types = F) %>% factorize_communities
fitness_transient <- read_csv(paste0(folder_data, "temp/15-fitness_transient.csv"), show_col_types = F) %>% factorize_communities
fitness_transient2 <- read_csv(paste0(folder_data, "temp/15-fitness_transient2.csv"), show_col_types = F) %>% factorize_communities
eq_freq_stable <- read_csv(paste0(folder_data, "temp/15-eq_freq_stable.csv"), show_col_types = F) %>% factorize_communities
eq_freq_transient <- read_csv(paste0(folder_data, "temp/15-eq_freq_transient.csv"), show_col_types = F) %>% factorize_communities
eq_freq_transient2 <- read_csv(paste0(folder_data, "temp/15-eq_freq_transient2.csv"), show_col_types = F) %>% factorize_communities


fitness_stable <- fitness_stable %>% left_join(eq_freq_stable) %>% replace_na(list(Significance = "p>=0.05"))
fitness_transient <- fitness_transient %>% left_join(eq_freq_transient) %>% replace_na(list(Significance = "p>=0.05"))

# 1. Check abundance vs. invasion fitness, ESVs present in stable communities, linear fit ----
p <- fitness_stable %>%
    ggplot() +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_point(aes(x = Relative_Abundance, y = Fitness), shape = 21) +
    geom_smooth(data = fitness_stable, aes(x = Relative_Abundance, y = Fitness, color = Significance), method = stats::lm, formula = y ~ x, se = F) +
    # Linear model predicted
    geom_vline(data = eq_freq_stable, aes(xintercept = PredictedEqAbundance), color = "green", linetype = 2) +
    # Mean of T9-12
    geom_vline(data = eq_freq_stable, aes(xintercept = EmpiricalEqAbundance), color = "navyblue", linetype = 2) +
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

ggsave(paste0(folder_data, "temp/15a-01-species_fitness_stable.png"), p, width = 12, height = 15)

# 2. Check abundance vs. invasion fitness, ESVs present in stable communities, Polynomial fit, T1-12 ----
p <- fitness_stable %>%
    ggplot() +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_point(aes(x = Relative_Abundance, y = Fitness), shape = 21) +
    geom_smooth(data = fitness_stable, aes(x = Relative_Abundance, y = Fitness, color = Significance),
                method = stats::loess, span = 0.9, formula = y ~ x, se = F) +
    # Linear model predicted
    geom_vline(data = eq_freq_stable, aes(xintercept = PredictedEqAbundance), color = "green", linetype = 2) +
    # Mean of T9-12
    geom_vline(data = eq_freq_stable, aes(xintercept = EmpiricalEqAbundance), color = "navyblue", linetype = 2) +
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

ggsave(paste0(folder_data, "temp/15a-02-species_fitness_stable_logfit.png"), p, width = 12, height = 15)


# 3. Check abundance vs. invasion fitness, extinct ESVs, linear fit T1-12 ----
p <- fitness_transient %>%
    ggplot() +
    geom_point(aes(x = Relative_Abundance, y = Fitness), shape = 21) +
    geom_smooth(aes(x = Relative_Abundance, y = Fitness, color = Significance),
                method = stats::lm, formula = y ~ x, se = F) +
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

ggsave(paste0(folder_data, "temp/15a-03-species_fitness_extinct.png"), p, width = 12, height = 9)


# 4. Check the ESVs in the 13 communities -----
plot_comm_ESV <- function (fitness_stable, comm) {
    n_facets = 9
    fitness_stable_comm <- fitness_stable %>%
        left_join(eq_freq_stable) %>%
        replace_na(list(Significance = "p>=0.05")) %>%
        filter(Community %in% comm)
    eq_freq_stable_comm <- eq_freq_stable %>%
        filter(Community %in% comm)


    n_ESVs <- length(unique(fitness_stable_comm$ESV_ID))

    p1 <- fitness_stable_comm %>%
        ggplot() +
        geom_hline(yintercept = 0, linetype = 2) +
        geom_point(aes(x = Relative_Abundance, y = Fitness), shape = 21) +
        geom_smooth(data = fitness_stable_comm, aes(x = Relative_Abundance, y = Fitness, color = Significance), method = stats::lm, formula = y ~ x, se = F) +
        # Linear model predicted
        geom_vline(data = eq_freq_stable_comm, aes(xintercept = PredictedEqAbundance), color = "green", linetype = 2) +
        # Mean of T9-12
        geom_vline(data = eq_freq_stable_comm, aes(xintercept = EmpiricalEqAbundance), color = "navyblue", linetype = 2) +
        scale_color_manual(values = c("p<0.05" = "maroon", "p>=0.05" = grey(0.8))) +
        scale_x_continuous(breaks = c(0.001, 0.01, 0.1, 1), trans = "log10") +
        scale_y_continuous(breaks = scales::pretty_breaks(n = 3)) +
        facet_wrap(~ESV_ID, nrow = 1) +
        theme_classic() +
        theme(
            axis.text = element_text(size = 8, angle = 30, hjust = 1),
            axis.title = element_text(size = 15),
            axis.title.x = element_blank(),
            strip.text = element_text(size = 6),
            panel.border = element_rect(color = 1, fill = NA),
            legend.position = "top"
        ) +
        guides(color = "none") +
        labs(x = expression(x[i]), y = expression(log(x[i+1]/x[i]))) +
        ggtitle(comm)

    plot_grid(p1, NULL, rel_widths = c(n_ESVs, n_facets - n_ESVs))

}

p1 <- plot_comm_ESV(fitness_stable, "C1R4")
p2 <- plot_comm_ESV(fitness_stable, "C2R6")
p3 <- plot_comm_ESV(fitness_stable, "C2R8")
p4 <- plot_comm_ESV(fitness_stable, "C8R4")

p <- plot_grid(p1,p2,p3,p4 + theme(axis.title.x = element_text(size = 5)), ncol = 1) + paint_white_background()

ggsave(paste0(folder_data, "temp/15a-04-species_fitness_stable_comm_ordered.png"), p, width = 10, height = 6)


# 5. Histogram: linear model predicted equilibrium frequency of stable vs. transient ESVs ----
eq_freq <- bind_rows(eq_freq_stable, eq_freq_transient)
eq_freq %>%
    filter(rho < 0) %>%
    pull(ESVType) %>%
    table() # For siginficant negaative correlation, 52 stable ESVs, 10 transient ESVs

p <- eq_freq %>%
    filter(rho < 0) %>%
    ggplot() +
    geom_histogram(aes(x = PredictedEqAbundance, fill = ESVType), color = 1, binwidth = 0.05) +
    scale_fill_manual(values = RColorBrewer::brewer.pal(3, "Set3"), label = c(stable = "stable ESV", transient = "transient ESV")) +
    scale_x_continuous(limits = c(-0.1,1.1), breaks = seq(0,1,.1), expand = c(0,0)) +
    scale_y_continuous(expand = c(0,0)) +
    theme_classic() +
    theme(
        legend.position = c(0.8, 0.8)
    ) +
    guides(fill = guide_legend(title = "")) +
    labs(x = "equilibirum frequency")
ggsave(paste0(folder_data, "temp/15a-05-lm_eq_stable_vs_transient.png"), p, width = 4, height = 3)


# 6. Correlation plot: for all stable ESVs, mean of T9-12 vs. lm predicted eq freq. color by siginficance ----
# x: the equilibrium frequency calculated as mean of T9-T12 vs.
# y: the equilibrium frequency predicted by the negative linear model

p <- eq_freq_stable %>%
    mutate(ShowNegFreqDep = case_when(
        rho < 0 ~ "neg freq dep",
        T ~ "no evidence of neg freq dep"
    )) %>%
    ggplot() +
    geom_abline(intercept = 0, slope = 1, linetype = 2, color = "red") +
    geom_point(aes(x = EmpiricalEqAbundance, y = PredictedEqAbundance, color = ShowNegFreqDep), shape = 21, size = 2, stroke = 1) +
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

ggsave(paste0(folder_data, "temp/15a-06-ESV_eq_freq_predicted.png"), p, width = 6, height = 4)

cor.test(eq_freq_stable$EmpiricalEqAbundance, eq_freq_stable$PredictedEqAbundance, method = "pearson") %>%
    tidy()


# 7. Correlation plot: for all stable ESVs, mean of T9-12 vs. lm predicted eq freq. color by communities used in this study ----
# x: the equilibrium frequency calculated as mean of T9-T12 vs.
# y: the equilibrium frequency predicted by the negative linear model

eq_freq_stable_comm <- eq_freq_stable %>%
    mutate(InThisStudy = case_when(
        Community %in% communities$Community ~ "four communities in current study",
        T ~ "other 22 communities"
    )) %>%
    arrange(desc(InThisStudy)) # 99 ESVs

eq_freq_stable_comm_filtered <- eq_freq_stable_comm %>%
    # Negative slope
    filter(Slope < 0)  # 95 ESVs

eq_freq_stable_comm_filtered %>%
    group_by(InThisStudy, Community) %>%
    count() %>%
    pull(InThisStudy) %>%
    table() # 4 communities in current study, 22 other communities

eq_freq_stable_comm_filtered %>%
    filter(InThisStudy == "four communities in current study") %>%
    count() %>%
    sum() # 16 red points. 18 ESVs in the current study and two of them have positive correlation


p <- eq_freq_stable_comm_filtered %>%
    ggplot() +
    geom_abline(intercept = 0, slope = 1, linetype = 2, color = "black") +
    geom_point(aes(x = EmpiricalEqAbundance, y = PredictedEqAbundance, color = InThisStudy), shape = 21, size = 2, stroke = 1) +
    scale_color_manual(values = c("four communities in current study" = "red", "other 22 communities" = grey(0.8))) +
    scale_x_continuous(limits = c(0,1)) +
    scale_y_continuous(limits = c(0,1)) +
    theme_classic() +
    theme(
        panel.border = element_rect(color = 1, fill = NA, linewidth = 1),
        axis.text = element_text(color = 1),
        legend.position = c(0.68, 0.15),
        legend.title = element_blank(),
        legend.margin = margin(0,3,1,1, unit = "mm"),
        legend.box.margin = margin(0,0,0,0, unit = "mm"),
        legend.text = element_text(size = 6),
        legend.key.size = unit(5, "mm"),
        legend.spacing.x = unit(1, "mm"),
        legend.spacing.y = unit(2, "mm"),
        legend.background = element_rect(color = 1, fill = NA)
    ) +
    guides(color = guide_legend()) +
    labs(x = "empirical equilibrium abundance", y = "equilibrium abundance\npredicted from assembly dynamics")
# There are 16 red points because two points have negative predicted ESV eq freq from the linear model
ggsave(paste0(folder_data, "temp/15a-07-ESV_eq_freq_predicted_comm.png"), p, width = 4, height = 4)

cor.test(eq_freq_stable_comm_filtered$EmpiricalEqAbundance, eq_freq_stable_comm_filtered$PredictedEqAbundance, method = "pearson") %>%
    tidy()

table(eq_freq_stable_comm$Slope < 0) # 95 ESVs have slope <0, 4 ESVs have slope > 0
table(eq_freq_stable_comm$PredictedEqAbundance > 0) # 94 ESVs have the predicted x intercept >0 , 5 ESVs have x intercept < 0

eq_freq_stable_comm %>%
    mutate(NegativeSlope = Slope < 0) %>%
    mutate(PositiveXintercept = PredictedEqAbundance > 0) %>%
    group_by(NegativeSlope, PositiveXintercept) %>%
    count()

# 8. Compare average fitness value vs. equilibirum freq ----

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

p <- fitness_eq_freq2_filtered %>%
    arrange(ESVType) %>%
    ggplot() +
    geom_vline(xintercept = 0, linewidth = 0.1, linetype = 2) +
    geom_hline(yintercept = 0, linewidth = 0.1, linetype = 2) +
    geom_point(aes(x = MeanFitness, y = PredictedEqAbundance, color = ESVType), shape = 21, size = 2, stroke = .5) +
    #geom_segment(aes(x = MeanFitness - 2*SdFitness, xend = MeanFitness + 2*SdFitness, y = Xintercept, yend = Xintercept, color = ESVType)) +
    scale_color_manual(values = c(stable = "firebrick1", transient = grey(0.7)), label = c(stable = "stable ESV", transient = "transient ESV")) +
    theme_classic() +
    theme(
        panel.border = element_rect(color = 1, fill = NA, linewidth = 1)
    ) +
    guides(color = guide_legend(title = "", override.aes = list(stroke = 1))) +
    labs(x = "average fitness value", y = "equilibrium frequency\npredicted from assembly dynamics")
ggsave(paste0(folder_data, "temp/15a-08-stable_vs_transient.png"), p, width = 5, height = 4)


# 9. The number of data points in each ESV panel ----
p <- fitness_eq_freq2_filtered %>%
    mutate(NumberPoint = factor(NumberPoint, 3:11)) %>%
    mutate(ESVType = factor(ESVType)) %>%
    group_by(ESVType, NumberPoint, .drop = F) %>%
    summarize(Count = n()) %>%
    ggplot() +
    geom_col(aes(x = NumberPoint, y = Count, fill = ESVType), position = position_dodge()) +
    scale_fill_manual(values = RColorBrewer::brewer.pal(3, "Set1"), label = c(stable = "stable ESV", transient = "transient ESV")) +
    theme_classic() +
    theme() +
    guides() +
    labs()
ggsave(paste0(folder_data, "temp/15a-09-stable_vs_transient_samplesize.png"), p, width = 5, height = 4)

# 10. Distribution of fitness value with error bar----
p1 <- fitness_eq_freq2_filtered %>%
    arrange(ESVType) %>%
    mutate(CommunityESV = factor(CommunityESV)) %>%
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
        axis.text.x = element_text(size = 10),
        axis.text.y = element_blank()
        #axis.text.x = element_text(size = 6, angle = 45, hjust = 0)
    ) +
    guides() +
    labs(x = "", y = "invasion fitness")

# predicted equilibirum fitness
p2 <- fitness_eq_freq2_filtered %>%
    arrange(ESVType) %>%
    mutate(CommunityESV = factor(CommunityESV)) %>%
    ggplot() +
    geom_point(aes(x = CommunityESV, y = PredictedEqAbundance, color = ESVType), shape = 21, size = 2, stroke = 1) +
    #geom_segment(aes(x = CommunityESV, xend = CommunityESV, y = MeanFitness - SdFitness, yend = MeanFitness + SdFitness, color = ESVType)) +
    geom_hline(yintercept = 0, linewidth = 0.5, linetype = 2) +
    scale_color_manual(values = c(stable = "firebrick1", transient = grey(0.7)), label = c(stable = "stable ESV", transient = "transient ESV")) +
    coord_flip() +
    scale_x_discrete(position = "top") +
    theme_classic() +
    theme(
        panel.border = element_rect(fill = NA, color = 1, linewidth = 1),
        panel.grid.major.y = element_line(color = grey(0.9)),
        legend.position = "top",
        legend.title = element_blank(),
        legend.margin = margin(0,3,1,1, unit = "mm"),
        legend.box.margin = margin(0,0,0,0, unit = "mm"),
        legend.text = element_text(size = 15),
        legend.key.size = unit(10, "mm"),
        legend.spacing.x = unit(1, "mm"),
        legend.spacing.y = unit(2, "mm"),
        legend.background = element_rect(color = NA, fill = NA),
        axis.title = element_text(size = 15),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 8, hjust = 0)
        #axis.text.x = element_text(size = 6, angle = 45, hjust = 0)
    ) +
    guides() +
    labs(x = "ESV", y = "predicted equilibrium frequency")

p <- plot_grid(p1, p2, nrow = 1, align = "h", axis = "tb", rel_widths = c(1,1.2))
ggsave(paste0(folder_data, "temp/15a-10-stable_vs_transient_fitness_value.png"), p, width = 10, height = 15)
























































