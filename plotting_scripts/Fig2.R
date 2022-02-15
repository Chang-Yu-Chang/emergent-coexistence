# Figure 2: pairs alike coexist more than pairs dissimilar
library(tidyverse)
library(tidymodels)
library(tidygraph)
library(cowplot)
library(ggsci)
library(ggpubr)
library(ggtree)
source(here::here("plotting_scripts/network_functions.R"))

isolates <- read_csv(here::here("data/output/isolates.csv"))
pairs <- read_csv(here::here("data/output/pairs.csv")) %>% mutate(InteractionType = ifelse(InteractionType == "neutrality", "coexistence", InteractionType))
pairs_meta <- read_csv(here::here("data/output/pairs_meta.csv")) %>% mutate(InteractionType = ifelse(InteractionType == "neutrality", "coexistence", InteractionType))
communities <- read_csv(here::here("data/output/communities.csv"))
load(here::here("data/output/network_community.Rdata"))
fermenter_color <- c("fermenter" = "#8A89C0", "respirator" = "#FFCB77")
dominant_color <- c("dominant" = "grey20", "subdominant" = "grey90")

# Figure 2A: diagram cartoon
p_A <- ggdraw() + draw_image(here::here("plots/cartoons/Fig2A.png")) + theme(plot.background = element_rect(fill = "white", color = NA))
#p_A <-  ggplot(mtcars, aes(x = wt, y = mpg)) + annotate("text", x = 0 , y = 0, label = "Cartoon for\nfermenters and respirators") + theme_void() + theme(plot.background = element_rect(fill = "white", color = NA))
ggsave(here::here("plots/Fig2A-functional_groups.png"), p_A, width = 5, height = 5)


# Figure 2B: Coexistence more likely in pairs alike
pairs_coexistence <- pairs_meta %>%
    filter(!is.na(PairFermenter)) %>%
    filter(Assembly == "self_assembly") %>%
    mutate(PairConspecific = ifelse(PairFermenter == "FF" | PairFermenter == "NN", "conspecific", ifelse(PairFermenter == "FN", "heterospecific", NA)))
pairs_count <- pairs_coexistence %>%
    group_by(PairConspecific) %>%
    summarize(Count = n())

p_B <- pairs_coexistence %>%
    group_by(InteractionType, PairConspecific) %>%
    summarize(Count = n()) %>%
    group_by(PairConspecific) %>%
    mutate(Fraction = Count / sum(Count)) %>%
    #filter(InteractionType == "coexistence") %>%
    ggplot() +
    geom_col(aes(x = PairConspecific, y = Fraction, fill = factor(InteractionType, c("exclusion", "coexistence"))), color = 1, width = .7, position = position_dodge(width = 0.8)) +
    # percentage
    geom_text(aes(x = PairConspecific, y = Fraction, label = paste0(round(Fraction, 3) * 100,"%"), group = factor(InteractionType, c("exclusion", "coexistence"))), size = 3, vjust = -1, position = position_dodge(width = 0.8)) +
    # sample size
    geom_text(data = pairs_count, aes(x = PairConspecific, y = 1, label = paste0("n=", Count)), vjust = 1) +
    scale_fill_manual(values = assign_interaction_color(level = "simple")) +
    scale_y_continuous(breaks = c(0,.5,1), limit = c(0, 1), expand = c(0,0), labels = c("0%", "50%", "100%")) +
    theme_classic() +
    theme(axis.title.x = element_blank(), axis.title.y = element_text(size = 12),
          axis.ticks.x = element_blank(),
          axis.text = element_text(size = 12, color = 1),
          axis.text.x = element_text(size = 12, color = 1, angle = 30, vjust = 1, hjust = 1)
    ) +
    guides(fill = "none") +
    labs(y = "Percentage") +
    draw_image(here::here("plots/cartoons/Fig2B_FF.png"), x = 1, y = 0, scale = .6, vjust = .4, hjust = .5) +
    draw_image(here::here("plots/cartoons/Fig2B_FR.png"), x = 2, y = 0, scale = .6, vjust = .4, hjust = .5) +
    draw_image(here::here("plots/cartoons/Fig2B_RR.png"), x = 1, y = 0, scale = .6, vjust = .25, hjust = .5)

ggsave(here::here("plots/Fig2B-pairs_alike.png"), p_B, width = 3, height = 3)


## Stat: whether FF and RR pairs coexist more often than FR pairs
### observation
observed_indep_statistic <- pairs_coexistence %>%
    select(PairConspecific, InteractionType) %>%
    specify(InteractionType ~ PairConspecific, success = "coexistence") %>%
    calculate(stat = "Chisq", order = c("conspecific", "heterospecific"))
### null
null_distribution_simulated <- pairs_coexistence %>%
    select(PairConspecific, InteractionType) %>%
    specify(InteractionType ~ PairConspecific, success = "coexistence") %>%
    hypothesize(null = "independence") %>%
    generate(reps = 1000, type = "permute") %>%
    calculate(stat = "Chisq", order = c("conspecific", "heterospecific"))
### p
null_distribution_simulated %>% get_p_value(obs_stat = observed_indep_statistic, direction = "greater")



# Figure 2C: r_glu_midhr per isolate
p_C <- pairs_coexistence %>%
    filter(PairConspecific == "conspecific") %>%
    filter(!is.na(r_glucose_midhr_d)) %>%
    ggplot(aes(x = InteractionType, y = r_glucose_midhr_d, fill = InteractionType)) +
    geom_hline(yintercept = 0, linetype = 2, color = "grey") +
    geom_boxplot(width = .5) +
    geom_point(aes(group = InteractionType), shape = 1, size = 1, position = position_jitterdodge(jitter.width = 0.3)) +
    # p value
    ggpubr::stat_compare_means(aes(group = InteractionType), label = "p.format", vjust = 0, method = "t.test") +
    scale_fill_manual(values = assign_interaction_color()) +
    scale_y_continuous(expand = expansion(mult = .1, add = 0), limits = c(-0.1, 0.14), breaks = scales::pretty_breaks(n = 3)) +
    theme_classic() +
    theme(legend.title = element_blank(), legend.position = "top",
          axis.text = element_text(size = 12, color = 1),
          axis.title = element_text(size = 15, color = 1),
          axis.text.x = element_text(size = 12, color = 1, angle = 30, vjust = 1, hjust = 1),
          panel.border = element_rect(color = NA, fill = NA, size = 1)) +
    guides(alpha = "none", fill = "none", color = "none") +
    labs(x = "", y = expression(r[A]-r[B]))
p_C
ggsave(here::here("plots/Fig2C-r_glu.png"), p_C, width = 4, height = 4)

## Stats: among conspecific, whether the dominant has a higher r_glu than subdominant
temp <- pairs_coexistence %>%
    filter(PairConspecific == "conspecific") %>%
    tidyr::drop_na(r_glucose_midhr_d) %>%
    select(InteractionType, r_glucose_midhr_d)
### One sample
temp %>%
    filter(InteractionType == "coexistence") %>%
    t_test(response = r_glucose_midhr_d)
temp %>%
    filter(InteractionType == "exclusion") %>%
    t_test(response = r_glucose_midhr_d)
### CI
x <- filter(temp, InteractionType == "exclusion") %>% pull(r_glucose_midhr_d)
t.test(x)
mean(x) + qt(0.975, length(x) - 1) * sd(x) / sqrt(length(x))
### Two sample
temp %>%
    t_test(r_glucose_midhr_d ~ InteractionType, order = c("coexistence", "exclusion"))



# Figure 2D: Amount of total acid secretion. X_sum_16hr
p_D <- pairs_coexistence %>%
    filter(PairConspecific == "conspecific") %>%
    tidyr::drop_na(X_sum_16hr_d) %>%
    ggplot(aes(x = InteractionType, y = X_sum_16hr_d, fill = InteractionType)) +
    geom_hline(yintercept = 0, linetype = 2, color = "grey") +
    geom_boxplot(width = 0.5) +
    geom_point(aes(group = InteractionType), shape = 1, size = 1, position = position_jitterdodge(jitter.width = 0.3)) +
    # p value
    ggpubr::stat_compare_means(aes(group = InteractionType), label = "p.format", vjust = 0, method = "t.test") +
    scale_fill_manual(values = assign_interaction_color()) +
    scale_y_continuous(expand = expansion(mult = .1, add = 0), limits = c(-20, 20), breaks = scales::pretty_breaks(n = 3)) +
    theme_classic() +
    theme(legend.title = element_blank(), legend.position = "top",
          axis.text = element_text(size = 12, color = 1),
          axis.title = element_text(size = 15, color = 1),
          axis.text.x = element_text(size = 12, color = 1, angle = 30, vjust = 1, hjust = 1),
          panel.border = element_rect(color = NA, fill = NA, size = 1)) +
    guides(alpha = "none", fill = "none", color = "none") +
    labs(x = "", y = expression(X[A]-X[B]))

ggsave(here::here("plots/Fig2D-secretion_total.png"), p_D, width = 4, height = 4)

## Stats: whether the dominant has a higher X_sum than subdominant
temp <- pairs_coexistence %>%
    filter(PairConspecific == "conspecific") %>%
    tidyr::drop_na(X_sum_16hr_d) %>%
    select(InteractionType, X_sum_16hr_d)
### One sample
temp %>%
    filter(InteractionType == "coexistence") %>%
    t_test(response = X_sum_16hr_d)
temp %>%
    filter(InteractionType == "exclusion") %>%
    t_test(response = X_sum_16hr_d)
### two sample
temp %>% t_test(X_sum_16hr_d ~ InteractionType, order = c("coexistence", "exclusion"))



p <- plot_grid(p_A, p_B, p_C, p_D, nrow = 1, labels = LETTERS[1:4], rel_widths = c(1.5, 1.5, 1, 1), axis = "tb", align = "h", scale = .9) + theme(plot.background = element_rect(fill = "white", color = NA))
ggsave(here::here("plots/Fig2.png"), p, width = 10, height = 3)





# Figure S8: growth curves, alpha by ranks
isolates_curves <- isolates %>%
    filter(str_detect(Community, "C\\d")) %>%
    mutate(Isolate = factor(Isolate)) %>%
    # Standardize rank
    group_by(Community) %>%
    mutate(Rank = Rank / n()) %>%
    left_join(read_csv(here::here("data/output/isolates_curves.csv"))) %>%
    filter(CS == "glucose") %>%
    group_by(ID, CS, Time) %>%
    arrange(ID, CS, Time)

## Replicate dummy variable
temp <- isolates_curves %>%
    distinct(Community, Isolate, Well) %>%
    group_by(ID, CS, Well) %>%
    summarize() %>%
    group_by(ID, CS) %>%
    mutate(Replicate = factor(1:n()))

p_S8 <- isolates_curves %>%
    left_join(temp) %>%
    filter(Replicate == 1) %>%
    mutate(Community = factor(Community, communities$Community)) %>%
    drop_na(Fermenter) %>%
    mutate(Fermenter = ifelse(Fermenter, "fermenter", "respirator")) %>%
    ggplot() +
    geom_line(aes(x = Time, y = OD620, color = Fermenter, group = interaction(Isolate,Replicate)), size = 1) +
    scale_color_manual(values = fermenter_color) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 3)) +
    facet_wrap(Community~., ncol = 4) +
    guides(color = guide_legend(title = "")) +
    theme_classic() +
    theme(panel.border = element_rect(fill = NA, color = 1), legend.position = "top") +
    labs(x = "Time (hr)", y = "OD (620 nm)")

ggsave(here::here("plots/FigS8-growth_curves.png"), p_S8, width = 8, height = 8)



# Figure S9: r_glu_midhr per isolate
## isolate r_glu
p1 <- isolates %>%
    filter(Assembly == "self_assembly") %>%
    drop_na(Fermenter, r_glucose_midhr) %>%
    mutate(Fermenter = ifelse(Fermenter, "fermenter", "respirator")) %>%
    ggplot(aes(x = Fermenter, y = r_glucose_midhr, fill = Fermenter)) +
    geom_boxplot() +
    geom_jitter(shape = 1, size = 1, stroke = 1) +
    ggpubr::stat_compare_means(comparisons = list(c("fermenter", "respirator")), label = "p.signif", method = "t.test", ) +
    scale_fill_manual(values = fermenter_color) +
    scale_y_continuous(limits = c(0, 0.23), breaks = scales::pretty_breaks(n = 3)) +
    theme_classic() +
    theme(axis.title.x = element_blank()) +
    guides(fill = "none") +
    labs(y = expression(r))

## interaction outcome as a function of  difference in r_glu
temp <- pairs_meta %>%
    filter(Assembly == "self_assembly") %>%
    drop_na(InteractionType, PairFermenter, r_glucose_midhr_d) %>%
    mutate(Coexistence = ifelse(InteractionType == "coexistence", 1, ifelse(InteractionType == "exclusion", 0, NA)))
p2 <- temp %>%
    ggplot(aes(x = r_glucose_midhr_d, y = Coexistence)) +
    geom_point(shape = 1, size = 2) +
    geom_smooth(formula = y ~ x, method = "glm") +
    ggpubr::stat_cor(label.x = -.1, label.y = 1.1, method = "pearson") +
    scale_y_continuous(breaks = c(0,1), labels = c("exclusion", "coexistence")) +
    theme_classic() +
    theme(axis.title.y = element_blank()) +
    labs(x = expression(r[A]-r[B]))
## Stats: correlation
cor.test(temp$r_glucose_midhr_d, temp$Coexistence) %>% broom::tidy()


## isolate r_glu in pairFermenter and InteractionType
p3 <- pairs_meta %>%
    filter(Assembly == "self_assembly") %>%
    select(PairFermenter, InteractionType, r_glucose_midhr1, r_glucose_midhr2) %>%
    pivot_longer(cols = starts_with("r_glucose_midhr"), names_pattern = "r_glucose_midhr(.)", names_to = "Isolate", values_to = "r_glucose_midhr") %>%
    drop_na(PairFermenter, r_glucose_midhr) %>%
    mutate(Isolate = ifelse(Isolate == 1, "dominant", "subdominant")) %>%
    ggplot(aes(x = InteractionType, y = r_glucose_midhr, fill = InteractionType, alpha = Isolate)) +
    geom_boxplot() +
    geom_point(aes(group = Isolate), alpha = 1, shape = 1, size = 1, stroke = 1, position = position_jitterdodge(jitter.width = 0.2)) +
    # p value
    ggpubr::stat_compare_means(aes(group = Isolate), label = "p.signif", method = "t.test") +
    scale_alpha_manual(values = c("dominant" = 1, "subdominant" = .2)) +
    scale_fill_manual(values = assign_interaction_color()) +
    scale_y_continuous(expand = expansion(mult = .1, add = 0)) +
    facet_grid(.~PairFermenter, labeller = labeller(PairFermenter = c(FF = "Fermenter-Fermenter", FN = "Fermenter-Respirator", NN = "Respirator-Respirator"))) +
    theme_classic() +
    theme(legend.title = element_blank(), legend.position = "right",
          strip.background = element_blank(),
          panel.border = element_rect(color = 1, fill = NA, size = .5)) +
    guides(fill = "none", alpha = "none") +
    labs(x = "", y = expression(r))

#
p4 <- pairs_meta %>%
    filter(Assembly == "self_assembly") %>%
    drop_na(PairFermenter, r_glucose_midhr_d) %>%
    ggplot(aes(x = InteractionType, y = r_glucose_midhr_d, fill = InteractionType)) +
    geom_hline(yintercept = 0, linetype = 2, color = "grey") +
    geom_boxplot(width = .5) +
    geom_point(aes(group = InteractionType), shape = 1, size = 1, stroke = 1, position = position_jitterdodge(jitter.width = 0.3)) +
    # p value
    ggpubr::stat_compare_means(aes(group = InteractionType), label = "p.signif", method = "t.test") +
    scale_fill_manual(values = assign_interaction_color()) +
    scale_y_continuous(expand = expansion(mult = .1, add = 0), breaks = scales::pretty_breaks(n = 3)) +
    facet_grid(.~PairFermenter, labeller = labeller(PairFermenter = c(FF = "Fermenter-Fermenter", FN = "Fermenter-Respirator", NN = "Respirator-Respirator"))) +
    theme_classic() +
    theme(legend.title = element_blank(), legend.position = "right",
          strip.background = element_blank(),
          panel.border = element_rect(color = 1, fill = NA, size = .5)) +
    guides(alpha = "none", fill = "none", color = "none") +
    labs(x = "", y = expression(r[A]-r[B]))

## Stats
### FF
pairs_meta %>%
    filter(Assembly == "self_assembly") %>%
    drop_na(PairFermenter, r_glucose_midhr_d) %>%
    filter(PairFermenter == "FN") %>%
    select(InteractionType, r_glucose_midhr_d) %>%
    t_test(r_glucose_midhr_d ~ InteractionType, order = c("coexistence", "exclusion"))
### FR
pairs_meta %>%
    filter(Assembly == "self_assembly") %>%
    drop_na(PairFermenter, r_glucose_midhr_d) %>%
    filter(PairFermenter == "FN") %>%
    select(InteractionType, r_glucose_midhr_d) %>%
    t_test(r_glucose_midhr_d ~ InteractionType, order = c("coexistence", "exclusion"))
### RR
pairs_meta %>%
    filter(Assembly == "self_assembly") %>%
    drop_na(PairFermenter, r_glucose_midhr_d) %>%
    filter(PairFermenter == "NN") %>%
    select(InteractionType, r_glucose_midhr_d) %>%
    t_test(r_glucose_midhr_d ~ InteractionType, order = c("coexistence", "exclusion"))


p_upper <- plot_grid(p1, p2, nrow = 1, labels = c(LETTERS[1:2], ""), rel_widths = c(1,1.5), axis = "tb", align = "h", scale = .9)
p_S9 <- plot_grid(p_upper, p3, p4, ncol = 1, labels = c("", LETTERS[3:4]),
                  rel_heights = c(1,1,1)) + theme(plot.background = element_rect(fill = "white", color = NA))
ggsave(here::here("plots/FigS9-r_glu.png"), p_S9, width = 6, height = 8)


## Stats: does difference in r_glu explain pairwise coexistence?
if (FALSE) {
pairs_meta %>%
    mutate_if(is.character, as.factor) %>%
    filter(!is.na(InteractionType)) %>%
    mutate(InteractionType = ifelse(InteractionType == "coexistence", 1, 0)) %>%
    glm(formula = InteractionType ~  r_glucose_midhr_d * PairFermenter, data = ., family = "binomial") %>%
    broom::tidy() %>%
    {.}

}
temp <- pairs_meta %>%
    filter(Assembly == "self_assembly") %>%
    tidyr::drop_na(r_glucose_midhr_d) %>%
    oelect(PairFermenter, InteractionType, r_glucose_midhr_d)
### one-sample
temp %>%
    filter(PairFermenter == "FF") %>%
    filter(InteractionType == "coexistence") %>%
    t_test(response = r_glucose_midhr_d)
temp %>%
    filter(PairFermenter == "FF") %>%
    filter(InteractionType == "exclusion") %>%
    t_test(response = r_glucose_midhr_d)
temp %>%
    filter(PairFermenter == "FN") %>%
    t_test(response = r_glucose_midhr_d)
temp %>%
    filter(PairFermenter == "NN") %>%
    t_test(response = r_glucose_midhr_d)

### two-sample
temp %>%
    filter(PairFermenter == "FF") %>%
    t_test(r_glucose_midhr_d ~ InteractionType, order = c("coexistence", "exclusion"))
temp %>%
    filter(PairFermenter == "FN") %>%
    t_test(r_glucose_midhr_d ~ InteractionType, order = c("coexistence", "exclusion"))
temp %>%
    filter(PairFermenter == "NN") %>%
    t_test(r_glucose_midhr_d ~ InteractionType, order = c("coexistence", "exclusion"))


# Figure S10: rmid_glu vs. number of wins
temp <- isolates %>%
    filter(Assembly == "self_assembly") %>%
    drop_na(r_glucose_midhr)
## raw
p1 <- temp %>%
    ggplot(aes(x = r_glucose_midhr, y = Score)) +
    geom_point(shape = 1, size = 2) +
    geom_smooth(method = "lm") +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 3)) +
    theme_classic() +
    labs(x = expression(r[glu]), y = "Competitive score") +
    ggpubr::stat_cor(label.x = .1, label.y = -10, method = "pearson")
## standardized
p2 <- temp %>%
    ggplot(aes(x = r_glucose_midhr, y = Score/Game)) +
    geom_point(shape = 1, size = 2) +
    geom_smooth(method = "lm") +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 3)) +
    theme_classic() +
    labs(x = expression(r[glu]), y = "Standardized score") +
    ggpubr::stat_cor(label.x = .1, label.y = -1, method = "pearson")

p_S10 <- plot_grid(p1, p2, labels = LETTERS[1:2], nrow = 1, scale = .9) + theme(plot.background = element_rect(fill = "white", color = NA))
ggsave(here::here("plots/FigS10-r_glu_score.png"), p_S10, width = 8, height = 3.5)

## Stats: correlation
cor.test(temp$Score, temp$r_glucose_midhr, alternative = "two.sided", method = "pearson") %>% broom::tidy()
cor.test(temp$Score/temp$Game, temp$r_glucose_midhr, alternative = "two.sided", method = "pearson") %>% broom::tidy()


# Figure S11: isolate byproduct acids measure
## Three acids
isolates_byproduct <- isolates %>%
    filter(Assembly == "self_assembly") %>%
    select(Community, Isolate, ID, Fermenter, ends_with("hr")) %>%
    pivot_longer(cols = ends_with("hr"), names_pattern = "(.*)_(.*)hr", names_to = c("Measure", "Time"), values_to = "Value") %>%
    filter(!is.na(Value))
p1 <- isolates_byproduct %>%
    # Only use acids
    filter(str_detect(Measure, "X_"), !str_detect(Measure, "sum")) %>%
    mutate(Measure = str_replace(Measure, "X_", "")) %>%
    ggplot(aes(x = Time, y = Value, group = ID, color = Fermenter, alpha = ID)) +
    geom_point() +
    geom_line() +
    scale_color_npg(labels = c("TRUE" = "Fermenter", "FALSE" = "Respirator"), breaks = c(T, F)) +
    facet_wrap(Measure~., nrow = 1, scales = "free_y") +
    theme_classic() +
    theme(legend.title = element_blank(), legend.position = "top",
          strip.background = element_rect(color = NA, fill = NA),
          panel.border = element_rect(fill = NA, color = 1)) +
    guides(alpha = "none") +
    labs(x = "", y = "Concentration (mM)")
p_upper <- get_legend(p1)
p1 <- p1 + theme(legend.position = "none")

## Total acids production
p2 <- isolates %>%
    select(ID, Fermenter, starts_with("X_sum")) %>%
    pivot_longer(cols = c(-ID, -Fermenter), names_to = c("Measure", "Time"), names_pattern = "(.*)_(.*)hr", values_to = "Value") %>%
    mutate(Measure = ifelse(Measure == "X_sum", "total acids", Measure)) %>%
    ggplot(aes(x = Time, y = Value, group = ID, color = Fermenter, alpha = ID)) +
    geom_point() +
    geom_line() +
    scale_color_npg(labels = c("TRUE" = "Fermenter", "FALSE" = "Respirator"), breaks = c(T, F)) +
    facet_wrap(Measure~., nrow = 1, scales = "free_y") +
    theme_classic() +
    theme(legend.title = element_blank(), legend.position = "none",
          strip.background = element_rect(color = NA, fill = NA),
          panel.border = element_rect(fill = NA, color = 1)) +
    guides(alpha = "none") +
    labs(x = "", y = "Concentration (mM)")

## Leakiness
p3 <- isolates %>%
    select(ID, Fermenter, starts_with("leakiness")) %>%
    pivot_longer(cols = c(-ID, -Fermenter), names_to = c("Measure", "Time"), names_pattern = "(.*)_(.*)hr", values_to = "Value") %>%
    ggplot(aes(x = Time, y = Value, group = ID, color = Fermenter, alpha = ID)) +
    geom_point() +
    geom_line() +
    scale_color_npg(labels = c("TRUE" = "Fermenter", "FALSE" = "Respirator"), breaks = c(T, F)) +
    facet_wrap(Measure~., nrow = 1, scales = "free_y") +
    theme_classic() +
    theme(legend.title = element_blank(), legend.position = "none",
          strip.background = element_rect(color = NA, fill = NA),
          panel.border = element_rect(fill = NA, color = 1)) +
    guides(alpha = "none") +
    labs(x = "Time (hr)", y = "leakiness")

## pH
p4 <- isolates_byproduct %>%
    filter(Measure == "pH") %>%
    ggplot(aes(x = Time, y = Value, group = ID, color = Fermenter, alpha = ID)) +
    geom_point() +
    geom_line() +
    scale_color_npg(labels = c("TRUE" = "Fermenter", "FALSE" = "Respirator"), breaks = c(T, F)) +
    facet_wrap(Measure~., nrow = 1, scales = "free_y") +
    theme_classic() +
    theme(legend.title = element_blank(), legend.position = "none",
          strip.background = element_rect(color = NA, fill = NA),
          panel.border = element_rect(fill = NA, color = 1)) +
    guides(alpha = "none") +
    labs(x = "", y = "pH")

p_lower <- plot_grid(p2, p3, p4, nrow = 1, axis = "tb", align = "h")
p_S11 <- plot_grid(p_upper, p1, p_lower, ncol = 1, rel_heights = c(1, 5, 5)) + theme(plot.background = element_rect(fill = "white", color = NA))
ggsave(here::here("plots/FigS11-byproduct.png"), p_S11, width = 9, height = 6)


# Figure S12: leakiness
temp <- isolates %>%
    filter(Assembly == "self_assembly") %>%
    drop_na(Fermenter, leakiness_16hr)
p_S12 <- temp %>%
    ggplot() +
    geom_boxplot(aes(x = Fermenter, y = leakiness_16hr, color = Fermenter), outlier.size = 2) +
    geom_jitter(aes(x = Fermenter, y = leakiness_16hr, color = Fermenter), shape = 1, size = 2, width = 0.3) +
    ggpubr::stat_compare_means(aes(group = Fermenter, x = Fermenter, y = leakiness_16hr), method = "t.test") +
    scale_x_discrete(labels = c("TRUE" = "Fermenter", "FALSE" = "Respirator")) +
    scale_y_continuous(limits = c(0, 0.65)) +
    scale_color_npg() +
    theme_classic() +
    theme(axis.title.x = element_blank(), legend.position = "none") +
    labs(y = "Leakiness")
ggsave(here::here("plots/FigS12-isolate_leakiness.png"), p_S12, width = 3, height = 3)

## Stats: are respirators less leakier than fermenters?
temp %>% t_test(leakiness_16hr ~ Fermenter, order = c("TRUE", "FALSE"))








# Figure S13: Amount of total acid secretion. X_sum_16hr
## isolate X_sum
p1 <- isolates %>%
    filter(Assembly == "self_assembly") %>%
    drop_na(Fermenter, X_sum_16hr) %>%
    mutate(Fermenter = ifelse(Fermenter, "fermenter", "respirator")) %>%
    ggplot(aes(x = Fermenter, y = X_sum_16hr, fill = Fermenter)) +
    geom_boxplot() +
    geom_jitter(shape = 1, size = 1, stroke = 1) +
    ggpubr::stat_compare_means(comparisons = list(c("fermenter", "respirator")), label = "p.signif", method = "t.test", ) +
    scale_fill_manual(values = fermenter_color) +
    scale_y_continuous(limits = c(0, 45), breaks = scales::pretty_breaks(n = 3)) +
    theme_classic() +
    theme(axis.title.x = element_blank()) +
    guides(fill = "none") +
    labs(y = expression(X))

## interaction outcome as a function of difference in X_sum
temp <- pairs_meta %>%
    filter(Assembly == "self_assembly") %>%
    drop_na(InteractionType, PairFermenter, X_sum_16hr_d) %>%
    mutate(Coexistence = ifelse(InteractionType == "coexistence", 1, ifelse(InteractionType == "exclusion", 0, NA)))
p2 <- temp %>%
    ggplot(aes(x = X_sum_16hr_d, y = Coexistence)) +
    geom_point(shape = 1, size = 2) +
    geom_smooth(formula = y ~ x, method = "glm") +
    ggpubr::stat_cor(label.x = -20, label.y = 1.1, method = "pearson") +
    scale_y_continuous(breaks = c(0,1), labels = c("exclusion", "coexistence")) +
    theme_classic() +
    theme(axis.title.y = element_blank()) +
    labs(x = expression(X[A]-X[B]))
## Stats: correlation
cor.test(temp$X_sum_16hr_d, temp$Coexistence) %>% broom::tidy()


## isolate r_glu in pairFermenter and InteractionType
p3 <- pairs_meta %>%
    filter(Assembly == "self_assembly") %>%
    select(PairFermenter, InteractionType, X_sum_16hr1, X_sum_16hr2) %>%
    pivot_longer(cols = starts_with("X_sum_16hr"), names_pattern = "X_sum_16hr(.)", names_to = "Isolate", values_to = "X_sum_16hr") %>%
    drop_na(PairFermenter, X_sum_16hr) %>%
    mutate(Isolate = ifelse(Isolate == 1, "dominant", "subdominant")) %>%
    ggplot(aes(x = InteractionType, y = X_sum_16hr, fill = InteractionType, alpha = Isolate)) +
    geom_boxplot() +
    geom_point(aes(group = Isolate), alpha = 1, shape = 1, size = 1, stroke = 1, position = position_jitterdodge(jitter.width = 0.2)) +
    # p value
    ggpubr::stat_compare_means(aes(group = Isolate), label = "p.signif", method = "t.test") +
    scale_alpha_manual(values = c("dominant" = 1, "subdominant" = .2)) +
    scale_fill_manual(values = assign_interaction_color()) +
    scale_y_continuous(expand = expansion(mult = .1, add = 0)) +
    facet_grid(.~PairFermenter, labeller = labeller(PairFermenter = c(FF = "Fermenter-Fermenter", FN = "Fermenter-Respirator", NN = "Respirator-Respirator"))) +
    theme_classic() +
    theme(legend.title = element_blank(), legend.position = "right",
          strip.background = element_blank(),
          panel.border = element_rect(color = 1, fill = NA, size = .5)) +
    guides(fill = "none") +
    labs(x = "", y = expression(X))

#
p4 <- pairs_meta %>%
    filter(Assembly == "self_assembly") %>%
    drop_na(PairFermenter, X_sum_16hr_d) %>%
    ggplot(aes(x = InteractionType, y = X_sum_16hr_d, fill = InteractionType)) +
    geom_hline(yintercept = 0, linetype = 2, color = "grey") +
    geom_boxplot(width = .5) +
    geom_point(aes(group = InteractionType), shape = 1, size = 1, stroke = 1, position = position_jitterdodge(jitter.width = 0.3)) +
    # p value
    ggpubr::stat_compare_means(aes(group = InteractionType), label = "p.signif", method = "t.test") +
    scale_fill_manual(values = assign_interaction_color()) +
    scale_y_continuous(expand = expansion(mult = .1, add = 0), breaks = scales::pretty_breaks(n = 3)) +
    facet_grid(.~PairFermenter, labeller = labeller(PairFermenter = c(FF = "Fermenter-Fermenter", FN = "Fermenter-Respirator", NN = "Respirator-Respirator"))) +
    theme_classic() +
    theme(legend.title = element_blank(), legend.position = "right",
          strip.background = element_blank(),
          panel.border = element_rect(color = 1, fill = NA, size = .5)) +
    guides(alpha = "none", fill = "none", color = "none") +
    labs(x = "", y = expression(X[A]-X[B]))

## Stats
### FF
pairs_meta %>%
    filter(Assembly == "self_assembly") %>%
    drop_na(PairFermenter, X_sum_16hr_d) %>%
    filter(PairFermenter == "FN") %>%
    select(InteractionType, X_sum_16hr_d) %>%
    t_test(X_sum_16hr_d ~ InteractionType, order = c("coexistence", "exclusion"))
### FR
pairs_meta %>%
    filter(Assembly == "self_assembly") %>%
    drop_na(PairFermenter, X_sum_16hr_d) %>%
    filter(PairFermenter == "FN") %>%
    select(InteractionType, X_sum_16hr_d) %>%
    t_test(X_sum_16hr_d ~ InteractionType, order = c("coexistence", "exclusion"))
### RR
pairs_meta %>%
    filter(Assembly == "self_assembly") %>%
    drop_na(PairFermenter, X_sum_16hr_d) %>%
    filter(PairFermenter == "NN") %>%
    select(InteractionType, X_sum_16hr_d) %>%
    t_test(X_sum_16hr_d ~ InteractionType, order = c("coexistence", "exclusion"))


p_upper <- plot_grid(p1, p2, nrow = 1, labels = c(LETTERS[1:2], ""), rel_widths = c(1,1.5), axis = "tb", align = "h", scale = .9)
p_S13 <- plot_grid(p_upper, p3, p4, ncol = 1, labels = c("", LETTERS[3:4]),
                  rel_heights = c(1,1,1), axis = "lr", align = "v") + theme(plot.background = element_rect(fill = "white", color = NA))
ggsave(here::here("plots/FigS13-secretion_total.png"), p_S13, width = 8, height = 8)

if (FALSE) {

temp <- pairs_meta %>%
    select(InteractionType, PairFermenter, starts_with("X_sum_16") & !ends_with("d")) %>%
    pivot_longer(cols = c(-InteractionType, -PairFermenter), names_to = c("Time", "Isolate"), names_pattern = "X_sum_(.*)hr(.)", names_transform = list(Time = as.numeric), values_to = "X_sum") %>%
    drop_na(PairFermenter, X_sum) %>%
    mutate(Isolate = ifelse(Isolate == 1, "dominant", "subdominant"))
p_S13 <- temp %>%
    ggplot(aes(x = InteractionType, y = X_sum, fill = Isolate)) +
    geom_boxplot() +
    geom_point(shape = 1, size = 1, position = position_jitterdodge(jitter.width = 0.2)) +
    ggpubr::stat_compare_means(aes(group = Isolate), label = "p.format", vjust = -.5, method = "t.test") +
    scale_fill_npg(labels = c("dominant", "subdominant"), name = "Isolate") +
    scale_y_continuous(expand = expansion(mult = .1, add = 0)) +
    facet_grid(Time~PairFermenter, scales = "free_y", labeller = labeller(PairFermenter = c(FF = "Fermenter-Fermenter", FN = "Fermenter-Respirator", NN = "Respirator-Respirator"))) +
    theme_classic() +
    theme(legend.title = element_blank(), legend.position = "top",
          strip.text.y = element_blank(),
          panel.border = element_rect(color = 1, fill = NA, size = .5)) +
    labs(x = "", y = expression(X[sum]))
ggsave(here::here("plots/FigS13-secretion_total.png"), p_S13, width = 6, height = 4)

### one-sample
temp %>%
    filter(PairFermenter == "FF") %>%
    filter(InteractionType == "coexistence") %>%
    t_test(X_sum ~ Isolate, order = c("dominant", "subdominant"))
temp %>%
    filter(PairFermenter == "FF") %>%
    filter(InteractionType == "exclusion") %>%
    t_test(X_sum ~ Isolate, order = c("dominant", "subdominant"))
temp %>%
    filter(PairFermenter == "FN") %>%
    t_test(response = r_glucose_midhr_d)
temp %>%
    filter(PairFermenter == "NN") %>%
    t_test(response = r_glucose_midhr_d)
}


# Figure S14: Scatter r_glu_d vs. sum_acids_d
p_S14 <- pairs_coexistence %>%
    drop_na(PairFermenter, r_glucose_midhr_d, X_sum_16hr_d) %>%
    ggplot() +
    geom_vline(xintercept = 0, linetype = 2) +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_point(aes(x = r_glucose_midhr_d, y = X_sum_16hr_d, color = InteractionType), shape = 21, size = 2, stroke = 1) +
    scale_color_manual(values = assign_interaction_color(level = "simple")) +
    facet_wrap(.~PairConspecific, scale = "free") +
    theme_classic() +
    theme(legend.title = element_blank(), legend.position = "top",
          axis.text = element_text(color = 1),
          strip.background = element_blank(),
          panel.background = element_rect(color = 1, fill = NA)) +
    labs(x = expression(r[A]-r[B]), y = expression(X[A]-X[B]))
ggsave(here::here("plots/FigS14-r_glu_secretion_d.png"), p_S14, width = 5, height = 3)

# Stats: does difference in r_glu and X_sum explain pairwise coexistence?
pairs_meta %>%
    mutate_if(is.character, as.factor) %>%
    filter(PairFermenter == "FF") %>%
    filter(!is.na(InteractionType)) %>%
    mutate(InteractionType = ifelse(InteractionType == "coexistence", 1, 0)) %>%
    glm(formula = InteractionType ~  r_glucose_midhr_d * X_sum_16hr_d, data = ., family = "binomial") %>%
    broom::tidy() %>%
    {.}
#p_top <- plot_grid(p_A, p_B, nrow = 1, labels = LETTERS[1:2], rel_widths = c(1.5,1), scale = .9)
#p_bottom <- plot_grid(p_C, p_D, nrow = 1, labels = LETTERS[3:5], axis = "tb", align = "h", scale = .9)
#p <- plot_grid(p_top, p_bottom, nrow = 1, rel_widths = c(3,2), axis = "tb", align = "h") + theme(plot.background = element_rect(fill = "white", color = NA))



## Stats: does difference in X_sum explain pairwise coexistence?
pairs_meta %>%
    mutate_if(is.character, as.factor) %>%
    filter(!is.na(InteractionType)) %>%
    #filter(PairFermenter == "FF") %>%
    mutate(InteractionType = ifelse(InteractionType == "coexistence", 1, 0)) %>%
    glm(formula = InteractionType ~  X_sum_16hr_d * PairFermenter, data = ., family = "binomial") %>%
    broom::tidy() %>%
    {.}

## Stats: does difference in r_glu*X_sum explain pairwise coexistence?
pairs_meta %>%
    mutate_if(is.character, as.factor) %>%
    filter(!is.na(InteractionType)) %>%
    #filter(PairFermenter == "FF") %>%
    mutate(InteractionType = ifelse(InteractionType == "coexistence", 1, 0)) %>%
    glm(formula = InteractionType ~  r_glucose_midhr_d*X_sum_16hr_d * PairFermenter, data = ., family = "binomial") %>%
    broom::tidy() %>%
    {.}



"
finish the rest of figures and see how they fit in the narrative
"






















# Figure 2D: Scatter r_glu_midhr vs. sum_acids
p_D <- isolates %>%
    filter(!is.na(Fermenter)) %>%
    ggplot() +
    geom_segment(data = pairs_meta, aes(x = r_glucose_midhr1, xend = r_glucose_midhr2, y = X_sum_16hr1, yend = X_sum_16hr2, linetype = InteractionType)) +
    geom_point(aes(x = r_glucose_midhr, y = X_sum_16hr, color = Fermenter), stroke = 1.5, size = 2, shape = 1) +
    scale_color_npg(labels = c(`TRUE` = "fermenter", `FALSE` = "respirator")) +
    scale_linetype_manual(values = c("coexistence" = 1, "exclusion" = 2)) +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 3)) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 3)) +
    theme_classic() +
    theme(legend.title = element_blank(), legend.position = "right") +
    labs(x = expression(r[glu]), y = expression(2*acetate + 3*lactate + 4*succinate(mM)))
ggsave(here::here("plots/Fig2D-r_glu_secretion.png"), p_D, width = 6, height = 4)




# Figure 2G: r_glu vs. leakiness
p_G <- isolates %>%
    filter(!is.na(Fermenter), !is.na(leakiness_16hr)) %>%
    ggplot(aes(x = r_glucose_midhr, y = leakiness_16hr, color = Fermenter)) +
    geom_point(shape = 1, size = 2) +
    geom_smooth(method = "lm") +
    scale_color_npg(labels = c("TRUE" = "Fermenter", "FALSE" = "Respirator")) +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 3)) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 3)) +
    theme_classic() +
    theme(legend.title = element_blank(), legend.position = "top") +
    labs(x = expression(r[glu]), y = "Leakiness")
ggsave(here::here("plots/Fig2G-r_glu_leakiness.png"), p_G, width = 3, height = 3)

## Stats: do the fermenters follow pareto front?
isolates %>%
    filter(!is.na(Fermenter), !is.na(leakiness_16hr)) %>%
    glm(formula = leakiness_16hr ~  r_glucose_midhr * Fermenter, data = .) %>%
    broom::tidy() %>%
    {.}


# Figure 2H: r_acids
p_H <- pairs_meta %>%
    filter(!is.na(PairFermenter)) %>%
    mutate(Pair = 1:n()) %>%
    select(Pair, PairFermenter, InteractionType, starts_with("r_") & contains("midhr") & !ends_with("d")) %>%
    pivot_longer(cols = contains("hr"), names_to = c("CarbonSource", "Isolate"), names_pattern = "r_(.*)_midhr(.)", values_drop_na = T, values_to = "r") %>%
    ggplot(aes(x = InteractionType, y = r, fill = Isolate)) +
    geom_boxplot() +
    geom_point(shape = 1, size = 1, position = position_jitterdodge(jitter.width = 0.2)) +
    # p value
    ggpubr::stat_compare_means(aes(group = Isolate), label = "p.signif", vjust = -.5, method = "t.test") +
    scale_fill_npg(labels = c("dominant", "subdominant"), name = "Isolate") +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 4), expand = expansion(mult = .2, add = 0)) +
    facet_grid(CarbonSource~PairFermenter, scales = "free_y", labeller = labeller(PairFermenter = c(FF = "Fermenter-Fermenter", FN = "Fermenter-Respirator", NN = "Respirator-Respirator"))) +
    theme_classic() +
    theme(legend.title = element_blank(), legend.position = "right",
          panel.border = element_rect(color = 1, fill = NA, size = 1)) +
    labs(x = "", y = expression(r))
ggsave(here::here("plots/Fig2H-r_acids.png"), p_H, width = 8, height = 12)



# Figure 2J: matrix with species name
append_meta_data <- function(graph, community_name, isolates, pairs) {
    graph %>%
        activate(nodes) %>%
        mutate(Community = community_name) %>%
        left_join(isolates) %>%
        activate(edges) %>%
        mutate(Community = community_name, From = from, To = to) %>%
        left_join(pairs)
}
plot_adjacent_matrix <- function(graph) {
    graph_ranked <- graph %>%
        activate(nodes) %>%
        select(Isolate, PlotRank) %>%
        activate(edges) %>%
        mutate(
            fromRank = .N()$PlotRank[match(from, .N()$Isolate)],
            toRank = .N()$PlotRank[match(to, .N()$Isolate)])
    isolates_comm <- graph %>%
        igraph::vertex.attributes() %>%
        as_tibble()
    n_nodes <- nrow(isolates_comm)

    axis_label <- isolates_comm %>%
        #mutate(Discription = glue("{Genus}\twin = {Win}, draw = {Draw}, lose = {Lose}")) %>%
        #mutate(Discription = paste0(Genus, "\nwin = ", Win,", draw = ", Draw,", lose = ", Lose)) %>%
        mutate(Discription = Genus) %>%
        select(Isolate, PlotRank, Discription) %>%
        arrange(PlotRank) %>%
        pull(Discription)

    # Color lable fermenter and respirator
    color_logic <- isolates_comm %>%
        arrange(PlotRank) %>%
        pull(Fermenter) %>%
        ifelse("#C3423F", "#2978A0")

    graph_ranked %>%
        filter(fromRank <= toRank) %>%
        bind_edges(tibble(from = 1:n_nodes, to = 1:n_nodes, fromRank = 1:n_nodes, toRank = 1:n_nodes, InteractionType = "self")) %>%
        as_tibble() %>%
        ggplot() +
        geom_tile(aes(x = toRank, y = fromRank, fill = InteractionType), width = 0.9, height = 0.9) +
        scale_x_continuous(breaks = 1:n_nodes, labels = axis_label, expand = c(0,0), position = "top") +
        scale_y_reverse(breaks = 1:n_nodes, labels = axis_label, position = "right", expand = c(0,0)) +
        scale_fill_manual(values = assign_interaction_color(level = "matrix")) +
        theme_classic() +
        theme(legend.position = "none",
              axis.title = element_blank(),
              axis.line = element_blank(), axis.ticks = element_blank(),
              axis.text.x.top = element_text(size = 15, angle = 90, face = "bold.italic", color = color_logic),
              axis.text.y = element_text(size = 15, face = "bold.italic", color = color_logic),
              plot.margin = margin(10,10,10,10, "mm")) +
        labs()

}

communities_network <- communities %>%
    filter(str_detect(Community, "C\\d")) %>%
    mutate(graph = net_list[str_detect(names(net_list), "C\\d")]) %>%
    rowwise() %>%
    mutate(graph_meta = list(append_meta_data(graph, Community, isolates, pairs))) %>%
    mutate(p_net_matrix_list = list(plot_adjacent_matrix(graph_meta)))
p_J <- plot_grid(plotlist = communities_network$p_net_matrix_list, nrow = 4, labels = communities_network$Community) +
    theme(plot.background = element_rect(fill = "white", color = NA))
ggsave(here::here("plots/Fig2J-matrix_taxa.png"), p_J, width = 20, height = 20)

# Figure 2K: isolate byproduct acids measure
isolates_byproduct <- isolates %>%
    filter(Assembly == "self_assembly") %>%
    select(Community, Isolate, ID, Fermenter, ends_with("hr")) %>%
    pivot_longer(cols = ends_with("hr"), names_pattern = "(.*)_(.*)hr", names_to = c("Measure", "Time"), values_to = "Value") %>%
    filter(!is.na(Value))
p1 <- isolates_byproduct %>%
    # Only use acids
    filter(str_detect(Measure, "X_"), !str_detect(Measure, "sum")) %>%
    mutate(Measure = str_replace(Measure, "X_", "")) %>%
    ggplot(aes(x = Time, y = Value, group = ID, color = Fermenter, alpha = ID)) +
    geom_point() +
    geom_line() +
    scale_color_npg(labels = c("TRUE" = "Fermenter", "FALSE" = "Respirator"), breaks = c(T, F)) +
    facet_wrap(Measure~., nrow = 1, scales = "free_y") +
    theme_classic() +
    theme(legend.title = element_blank(), legend.position = "top",
          strip.background = element_rect(color = NA, fill = NA),
          panel.border = element_rect(fill = NA, color = 1)) +
    guides(alpha = "none") +
    labs(x = "Time (hr)", y = "Concentration (mM)")
p_upper <- get_legend(p1)
p1 <- p1 + theme(legend.position = "none")

p2 <- isolates_byproduct %>%
    filter(Measure == "pH") %>%
    ggplot(aes(x = Time, y = Value, group = ID, color = Fermenter, alpha = ID)) +
    geom_point() +
    geom_line() +
    scale_color_npg(labels = c("TRUE" = "Fermenter", "FALSE" = "Respirator"), breaks = c(T, F)) +
    facet_wrap(Measure~., nrow = 1, scales = "free_y") +
    theme_classic() +
    theme(legend.title = element_blank(), legend.position = "none",
          strip.background = element_rect(color = NA, fill = NA),
          panel.border = element_rect(fill = NA, color = 1)) +
    guides(alpha = "none") +
    labs(x = "Time (hr)", y = "pH")

p_lower <- plot_grid(p1, p2, nrow = 1, rel_widths = c(3,1), axis = "tb", align = "h")
p_K <- plot_grid(p_upper, p_lower, ncol = 1, rel_heights = c(1, 5)) + theme(plot.background = element_rect(fill = "white", color = NA))
ggsave(here::here("plots/Fig2K-byproduct.png"), p_K, width = 9, height = 3)


# Figure 2L: tree, all
load(here::here("data/output/tree.Rdata"))
p_L <- tree_meta %>%
    #ggtree(branch.length = "none", layout = "circular") +
    ggtree(branch.length = "none") +
    geom_tippoint(aes(colour = Fermenter), size = 3) +
    scale_color_npg(name = "", label = c("TRUE" = "fermenter", "FALSE" = "respirator")) +
    ggnewscale::new_scale_color() +
    geom_tiplab(aes(label = paste0(Community, " ", Genus, " sp.")), offset = 1) +
    scale_x_continuous(limits = c(0, 30)) +
    theme_void() +
    theme(plot.background = element_rect(fill = "white", color = NA)) +
    labs()
ggsave(here::here("plots/Fig2L-tree.png"), p_L, width = 10, height = 10)

# Figure 2M: tree for each community
plot_tree <- function(tree) {
    tree %>%
        #ggtree(branch.length = "none") +
        ggtree() +
        geom_tippoint(aes(colour = Fermenter), size = 5) +
        scale_color_npg(name = "", label = c("TRUE" = "Fermenter", "FALSE" = "Respirator")) +
        new_scale_color() +
        geom_tiplab(aes(label = paste0(Genus, " sp.")), offset = 0.05) +
        scale_x_continuous(expand = expansion(.5, .1)) +
        theme_void() +
        theme(legend.position = "none")
}
here::here("data/output/tree.Rdata")
## Get legend from fig 2L
l <- get_legend(p_L  + theme(legend.text = element_text(size = 20)))
p_M <- plot_grid(plotlist = c(communities_tree$tree_plot, list(l)), labels = communities_tree$Community, ncol = 2) +
    theme(plot.background = element_rect(fill = "white", color = NA))
ggsave(here::here("plots/Fig2M-tree_community.png"), p_M, width = 15, height = 15)


# Figure 2N: growth trait differences in a community
isolates_plot <- isolates %>% filter(Assembly == "self_assembly") %>% filter(!is.na(Fermenter))
pairs_plot <- pairs_meta %>% filter(Assembly == "self_assembly")
p_N <- isolates_plot %>%
    ggplot(aes(x = Rank, y = r_glucose_midhr)) +
    geom_segment(data = pairs_plot, aes(x = Rank1, xend = Rank2, y = r_glucose_midhr1, yend = r_glucose_midhr2, linetype = InteractionType)) +
    geom_point(aes(color = Fermenter), shape = 1, size = 3, stroke = 2) +
    scale_color_npg(name = "", label = c("TRUE" = "fermenter", "FALSE" = "respirator")) +
    scale_linetype_manual(name = "", values = c("coexistence" = 1, "exclusion" = 0)) +
    scale_x_continuous(breaks = 1:12) +
    facet_wrap(Community ~., scales = "free") +
    theme_classic() +
    theme(panel.border = element_rect(color = 1, fill = NA, size = 1)) +
    labs(x = "Rank", y = expression(r[glu]))
ggsave(here::here("plots/Fig2N-r_glu_rank.png"), p_N, width = 10, height = 8)



#



# Figure 2D: r_glu versus rank for C11R2
isolates_plot <- isolates %>% filter(Assembly == "self_assembly", Community == "C11R2") %>% filter(!is.na(Fermenter))
pairs_plot <- pairs_meta %>% filter(Assembly == "self_assembly", Community == "C11R2")
p_D <- isolates_plot %>%
    ggplot(aes(x = Rank, y = r_glucose_midhr)) +
    geom_segment(data = pairs_plot, aes(x = Rank1, xend = Rank2, y = r_glucose_midhr1, yend = r_glucose_midhr2, linetype = InteractionType)) +
    geom_point(aes(color = Fermenter), shape = 1, size = 3, stroke = 2) +
    scale_color_npg(name = "", label = c("TRUE" = "fermenter", "FALSE" = "respirator")) +
    scale_linetype_manual(name = "", values = c("coexistence" = 1, "exclusion" = 0)) +
    scale_x_continuous(breaks = 1:12) +
    theme_classic() +
    theme(panel.border = element_rect(color = 1, fill = NA, size = 1)) +
    labs(x = "Rank", y = expression(r[glu]))
ggsave(here::here("plots/Fig2D-r_glu_rank_C11R2.png"), p_D, width = 4, height = 3)


















