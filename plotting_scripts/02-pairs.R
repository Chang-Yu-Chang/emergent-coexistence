# Figure 2: pairs alike coexist more than pairs dissimilar

library(tidyverse)
library(tidymodels)
library(cowplot)
library(ggsci)
library(ggpubr)
library(ggtree)
library(tidygraph)
library(ggraph)
source(here::here("plotting_scripts/network_functions.R"))

isolates <- read_csv(here::here("data/output/isolates.csv"), col_types = cols())
communities <- read_csv(here::here("data/output/communities.csv"), col_types = cols())
pairs <- read_csv(here::here("data/output/pairs.csv"), col_types = cols()) %>% mutate(InteractionType = ifelse(InteractionType == "neutrality", "coexistence", InteractionType))
pairs_meta <- read_csv(here::here("data/output/pairs_meta.csv"), col_types = cols()) %>% mutate(InteractionType = ifelse(InteractionType == "neutrality", "coexistence", InteractionType))
load(here::here("data/output/network_community.Rdata"))
load("~/Dropbox/lab/invasion-network/data/output/network_randomized.Rdata") # Randomized networks


# Figure 3A: one example community of crossfeeding netorks
## growth rate and secretion
temp <- isolates %>%
    filter(Assembly == "self_assembly") %>%
    select(ID, Fermenter, starts_with("r_"), starts_with("X_")) %>%
    select(ID, Fermenter, ends_with("midhr"), ends_with("16hr")) %>%
    drop_na(r_glucose_16hr, X_acetate_16hr)
# Fermenter
temp2 <- isolates %>%
    mutate(Node = as.character(ID), Fermenter = ifelse(Fermenter, "fermenter", ifelse(Fermenter == FALSE, "respirator", NA))) %>%
    select(Node, Community, Fermenter)


isolates_in <- temp %>%
    select(ID, Fermenter, ends_with("midhr")) %>%
    pivot_longer(cols = starts_with("r_"), names_to = "CarbonSource", names_pattern = "r_(.+)_", values_to = "GrowthRate") %>%
    filter(!CarbonSource %in% c("ketogluconate", "gluconate")) %>%
    mutate(ID = as.character(ID)) %>%
    mutate(CarbonSource = factor(CarbonSource, c("glucose", "acetate", "lactate", "succinate"))) %>%
    arrange(Fermenter, ID, CarbonSource) %>%
    mutate(FermenterID = rep(1:(n()/4), each = 4))
isolates_out <- temp %>%
    select(ID, Fermenter, ends_with("16hr") & starts_with("X")) %>%
    pivot_longer(cols = starts_with("X_"), names_to = "CarbonSource", names_pattern = "X_(.+)_", values_to = "Secretion") %>%
    filter(CarbonSource != "sum", !CarbonSource %in% c("ketogluconate", "gluconate")) %>%
    mutate(CarbonSource = factor(CarbonSource, c("glucose", "acetate", "lactate", "succinate"))) %>%
    mutate(ID = as.character(ID)) %>%
    arrange(Fermenter, ID, CarbonSource) %>%
    mutate(FermenterID = rep(1:(n()/3), each = 3))

## Make the big bipartite network
nodes <- bind_rows(
    tibble(Node = unique(isolates_in$ID)) %>% mutate(Type = "consumer"),
    tibble(Node = unique(isolates_in$CarbonSource)) %>% mutate(Type = "resource")
) %>%
    left_join(temp2)
edges <- bind_rows(
    isolates_in %>% mutate(from = CarbonSource, to = ID, Strength = GrowthRate, Direction = "consumed", .keep = "unused"),
    isolates_out %>% mutate(from = ID, to = CarbonSource, Strength = Secretion, Direction = "secreted", .keep = "unused")
) %>%
    mutate(BinaryStrength = ifelse(Strength == 0, 0, 1)) %>%
    select(from, to, everything())

example_comm <- "C1R2"
nodes_example <- nodes %>% filter(Community == example_comm | Type == "resource")
edges_example <- edges %>% filter(from %in% nodes_example$Node & to %in% nodes_example$Node) %>%
    mutate(from = match(from, nodes_example$Node), to = match(to, nodes_example$Node))

g <- tbl_graph(nodes_example, edges_example) %>%
    activate(nodes) %>%
    mutate(showResourceID = ifelse(Type == "resource", Node, ""),
           showNode = ifelse(Type == "resource", F, T))

node_size = 5
pA <- g %>%
    activate(nodes) %>%
    group_by(Type) %>%
    mutate(x = 1:n(), y = ifelse(Type == "consumer", 1, 0)) %>%
    activate(edges) %>%
    filter(Strength != 0) %>%
    ggraph(layout = "nicely") +
    geom_node_point(size = node_size, shape = 21) +
    geom_node_text(aes(label = showResourceID), size = node_size/2) +
    #geom_edge_arc(aes(color = Direction), strength = 0.02) +
    geom_edge_arc(aes(color = Direction), strength = 0.02, width = 1, start_cap = circle(node_size/2+1, "mm"), end_cap = circle(node_size/2+1, "mm")) +
    #scale_color_manual(values = fermenter_color) +
    scale_fill_manual(values = c("consumer" = "grey50", "resource" = "grey90")) +
    scale_edge_color_manual(values = c("consumed" = "#557BAA", "secreted" = "#DB7469")) +
    scale_y_continuous(limits = c(0, 1.1)) +
    theme_void() +
    theme(legend.title = element_blank(), legend.position = "top") +
    #guides(width = "none") +
    labs() +
    paint_white_background()

ggsave(here::here("plots/Fig3A-example_crossfeeding.png"), pA, width = 5, height = 5)


# Figure 3B: d_r explains ranking
pairs_coexistence <- pairs_meta %>%
    filter(!is.na(PairFermenter)) %>%
    filter(Assembly == "self_assembly") %>%
    mutate(PairConspecific = ifelse(PairFermenter == "FF" | PairFermenter == "NN", "conspecific", ifelse(PairFermenter == "FN", "heterospecific", NA))) %>%
    mutate(Rank_d = Rank1 - Rank2, Score_d = Score1 - Score2)

pB <- pairs_coexistence %>%
    ggplot(aes(x = abs(r_glucose_midhr_d), y = abs(Rank_d))) +
    geom_point(size = 2, shape = 21, stroke = 1) +
    geom_smooth(method = "lm") +
    theme_classic() +
    labs(x = expression(r[1]-r[2]), y = expression(Rank[1]-Rank[2]))

# Figure 3C: d_X explaos ranking
pC <- pairs_coexistence %>%
    ggplot(aes(x = X_sum_16hr_d, y = Rank_d)) +
    geom_point(size = 2, shape = 21, stroke = 1) +
    geom_smooth(method = "lm") +
    theme_classic() +
    labs(x = expression(X[1]-X[2]), y = expression(Rank[1]-Rank[2]))

# Figure 3D: together
pD <- pairs_coexistence %>%
    drop_na(PairFermenter, r_glucose_midhr_d, X_sum_16hr_d) %>%
    ggplot() +
    geom_vline(xintercept = 0, linetype = 2) +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_point(aes(x = r_glucose_midhr_d, y = X_sum_16hr_d, color = InteractionType), shape = 21, size = 2, stroke = 1) +
    scale_color_manual(values = assign_interaction_color(level = "simple")) +
    #facet_wrap(.~PairConspecific, scale = "free") +
    theme_classic() +
    theme(legend.title = element_blank(), legend.position = "top",
          axis.text = element_text(color = 1),
          strip.background = element_blank(),
          panel.background = element_rect(color = 1, fill = NA)) +
    labs(x = expression(r[1]-r[2]), y = expression(X[1]-X[2]))

# Figure 3E: simulation
pE <- pairs_pool_meta %>%
    #filter(InteractionType != "no-growth") %>%
    ggplot() +
    geom_vline(xintercept = 0, linetype = 2) +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_point(aes(x = d_ConsumptionRate, y = d_CrossFeedingPotential, color = InteractionType), shape = 21, size = 2, stroke = .5) +
    scale_color_manual(values = c(assign_interaction_color())) +
    # facet_grid(.~PairConspecific) +
    theme_classic() +
    theme(legend.position = "top", strip.background = element_blank(), panel.background = element_rect(color = 1)) +
    guides(color = "none") +
    labs()


p <- plot_grid(pA, pB, pC, pD, pE, nrow = 2, labels = LETTERS[1:5], scale = c(.8, .9, .9, .9, .9),
               axis = "lrt", align = "h") + paint_white_background()


ggsave(here::here("plots/Fig3.png"), p, width = 9, height = 6)













#---------------------------------------------

# Figure 2S1: fermenter and respirator cartoon
pS1 <- ggdraw() + draw_image(here::here("plots/cartoons/Fig2A.png")) + theme(plot.background = element_rect(fill = "white", color = NA))
ggsave(here::here("plots/Fig2S1-functional_groups.png"), pS1, width = 5, height = 5)

# Figure 2S2:
pS2 <- pairs_coexistence %>%
    group_by(InteractionType, PairConspecific) %>%
    summarize(Count = n()) %>%
    group_by(PairConspecific) %>%
    mutate(Fraction = Count / sum(Count)) %>%
    #filter(InteractionType == "coexistence") %>%
    ggplot() +
    geom_col(aes(x = PairConspecific, y = Fraction, fill = factor(InteractionType, c("exclusion", "coexistence"))), color = 1, width = .7, position = "fill") +
    # percentage
    #geom_text(aes(x = PairConspecific, y = Fraction, label = paste0(round(Fraction, 3) * 100,"%"), group = factor(InteractionType, c("exclusion", "coexistence"))), size = 3, vjust = -1, position = position_dodge(width = 0.8)) +
    # sample size
    geom_text(data = pairs_count, aes(x = PairConspecific, y = 1, label = paste0("n=", Count)), vjust = 2) +
    scale_fill_manual(values = assign_interaction_color(level = "simple")) +
    scale_y_continuous(breaks = c(0,.5,1), limit = c(0, 1), expand = c(0,0), labels = c("0%", "50%", "100%")) +
    theme_classic() +
    theme(axis.title.x = element_blank(), axis.title.y = element_text(size = 12),
          axis.ticks.x = element_blank(),
          axis.text = element_text(size = 12, color = 1),
          axis.text.x = element_text(size = 12, color = 1, angle = 30, vjust = 1, hjust = 1)
    ) +
    guides(fill = "none") +
    labs(y = "Fracrion of coexistence") +
    draw_image(here::here("plots/cartoons/Fig2B_FF.png"), x = 1, y = 0, scale = .6, vjust = .4, hjust = .5) +
    draw_image(here::here("plots/cartoons/Fig2B_FR.png"), x = 2, y = 0, scale = .6, vjust = .4, hjust = .5) +
    draw_image(here::here("plots/cartoons/Fig2B_RR.png"), x = 1, y = 0, scale = .6, vjust = .25, hjust = .5)

ggsave(here::here("plots/Fig2S2-pairs_alike.png"), pS2, width = 3, height = 3)

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


# Figure 2S3: r_glu_midhr per isolate
pS3 <- pairs_coexistence %>%
    filter(PairConspecific == "conspecific") %>%
    filter(!is.na(r_glucose_midhr_d)) %>%
    mutate(InteractionType = factor(InteractionType, c("exclusion", "coexistence"))) %>%
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

ggsave(here::here("plots/Fig2S3-r_glu.png"), pS3, width = 4, height = 4)

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



# Figure 2S4: Amount of total acid secretion. X_sum_16hr
pS4 <- pairs_coexistence %>%
    filter(PairConspecific == "conspecific") %>%
    drop_na(X_sum_16hr_d) %>%
    mutate(InteractionType = factor(InteractionType, c("exclusion", "coexistence"))) %>%
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

ggsave(here::here("plots/Fig2S4-secretion_total.png"), pS4, width = 4, height = 4)

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



# Figure 2S5: growth curves, alpha by ranks
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

pS5 <- isolates_curves %>%
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

ggsave(here::here("plots/Fig2S5-growth_curves.png"), pS5, width = 8, height = 8)



# Figure 2S6: r_glu_midhr per isolate
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
pS6 <- plot_grid(p_upper, p3, p4, ncol = 1, labels = c("", LETTERS[3:4]),
                  rel_heights = c(1,1,1)) + theme(plot.background = element_rect(fill = "white", color = NA))
ggsave(here::here("plots/Fig2S6-r_glu.png"), pS6, width = 6, height = 8)


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


# Figure 2S7: rmid_glu vs. number of wins
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
## r1-r2 versus rank1 - rank2
p3 <- pairs_meta %>%
    mutate(Rank_d = Rank1 - Rank2) %>%
    ggplot(aes(x = r_glucose_midhr_d, y = Rank_d)) +
    geom_vline(xintercept = 0, linetype = 2, color = 1) +
    geom_point(aes(color = InteractionType), shape = 21, size = 2, stroke = 1) +
    geom_smooth(method = "lm") +
    scale_color_manual(values = assign_interaction_color()) +
    theme_classic() +
    theme(legend.position = "top", legend.title = element_blank()) +
    labs(x = expression(r[1]-r[2]), y = expression(Rank[1]-Rank[2]))

## rank1 - rank2 vs. r1 - r2
p4 <- pairs_meta %>%
    mutate(Rank_d = Rank1 - Rank2) %>%
    ggplot(aes(x = Rank_d, y = r_glucose_midhr_d)) +
    geom_vline(xintercept = 0, linetype = 2, color = 1) +
    geom_point(aes(color = InteractionType), shape = 21, size = 2, stroke = 1) +
    geom_smooth(method = "lm") +
    scale_color_manual(values = assign_interaction_color()) +
    #scale_x_continuous() +
    theme_classic() +
    theme(legend.position = "top", legend.title = element_blank()) +
    labs(x = expression(Rank[1]-Rank[2]), y = expression(r[1]-r[2]))



pS7 <- plot_grid(p1, p2, p3, p4, labels = LETTERS[1:4], nrow = 2, scale = .9) + theme(plot.background = element_rect(fill = "white", color = NA))
ggsave(here::here("plots/Fig2S7-r_glu_score.png"), pS7, width = 8, height = 8)

## Stats: correlation
cor.test(temp$Score, temp$r_glucose_midhr, alternative = "two.sided", method = "pearson") %>% broom::tidy()
cor.test(temp$Score/temp$Game, temp$r_glucose_midhr, alternative = "two.sided", method = "pearson") %>% broom::tidy()


# Figure 2S8: isolate byproduct acids measure
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
pS8 <- plot_grid(p_upper, p1, p_lower, ncol = 1, rel_heights = c(1, 5, 5)) + theme(plot.background = element_rect(fill = "white", color = NA))
ggsave(here::here("plots/Fig2S8-byproduct.png"), pS8, width = 9, height = 6)


# Figure 2S9: leakiness
temp <- isolates %>%
    filter(Assembly == "self_assembly") %>%
    drop_na(Fermenter, leakiness_16hr)
pS9 <- temp %>%
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
ggsave(here::here("plots/Fig2S9-isolate_leakiness.png"), pS9, width = 3, height = 3)

## Stats: are respirators less leakier than fermenters?
temp %>% t_test(leakiness_16hr ~ Fermenter, order = c("TRUE", "FALSE"))


# Figure 2S10: Amount of total acid secretion. X_sum_16hr
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
pS10 <- plot_grid(p_upper, p3, p4, ncol = 1, labels = c("", LETTERS[3:4]),
                  rel_heights = c(1,1,1), axis = "lr", align = "v") + theme(plot.background = element_rect(fill = "white", color = NA))
ggsave(here::here("plots/Fig2S10-secretion_total.png"), pS10, width = 8, height = 8)


# Figure 2S11: Scatter r_glu_d vs. sum_acids_d
pS11 <- pairs_coexistence %>%
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
ggsave(here::here("plots/Fig2S11-r_glu_secretion_d.png"), pS11, width = 5, height = 3)

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
    filter(PairFermenter == "FF") %>%
    mutate(InteractionType = ifelse(InteractionType == "coexistence", 1, 0)) %>%
    glm(formula = InteractionType ~  r_glucose_midhr_d*X_sum_16hr_d * PairFermenter, data = ., family = "binomial") %>%
    broom::tidy() %>%
    {.}


# Figure 2S12: Scatter r_glu_midhr vs. sum_acids
pS12 <- isolates %>%
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
ggsave(here::here("plots/Fig2S12-r_glu_secretion.png"), pS12, width = 6, height = 4)




# Figure S13: r_glu vs. leakiness
pS13 <- isolates %>%
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
ggsave(here::here("plots/Fig2S13-r_glu_leakiness.png"), pS13, width = 3, height = 3)

## Stats: do the fermenters follow pareto front?
isolates %>%
    filter(!is.na(Fermenter), !is.na(leakiness_16hr)) %>%
    glm(formula = leakiness_16hr ~  r_glucose_midhr * Fermenter, data = .) %>%
    broom::tidy() %>%
    {.}


# Figure 2S14: r_acids
pS14 <- pairs_meta %>%
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
ggsave(here::here("plots/Fig2S14-r_acids.png"), pS14, width = 8, height = 12)



# Figure 2S15: matrix with species name
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
pS15 <- plot_grid(plotlist = communities_network$p_net_matrix_list, nrow = 4, labels = communities_network$Community) +
    theme(plot.background = element_rect(fill = "white", color = NA))
ggsave(here::here("plots/Fig2S15-matrix_taxa.png"), pS15, width = 20, height = 20)


# Figure 2S16: isolate byproduct acids measure
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
pS16 <- plot_grid(p_upper, p_lower, ncol = 1, rel_heights = c(1, 5)) + theme(plot.background = element_rect(fill = "white", color = NA))
ggsave(here::here("plots/Fig2S16-byproduct.png"), pS16, width = 9, height = 3)


# Figure 2S17: tree, all
"deprecated. fix it"
load(here::here("data/output/tree.Rdata"))
pS17 <- tree_meta %>%
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
ggsave(here::here("plots/Fig2S17-tree.png"), pS17, width = 10, height = 10)

# Figure 2S18: tree for each community
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
## Get legend
l <- get_legend(p_L  + theme(legend.text = element_text(size = 20)))
pS18 <- plot_grid(plotlist = c(communities_tree$tree_plot, list(l)), labels = communities_tree$Community, ncol = 2) +
    theme(plot.background = element_rect(fill = "white", color = NA))
ggsave(here::here("plots/Fig2S18-tree_community.png"), pS18, width = 15, height = 15)


# Figure 2S19: growth trait differences in a community
isolates_plot <- isolates %>% filter(Assembly == "self_assembly") %>% filter(!is.na(Fermenter))
pairs_plot <- pairs_meta %>% filter(Assembly == "self_assembly")
pS19 <- isolates_plot %>%
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
ggsave(here::here("plots/Fig2S19-r_glu_rank.png"), pS19, width = 10, height = 8)


# Figure 2S20: r_glu versus rank for C11R2
isolates_plot <- isolates %>% filter(Assembly == "self_assembly", Community == "C11R2") %>% filter(!is.na(Fermenter))
pairs_plot <- pairs_meta %>% filter(Assembly == "self_assembly", Community == "C11R2")
pS20 <- isolates_plot %>%
    ggplot(aes(x = Rank, y = r_glucose_midhr)) +
    geom_segment(data = pairs_plot, aes(x = Rank1, xend = Rank2, y = r_glucose_midhr1, yend = r_glucose_midhr2, linetype = InteractionType)) +
    geom_point(aes(color = Fermenter), shape = 1, size = 3, stroke = 2) +
    scale_color_npg(name = "", label = c("TRUE" = "fermenter", "FALSE" = "respirator")) +
    scale_linetype_manual(name = "", values = c("coexistence" = 1, "exclusion" = 0)) +
    scale_x_continuous(breaks = 1:12) +
    theme_classic() +
    theme(panel.border = element_rect(color = 1, fill = NA, size = 1)) +
    labs(x = "Rank", y = expression(r[glu]))
ggsave(here::here("plots/Fig2S20-r_glu_rank_C11R2.png"), pS20, width = 4, height = 3)


# Figure 2S21. Genus level metabolic traits
temp <- isolates %>%
    drop_na(Genus) %>%
    filter(Family %in% c("Enterobacteriaceae", "Pseudomonadaceae")) %>%
    arrange(Family) %>%
    mutate(Genus = factor(Genus, levels = unique(Genus)))

p1 <- temp %>%
    ggplot(aes(x = Genus, y = r_glucose_midhr, color = Family)) +
    geom_boxplot() +
    geom_jitter() +
    theme_classic() +
    labs()

p2 <- temp %>%
    ggplot(aes(x = Genus, y = X_sum_16hr, color = Family)) +
    geom_boxplot() +
    geom_jitter() +
    theme_classic() +
    labs()

pS21 <- plot_grid(p1, p2, ncol = 1, axis = "lr", align = "v", labels = c("A", "B"))
ggsave(here::here("plots/Fig2S21-isolate_trait.png"), pS21, width = 8, height = 6)


# Figure 2S22. Pairwise coexistence at genus resolution
temp <- pairs_meta %>%
    filter(Assembly == "self_assembly") %>%
    filter(PairFermenter == "FF") %>%
    filter(Family1 == Family2) %>%
    #select(InteractionType, Family1, Family2, PairFermenter, Genus1, Genus2) %>%
    drop_na(Genus1, Genus2) %>%
    mutate(InteractionType = factor(InteractionType, c("exclusion", "coexistence"))) %>%
    mutate(PairGenusConspecific = ifelse(Genus1 == Genus2, "conspecific", "heterospecific")) %>%
    unite(col = "PairGenus", Genus1, Genus2)

pS22 <- temp %>%
    group_by(PairGenusConspecific, InteractionType) %>%
    summarize(Count = n()) %>%
    mutate(TotalCount = sum(Count)) %>%
    ggplot() +
    #geom_col(aes(x = PairGenusConspecific, y = Count, fill = InteractionType)) +
    geom_col(aes(x = PairGenusConspecific, y = Count, fill = InteractionType), position = "fill") +
    geom_text(aes(x = PairGenusConspecific, label = paste0("n=", TotalCount)), y = 1, vjust = 2) +
    scale_fill_manual(values = assign_interaction_color()) +
    theme_classic() +
    theme(legend.position = "top", legend.title = element_blank()) +
    labs(x = "Within-fermenter pairs", y = "Fraction")

ggsave(here::here("plots/Fig2S22-isolate_trait.png"), pS22, width = 3, height = 3)


# Figure 2S23: cross-feeding networks
## Using the cross-feeding assay from Sylvie
## Using the secretion profile from Sylvie
temp <- isolates %>%
    filter(Assembly == "self_assembly") %>%
    select(ID, Fermenter, starts_with("r_"), starts_with("X_")) %>%
    select(ID, Fermenter, ends_with("midhr"), ends_with("16hr")) %>%
    drop_na(r_glucose_16hr, X_acetate_16hr)
temp2 <- isolates %>%
    mutate(Node = as.character(ID), Fermenter = ifelse(Fermenter, "fermenter", ifelse(Fermenter == FALSE, "respirator", NA)), .keep = "used") %>%
    select(Node, Fermenter)

isolates_in <- temp %>%
    select(ID, Fermenter, ends_with("midhr")) %>%
    pivot_longer(cols = starts_with("r_"), names_to = "CarbonSource", names_pattern = "r_(.+)_", values_to = "GrowthRate") %>%
    filter(!CarbonSource %in% c("ketogluconate", "gluconate")) %>%
    mutate(ID = as.character(ID)) %>%
    mutate(CarbonSource = factor(CarbonSource, c("glucose", "acetate", "lactate", "succinate"))) %>%
    arrange(Fermenter, CarbonSource)
isolates_out <- temp %>%
    select(ID, Fermenter, ends_with("16hr") & starts_with("X")) %>%
    pivot_longer(cols = starts_with("X_"), names_to = "CarbonSource", names_pattern = "X_(.+)_", values_to = "Secretion") %>%
    filter(CarbonSource != "sum", !CarbonSource %in% c("ketogluconate", "gluconate")) %>%
    mutate(CarbonSource = factor(CarbonSource, c("glucose", "acetate", "lactate", "succinate"))) %>%
    mutate(ID = as.character(ID)) %>%
    arrange(Fermenter, CarbonSource)
nodes <- bind_rows(tibble(Node = unique(isolates_in$ID)) %>% mutate(Type = "consumer", y = 1, x = 1:n()),
                   tibble(Node = unique(isolates_in$CarbonSource)) %>% mutate(Type = "resource", y = 0, x = seq(10, 50, length = n()))) %>%
    left_join(temp2)
edges <- bind_rows(isolates_in %>% mutate(from = CarbonSource, to = ID, Strength = GrowthRate, Direction = "consumed", .keep = "unused"),
                   isolates_out %>% mutate(from = ID, to = CarbonSource, Strength = Secretion, Direction = "secreted", .keep = "unused")) %>%
    mutate(BinaryStrength = ifelse(Strength == 0, 0, 1))

p1 <- edges %>%
    ggplot() +
    geom_histogram(aes(x = Strength, fill = Direction), alpha = .5, position = "dodge") +
    facet_grid(.~Direction, scales = "free_x") +
    theme_classic() +
    theme(panel.border = element_rect(color = 1, fill = NA)) +
    guides(fill = "none") +
    labs()

g <- tbl_graph(nodes, edges) %>%
    activate(nodes) %>%
    mutate(showResourceID = ifelse(Type == "resource", Node, "o"))
node_size = 2
p2 <- g %>%
    activate(edges) %>%
    filter(Direction == "consumed") %>%
    filter(Strength != 0) %>%
    ggraph(layout = "nicely") +
    geom_node_text(aes(label = showResourceID, color = Fermenter), size = 5) +
    geom_edge_link(aes(color = Direction), width = .2, start_cap = circle(node_size/2+1, "mm"), end_cap = circle(node_size/2+1, "mm")) +
    annotate(geom = "text", label = "respirator", x = 10, y = 1, color = "#FFCB77", vjust = -1, size = 5) +
    annotate(geom = "text", label = "fermenter", x = 40, y = 1, color = "#8A89C0", vjust = -1, size = 5) +
    scale_color_manual(values = fermenter_color) +
    scale_fill_manual(values = c("consumer" = "grey50", "resource" = "grey90")) +
    scale_edge_color_manual(values = c("consumed" = "#557BAA", "secreted" = "#DB7469")) +
    scale_y_continuous(limits = c(0, 1.1)) +
    theme_void() +
    theme(legend.title = element_blank(), legend.position = "none") +
    guides(width = "none") +
    labs() +
    ggtitle("consumption network")

p3 <- g %>%
    activate(edges) %>%
    filter(Direction == "secreted") %>%
    filter(Strength != 0) %>%
    ggraph(layout = "nicely") +
    geom_node_text(aes(label = showResourceID, color = Fermenter), size = 5) +
    geom_edge_link(aes(color = Direction), width = .2, start_cap = circle(node_size/2+1, "mm"), end_cap = circle(node_size/2+1, "mm")) +
    annotate(geom = "text", label = "respirator", x = 10, y = 1, color = "#FFCB77", vjust = -1, size = 5) +
    annotate(geom = "text", label = "fermenter", x = 40, y = 1, color = "#8A89C0", vjust = -1, size = 5) +
    scale_color_manual(values = fermenter_color) +
    scale_fill_manual(values = c("consumer" = "grey50", "resource" = "grey90")) +
    scale_edge_color_manual(values = c("consumed" = "#557BAA", "secreted" = "#DB7469")) +
    scale_y_continuous(limits = c(0, 1.1)) +
    theme_void() +
    theme(legend.title = element_blank(), legend.position = "none") +
    guides(width = "none") +
    labs() +
    ggtitle("secretion network")


pS23 <- plot_grid(p1, p2, p3, ncol = 1, scale = .9, labels = c("A", "", "")) + paint_white_background()
ggsave(here::here("plots/Fig2S23-crossfeeding_network.png"), pS23, width = 8, height = 10)



# Figure 2S24: cross-feeding networks in heatmaps
temp <- isolates %>%
    filter(Assembly == "self_assembly") %>%
    select(ID, Fermenter, starts_with("r_"), starts_with("X_")) %>%
    select(ID, Fermenter, ends_with("midhr"), ends_with("16hr")) %>%
    drop_na(r_glucose_16hr, X_acetate_16hr)
temp2 <- isolates %>%
    mutate(Node = as.character(ID), Fermenter = ifelse(Fermenter, "fermenter", ifelse(Fermenter == FALSE, "respirator", NA)), .keep = "used") %>%
    select(Node, Fermenter)

isolates_in <- temp %>%
    select(ID, Fermenter, ends_with("midhr")) %>%
    pivot_longer(cols = starts_with("r_"), names_to = "CarbonSource", names_pattern = "r_(.+)_", values_to = "GrowthRate") %>%
    filter(!CarbonSource %in% c("ketogluconate", "gluconate")) %>%
    mutate(ID = as.character(ID)) %>%
    mutate(CarbonSource = factor(CarbonSource, c("glucose", "acetate", "lactate", "succinate"))) %>%
    arrange(Fermenter, ID, CarbonSource) %>%
    mutate(FermenterID = rep(1:(n()/4), each = 4))
isolates_out <- temp %>%
    select(ID, Fermenter, ends_with("16hr") & starts_with("X")) %>%
    pivot_longer(cols = starts_with("X_"), names_to = "CarbonSource", names_pattern = "X_(.+)_", values_to = "Secretion") %>%
    filter(CarbonSource != "sum", !CarbonSource %in% c("ketogluconate", "gluconate")) %>%
    mutate(CarbonSource = factor(CarbonSource, c("glucose", "acetate", "lactate", "succinate"))) %>%
    mutate(ID = as.character(ID)) %>%
    arrange(Fermenter, ID, CarbonSource) %>%
    mutate(FermenterID = rep(1:(n()/3), each = 3))

# Growth rate
p1 <- isolates_in %>%
    ggplot() +
    geom_tile(aes(x = CarbonSource, y = FermenterID, fill = GrowthRate), color = 1) +
    scale_fill_gradient2() +
    scale_y_continuous(limits = c(0, 60), breaks = 1:59) +
    theme_classic() +
    theme(legend.position = "top") +
    labs()

# Secretion
p2 <- isolates_out %>%
    ggplot() +
    geom_tile(aes(x = CarbonSource, y = FermenterID, fill = Secretion), color = 1) +
    scale_fill_gradient2() +
    scale_y_continuous(limits = c(0, 60), breaks = 1:59) +
    theme_classic() +
    theme(legend.position = "top", axis.title.y = element_blank()) +
    labs()

pS24 <- plot_grid(p1, p2, ncol = 2, scale = .9, labels = c("", "")) + paint_white_background()
ggsave(here::here("plots/Fig2S24-cross_feeding.png"), pS24, width = 8, height = 10)





