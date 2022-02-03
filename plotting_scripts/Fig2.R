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

# Figure 2A: diagram cartoon
#p_A <- ggdraw() + draw_image(here::here("plots/cartoons/Fig1A.png")) + theme(plot.background = element_rect(fill = "white", color = NA))
p_A <-  ggplot(mtcars, aes(x = wt, y = mpg)) + annotate("text", x = 0 , y = 0, label = "Cartoon for\nfermenters and respirators") + theme_void() + theme(plot.background = element_rect(fill = "white", color = NA))
ggsave(here::here("plots/Fig2A-functional_groups.png"), p_A, width = 2, height = 2)


# Figure 2B: Coexistence more likely in pairs alike
pairs_coexistence <- pairs %>%
    filter(!is.na(PairFermenter)) %>%
    filter(Assembly == "self_assembly")
pairs_count <- pairs_coexistence %>%
    group_by(PairFermenter) %>%
    summarize(Count = n())

p_B <- pairs_coexistence %>%
    ggplot() +
    geom_bar(aes(x = PairFermenter, fill = InteractionType), color = 1, position = position_fill()) +
    geom_text(data = pairs_count, aes(x = PairFermenter, y = 1, label = paste0("n=", Count)), vjust = -.5) +
    geom_rect(xmin = -Inf, xmax = Inf, ymin = 0, ymax = 1, color = grey(0.1), fill = NA, size = .5) +
    scale_fill_manual(values = assign_interaction_color(level = "simple")) +
    #scale_x_discrete(labels = c(FF = "Fermenter-\nFermenter", FN = "Fermenter-\nRespirator", NN = "Respirator-\nRespirator")) +
    scale_y_continuous(breaks = c(0,.5,1), limit = c(0, 1.15), expand = c(0,0)) +
    theme_classic() +
    theme(axis.title.x = element_blank(), axis.line.y = element_blank(),
          axis.text = element_text(size = 10, color = 1),
          legend.title = element_blank(), legend.position = "top") +
    guides(fill = F) +
    labs(y = "Fraction")
ggsave(here::here("plots/Fig2B-pairs_alike.png"), p_B, width = 3, height = 3)


# Figure 2C: Scatter r_glu_d vs. sum_acids_d
p_C <- pairs_meta %>%
    filter(!is.na(PairFermenter)) %>%
    # Use only FF pairs
    filter(PairFermenter == "FF") %>%
    ggplot() +
    geom_vline(xintercept = 0, linetype = 2) +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_point(aes(x = r_glucose_midhr_d, y = X_sum_16hr_d, color = InteractionType),
               shape = 21, size = 2, stroke = 1) +
    #scale_shape_manual(values = c(coexistence = 16, exclusion = 1)) +
    #scale_color_npg(labels = c("Fermenter-Fermenter", "Fermenter-Respirator", "Respirator-Respirator")) +
    scale_color_manual(values = assign_interaction_color(level = "simple")) +
    theme_classic() +
    theme(legend.title = element_blank(), legend.position = "top", axis.text = element_text(color = 1)) +
    labs(x = expression(r[A]-r[B]), y = expression(X[A]-X[B]))
ggsave(here::here("plots/Fig2C-r_glu_secretion_d.png"), p_C, width = 4, height = 3)

# Stats: does difference in r_glu and X_sum explain pairwise coexistence?
pairs_meta %>%
    mutate_if(is.character, as.factor) %>%
    filter(PairFermenter == "FF") %>%
    filter(!is.na(InteractionType)) %>%
    mutate(InteractionType = ifelse(InteractionType == "coexistence", 1, 0)) %>%
    glm(formula = InteractionType ~  r_glucose_midhr_d * X_sum_16hr_d, data = ., family = "binomial") %>%
    broom::tidy() %>%
    {.}


p <- plot_grid(p_A, p_B, p_C, nrow = 1, axis = "lrbt", align = "hv", rel_widths = c(1.5,2,2),
               labels = LETTERS[1:3], scale = .9) + theme(plot.background = element_rect(fill = "white", color = NA))
ggsave(here::here("plots/Fig2.png"), p, width = 8, height = 3)





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
    ggplot() +
    #geom_point(aes(x = Time, y = OD620, color = Replicate, group = Isolate), shape = 1) +
    geom_line(aes(x = Time, y = OD620, group = interaction(Isolate,Replicate), alpha = Rank), size = 1) +
    #scale_alpha_continuous(range = c(1,0.2), labels = c("")) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 3)) +
    facet_wrap(Community~., ncol = 4) +
    theme_classic() +
    theme(panel.border = element_rect(fill = NA, color = 1), legend.position = "none") +
    labs(x = "Time (hr)", y = "OD (620 nm)")

ggsave(here::here("plots/FigS8-growth_curves.png"), p_S8, width = 8, height = 8)



# Figure S9: r_glu_midhr per isolate
p_S9 <- pairs_meta %>%
    filter(!is.na(PairFermenter)) %>%
    select(PairFermenter, InteractionType, r_glucose_midhr1, r_glucose_midhr2) %>%
    pivot_longer(cols = starts_with("r_glucose_midhr"), names_to = "Isolate", values_to = "r_glucose_midhr") %>%
    filter(!is.na(r_glucose_midhr)) %>%
    ggplot(aes(x = InteractionType, y = r_glucose_midhr, fill = Isolate)) +
    geom_boxplot() +
    geom_point(shape = 1, size = 1, position = position_jitterdodge(jitter.width = 0.2)) +
    # p value
    ggpubr::stat_compare_means(aes(group = Isolate), label = "p.format", vjust = -.5, method = "t.test") +
    scale_fill_npg(labels = c("dominant", "subdominant"), name = "Isolate") +
    scale_y_continuous(expand = expansion(mult = .1, add = 0)) +
    facet_grid(.~PairFermenter, scales = "free_y", labeller = labeller(PairFermenter = c(FF = "Fermenter-Fermenter", FN = "Fermenter-Respirator", NN = "Respirator-Respirator"))) +
    theme_classic() +
    theme(legend.title = element_blank(), legend.position = "top",
          panel.border = element_rect(color = 1, fill = NA, size = .5)) +
    labs(x = "", y = expression(r[glu]))
ggsave(here::here("plots/FigS9-r_glu.png"), p_S9, width = 6, height = 4)


## Stats: does difference in r_glu explain pairwise coexistence?
pairs_meta %>%
    mutate_if(is.character, as.factor) %>%
    filter(!is.na(InteractionType)) %>%
    mutate(InteractionType = ifelse(InteractionType == "coexistence", 1, 0)) %>%
    glm(formula = InteractionType ~  r_glucose_midhr_d * PairFermenter, data = ., family = "binomial") %>%
    broom::tidy() %>%
    {.}


# Figure S10: Amount of total acid secretion. X_sum_16hr
p_S10 <- pairs_meta %>%
    select(InteractionType, PairFermenter, starts_with("X_sum") & !ends_with("d")) %>%
    pivot_longer(cols = c(-InteractionType, -PairFermenter), names_to = c("Time", "Isolate"), names_pattern = "X_sum_(.*)hr(.)", names_transform = list(Time = as.numeric), values_to = "X_sum") %>%
    filter(!is.na(PairFermenter)) %>%
    filter(Time == 16) %>%
    filter(!is.na(X_sum)) %>%
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
ggsave(here::here("plots/FigS10-secretion_total.png"), p_S10, width = 6, height = 4)

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
names(isolates)
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
p_S12 <- isolates %>%
    filter(!is.na(Fermenter), !is.na(leakiness_16hr)) %>%
    ggplot() +
    geom_boxplot(aes(x = Fermenter, y = leakiness_16hr, color = Fermenter), outlier.size = 2) +
    geom_jitter(aes(x = Fermenter, y = leakiness_16hr, color = Fermenter), shape = 1, size = 2, width = 0.3) +
    ggpubr::stat_compare_means(aes(group = Fermenter, x = Fermenter, y = leakiness_16hr)) +
    scale_x_discrete(labels = c("TRUE" = "Fermenter", "FALSE" = "Respirator")) +
    scale_y_continuous(limits = c(0, 0.65)) +
    scale_color_npg() +
    theme_classic() +
    theme(axis.title.x = element_blank(), legend.position = "none") +
    labs(y = "Leakiness")
ggsave(here::here("plots/FigS12-isolate_leakiness.png"), p_S12, width = 3, height = 3)




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


# Figure 2I: rmax_glu vs. number of wins
## raw
p1 <- isolates %>%
    ggplot(aes(x = r_glucose_maxhr, y = Score)) +
    geom_point(shape = 1, size = 2) +
    geom_smooth(method = "lm") +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 3)) +
    theme_classic() +
    labs(x = expression(r[glu]), y = "Competitive score") +
    ggpubr::stat_cor(label.x = .7, label.y = -10, method = "pearson")
## standardized
p2 <- isolates %>%
    ggplot(aes(x = r_glucose_maxhr, y = Score/Game)) +
    geom_point(shape = 1, size = 2) +
    geom_smooth(method = "lm") +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 3)) +
    theme_classic() +
    labs(x = expression(r[glu]), y = "Standardized score") +
    ggpubr::stat_cor(label.x = .7, label.y = -1, method = "pearson")

p_I <- plot_grid(p1, p2, labels = LETTERS[1:2], nrow = 1)
ggsave(here::here("plots/Fig2I-r_glu_score.png"), p_I, width = 8, height = 3.5)

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


















