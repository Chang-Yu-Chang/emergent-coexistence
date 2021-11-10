# Figure 3: pairs alike coexist more than pairs dissimilar
library(tidyverse)
library(tidymodels)
library(cowplot)
library(ggsci)

isolates <- read_csv(here::here("data/output/isolates.csv"))
pairs <- read_csv(here::here("data/output/pairs.csv")) %>% mutate(InteractionType = ifelse(InteractionType == "neutrality", "coexistence", InteractionType))
pairs_meta <- read_csv(here::here("data/output/pairs_meta.csv")) %>% mutate(InteractionType = ifelse(InteractionType == "neutrality", "coexistence", InteractionType))

interaction_type <- c("exclusion", "coexistence")
interaction_color = c("#DB7469", "#557BAA")
names(interaction_color) <- interaction_type

# Coexistence more likely in pairs alike

#
isolates %>%
    ggplot() +
    geom_histogram(aes(x = X_sum_16hr, fill = Fermenter), color = 1) +
    theme_classic()



# r_glu
p1 <- pairs_meta %>%
    filter(!is.na(PairFermenter)) %>%
    select(PairFermenter, InteractionType, r_glucose1, r_glucose2) %>%
    pivot_longer(cols = starts_with("r_glucose"), names_to = "Isolate", values_to = "r_glucose") %>%
    ggplot(aes(x = InteractionType, y = r_glucose, fill = Isolate)) +
    geom_boxplot() +
    geom_point(shape = 1, size = 2, position = position_jitterdodge(jitter.width = 0.2)) +
    scale_fill_npg(labels = c("dominant", "subdominant"), name = "Isolate") +
    facet_grid(.~PairFermenter, scales = "free_y", labeller = labeller(PairFermenter = c(FF = "Fermenter-Fermenter", FN = "Fermenter-Respirator", NN = "Respirator-Respirator"))) +
    theme_classic() +
    theme(legend.title = element_blank(), legend.position = "top") +
    labs(x = "", y = expression(r[glu]))
p1

# Amount of total acid secretion. sum_A, sum_B
p2 <- pairs_meta %>%
    select(InteractionType, PairFermenter, starts_with("X_sum") & !ends_with("d")) %>%
    pivot_longer(cols = c(-InteractionType, -PairFermenter), names_to = c("Time", "Isolate"), names_pattern = "X_sum_(.*)hr(.)", names_transform = list(Time = as.numeric), values_to = "X_sum") %>%
    filter(!is.na(PairFermenter)) %>%
    filter(Time == 16) %>%
    ggplot(aes(x = InteractionType, y = X_sum, fill = Isolate)) +
    geom_boxplot() +
    geom_point(shape = 1, size = 2, position = position_jitterdodge(jitter.width = 0.2)) +
    scale_fill_npg(labels = c("dominant", "subdominant"), name = "Isolate") +
    facet_grid(Time~PairFermenter, scales = "free_y", labeller = labeller(PairFermenter = c(FF = "Fermenter-Fermenter", FN = "Fermenter-Respirator", NN = "Respirator-Respirator"))) +
    theme_classic() +
    theme(legend.title = element_blank(), legend.position = "none") +
    labs(x = "", y = expression(X[sum]))
p2
# Scatter r_glu vs. sum_acids
p3 <- isolates %>%
    filter(!is.na(Fermenter)) %>%
    ggplot() +
    geom_segment(data = pairs_meta, aes(x = r_glucose1, xend = r_glucose2, y = X_sum_16hr1, yend = X_sum_16hr2, color = InteractionType, linetype = InteractionType)) +
    geom_point(aes(x = r_glucose, y = X_sum_16hr, shape = Fermenter), stroke = 1, size = 2) +
    scale_color_manual(values = c("coexistence" = "black", "exclusion" = "grey")) +
    scale_linetype_manual(values = c("coexistence" = 1, "exclusion" = 3)) +
    scale_shape_manual(values = c(`TRUE` = 1, `FALSE` = 2), labels = c("Fermenter", "Respirator")) +
    theme_classic() +
    theme(legend.title = element_blank(), legend.position = "right") +
    labs(x = expression(r[glu]), y = expression(2*acetate + 3*lactate + 4*succinate(mM)))
p3
#ggsave(here::here("plots/Fig_pairs-glu_acids.png"), p3, width = 6, height = 4)


# Scatter r_glu_d vs. sum_acids_d
p4 <- pairs_meta %>%
    filter(!is.na(PairFermenter)) %>%
    ggplot() +
    geom_vline(xintercept = 0, linetype = 2) +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_point(aes(x = r_glu_d, y = X_sum_16hr_d, shape = InteractionType, color = PairFermenter), size = 2, stroke = .5) +
    #scale_color_manual(values = interaction_color) +
    scale_shape_manual(values = c(coexistence = 16, exclusion = 1)) +
    scale_color_npg(labels = c("Fermenter-Fermenter", "Fermenter-Respirator", "Respirator-Respirator")) +
    theme_classic() +
    theme(legend.title = element_blank()) +
    labs(x = expression(r[A]-r[B]), y = expression(X[A]-X[B]))
p4
#ggsave(here::here("plots/Fig_pairs-glu_acids_d_scatter.png"), p4, width = 6, height = 4)


#
p <- plot_grid(p1, p2, align = "v", axis = "lr", ncol = 1)
#ggsave(here::here("plots/Fig3.png"), p, width = 8, height = 8)


# Stats
pairs_meta <- pairs_meta %>%
    mutate_if(is.character, as.factor)


pairs_meta %>%

    #filter(PairFermenter == "FF") %>% filter(!is.na(InteractionType)) %>%
    mutate(InteractionType = ifelse(InteractionType == "coexistence", 1, 0)) %>%
    glm(formula = InteractionType ~ r_glu_d * X_sum_16hr_d, data = ., family = "binomial") %>%
    broom::tidy() %>%
    #select(term, estimate, p.value) %>%
    {.}

logistic_reg() %>%
    set_engine("glm") %>%
    fit(InteractionType ~ r_glu_d * X_sum_16hr_d, data = pairs_meta) %>%
    tidy() %>%










if (FALSE) {
    ## Boxplot difference in the total amount of secretion
    p2 <- pairs_meta %>%
        select(Community, ID1, ID2, r_glu_d, InteractionType, PairFermenter, contains("_sum_")) %>%
        pivot_longer(cols = ends_with("hr_d"), names_to = "Time", values_to = "byproduct_sum_d") %>%
        mutate(Time = sub("byproduct_sum_", "", sub("hr_d", "", Time))) %>%
        filter(!is.na(Time), !is.na(PairFermenter)) %>%
        ggplot(aes(x = InteractionType, y = byproduct_sum_d, fill = InteractionType)) +
        geom_boxplot() +
        geom_point(aes(color = Isolate), shape = 1, size = 2, position = position_jitterdodge(jitter.width = 0.2)) +
        #geom_jitter(shape = 1, size = 2) +
        geom_hline(yintercept = 0, color = "red") +
        scale_fill_manual(values = c("coexistence"=grey(0.5), "exclusion"=grey(0.9))) +
        facet_grid(Time~PairFermenter, scales = "free_y",
                   labeller = labeller(PairFermenter = c(FF = "Fermenter-Fermenter", FN = "Fermenter-Respirator", NN = "Respirator-Respirator"))) +
        theme_classic() +
        theme(legend.title = element_blank()) +
        labs(x = "", y = expression(sum[A] - sum[B]))
    ggsave(here::here("plots/Fig_pairs-total_acids_boxplot.png"), p, width = 8, height = 8)

    ## Scatterplot, r_glu_d vs. sum_d
    p <- pairs_meta %>%
        filter(!is.na(PairFermenter)) %>%
        select(Community, ID1, ID2, r_glu_d, InteractionType, PairFermenter, contains("_sum_")) %>%
        pivot_longer(cols = ends_with("hr_d"), names_to = "Time", values_to = "byproduct_sum_d") %>%
        mutate(Time = sub("byproduct_sum_", "", sub("hr_d", "", Time))) %>%
        ggplot() +
        geom_vline(xintercept = 0, linetype = 2) +
        geom_hline(yintercept = 0, linetype = 2) +
        geom_point(aes(x = r_glu_d, y = byproduct_sum_d, shape = InteractionType, color = PairFermenter), size = 2) +
        scale_shape_manual(values = c("coexistence" = 16, "exclusion" = 1)) +
        scale_color_discrete(labels = c(FF = "Fermenter-Fermenter", FN = "Fermenter-Respirator", NN = "Respirator-Respirator")) +
        facet_grid(.~Time, labeller = labeller(Time = c(`16` = "16hr", `28` = "28hr", `48` = "48hr"))) +
        theme_classic() +
        theme(legend.title = element_blank(), legend.position = "top", legend.direction = "horizontal") +
        labs(x = expression(r[A_glu] - r[B_glu]), y = expression(sum[A]-sum[B]))
    ggsave(here::here("plots/Fig_pairs-total_acids_scatter.png"), p, width = 10, height = 4)


}

