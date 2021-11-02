# Figure 3: pairs alike coexist more than pairs dissimilar
library(tidyverse)
library(cowplot)
library(ggsci)

isolates <- read_csv(here::here("data/output/isolates.csv"))
pairs_meta <- read_csv(here::here("data/output/pairs_meta.csv"))

interaction_type <- c("exclusion", "coexistence")
interaction_color = c("#DB7469", "#557BAA")
names(interaction_color) <- interaction_type

# r_glu vs. sum_acids
pairs_glu <- pairs_meta %>%
    select(Community, ID1, ID2, PairFermenter, InteractionType, `D-Glucose1`, `D-Glucose2`, Fermenter1, Fermenter2) %>%
    pivot_longer(cols = starts_with("D-"), names_to = "Isolate", values_to = "r_glu") %>%
    mutate(Isolate = sub("D-Glucose", "", Isolate)) %>%
    mutate(Fermenter = ifelse(Isolate == 1, Fermenter1, Fermenter2)) %>%
    select(-Fermenter1, -Fermenter2)

pairs_meta_sum <- pairs_meta %>%
    select(Community, ID1, ID2, r_glu_d, InteractionType, PairFermenter, contains("ByproductSum")) %>%
    pivot_longer(cols = contains("Sum"), names_to = "Time", values_to = "byproduct_sum") %>%
    separate(col = Time, into = c("Time", "Isolate"), sep = "hr") %>%
    mutate(Time = sub("ByproductSum_", "", sub("", "", Time))) %>%
    filter(!is.na(PairFermenter), Time != 0)
pairs_sum <- pairs_meta_sum %>%  filter(Time == 16) %>% select(-r_glu_d)
pairs_glu_sum <- pairs_glu %>% left_join(pairs_sum)
pairs_glu_sum_w <- pairs_glu_sum %>% select(-Fermenter) %>%
    pivot_wider(names_from = Isolate, values_from = c(r_glu, byproduct_sum))


# r_glu
p1 <- pairs_meta %>%
    filter(!is.na(PairFermenter)) %>%
    select(PairFermenter, InteractionType, `D-Glucose1`, `D-Glucose2`) %>%
    pivot_longer(cols = starts_with("D-"), names_to = "Isolate", values_to = "r_glu") %>%
    ggplot(aes(x = InteractionType, y = r_glu, fill = Isolate)) +
    geom_boxplot() +
    geom_point(shape = 1, size = 2, position = position_jitterdodge(jitter.width = 0.2)) +
    scale_fill_npg(labels = c("dominant", "subdominant"), name = "Isolate") +
    facet_grid(.~PairFermenter, scales = "free_y", labeller = labeller(PairFermenter = c(FF = "Fermenter-Fermenter", FN = "Fermenter-Respirator", NN = "Respirator-Respirator"))) +
    theme_classic() +
    theme(legend.title = element_blank(), legend.position = "top") +
    labs(x = "", y = expression(r[glu]))

# Amount of total acid secretion. sum_A, sum_B
p2 <- pairs_meta_sum %>%
    filter(Time == 16) %>%
    ggplot(aes(x = InteractionType, y = byproduct_sum, fill = Isolate)) +
    geom_boxplot() +
    geom_point(shape = 1, size = 2, position = position_jitterdodge(jitter.width = 0.2)) +
    scale_fill_npg(labels = c("dominant", "subdominant"), name = "Isolate") +
    facet_grid(Time~PairFermenter, scales = "free_y", labeller = labeller(PairFermenter = c(FF = "Fermenter-Fermenter", FN = "Fermenter-Respirator", NN = "Respirator-Respirator"))) +
    theme_classic() +
    theme(legend.title = element_blank(), legend.position = "none") +
    labs(x = "", y = expression(2*acetate + 3*lactate + 4*succinate(mM)))

# Scatter r_glu vs. sum_acids ----
p3 <- pairs_glu_sum %>%
    filter(!is.na(Fermenter)) %>%
    ggplot() +
    geom_segment(data = pairs_glu_sum_w, aes(x = r_glu_1, xend = r_glu_2, y = byproduct_sum_1, yend = byproduct_sum_2, color = InteractionType, linetype = InteractionType)) +
    geom_point(aes(x = r_glu, y = byproduct_sum, shape = Fermenter), stroke = 1, size = 2) +
    #scale_color_manual(values = interaction_color) +
    scale_color_manual(values = c("coexistence" = "black", "exclusion" = "grey")) +
    scale_linetype_manual(values = c("coexistence" = 1, "exclusion" = 3)) +
    scale_shape_manual(values = c(`TRUE` = 1, `FALSE` = 2), labels = c("Fermenter", "Respirator")) +
    theme_classic() +
    theme(legend.title = element_blank(), legend.position = "right") +
    labs(x = expression(r[glu]), y = expression(2*acetate + 3*lactate + 4*succinate(mM)))

ggsave(here::here("plots/Fig_pairs-glu_acids.png"), p3, width = 6, height = 4)

p <- plot_grid(p1, p2, align = "v", axis = "lr", ncol = 1)
ggsave(here::here("plots/Fig3.png"), p, width = 8, height = 8)



#
pairs_glu_sum_w %>%
    filter(PairFermenter == "FN") %>%
    filter(!is.na(InteractionType)) %>%
    mutate(InteractionType = ifelse(InteractionType == "coexistence", 1, 0)) %>%
    mutate(r_glu_d = r_glu_1 - r_glu_2, byproduct_sum_d = byproduct_sum_1 - byproduct_sum_2) %>%
    glm(formula = InteractionType ~ r_glu_d * byproduct_sum_d, data = ., family = "binomial") %>%
    broom::tidy() %>%
    select(term, estimate, p.value)







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

