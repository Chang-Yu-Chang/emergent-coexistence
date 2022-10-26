#' Script for supplement figures

library(tidyverse)
library(cowplot)
library(RColorBrewer)
library(ggsci)
library(officer)
library(flextable)
source(here::here("analysis/00-metadata.R"))

isolates <- read_csv(paste0(folder_data, "output/isolates.csv"), show_col_types = F)
pairs <- read_csv(paste0(folder_data, "output/pairs.csv"), show_col_types = F)
pairs_example_outcomes_finer <- read_csv(paste0(folder_data, "output/pairs_example_outcomes_finer.csv"), show_col_types = F)
pairs_example_outcomes <- read_csv(paste0(folder_data, "output/pairs_example_outcomes.csv"), show_col_types = F)
communities <- read_csv(paste0(folder_data, "temp/00c-communities.csv"), show_col_types = F)
accuracy <- read_csv(paste0(folder_data, "temp/91-accuracy.csv"), show_col_types = F)
pairs_freq <- read_csv(paste0(folder_data, "temp/93-pairs_freq.csv"), show_col_types = F)
load(paste0(folder_data, "temp/94-communities_network.Rdata"))
communities_hierarchy <- read_csv(paste0(folder_data, "temp/94-communities_hierarchy.csv"), show_col_types = F)


# Figure SXX isolate abundance in community ----
## Panel A. cartoon of self-assembly experiment and isolate abundance
p1 <- ggdraw() + draw_image(here::here("plots/cartoons/FigS1.png")) + paint_white_background()

## Panel B. isolate abundance
color_sets <- tibble(Color = c("yellow", "deepskyblue3", "blue", "darkorchid2", "firebrick", "orange2", "grey"),
                     Family = c("Aeromonadaceae", "Enterobacteriaceae", "Moraxellaceae", "Pseudomonadaceae","Comamonadaceae","Alcaligenaceae", "Sphingobacteriaceae"))
p2 <- isolates %>%
    mutate(Community = factor(Community, communities$Community)) %>%
    filter(!is.na(RelativeAbundance)) %>%
    arrange(Community) %>%
    ggplot() +
    geom_bar(aes(x = Community, y = RelativeAbundance, fill = Family), size = .3, color = "grey30", position = "stack", stat = "identity") +
    theme_bw() +
    scale_fill_manual(values = setNames(color_sets$Color, color_sets$Family)) +
    scale_x_discrete(labels = 1:13) +
    scale_y_continuous(breaks = c(0, .5, 1), expand = c(0,0), limits = c(0, 1)) +
    theme(
        panel.grid = element_blank(),
        axis.title = element_text(size = 12),
        panel.border = element_rect(color = 1, size = 1)) +
    labs(y = "Relative abundance")

## Stats
isolates %>%
    group_by(Community) %>%
    summarize(Total = sum(RelativeAbundance, na.rm = T)) %>%
    summarize(Mean = mean(Total))

p <- plot_grid(p1, p2, nrow = 1, scale = c(.8, .9), rel_widths = c(1, 1.5), labels = c("A", "B")) + theme(plot.background = element_rect(color = NA, fill = "white"))
ggsave(here::here("plots/FigS-isolate_abundance.png"), p, width = 10, height = 3)


# Figure SXX machine vs. human ----
pairs_T8_combined <- read_csv(paste0(folder_data, "temp/92-pairs_T8_combined.csv"), show_col_types = F)

p1 <- pairs_T8_combined %>%
    ggplot() +
    geom_abline(slope = 1, intercept = 0, color = "red", linetype = 2) +
    geom_point(aes(x = TotalCount_human, y = TotalCount_machine), shape = 21, size = 2) +
    scale_x_log10() +
    scale_y_log10() +
    theme_classic() +
    labs(x = "Segmentation CFU count", y = "Manual CFU count")

p2 <- pairs_T8_combined %>%
    filter(!is.na(Isolate1Count_human)) %>%
    ggplot() +
    geom_abline(slope = 1, intercept = 0, color = "red", linetype = 2) +
    geom_hline(yintercept = c(0,1), color = gray(.8), linetype = 2) +
    geom_vline(xintercept = c(0,1), color = gray(.8), linetype = 2) +
    geom_point(aes(x = Isolate1CFUFreq_human, y = Isolate1CFUFreq_machine),
               shape = 21, size = 2, stroke = .4) +
    scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 21)) +
    theme_classic() +
    labs(x = "Random Forest CFU frequency", y = "Manual CFU frequency")

p <- plot_grid(p1, p2, nrow = 1, axis = "tblr", align = "h", scale = .9, labels = c("A", "B")) +
    paint_white_background()
ggsave(here::here("plots/FigS-human_machine_comparison.png"), p, width = 8, height = 4)

## Stats
pairs_T8_combined %>%
    filter(!is.na(Isolate1Count_human)) %>%
    lm(TotalCount_human ~ TotalCount_machine, data = .) %>%
    summary()
pairs_T8_combined %>%
    filter(!is.na(Isolate1Count_human)) %>%
    lm(Isolate1CFUFreq_human ~ Isolate1CFUFreq_machine, data = .) %>%
    summary()


# Table SXX image object features ----
features_example <- read_csv(paste0(folder_pipeline, "images/D-07-feature/merged/D_T8_C1R2_1.csv"), show_col_types = F)


features <- tibble(Feature = names(features_example)) %>%
    filter(!Feature == "ObjectID") %>%
    mutate(`Feature type` = case_when(
        str_sub(Feature, 1, 1) == "s" ~ "shape",
        str_sub(Feature, 1, 1) == "m" ~ "moment",
        str_sub(Feature, 1, 1) == "b" & str_sub(Feature, 1, 6) != "b.tran" ~ "intensity",
        str_sub(Feature, 1, 6) == "b.tran" ~ "transect"
    )) %>%
    mutate(Description =
               c("area size (in pixels)",
                 "perimeter (in pixels)",
                 "mean radius (in pixels)",
                 "standard deviation of the mean radius (in pixels)",
                 "min radius (in pixels)",
                 "max radius (in pixels)",
                 "center of mass x (in pixels)",
                 "center of mass y (in pixels)",
                 "elliptical fit major axis (in pixels)",
                 "elliptical eccentricity defined by sqrt(1-minoraxis^2/majoraxis^2). Circle eccentricity is 0 and straight line eccentricity is 1.",
                 "object angle (in radians)",
                 "green channel mean intensity",
                 "green channel standard deviation intensity",
                 "green channel mad intensity",
                 "green channel 1st quantile intensity",
                 "green channel 5th quantile intensity",
                 "green channel 10th quantile intensity",
                 "green channel 20th quantile intensity",
                 "green channel 50th quantile intensity",
                 "green channel 80th quantile intensity",
                 "green channel 90th quantile intensity",
                 "green channel 95th quantile intensity",
                 "green channel 99th quantile intensity",
                 "green channel transect mean intensity",
                 "green channel transect standard deviation intensity",
                 "green channel transect mad intensity",
                 "green channel the inmost transect pixel intensity",
                 "green channel the outmost transect pixel intensity",
                 "green channel difference between b.center and b.periphery",
                 "green channel pixel intensity at the 5% on the scaled transect",
                 "green channel pixel intensity at the 10% on the scaled transect",
                 "green channel pixel intensity at the 50% on the scaled transect",
                 "green channel pixel intensity at the 90% on the scaled transect",
                 "green channel pixel intensity at the 95% on the scaled transect",
                 "red channel mean intensity",
                 "red channel standard deviation intensity",
                 "red channel mad intensity",
                 "blue channel mean intensity",
                 "blue channel standard deviation intensity",
                 "blue channel mad intensity"
               ))

ft <- features %>%
    flextable() %>%
    width(j = 1:2, width = 2) %>%
    width(j = 3, width = 5)

save_as_image(ft, here::here("plots/TableS-features.png"), webshot = "webshot2")


# Table SXX List of isolates and images used for monocultures ----
isolates_epsilon <- read_csv(paste0(folder_data, "temp/06-isolates_epsilon.csv"), show_col_types = F) %>%
    select(Batch, Community, Isolate, Time, image_name, ColonyCount) %>%
    mutate(Community = factor(Community, communities$Community)) %>%
    arrange(Batch, Community) %>%
    rename(`Image name` = image_name, `Colony count` = ColonyCount)
t1 <- isolates_epsilon %>% slice(1:35)
t2 <- isolates_epsilon %>% slice(36:n())

ft1 <- t1 %>%
    flextable() %>%
    width(j = c(1:4, 6), width = .8) %>%
    width(j = 5, width = 1) %>%
    width(j = 6, width = 1.5) %>%
    highlight(i = which(t1$Time %in% c("T0", "T1")), j = 4, color = "yellow")
ft2 <- t2 %>%
    slice(1:35) %>%
    flextable() %>%
    width(j = c(1:4, 6), width = .8) %>%
    width(j = 5, width = 1) %>%
    width(j = 6, width = 1.5) %>%
    highlight(i = which(t2$Time %in% c("T0", "T1")), j = 4, color = "yellow")

save_as_image(ft1, here::here("plots/TableS-monoculture_1.png"), webshot = "webshot2")
save_as_image(ft2, here::here("plots/TableS-monoculture_2.png"), webshot = "webshot2")



# Table SXX fitness function ----
make_interaction_type <- function () {
    #' This function generates the fitness function table.
    #' There are a total of 27 possibilities
    interaction_type <- tibble(
        FromRare = rep(c(1, -1, 0), each = 9),
        FromMedium = rep(rep(c(1, -1, 0), each = 3), 3),
        FromAbundant = rep(c(1, -1, 0), 9),
        InteractionType = NA,
        InteractionTypeFiner = NA
    )
    ## Assign interaction types to combinations of frequency changes signs
    interaction_type$InteractionType[c(1,10,13,14)] <- "exclusion"
    interaction_type$InteractionType[c(2:6,8,11,20,23, 9,18,21,24,25,26,27)] <- "coexistence"

    ## Assign finer interaction types to combinations of frequency changes signs
    interaction_type$InteractionTypeFiner[c(1,14)] <- "competitive exclusion"
    interaction_type$InteractionTypeFiner[c(10,13)] <- "mutual exclusion"
    interaction_type$InteractionTypeFiner[c(2,5,8)] <- "stable coexistence"
    interaction_type$InteractionTypeFiner[c(4,6,11,20)] <- "frequency-dependent coexistence"
    interaction_type$InteractionTypeFiner[c(3)] <- "coexistence at 95%"
    interaction_type$InteractionTypeFiner[c(23)] <- "coexistence at 5%"
    interaction_type$InteractionTypeFiner[c(9,18,21,24:26)] <- "neutrality"
    interaction_type$InteractionTypeFiner[c(9,18,21,24:26)] <- "2-freq neutrality"
    interaction_type$InteractionTypeFiner[c(27)] <- "3-freq neutrality"
    interaction_type <- interaction_type %>%  mutate(FitnessFunction = paste(FromRare, FromMedium, FromAbundant, sep = "_"))
}
reformat <- function(x) {
    x <- ifelse(x == "1", "+", x)
    x <- ifelse(x == "-1", "-", x)
    x <- factor(x, c("+", "-", "0"))
    return(x)
}
interaction_type <- make_interaction_type()

pairs_interaction <- pairs %>%
    group_by(InteractionType, InteractionTypeFiner, FitnessFunction) %>%
    count(name = "Count")

ft <- interaction_type %>%
    left_join(pairs_interaction) %>%
    replace_na(list(Count = 0)) %>%
    mutate(across(starts_with("From"), reformat)) %>%
    select(-FitnessFunction) %>%
    rename(`Outcome` = InteractionType, `Finer outcome` = InteractionTypeFiner, `From rare` = FromRare, `From medium` = FromMedium, `From abundant` = FromAbundant) %>%
    replace_na(list(Outcome = "unknown")) %>%
    mutate(` ` = 1:n()) %>%
    select(` `, everything()) %>%
    flextable() %>%
    width(j = 1, width = .5) %>%
    width(j = 2:5, width = 1.3) %>%
    width(j = 6, width = 3) %>%
    vline(j = 1, border = NULL, part = "all")

save_as_image(ft, here::here("plots/TableS-fitness_function.png"), webshot = "webshot2")
















# Exploratory ----

## mismatch in different groups

p1 <- pairs %>%
    filter(!is.na(InteractionType)) %>%
    ggplot(aes(x = PairFermenter, y = Mismatch, color = PairFermenter)) +
    geom_boxplot(shape = 21, size = .5, position = position_dodge(width = 0.8)) +
    geom_point(shape = 21, size = 1, stroke = .5, position = position_jitterdodge(dodge.width = 0.8, jitter.width = .2)) +
    scale_color_npg() +
    theme_classic()


p2 <- pairs %>%
    filter(!is.na(InteractionType)) %>%
    ggplot(aes(x = InteractionType, y = Mismatch, color = PairFermenter)) +
    geom_boxplot(shape = 21, size = .5, position = position_dodge(width = 0.8)) +
    geom_point(shape = 21, size = 1, stroke = .5, position = position_jitterdodge(dodge.width = 0.8, jitter.width = .2)) +
    scale_color_npg() +
    theme_classic()


p3 <- pairs %>%
    filter(!is.na(InteractionType)) %>%
    group_by(PairFermenter, InteractionType) %>%
    count(name = "Count") %>%
    group_by(PairFermenter) %>%
    mutate(Fraction = Count / sum(Count), TotalCount = sum(Count)) %>%
    ggplot() +
    geom_col(aes(x = PairFermenter, y = Fraction, fill = InteractionType), color = 1) +
    geom_text(aes(x = PairFermenter, label = paste0("n=", TotalCount)), y = 0.9) +
    scale_fill_manual(values = assign_interaction_color()) +
    theme_classic() +
    theme(legend.title = element_blank()) +
    guides()

## Mismatch vs. coexistence
p4 <- pairs %>%
    filter(!is.na(InteractionType)) %>%
    mutate(InteractionType = ifelse(InteractionType == "coexistence", 1, 0)) %>%
    ggplot() +
    geom_point(aes(x = Mismatch, y = InteractionType), shape = 21, size = 2) +
    geom_smooth(aes(x = Mismatch, y = InteractionType), method = "glm", method.args = list(family = binomial)) +
    theme_classic()


model_logit <- pairs %>%
    filter(!is.na(InteractionType)) %>%
    mutate(InteractionType = ifelse(InteractionType == "exclusion", 1, 0)) %>%
    glm(InteractionType ~ Mismatch, family = "binomial", data = .)
summary(model_logit)





