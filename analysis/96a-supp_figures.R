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
communities <- read_csv(paste0(folder_data, "temp/00c-communities.csv"), show_col_types = F)
accuracy <- read_csv(paste0(folder_data, "temp/91-accuracy.csv"), show_col_types = F)
pairs_freq <- read_csv(paste0(folder_data, "temp/93-pairs_freq.csv"), show_col_types = F)
load(paste0(folder_data, "temp/95-communities_network.Rdata"))
communities_hierarchy <- read_csv(paste0(folder_data, "temp/95-communities_hierarchy.csv"), show_col_types = F)

# 0. Clean up column factors ----
# Arrange communities by size
communities <- communities %>%
    arrange(CommunitySize) %>%
    mutate(Community = factor(Community, Community))

# Clean up the pairs data ----
pairs <- pairs %>%
    # Remove no-colony pairs
    unite(col = "Pair", Community, Isolate1, Isolate2, sep = "_", remove = F) %>%
    filter(!(Pair %in% pairs_no_colony)) %>%
    # Remove low-accuracy model pairs
    filter(AccuracyMean > 0.9)

# Figure S4 model accuracy----
accuracy_to_plot <- accuracy %>%
    # Remove pairs that have cocultures with no colony
    unite(col = "Pair", Community, Isolate1, Isolate2, sep = "_", remove = F) %>%
    filter(!(Pair %in% pairs_no_colony))

p1 <- accuracy_to_plot %>%
    ggplot() +
    geom_histogram(aes(x = Accuracy), color = 1, binwidth = 0.01, breaks = seq(0.6,1,.01), fill = NA) +
    geom_text(x = -Inf, y = Inf, label = paste0("N=", nrow(accuracy_to_plot)), vjust = 2, hjust = -1) +
    geom_vline(xintercept = 0.9, color = "red", linetype = 2) +
    theme_classic() +
    labs(x = "Accuracy", y = "Count")


accuracy_to_plot_count <- accuracy_to_plot %>%
    group_by(AccuracyPassThreshold) %>%
    count(name = "Count")
p2 <- accuracy_to_plot_count %>%
    ggplot() +
    geom_col(aes(x = AccuracyPassThreshold, y = Count), color = 1, fill = NA, width = .8) +
    geom_text(aes(x = AccuracyPassThreshold, y = Count, label = paste0("n=", Count)), vjust = -1) +
    geom_text(x = -Inf, y = Inf, label = paste0("N=", nrow(accuracy_to_plot)), vjust = 2, hjust = -1) +
    scale_x_discrete(labels = c("FALSE" = "Accuracy<0.9", "TRUE" = "Accuracy>0.9")) +
    scale_y_continuous(limits = c(0, 560)) +
    theme_classic() +
    labs(x = "")
p <- plot_grid(p1, p2, nrow = 1, labels = LETTERS[1:2], scale = 0.9) + paint_white_background()
ggsave(here::here("plots/FigS4-random_forest_accuracy.png"), p, width = 8, height = 4)



# Figure S5 machine vs. human ----
pairs_to_include <- pairs %>%
    select(Community, Isolate1, Isolate2) %>%
    mutate(Include = T)
pairs_T8_combined <- read_csv(paste0(folder_data, "temp/92-pairs_T8_combined.csv"), show_col_types = F) %>%
    left_join(pairs_to_include) %>%
    filter(!is.na(Include))

p1 <- pairs_T8_combined %>%
    ggplot() +
    geom_abline(slope = 1, intercept = 0, color = "red", linetype = 2) +
    geom_point(aes(x = TotalCount_human, y = TotalCount_machine), shape = 21, size = 2) +
    geom_text(x = -Inf, y = Inf, label = paste0("N=", nrow(pairs_T8_combined)), vjust = 2, hjust = -1) +
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
    geom_point(aes(x = Isolate1CFUFreq_human, y = Isolate1CFUFreq_machine), shape = 21, size = 2, stroke = .4) +
    geom_text(x = 0.5, y = 0.9, label = paste0("N=", pairs_T8_combined %>% filter(!is.na(Isolate1Count_human)) %>% nrow())) +
    scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 21)) +
    theme_classic() +
    labs(x = "Random Forest CFU frequency", y = "Manual CFU frequency")

p <- plot_grid(p1, p2, nrow = 1, axis = "tblr", align = "h", scale = .9, labels = c("A", "B")) +
    paint_white_background()

ggsave(here::here("plots/FigS5-human_machine_comparison.png"), p, width = 8, height = 4)

## cocultures without human results
pairs_T8_combined %>%
    filter(is.na(Isolate1Count_human)) %>%
    nrow()
## R-squared
pairs_T8_combined %>%
    filter(!is.na(Isolate1Count_human)) %>%
    lm(TotalCount_human ~ TotalCount_machine, data = .) %>%
    summary()
pairs_T8_combined %>%
    filter(!is.na(Isolate1Count_human)) %>%
    lm(Isolate1CFUFreq_human ~ Isolate1CFUFreq_machine, data = .) %>%
    summary()


# Figure S6 pairwise competition between highly abundant species ----
# This csv keeps all Sangers that match to one ESV with up to 0-2 mismatch
isolates_abundance_loose <- read_csv(paste0(folder_data, "temp/13-isolates_abundance_loose.csv"), show_col_types = F)

# Highly abundant isolate species have >0.05 relative ESV abundance
isolates_dominant <- isolates_abundance_loose %>%
    select(Community, Isolate, RelativeAbundance) %>%
    filter(RelativeAbundance > 0.05)

#
pairs_dominant <- pairs %>%
    select(PairID, Community, Isolate1, Isolate2, InteractionType) %>%
    left_join(rename(isolates_dominant, Isolate1 = Isolate, RelativeAbundance1 = RelativeAbundance)) %>%
    left_join(rename(isolates_dominant, Isolate2 = Isolate, RelativeAbundance2 = RelativeAbundance)) %>%
    filter(!is.na(RelativeAbundance1), !is.na(RelativeAbundance2)) %>%
    left_join(pairs)

# Stats
pairs_dominant %>%
    group_by(InteractionType) %>%
    count(name = "Count") %>%
    ungroup() %>%
    mutate(Fraction = Count / sum(Count))

#
p <- pairs_dominant %>%
    filter(!is.na(FitnessFunction)) %>%
    group_by(Community, InteractionType) %>%
    count(name = "Count") %>%
    group_by(Community) %>% mutate(Fraction = Count / sum(Count), TotalCount = sum(Count)) %>%
    ungroup() %>%
    mutate(Community = factor(Community, communities$Community)) %>%
    arrange(Community) %>%
    left_join(communities, by = "Community") %>%
    mutate(CommunityLabel = factor(CommunityLabel)) %>%
    replace_na(list(InteractionType = "unknown")) %>%
    ggplot() +
    geom_col(aes(x = Community, fill = InteractionType, y = Fraction), color = 1, width = .8, size = .5) +
    #geom_text(aes(x = CommunityLabel, y = .9, label = paste0("n=", CommunityPairSize)), vjust = -.5, size = 3) +
    geom_text(aes(x = Community, y = .9, label = paste0("n=", TotalCount))) +
    scale_fill_manual(values = assign_interaction_color(), breaks = c("coexistence", "exclusion", "unknown")) +
    scale_y_continuous(breaks = c(0,.5,1), limit = c(0, 1), expand = c(0,0)) +
    theme_classic() +
    theme(legend.text = element_text(size = 12),
          axis.text = element_text(color = 1, size = 12),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          axis.title = element_text(color = 1, size = 12),
          legend.title = element_blank(),
          legend.position = "top") +
    labs(x = "Community", y = "Fraction")

ggsave(here::here("plots/FigS6-pairwise_competition_dominant.png"), p, width = 8, height = 4)



# Table S1 image object features ----
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
               )) %>%
    mutate(` ` = 1:n()) %>%
    select(` `, everything())

ft <- features %>%
    flextable() %>%
    width(j = 1, width = 1) %>%
    width(j = 2:3, width = 2) %>%
    width(j = 4, width = 5) %>%
    vline(j = 1, border = NULL, part = "all")

save_as_image(ft, here::here("plots/TableS1-features.png"), webshot = "webshot2")


# Table S2 List of isolates and images used for monocultures ----
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

save_as_image(ft1, here::here("plots/TableS2-monoculture_1.png"), webshot = "webshot2")
save_as_image(ft2, here::here("plots/TableS2-monoculture_2.png"), webshot = "webshot2")



# Table S3 fitness function ----
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

save_as_image(ft, here::here("plots/TableS3-fitness_function.png"), webshot = "webshot2")















