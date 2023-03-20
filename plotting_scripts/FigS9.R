library(tidyverse)
library(cowplot)
library(broom)
library(infer) # For tidyverse statistics
source(here::here("analysis/00-metadata.R"))

pairs <- read_csv(paste0(folder_data, "output/pairs.csv"), show_col_types = F)
communities <- read_csv(paste0(folder_data, "temp/00c-communities.csv"), show_col_types = F)

# Clean up the pairs data ----
pairs <- pairs %>%
    # Remove no-colony pairs
    unite(col = "Pair", Community, Isolate1, Isolate2, sep = "_", remove = F) %>%
    filter(!(Pair %in% pairs_no_colony)) %>%
    # Remove low-accuracy model pairs
    filter(AccuracyMean > 0.9)

pairs_fermenter_group <- pairs %>%
    filter(!is.na(PairFermenter)) %>%
    filter(!is.na(InteractionType), InteractionType != "unknown") %>%
    mutate(InteractionType = factor(InteractionType, c("coexistence", "exclusion"))) %>%
    group_by(PairFermenter, InteractionType) %>%
    summarize(Count = n()) %>%
    mutate(Fraction = Count / sum(Count), TotalCount = sum(Count))

p <- pairs_fermenter_group %>%
    ggplot() +
    geom_col(aes(x = InteractionType, y = Fraction, fill = InteractionType), color = 1, width = 0.7, position = position_dodge()) +
    geom_text(aes(x = InteractionType, y = Fraction, label = Count), vjust = 2, position = position_dodge(width = 0.9)) +
    scale_fill_manual(values = interaction_color) +
    facet_grid(.~PairFermenter, labeller = as_labeller(c(FF = "fermenter-fermenter", FR = "fermenter-respirator", RR = "respirator-respirator"))) +

    # geom_col(aes(x = PairFermenter, y = Fraction, fill = InteractionType), color = 1, width = 0.7) +
    # geom_text(aes(x = PairFermenter, label = paste0("n=", TotalCount)), y = 0.9) +
    # scale_fill_manual(values = interaction_color) +
    # scale_x_discrete(labels = c(FF = "fermenter-fermenter", FR = "fermenter-respirator", RR = "respirator-respirator")) +
    # scale_y_continuous(expand = c(0,0), breaks = seq(0, 1, 0.2)) +
    theme_classic() +
    theme(panel.border = element_rect(color = 1, fill = NA),
          axis.text.x = element_text(angle = 20, hjust = 1)) +
    guides(fill = "none") +
    labs(x = "")

matrix(
    c(pairs_fermenter_group %>% filter(PairFermenter == "FF") %>% pull(Count),
    pairs_fermenter_group %>% filter(PairFermenter == "FR") %>% pull(Count),
    pairs_fermenter_group %>% filter(PairFermenter == "RR") %>% pull(Count)),
    nrow = 2, byrow = F
) %>%
    chisq.test

ggsave(here::here("plots/FigS9-fermenter_vs_coexistence.png"), p, width = 5, height = 4)

