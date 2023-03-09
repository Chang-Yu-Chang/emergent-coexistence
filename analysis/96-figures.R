#' Script for figures

library(tidyverse)
library(cowplot)
library(tidygraph)
library(ggraph)
library(gridExtra)
source(here::here("analysis/00-metadata.R"))

isolates <- read_csv(paste0(folder_data, "output/isolates.csv"), show_col_types = F)
pairs <- read_csv(paste0(folder_data, "output/pairs.csv"), show_col_types = F)
communities <- read_csv(paste0(folder_data, "temp/00c-communities.csv"), show_col_types = F)
pairs_freq <- read_csv(paste0(folder_data, "temp/93-pairs_freq.csv"), show_col_types = F)
load(paste0(folder_data, "temp/95-communities_network.Rdata"))
communities_hierarchy <- read_csv(paste0(folder_data, "temp/95-communities_hierarchy.csv"), show_col_types = F)

# 0. Clean up column factors ----
# Arrange communities by size
communities <- communities %>%
    mutate(Community = factor(Community, Community))

# Count pairwise competition outcomes ----
pairs %>%
    group_by(InteractionType) %>%
    count(name = "Count") %>%
    ungroup() %>%
    mutate(Fraction = Count / sum(Count))
pairs %>%
    group_by(InteractionType, InteractionTypeFiner) %>%
    count(name = "Count") %>%
    ungroup() %>%
    mutate(Fraction = Count / sum(Count)) %>%
    mutate(InteractionTypeFiner = factor(InteractionTypeFiner, interaction_type_finer)) %>%
    arrange(InteractionTypeFiner)







