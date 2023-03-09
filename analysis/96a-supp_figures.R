#' Script for supplement figures

library(tidyverse)
library(cowplot)
library(broom)
library(infer) # For tidyverse statistics
library(ggsignif) # For adding asterisks to boxplots
library(RColorBrewer) # For color
library(ggsci) # For color
library(ggraph)
library(tidygraph)
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




## Analysis
lm(meanIsolate1CFUFreqMean ~ RelativeAbundance1, data = pairs_freq_ESV) %>% summary()



