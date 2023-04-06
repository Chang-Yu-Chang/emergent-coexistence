#' This script checks the community abundance data
library(tidyverse)
source(here::here("analysis/00-metadata.R"))

# communities <- read_csv(paste0(folder_data, "temp/00c-communities.csv"), show_col_types = F) %>%
#     mutate(Community = factor(Community, Community))
emergent_communities <- read_csv(paste0(folder_data, "temp/13-emergent_communities.csv"), show_col_types = F)

communities_abundance <- emergent_communities %>%
    # Remove Leucine and Citrate communities
    filter(Carbon_Source == "Glucose" | Carbon_Source == "Original") %>%
    mutate(Community = factor(paste0("C", Inoculum, "R", Replicate), paste0("C", rep(1:12, each = 8), "R", rep(1:8, 12))))%>%
    arrange(Community, Family, Transfer, ESV)

write_csv(communities_abundance, paste0(folder_data, 'temp/14-communities_abundance.csv'))
