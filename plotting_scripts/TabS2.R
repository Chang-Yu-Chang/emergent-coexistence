library(tidyverse)
library(cowplot)
library(officer)
library(flextable)
source(here::here("analysis/00-metadata.R"))

# Table S2 community overview and labels ----
isolates <- read_csv(paste0(folder_data, "output/isolates.csv"), show_col_types = F)
pairs <- read_csv(paste0(folder_data, "output/pairs.csv"), show_col_types = F)
communities <- read_csv(paste0(folder_data, "temp/00c-communities.csv"), show_col_types = F)

temp1 <- isolates %>%
    distinct(Community, Batch) %>%
    filter(!(Community == "C11R1" & Batch == "B2")) %>%
    rows_update(tibble(Community = "C11R1", Batch = "B2 and C"), by = "Community")
temp2 <- pairs %>%
    group_by(Community) %>%
    filter(!is.na(InteractionType)) %>%
    filter(AccuracyMean > 0.9) %>%
    count(name = "ActualPairs") %>%
    ungroup()
ft <- communities %>%
    left_join(temp1, by = "Community") %>%
    left_join(temp2, by = "Community") %>%
    select(CommunityLabel, Community, Batch, CommunitySize, CommunityPairSize, ActualPairs) %>%
    mutate(CommunityLabel = as.character(CommunityLabel)) %>%
    rename(`Community` = CommunityLabel, `Internal label` = Community, `Number of isolates` = CommunitySize,
           `Number of tested pairs` = CommunityPairSize, `Number of applicable pairs` = ActualPairs) %>%
    janitor::adorn_totals() %>%
    flextable() %>%
    width(j = 1:2, width = 1.1) %>%
    width(j = 3, width = 0.8) %>%
    width(j = 4:6, width = 1.3) %>%
    hline(i = 13, border = fp_border(color = "black", style = "solid", width = 2))

save_as_image(ft, here::here("plots/TableS2-communities.png"), webshot = "webshot2")
