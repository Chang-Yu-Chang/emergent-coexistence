library(tidyverse)
library(cowplot)
library(officer)
library(flextable)
source(here::here("processing_scripts/00-metadata.R"))

# Table S2 community overview and labels ----
isolates <- read_csv(paste0(folder_data, "output/isolates_remained.csv"), show_col_types = F)
pairs <- read_csv(paste0(folder_data, "output/pairs.csv"), show_col_types = F)
pairs_remained <- read_csv(paste0(folder_data, "output/pairs_remained.csv"), show_col_types = F)
communities <- read_csv(paste0(folder_data, "output/communities_remained.csv"), show_col_types = F)

temp1 <- isolates %>%
    distinct(Community, Batch) %>%
    filter(!(Community == "C11R1" & Batch == "B2")) %>%
    rows_update(tibble(Community = "C11R1", Batch = "B2 and C"), by = "Community")
temp2 <- pairs %>%
    group_by(Community) %>%
    count(name = "TotalPairs") %>%
    ungroup()
temp3 <- pairs_remained %>%
    group_by(Community) %>%
    count(name = "TestedPairs") %>%
    ungroup()

ft <- communities %>%
    left_join(temp1, by = "Community") %>%
    left_join(temp2, by = "Community") %>%
    left_join(temp3, by = "Community") %>%
    select(CommunityLabel, Community, Batch, CommunitySize, TotalPairs, TestedPairs) %>%
    mutate(CommunityLabel = as.character(CommunityLabel)) %>%
    rename(` ` = CommunityLabel, `Community` = Community,
           `Number of isolates` = CommunitySize,
           `Number of total pairs` = TotalPairs,
           `Number of tested pairs` = TestedPairs) %>%
    janitor::adorn_totals() %>%
    flextable() %>%
    width(j = 1:2, width = 1.1) %>%
    width(j = 3, width = 0.8) %>%
    width(j = 4:6, width = 1.3) %>%
    hline(i = 13, border = fp_border(color = "black", style = "solid", width = 2))

save_as_image(ft, here::here("plots/TableS2-communities.png"), webshot = "webshot2")
