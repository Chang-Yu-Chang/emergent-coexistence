library(tidyverse)
library(cowplot)
library(officer)
library(flextable)
source(here::here("processing_scripts/00-metadata.R"))

communities <- read_csv(paste0(folder_data, "output/communities_remained.csv"), show_col_types = F)

# Table S3 list of isolates, images used for monocultures, CFU, OD ----
isolates_epsilon <- read_csv(paste0(folder_data, "temp/06-isolates_epsilon.csv"), show_col_types = F) %>%
    select(Batch, Community, Isolate, Time, image_name, ColonyCount, OD620, Epsilon) %>%
    mutate(Community = factor(Community, communities$Community)) %>%
    mutate(OD620 = round(OD620, 3)) %>%
    mutate(Epsilon = format(Epsilon, scientific = T, digits = 2)) %>%
    arrange(Batch, Community) %>%
    rename(`Image name` = image_name, `Colony count` = ColonyCount) %>%
    left_join(communities) %>%
    select(Batch, ` ` = CommunityLabel, everything(), -CommunitySize, -CommunityPairSize, -ESVRichness)


t1 <- isolates_epsilon %>% slice(1:35)
t2 <- isolates_epsilon %>% slice(36:n())

ft1 <- t1 %>%
    flextable() %>%
    width(j = 2, width = .4) %>%
    width(j = c(1,3:5), width = .8) %>%
    width(j = 6, width = 1) %>%
    width(j = 7, width = 1.2) %>%
    highlight(i = which(t1$Time %in% c("T0", "T1")), j = 5, color = "yellow")
ft2 <- t2 %>%
    slice(1:35) %>%
    flextable() %>%
    width(j = 2, width = .4) %>%
    width(j = c(1,3:5), width = .8) %>%
    width(j = 6, width = 1) %>%
    width(j = 7, width = 1.2) %>%
    highlight(i = which(t2$Time %in% c("T0", "T1")), j = 5, color = "yellow")

save_as_image(ft1, here::here("plots/TableS3-monoculture_1.png"), webshot = "webshot2")
save_as_image(ft2, here::here("plots/TableS3-monoculture_2.png"), webshot = "webshot2")


