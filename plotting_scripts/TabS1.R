library(tidyverse)
library(cowplot)
library(flextable)
source(here::here("analysis/00-metadata.R"))

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


