#' This script reads the colony position from feature/ and crop one colony out from the original images
#' The cropped images are later included in Fig S6 using Illustrator

library(EBImage)
library(tidyverse)
source(here::here("processing_scripts/00-metadata.R"))

paste0(folder_pipeline, "images/")
communities <- read_csv(paste0(folder_data, "temp/00c-communities.csv"), show_col_types = F)
list_image_mapping_folder_master <- read_csv(here::here("image_scripts/mapping_files/00-list_image_mapping_folder_master.csv"), show_col_types = F)

# Total number of isolate monoculture images should be 62+8=70
# 62 isolates including C C11R1 1-9
# 8 isolates B2 C11R1 2-9
length(unique(c(list_image_mapping_folder_master$image_name_isolate1, list_image_mapping_folder_master$image_name_isolate2)))
list_image_mapping_folder_master %>% distinct(Batch, Community, image_name_isolate1, .keep_all = T)


if (!dir.exists(paste0(folder_pipeline, "images/isolate_colony/"))) dir.create(paste0(folder_pipeline, "images/isolate_colony/"))
if (!dir.exists(paste0(folder_pipeline, "images/isolate_plate/"))) dir.create(paste0(folder_pipeline, "images/isolate_plate/"))

list_image_isolates <- list_image_mapping_folder_master %>%
    distinct(image_name_isolate1, .keep_all = T) %>%
    mutate(Community = factor(Community, paste0("C", rep(1:12, each = 8), "R", rep(1:8, 12)))) %>%
    arrange(Community, Batch, Isolate1)

window_length <- 100
list_image_crop <- rep(list(NA), nrow(list_image_isolates))

for (i in 1:length(list_image_crop)) {
    img <- readImage(paste0(list_image_isolates$folder_original[i], list_image_isolates$image_name_isolate1[i], ".tiff"))
    object_feature <- read_csv(paste0(list_image_isolates$folder_feature[i], "green/", list_image_isolates$image_name_isolate1[i], ".csv"), show_col_types = F)
    colony_center_x <- round(object_feature$m.cx[2],0)
    colony_center_y <- round(object_feature$m.cy[2],0)
    img_crop <- img[((colony_center_x-window_length):(colony_center_x+window_length)), ((colony_center_y-window_length):(colony_center_y+window_length)), 1:3]
    list_image_crop[[i]] <- img_crop
    cat(" ", i)

    writeImage(img_crop, paste0(folder_pipeline, "images/isolate_colony/", list_image_isolates$image_name_isolate1[i], ".tiff"))
    writeImage(img, paste0(folder_pipeline, "images/isolate_plate/", list_image_isolates$image_name_isolate1[i], ".tiff"))
}


