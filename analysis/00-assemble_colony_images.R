library(EBImage)
library(tidyverse)

paste0(folder_pipeline, "images/")
communities <- read_csv(paste0(folder_data, "temp/00c-communities.csv"), show_col_types = F)
list_image_mapping_folder_master <- read_csv(here::here("analysis/mapping_files/00-list_image_mapping_folder_master.csv"), show_col_types = F)


# Total number of isolate monoculture images should be 78
list_image_mapping_folder_master %>%
    distinct(image_name_isolate1) %>%
    nrow()
#' C11R2 isolate 13 is a contaminant
#' B2 C11R1 isolate 1 is contaminated so Batch C C11R1 isolates 1-9 is to replace all B2 pairs with isolate 1
#' 78 - 1 - 9 = 68 isolates

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


