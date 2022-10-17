# Read and write all immediate images/csv/figures/in the test example for examination

library(tidyverse)
library(EBImage)
library(cowplot)
source(here::here("analysis/00-list_images.R"))

list_images <- read_csv("~/Desktop/Lab/emergent-coexistence/analysis/00-list_images-D-green.csv", show_col_types = F)
list_image_mapping <- read_csv(paste0(folder_script, "00-list_image_mapping-D.csv"), show_col_types = F)

list_image_mapping_folder <- list_image_mapping %>%
    left_join(rename(list_images, image_name_pair = image_name), by = "image_name_pair") %>%
    left_join(select(list_images, image_name_isolate1 = image_name), by = "image_name_isolate1") %>%
    left_join(select(list_images, image_name_isolate2 = image_name), by = "image_name_isolate2")


# To check colonies on plates
images_tocheck <- c("D_T8_C1R2_2", "D_T8_C1R2_5-95_1_2")
images_tocheck_index <- which(list_images$image_name %in% images_tocheck)
folder_examples <- paste0(folder_main, "examples/")

for (i in images_tocheck_index) {
    image_name <- list_images$image_name[i]
    # 00. original
    image_original <- readImage(paste0(list_images$folder_original[i], image_name, ".tiff"))
    writeImage(image_original, paste0(folder_examples, image_name, "-00-original.tiff"))
    # 01. green channel
    image_green <- readImage(paste0(list_images$folder_channel[i], list_images$color_channel[i], "/", image_name, ".tiff"))
    writeImage(image_green, paste0(folder_examples, image_name, "-01-green_channel.tiff"))
    # 02. background subtraction
    image_green_rolled <- readImage(paste0(list_images$folder_rolled[i], list_images$color_channel[i], "/", image_name, ".tiff"))
    writeImage(image_green_rolled, paste0(folder_examples, image_name, "-02-green_rolled.tiff"))
    # 03. threshold
    image_green_threshold <- readImage(paste0(list_images$folder_threshold[i], list_images$color_channel[i], "/", image_name, ".tiff"))
    writeImage(image_green_threshold, paste0(folder_examples, image_name, "-03-green_threshold.tiff"))
    # 04. round objects
    image_green_round <- readImage(paste0(list_images$folder_round[i], list_images$color_channel[i], "/", image_name, ".tiff"))
    writeImage(image_green_round, paste0(folder_examples, image_name, "-04-green_round.tiff"))
    # 05. watershed
    image_green_watershed <- readImage(paste0(list_images$folder_watershed[i], list_images$color_channel[i], "/", image_name, ".tiff"))
    writeImage(image_green_watershed, paste0(folder_examples, image_name, "-05-green_watershed.tiff"))
    # 06. transect
    image_green_transect <- readImage(paste0(list_images$folder_transect[i], list_images$color_channel[i], "/", image_name, ".tiff"))
    writeImage(image_green_transect, paste0(folder_examples, image_name, "-06-green_transect.tiff"))
    # 07. features
    table_feature <- read_csv(paste0(list_images$folder_feature[i], "merged/", image_name, ".csv"), show_col_types = F)
    write_csv(table_feature, paste0(folder_examples, image_name, "-07-green_feature.csv"))
    # cluster figures. Only if it's a pair image
    # monoculture has 4, where coculture has 6 segments
    if (length(unlist(str_split(image_name, "_"))) == 6) {
        # 08. logistic regression
        table_logit <- read_csv(paste0(list_images$folder_logit[i], image_name, ".csv"), show_col_types = F)
        write_csv(table_feature, paste0(folder_examples, image_name, "-07-green_feature.csv"))
        # 09. random forest
        p_cluster <- ggdraw() + draw_image(paste0(list_images$folder_random_forest[i], image_name, ".png"))
        ggsave(paste0(folder_examples, image_name, "-09-green_cluster.png"), plot = p_cluster, width = 15, height = 8)
    }
    cat("\n ", image_name)
}










