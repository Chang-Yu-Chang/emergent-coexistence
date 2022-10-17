library(tidyverse)
library(EBImage)

list_images <- read_csv(commandArgs(trailingOnly = T)[1], show_col_types = F) # 00-list_images-D.csv

#list_images <- read_csv("~/Desktop/lab/emergent-coexistence/analysis/00-list_images-D-green.csv", show_col_types = F)

paste_folder_name <- function (image_type = "channel", channel = "green") {
    paste0(list_images[i,paste0("folder_", image_type)], channel, "/")
}

# i = which(list_images$image_name == "B2_T1_C11R1_3")

for (i in 1:nrow(list_images)) {
    image_name <- list_images$image_name[i]
    color_channel <- list_images$color_channel[i]

    # 0. original image
    image_original <- readImage(paste0(list_images$folder_original[i], image_name, ".tiff"))

    # 1. channel
    temp <- image_original
    colorMode(temp) = Grayscale
    image_channel <- temp[,,match(color_channel, c("red", "green", "blue"))]
    writeImage(image_channel, paste0(paste_folder_name("channel", color_channel), image_name, ".tiff"))
    cat("\n", color_channel, " channel\t", i, "/", nrow(list_images), "\t", list_images$image_name[i])
}


