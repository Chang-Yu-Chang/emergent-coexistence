#' This script extracts the three color channels from the oroginal color images
#'
#' To use this Rscript, in bash environment:
#' Rscript 01-channel.R list_images.csv
#'
#' For example:
#' Rscript 01-channel.R mapping_files/00-list_images-B2-green.csv

library(tidyverse)
library(EBImage)

list_images <- read_csv(commandArgs(trailingOnly = T)[1], show_col_types = F) # 00-list_images-D.csv

paste_folder_name <- function (image_type = "channel", channel = "green") {
    paste0(list_images[i,paste0("folder_", image_type)], channel, "/")
}

for (i in 1:nrow(list_images)) {
    image_name <- list_images$image_name[i]
    color_channel <- list_images$color_channel[i]

    # 0. original image
    image_original <- readImage(paste0(list_images$folder_original[i], image_name, ".tiff"))

    # 1. channel
    temp <- image_original
    image_channel <- temp[,,match(color_channel, c("red", "green", "blue"))]
    colorMode(image_channel) = Grayscale
    writeImage(image_channel, paste0(paste_folder_name("channel", color_channel), image_name, ".tiff"))
    cat("\n", color_channel, " channel\t", i, "/", nrow(list_images), "\t", list_images$image_name[i])
}


