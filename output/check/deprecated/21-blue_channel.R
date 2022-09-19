library(tidyverse)
library(EBImage)

list_images <- read_csv(commandArgs(trailingOnly = T)[1])
#list_images <- read_csv("~/Desktop/Lab/emergent-coexistence/output/check/00-list_images-D.csv")

for (i in 1:nrow(list_images)) {
    image_name <- list_images$image_name[i]
    folder_original <- list_images$folder_original[i]
    folder_blue <- list_images$folder_blue[i]

    # 0. original image
    image_original <- readImage(paste0(folder_original, image_name, ".tiff"))

    # 1. blue channel
    temp <- image_original
    colorMode(temp) = Grayscale
    image_blue <- temp[,,2]
    writeImage(image_blue, paste0(folder_blue, image_name, ".tiff"))
    cat("\nblue channel\t", i, "/", nrow(list_images), "\t", list_images$image_name[i])
}


