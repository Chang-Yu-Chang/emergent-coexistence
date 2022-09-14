library(tidyverse)
library(EBImage)

list_images <- read_csv(commandArgs(trailingOnly = T)[1])
#list_images <- read_csv("~/Desktop/Lab/emergent-coexistence/output/check/00-list_images-D.csv", show_col_types = F)
#list_images <- read_csv("~/Desktop/Lab/emergent-coexistence/output/check/00-list_images-C2.csv", show_col_types = F)
#list_images <- read_csv("~/Desktop/Lab/emergent-coexistence/output/check/00-list_images-B2.csv", show_col_types = F)
#list_images <- read_csv("~/Desktop/Lab/emergent-coexistence/output/check/00-list_images-C.csv", show_col_types = F)
i = which(list_images$image_name == "C_stock_C11R1_1")
#which(str_detect(list_images$image_name, "T0"))
i=3
#for (i in 8:9) {
for (i in 1:nrow(list_images)) {
    image_name <- list_images$image_name[i]
    folder_original <- list_images$folder_original[i]
    folder_green <- list_images$folder_green[i]

    # 0. original image
    image_original <- readImage(paste0(folder_original, image_name, ".tiff"))

    # 1. Green channel
    temp <- image_original
    colorMode(temp) = Grayscale
    image_green <- temp[,,2]
    writeImage(image_green, paste0(folder_green, image_name, ".tiff"))
    cat("\ngreen channel\t", i, "/", nrow(list_images), "\t", list_images$image_name[i])
}


