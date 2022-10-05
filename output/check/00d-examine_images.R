# Read and write all immediate images/csv/figures/in the test example for examination
library(tidyverse)
library(EBImage)
library(cowplot)
library(magick)

folder_main <- "~/Dropbox/lab/emergent-coexistence/data/raw/plate_scan/emergent_coexistence_plate_scan_check/"
folder_script <- "~/Desktop/Lab/emergent-coexistence/output/check/"

list_images <- read_csv("~/Desktop/Lab/emergent-coexistence/output/check/00-list_images-D.csv", show_col_types = F)
list_image_mapping <- read_csv(paste0(folder_script, "00-list_image_mapping-D.csv"), show_col_types = F)
list_image_mapping_folder <- list_image_mapping %>%
    left_join(select(list_images, image_name_pair = image_name, folder_feature_pair = folder_green_feature, folder_green_cluster), by = "image_name_pair") %>%
    left_join(select(list_images, image_name_isolate1 = image_name, folder_feature_isolate1 = folder_green_feature), by = "image_name_isolate1") %>%
    left_join(select(list_images, image_name_isolate2 = image_name, folder_feature_isolate2 = folder_green_feature), by = "image_name_isolate2")


# to check many colonies on plates
#images_tocheck <- c("D_T8_C1R6_5-95_2_1", "D_T8_C1R5_5-95_3_1", "D_T8_C1R6_50-50_2_5", "D_T8_C1R7_5-95_1_6", "D_T8_C1R7_5-95_5_2")
images_tocheck <- c("D_T8_C1R2_2")
images_tocheck_index <- which(list_images$image_name %in% images_tocheck)

folder_examination <- paste0(folder_main, "examination/")

for (i in images_tocheck_index) {
    image_name <- list_images$image_name[i]
    # 00. original
    image_original <- readImage(paste0(list_images$folder_original[i], image_name, ".tiff"))
    writeImage(image_original, paste0(folder_examination, image_name, "-00-original.tiff"))
    # 01. green channel
    image_green <- readImage(paste0(list_images$folder_green[i], image_name, ".tiff"))
    writeImage(image_green, paste0(folder_examination, image_name, "-01-green_channel.tiff"))
    # 02. background subtraction
    image_green_rolled <- readImage(paste0(list_images$folder_green_rolled[i], image_name, ".tiff"))
    writeImage(image_green_rolled, paste0(folder_examination, image_name, "-02-green_rolled.tiff"))
    # 06. watershed
    image_green_watershed <- readImage(paste0(list_images$folder_green_watershed[i], image_name, ".tiff"))
    writeImage(image_green_watershed, paste0(folder_examination, image_name, "-06-green_watershed.tiff"))
    # 07 feature table
    table_green_feature <- read_csv(paste0(list_images$folder_green_feature[i], image_name, ".csv"), show_col_types = F)
    write_csv(table_green_feature, paste0(folder_examination, image_name, "-07-green_feature.csv"))
    # 08. transection
    image_green_transection <- readImage(paste0(list_images$folder_green_transection[i], image_name, ".tiff"))
    writeImage(image_green_transection, paste0(folder_examination, image_name, "-08-green_transection.tiff"))
    # 09 cluster figures. Only if it's a pair image
    if (length(unlist(str_split(image_name, "_"))) == 6) { # pair has four
        p_cluster <- ggdraw() +
            draw_image(paste0(list_images$folder_green_cluster[i], image_name, ".png"))
        ggsave(paste0(folder_examination, image_name, "-09-green_cluster.png"), plot = p_cluster, width = 15, height = 8)
    }
    cat("\n ", image_name)
}


# morphology: Perform morphological operations on images
# x = readImage(system.file("images", "shapes.png", package="EBImage"))
# kern = makeBrush(5, shape='diamond')
#
# display(x)
# display(kern, title='Structuring element')
# display(erode(x, kern), title='Erosion of x') # cut each pixel using the pattern kern
# display(dilate(x, kern), title='Dilatation of x') # dilate each pixel using the paterrn kern
# display(opening(x,kern)) # erosion followed by dilation
# display(makeBrush(size = , shape='diamond'))
# display(makeBrush(size = 101, shape='disc', step = T)) # THis may be use to identify the objects in the plates
# display(makeBrush(size = 99, shape='disc', step=FALSE))
# display(2000*makeBrush(501, shape='Gaussian', sigma=10))


# ## load images
# nuc = readImage(system.file('images', 'nuclei.tif', package='EBImage'))
# display(nuc, title='nucleus', method = "raster", frame = 1)
# cel = readImage(system.file('images', 'cells.tif', package='EBImage'))
# display(cel, title='cell', method = "raster", frame = 1)
# img = rgbImage(green=cel, blue=nuc)
# display(img, title='Cells', method = "raster", frame = 1)
#
# ## segment nuclei
# nuc = readImage(system.file('images', 'nuclei.tif', package='EBImage'))
# display(nuc, title='nucleus', method = "raster", frame = 1)
# nmask = thresh(nuc, 10, 10, 0.05) # this may be the solution to replace otsu
# display(nmask, method = "raster", frame = 1)
# # nmask <- nuc > otsu(nuc)
# # display(nmask, method = "raster", frame = 1)
# nmask = opening(nmask, makeBrush(5, shape='disc')) # erosion and dilation to remove wierd shits at the object border
# display(nmask, method = "raster", frame = 1)
#
# nmask = fillHull(nmask) # fill hollow center
# display(nmask, method = "raster", frame = 1)
# nmask = bwlabel(nmask) # label
# display(nmask, method = "raster", frame = 1)
# #display(normalize(nmask), title='Cell nuclei mask')
#
# ## segment cells, using propagate and nuclei as 'seeds'
# ctmask = opening(cel>0.1, makeBrush(5, shape='disc')) # morphology: Perform morphological operations on images
# display(ctmask, method = "raster", frame = 1)
# #in EBImage: Image processing and analysis toolbox for R
# #Functions to perform morphological operations on binary and grayscale images.
# display(ctmask, title='Cell mask', method = "raster", frame = 1)
# cmask = propagate(cel, nmask, ctmask)
# display(cmask, title='Cell mask', method = "raster", frame = 1)
# display(normalize(cmask), title='Cell mask', method = "raster", frame = 1)
#
# ## using paintObjects to highlight objects
# #res = paintObjects(cmask, img, col='#ff00ff')
# res = paintObjects(nmask, img, col='#ffff00')
# #display(img, title='Segmented cells', method = "raster", frame = 1)
# display(res, title='Segmented cells', method = "raster", frame = 1)
#
#
#
#
#
#
#
#
#
#
#
#
