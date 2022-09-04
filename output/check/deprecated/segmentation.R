#' Segementation

library(tidyverse)
library(EBImage)

folder_main <- "~/Dropbox/lab/emergent-coexistence/data/raw/plate_scan/emergent_coexistence_plate_scan_check/"

# Read image
#' ID, original image folder, background subtraction image folder, output folder
list_individuals <- list(
    #c("check/C_green/", "C_T8_C11R1_1", "example/"),
    c("D_T8_C1R2_1", "transtivity-D/pairs/", "check/D_green/", "example/"),
    c("D_T8_C1R2_2", "transtivity-D/pairs/", "check/D_green/", "example/"),
    c("D_T8_C1R2_5-95_1_2", "transtivity-D/pairs/", "check/D_green/", "example/")
)
i = 1

folder_source <- list_individuals[[i]][1]
img_name <- list_individuals[[i]][2]
folder_output <- list_individuals[[i]][3]

image <- readImage(paste0(folder_main, folder_source, img_name, ".tiff"))
#display(image)


# 1. Thresholding
#' Here because the images have undergone greyscale and background subtraction
#' so I apply a global threshold to the image
threshold <- otsu(image_grey)
image_thresholded <- image < threshold
#display(image_thresholded)

# 2. Detect round shaped object and remove super small size
#' To select potential colonies: area, perimeter, and circularity
#' This step is implemented here to prevent over segmentation and reduce segmentation load
image_object <- bwlabel(image_thresholded)
object_shape <- computeFeatures.shape(image_object) %>% as_tibble(rownames = "ObjectID")
#display(colorLabels(image_object))

object_shape_round <- object_shape %>%
    # Area
    filter(s.area > 500 & s.area < 500000) %>%
    # Roundness = 1 means a perfect circle
    mutate(Roundness = (s.radius.max - s.radius.min)/2) %>%
    filter(Roundness > 0.1 & Roundness < 50) %>%
    # Circularity = 1 means a perfect circle and goes down to 0 for non-cicular shapes
    mutate(Circularity = 4 * pi * s.area / s.perimeter^2) %>%
    filter(Circularity > 0.3) %>%
    arrange(desc(s.area))
object_ID_nonround <- which(!(object_shape$ObjectID %in% object_shape_round$ObjectID))

image_round <- rmObjects(image_object, object_ID_nonround, reenumerate = F)
#display(colorLabels(image_round))

# 3. Distance map
# The distance map contains for each pixel the distance to the nearest background pixe
image_distancemap <- distmap(image_round)
#display(normalize(image_distancemap), title = 'Distance map')

# 4. Watershed
#' This is the actual segmentation step
image_watershed <- watershed(image_distancemap, tolerance = 1)
display(colorLabels(image_watershed))

# 5. Use the segmented objects to locate the intensity from the greyscale image
object_intensity <- computeFeatures.basic(image_watershed, image) %>% as_tibble(rownames = "ObjectID")

object_intensity %>%
    pivot_longer(-ObjectID, names_to = "Property", values_to = "Value") %>%
    ggplot() +
    geom_histogram(aes(x = Value)) +
    facet_wrap(.~Property, scales = "free_x") +
    theme_classic()


# Output the example
image <- readImage(paste0(folder_main, folder_source, img_name, ".tiff"))

writeImage(image_object, paste0(folder_main, folder_output, sprintf("%02d", i), "-example_04-object.tiff"))
writeImage(image_round, paste0(folder_main, folder_output, sprintf("%02d", i), "-example_05-round_object.tiff"))
writeImage(normalize(image_distancemap), paste0(folder_main, folder_output, sprintf("%02d", i), "-example_06-distance_map.tiff"))
writeImage(colorLabels(image_watershed), paste0(folder_main, folder_output, sprintf("%02d", i), "-example_07-watershed.tiff"))











