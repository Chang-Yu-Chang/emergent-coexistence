#' This is the step 3 in couting the colony features
#' TSA plate image, green channel, grey scale, subtract background
#'  The processed image files are stored in folders "Batch"_green/
#'  - check/B2_green/
#'  - check/D_green/
#'  - etc..

library(tidyverse)
library(EBImage)

# Read image
list_individuals <- list(
    c("check/C_green/", "C_T8_C11R1_1", "test/"),
    c("check/C_green/", "C_T8_C11R1_2", "test/"),
    c("check/C_green/", "C_T8_C11R1_4", "test/")
)
i=3
folder_source <- list_individuals[[i]][1]
img_name <- list_individuals[[i]][2]
folder_output <- list_individuals[[i]][3]
output_file <- paste0("~/Dropbox/lab/emergent-coexistence/data/raw/plate_scan/emergent_coexistence_plate_scan_check/", folder_output, img_name, ".tiff")
cat("\n", i, "/", length(list_individuals), list_individuals[[i]][2])

image <- readImage(paste0("~/Dropbox/lab/emergent-coexistence/data/raw/plate_scan/emergent_coexistence_plate_scan_check/", folder_source, img_name, ".tiff"))
#image <- readImage(system.file("images", "nuclei.tif", package="EBImage"))[,,1]
display(image)

# 1. Thresholding
#' Here because the images have undergone greyscale and background subtraction
#' so I use a global threshold
threshold <- otsu(image)
image_thresholded <- image < threshold
display(image_thresholded)

# 2. Distance map
# The distance map contains for each pixel the distance to the nearest background pixe
image_distancemap <- distmap(image_thresholded)
display(normalize(image_distancemap), title='Distance map')

# 3. Watershed
image_watershed <- watershed(image_distancemap, tolerance = .1)
display(normalize(image_watershed), title='watershed')

# 3. Image segmentation
# performs partitioning of an image, and is typically used to identify objects in an image.
image_segmented <- bwlabel(image_watershed)
display(colorLabels(image_segmented))

# 4. Clean up.
# For example, remove the apparely too large objects (the plate boundary) and too small ones

# 5. compute object properties


# x = thresh(y, 10, 10, 0.05)
# x = opening(x, makeBrush(5, shape='disc'))
# x = bwlabel(x)
# display(y, title="Cell nuclei")
# display(x, title="Segmented nuclei")













