library(tidyverse)
library(EBImage)

list_images <- read_csv(commandArgs(trailingOnly = T)[1], show_col_types = F)
#list_images <- read_csv("~/Desktop/Lab/emergent-coexistence/output/check/00-list_images-D.csv", show_col_types = F)

detect_nonround_object <- function (image_object, watershed = F) {
    # Reomve too large or too small objects before watershed to reduce computational load
    if (!watershed) {
        # Check if the are away from the image border (use 100 pixel)
        oc <- ocontour(image_object)
        inside <- sapply(oc, function (x) {
            if (all(x[,1] > 100 & x[,1] < 2120 & x[,2] > 100 & x[,2] < 2120)) {
                return (T)
            } else return(F)
        })
        object_shape <- computeFeatures.shape(image_object) %>% as_tibble(rownames = "ObjectID")
        object_shape_round <- object_shape %>%
            mutate(inside = inside) %>%
            # Area
            filter(s.area > 300 & s.area < 20000) %>%
            #
            filter(inside)

    }

    # Filter for circularity only after watershed segmentation
    if (watershed) {
        object_shape <- computeFeatures.shape(image_object) %>% as_tibble(rownames = "ObjectID")
        object_moment <- computeFeatures.moment(image_object) %>% as_tibble(rownames = "ObjectID")

        object_shape_round <- object_shape %>%
            left_join(object_moment, by = "ObjectID") %>%
            # Area. Remove super small object after segementation
            filter(s.area > 300 & s.area < 20000) %>%
            # Circularity = 1 means a perfect circle and goes down to 0 for non-circular shapes
            mutate(Circularity = 4 * pi * s.area / s.perimeter^2) %>%
            filter(Circularity > 0.7) %>%
            # Remove tape and label that has really large variation in radius
            filter(s.radius.sd/s.radius.mean < 0.2) %>%
            filter(m.eccentricity < 0.8) # Circle eccentricity=0, straight line eccentricity=1
    }

    # Arrange by area size
    object_shape_round <- object_shape_round %>% arrange(desc(s.area))
    object_ID_nonround <- object_shape$ObjectID[!(object_shape$ObjectID %in% object_shape_round$ObjectID)]
    return(object_ID_nonround)
}
#i = which(list_images$image_name %in% c("D_T8_C1R2_5-95_2_1"))
#i = which(list_images$image_name %in% c("D_T8_C1R2_5-95_2_4"))
#i = which(list_images$image_name %in% c("D_T8_C11R5_50-50_1_4"))
for (i in 1:nrow(list_images)) {
    image_name <- list_images$image_name[i]
    image_rolled <- readImage(paste0(list_images$folder_green_rolled[i], image_name, ".tiff"))

    # 3. Thresholding
    threshold <- otsu(image_rolled)
    image_thresholded <- image_rolled < threshold
    cat("\nthreshold")

    # 4. Detect round shaped object and remove super small size
    image_object <- bwlabel(image_thresholded)
    object_ID_nonround <- detect_nonround_object(image_object, watershed = F)
    image_round <- rmObjects(image_object, object_ID_nonround, reenumerate = T)
    cat("\tround object")

    # 5-6. Watershed
    image_distancemap <- distmap(image_round)
    cat("\tdistance map")
    image_watershed <- watershed(image_distancemap, tolerance = 1)
    ## After watershed, apply a second filter removing objects that are too small to be colonies
    object_ID_nonround2 <- detect_nonround_object(image_watershed, watershed = T)
    image_watershed2 <- rmObjects(image_watershed, object_ID_nonround2, reenumerate = T)
    save(image_watershed2, file = paste0(list_images$folder_green_watershed_file[i], image_name, ".RData")) # save watersed image object
    writeImage(colorLabels(image_watershed2), paste0(list_images$folder_green_watershed[i], image_name, ".tiff")) # save
    cat("\twatershed\t", i, "/", nrow(list_images), "\t", list_images$image_name[i])

}



























