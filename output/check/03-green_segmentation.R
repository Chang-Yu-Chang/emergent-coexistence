library(tidyverse)
library(EBImage)

list_images <- read_csv(commandArgs(trailingOnly = T)[1], show_col_types = F)
#list_images <- read_csv("~/Desktop/Lab/emergent-coexistence/output/check/00-list_images-D.csv", show_col_types = F)
compute_feature <- function (image_object, image_intensity) {
    computeFeatures(
        image_object, image_intensity,
        methods.noref = c("computeFeatures.shape"),
        methods.ref = c("computeFeatures.basic", "computeFeatures.moment"),
        basic.quantiles = c(0.01)) %>%
        as_tibble(rownames = "ObjectID") %>%
        # Remove duplicatedly calculated properties
        select(ObjectID, starts_with("x.0"), starts_with("x.Ba")) %>%
        # Remove the redundant prefix
        rename_with(function(x) str_replace(x,"x.0.", ""), starts_with("x.0")) %>%
        rename_with(function(x) str_replace(x,"x.Ba.", ""), starts_with("x.Ba")) %>%
        select(ObjectID, starts_with("b."), starts_with("s."), starts_with("m."))
}
detect_nonround_object <- function (image_object, image_intensity = NULL, watershed = F) {
    # Reomve too large or too small objects before watershed to reduce computational load
    if (!watershed) {
        # Check if the are away from the image border (use 100 pixel)
        oc <- ocontour(image_object)
        inside <- sapply(oc, function (x) {
            if (all(x[,1] > 100 & x[,1] < (nrow(image_rolled)-100) & x[,2] > 100 & x[,2] < (ncol(image_rolled)-100))) {
                return (T)
            } else return(F)
        })
        object_feature <- computeFeatures.shape(image_object) %>% as_tibble(rownames = "ObjectID")
        object_shape_round <- object_feature %>%
            mutate(inside = inside) %>%
            # Area
            filter(s.area > 300 & s.area < 20000) %>%
            # Contour located away from the edges
            filter(inside)

    }

    # Filter for circularity only after watershed segmentation
    if (watershed) {
        object_feature <- compute_feature(image_object, image_intensity)
        #object_feature <- compute_feature(image_watershed, image_rolled)

        # Remove segmented objects based on shape
        object_shape_round <- object_feature %>%
            # Area. Remove super small object after segementation
            filter(s.area > 300 & s.area < 20000) %>%
            # Circularity = 1 means a perfect circle and goes down to 0 for non-circular shapes
            mutate(Circularity = 4 * pi * s.area / s.perimeter^2) %>%
            filter(Circularity > 0.7) %>%
            # Remove tape and label that has really large variation in radius
            filter(s.radius.sd/s.radius.mean < 0.2) %>%
            filter(m.eccentricity < 0.8 & m.eccentricity != 0) # Circle eccentricity=0, straight line eccentricity=1


        # Remove outliers by b.sd/b.mean ratio
        object_shape_round <- object_shape_round %>%
            ungroup() %>%
            mutate(b.sd_over_mean = b.sd/b.mean) %>%
            mutate(b.sd_over_mean.up = quantile(b.sd_over_mean, .75) + 5 * IQR(b.sd_over_mean),
                   b.sd_over_mean.low = quantile(b.sd_over_mean, .25) - 5 * IQR(b.sd_over_mean)) %>%
            filter(b.sd_over_mean < b.sd_over_mean.up & b.sd_over_mean > b.sd_over_mean.low)
    }

    # Arrange by area size
    object_shape_round <- object_shape_round %>% arrange(desc(s.area))
    object_ID_nonround <- object_feature$ObjectID[!(object_feature$ObjectID %in% object_shape_round$ObjectID)]
    return(object_ID_nonround)
}

i = which(list_images$image_name %in% c("D_T8_C11R5_50-50_1_2"))

#i=1

for (i in 1:nrow(list_images)) {
    image_name <- list_images$image_name[i]

    # Rolled image
    image_rolled <- readImage(paste0(list_images$folder_green_rolled[i], image_name, ".tiff"))

    # 3. Thresholding
    #' thresh() applies a sliding window to calculate the local threshold
    #' opening() brushes the image to remove very tiny objects
    image_threshold <- thresh(-image_rolled, w = 150, h = 150, offset = 0.01) %>%
        opening(makeBrush(11, shape='disc'))
    writeImage(image_threshold, paste0(list_images$folder_green_threshold[i], image_name, ".tiff"))
    cat("\nthreshold")

    # 4. Detect round shaped object and remove super small size
    image_object <- bwlabel(image_threshold)
    object_ID_nonround <- detect_nonround_object(image_object, watershed = F)
    image_round <- rmObjects(image_object, object_ID_nonround, reenumerate = T)
    writeImage(image_round, paste0(list_images$folder_green_round[i], image_name, ".tiff"))
    cat("\tround object")

    # 5. Watershed
    image_distancemap <- distmap(image_round)
    cat("\tdistance map")
    image_watershed <- watershed(image_distancemap, tolerance = 1)
    ## After watershed, apply a second filter removing objects that are too small to be colonies
    object_ID_nonround2 <- detect_nonround_object(image_watershed, image_rolled, watershed = T)
    image_watershed2 <- rmObjects(image_watershed, object_ID_nonround2, reenumerate = T)
    save(image_watershed2, file = paste0(list_images$folder_green_watershed[i], image_name, ".RData")) # save watersed image object
    writeImage(colorLabels(image_watershed2), paste0(list_images$folder_green_watershed[i], image_name, ".tiff"))
    cat("\twatershed\t", i, "/", nrow(list_images), "\t", list_images$image_name[i])
}

