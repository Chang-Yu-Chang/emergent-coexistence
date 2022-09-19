library(tidyverse)
library(EBImage)

list_images <- read_csv(commandArgs(trailingOnly = T)[1])
#list_images <- read_csv("~/Desktop/Lab/emergent-coexistence/output/check/00-list_images-D.csv")

detect_nonround_object <- function (image_object, watershed = F) {
    # Reomve too large or too small objects before watershed to reduce computational load
    if (!watershed) {
        object_shape <- computeFeatures.shape(image_object) %>% as_tibble(rownames = "ObjectID")
        object_shape_round <- object_shape %>%
            # Area
            filter(s.area > 300 & s.area < 20000)
    }

    # Filter for circularity only after watershed segmentation
    if (watershed) {
        object_shape <- computeFeatures.shape(image_object) %>% as_tibble(rownames = "ObjectID")
        object_moment <- computeFeatures.moment(image_object) %>% as_tibble(rownames = "ObjectID")

        object_shape_round <- object_shape %>%
            left_join(object_moment, by = "ObjectID") %>%
            # Circularity = 1 means a perfect circle and goes down to 0 for non-circular shapes
            mutate(Circularity = 4 * pi * s.area / s.perimeter^2) %>%
            filter(Circularity > 0.3) %>%
            # Remove tape and label that has really large variation in radius
            filter(s.radius.sd/s.radius.mean < 1/2) %>%
            filter(m.eccentricity < 0.8) # Circle eccentricity=0, straight line eccentricity=1
    }

    # Arrange by area size
    object_shape_round <- object_shape_round %>% arrange(desc(s.area))
    object_ID_nonround <- object_shape$ObjectID[!(object_shape$ObjectID %in% object_shape_round$ObjectID)]
    return(object_ID_nonround)
}

for (i in 1:nrow(list_images)) {
    image_name <- list_images$image_name[i]
    image_rolled <- readImage(paste0(list_images$folder_red_rolled[i], image_name, ".tiff"))

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
    save(image_watershed2, file = paste0(list_images$folder_red_watershed_file[i], image_name, ".RData")) # save watersed image object
    writeImage(colorLabels(image_watershed2), paste0(list_images$folder_red_watershed[i], image_name, ".tiff")) # save
    cat("\twatershed\t", i, "/", nrow(list_images), "\t", list_images$image_name[i])

    # # 7. Calculate feature
    # ## Compute feature. It is NULL if no object (no colony)
    # object_feature <- computeFeatures(
    #     image_watershed2, image_rolled,
    #     methods.noref = c("computeFeatures.shape"),
    #     methods.ref = c("computeFeatures.basic", "computeFeatures.moment"),
    #     basic.quantiles = c(0.01, 0.05, seq(0.1, 0.9, by = .1), 0.95, 0.99)
    # )
    #
    # ## Execute the name cleanup only if there is at least 1 object
    # if (is_null(object_feature)) {
    #     #write_csv(object_feature, paste0(list_images$folder_red_feature[i], image_name, ".csv"))
    #     cat("\tno object\t", i, "/", nrow(list_images), "\t", list_images$image_name[i])
    # }
    #
    # if (!is_null(object_feature)) {
    #     object_feature <- object_feature %>%
    #         as_tibble(rownames = "ObjectID") %>%
    #         # Remove duplicatedly calculated properties
    #         select(ObjectID, starts_with("x.0"), starts_with("x.Ba")) %>%
    #         # Remove the redundant prefix
    #         rename_with(function(x) str_replace(x,"x.0.", ""), starts_with("x.0")) %>%
    #         rename_with(function(x) str_replace(x,"x.Ba.", ""), starts_with("x.Ba")) %>%
    #         #
    #         select(ObjectID, starts_with("b."), starts_with("s."), starts_with("m.")) %>%
    #         # Remove the round plate lid crack that has high intensity variability
    #         filter(b.sd < 0.3)
    #
    #     write_csv(object_feature, paste0(list_images$folder_red_feature[i], image_name, ".csv"))
    #     cat("\tfeature\t", i, "/", nrow(list_images), "\t", list_images$image_name[i])
    # }

}



























