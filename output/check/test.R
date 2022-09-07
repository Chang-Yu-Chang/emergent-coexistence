library(EBImage)
library(tidyverse)

nuc = readImage(system.file('images', 'nuclei.tif', package='EBImage'))
cel = readImage(system.file('images', 'cells.tif', package='EBImage'))
img = rgbImage(green=cel, blue=nuc)
display(img, title='Cells')

## segment nuclei
nmask = thresh(nuc, 10, 10, 0.05)
nmask = opening(nmask, makeBrush(5, shape='disc'))
nmask = fillHull(nmask)
nmask = bwlabel(nmask)
display(normalize(nmask), title='Cell nuclei mask')

## segment cells, using propagate and nuclei as 'seeds'
ctmask = opening(cel>0.1, makeBrush(5, shape='disc'))
cmask = propagate(cel, nmask, ctmask)
display(normalize(cmask), title='Cell mask')


## using paintObjects to highlight objects
res = paintObjects(cmask, img, col='#ff00ff')
#res = paintObjects(nmask, res, col='#ffff00')
display(res, title='Segmented cells')

list_images <- read_csv("~/Desktop/Lab/emergent-coexistence/output/check/00-list_images-D.csv")

#i = which(list_images$image_name=="D_T8_C1R2_2")
i = which(list_images$image_name=="D_T8_C1R2_1")
image_name <- list_images$image_name[i]
image_rolled <- readImage(paste0(list_images$folder_green_rolled[i], image_name, ".tiff"))
load(paste0(list_images$folder_green_watershed_file[i], image_name, ".RData")) # this should contain one R object image_watershed2

range(image_rolled)
threshold <- otsu(image_rolled, range = c(0.1, .9))
image_thresholded <- image_rolled < 0.9
display(image_thresholded, method = "raster")

paintObjects(image_watershed2, image_thresholded, col = "") %>%
    display(method = "raster")

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
image_object <- bwlabel(image_thresholded)
object_ID_nonround <- detect_nonround_object(image_object, watershed = F)
image_round <- rmObjects(image_object, object_ID_nonround, reenumerate = T)
image_distancemap <- distmap(image_round)
image_watershed <- watershed(image_distancemap, tolerance = 1)
paintObjects(image_watershed, image_rolled, col = "red") %>%
    display(method = "raster")















