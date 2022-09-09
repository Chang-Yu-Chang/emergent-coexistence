library(tidyverse)
library(EBImage)
library(EBImageExtra)
library(purrr)

list_images <- read_csv(commandArgs(trailingOnly = T)[1], show_col_types = F)
#list_images <- read_csv("~/Desktop/Lab/emergent-coexistence/output/check/00-list_images-D.csv", show_col_types = F)
extract_transection <- function (watershed, ref) {
    #' This function searches for all objects on an image and return the pixel intensity along the transection of each object
    #' Arguments:
    #'  watershed is the watershed image with object labels
    #'  ref is the original image with intensity data. This has to be single-channel
    oc <- ocontour(watershed) # all contour points
    transection <- NULL
    #max_radius <- NULL # A test for object ID -> good. As long I use the oc output, they should be good
    for (j in 1:length(oc)) {
        z <- oc[[j]] # coordinate of contour pixels
        cz <- apply(z, 2, mean) # object center
        radius <- sqrt(rowSums((z - rep(cz, each=nrow(z)))^2)) # all radius from the center pixel to the contour pixels
        mz <- z[which.min(abs(radius - median(radius))),] # the coordinate of contour pixels with the median radius (or closest to the median)
        # The coordinates of all along the line to the maximum line: this is arranged such that it's always from the center to the periphery
        lz <- bresenham(c(cz[1], mz[1]), c(cz[2], mz[2]))
        radial_gradient <- NULL
        for (k in 1:length(lz$x)) radial_gradient[k] <- ref@.Data[lz$x[k], lz$y[k]]

        transection[[j]] <- tibble(x = lz$x, y = lz$y, Intensity = radial_gradient)
    }

    return(transection)
}
draw_pixels <- function (img, pixel.x, pixel.y) {
    #' This function takes an image and draw red on the assigned pixels
    stopifnot(length(pixel.x) == length(pixel.y))

    # Render the gray-scale image back into color mode
    img <- abind(img, img, img, along = 3) %>% Image(colormode = "Color")

    # Create a logical mask
    m <- Image(T, dim(img)[1:2])
    for (i in 1:length(pixel.x)) m[pixel.x[i], pixel.y[i]] <- F
    M <- abind(m, m, m, along = 3)

    # Combine with solid colored image, replace appropriate pixels in image
    mask <- Image("red", dim = dim(img)[1:2], colormode = "Color")
    ans <- img * M + mask * !M

    return(ans)
}
i = which(list_images$image_name == "D_T8_C4R1_50-50_1_3_-4")

#i = which(list_images$image_name %in% c("D_T8_C1R2_5-95_2_1")
images_tocheck_index = which(list_images$image_name %in% c("D_T8_C1R2_5-95_2_1",
                                        "D_T8_C1R2_5-95_2_4",
                                        "D_T8_C1R6_5",
                                        "D_T8_C1R7_50-50_3_4",
                                        "D_T8_C11R5_1",
                                        "D_T8_C11R5_50-50_1_4"))

#for (i in images_tocheck_index) {
for (i in 1:nrow(list_images)) {
    image_name <- list_images$image_name[i]
    image_rolled <- readImage(paste0(list_images$folder_green_rolled[i], image_name, ".tiff"))
    load(paste0(list_images$folder_green_watershed_file[i], image_name, ".RData")) # this should contain one R object image_watershed2

    # 7.1 Extract transection data
    transection <- extract_transection(image_watershed2, image_rolled) %>%
        lapply(function(x) mutate(x, DistanceToCenter = 1:length(x))) %>%
        bind_rows(.id = "ObjectID")
    write_csv(transection, file = paste0(list_images$folder_green_transection[i], image_name, ".csv"))
    cat("\ntransection")

    # Save the transection image
    image_transect <- draw_pixels(image_rolled, transection$x, transection$y)
    writeImage(image_transect, paste0(list_images$folder_green_transection[i], image_name, ".tiff"))
    cat("\tdraw transection")

    # 7.2 Smooth the curves using loess
    transection_smooth <- transection %>%
        nest(data = -ObjectID) %>%
        mutate(mod = map(data, loess, formula = Intensity ~ DistanceToCenter),
               fitted = map(mod, `[[`, "fitted")) %>%
        select(-mod) %>%
        unnest(cols = c(data, fitted))
    cat("\tsmooth curve")

    # 7.3 calculate features using the transection gradient
    transection_feature <- transection %>%
        group_by(ObjectID) %>%
        summarize(b.tran.mean = mean(Intensity), # summary statistics
                  b.tran.sd = sd(Intensity),
                  b.tran.mad = mad(Intensity),
                  b.center = Intensity[1], # center intensity
                  b.periphery = last(Intensity) ) %>% # periphery intensity
        mutate(b.diff.cp = b.periphery - b.center)

    # 7.4 Compute feature. The output table is NULL if no object (no colony)
    object_feature <- computeFeatures(
        image_watershed2, image_rolled,
        methods.noref = c("computeFeatures.shape"),
        methods.ref = c("computeFeatures.basic", "computeFeatures.moment"),
        basic.quantiles = c(0.01, 0.05, seq(0.1, 0.9, by = .1), 0.95, 0.99)
    )

    ## Execute the name cleanup only if there is at least 1 object
    if (is_null(object_feature)) {
        cat("\tno object\t", i, "/", nrow(list_images), "\t", image_name)
    } else if (!is_null(object_feature)) {
        object_feature <- object_feature %>%
            as_tibble(rownames = "ObjectID") %>%
            # Remove duplicatedly calculated properties
            select(ObjectID, starts_with("x.0"), starts_with("x.Ba")) %>%
            # Remove the redundant prefix
            rename_with(function(x) str_replace(x,"x.0.", ""), starts_with("x.0")) %>%
            rename_with(function(x) str_replace(x,"x.Ba.", ""), starts_with("x.Ba")) %>%
            select(ObjectID, starts_with("b."), starts_with("s."), starts_with("m.")) %>%
            # # Remove the round plate lid crack that has high intensity variability
            # filter(b.sd < 0.3)
            # Join the transection summary statistic features
            left_join(transection_feature, by = "ObjectID")

        write_csv(object_feature, paste0(list_images$folder_green_feature[i], image_name, ".csv"))
        cat("\tfeature\t", i, "/", nrow(list_images), "\t", image_name)
    }
}
