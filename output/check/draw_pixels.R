
image_rolled


extract_transection <- function (watershed, ref) {
    #' This function search all objects on an image and return the pixel intensity along the transection of each object
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
        # max_radius[j] <- max(radius)
        # The coordinates of all along the line to the maximum line: this is arranged such that it's always from the center to the periphery
        lz <- bresenham(c(cz[1], mz[1]), c(cz[2], mz[2]))
        radial_gradient <- NULL
        for (k in 1:length(lz$x)) radial_gradient[k] <- ref@.Data[lz$x[k], lz$y[k]]

        transection[[j]] <- tibble(x = lz$x, y = lz$y, Intensity = radial_gradient)
    }

    return(transection)
}


transections <- NULL
for (i in 1:nrow(list_images)) {
    image_name <- list_images$image_name[i]
    image_rolled <- readImage(paste0(list_images$folder_green_rolled[i], image_name, ".tiff"))
    load(paste0(list_images$folder_green_watershed_file[i], image_name, ".RData")) # this should contain one R object image_watershed2

    transections[[i]] <- extract_transection(image_watershed2, image_rolled) %>%
        lapply(function(x) mutate(x, DistanceToCenter = 1:length(x))) %>%
        bind_rows(.id = "ObjectID")
}

transections[[1]]

i=1
image_name <- list_images$image_name[i]
image_rolled <- readImage(paste0(list_images$folder_green_rolled[i], image_name, ".tiff"))

#m[transections[[i]]$x[j], transections[[i]]$y[j]
draw_pixels <- function (img, pixel.x, pixel.y) {

    # Check argument length
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

i=1
image_name <- list_images$image_name[i]
image_rolled <- readImage(paste0(list_images$folder_green_rolled[i], image_name, ".tiff"))
image_transect <- draw_pixels(image_rolled, transections[[i]]$x, transections[[i]]$y)
display(img, method = "raster")
display(image_transect, method = "raster")


