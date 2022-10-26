#' This script computes the object features using 1) a roadmap image and 2) a reference image
#'
#' To use this Rscript, in bash environment:
#' Rscript 04-feature.R list_images.csv
#'
#' For example:
#' Rscript 04-feature.R mapping_files/00-list_images-B2-green.csv

library(tidyverse)
library(EBImage)
library(EBImageExtra) # for the bresenham algorithm. `devtools::install_github("ornelles/EBImageExtra")`
library(purrr) # for applying functional programming to transect curve smoothing
source(here::here("analysis/00-metadata.R"))

list_images <- read_csv(commandArgs(trailingOnly = T)[1], show_col_types = F)
#list_images <- read_csv("~/Desktop/lab/emergent-coexistence/analysis/00-list_images-D-green.csv", show_col_types = F)

paste_folder_name <- function (image_type = "channel", channel = "green") {
    paste0(list_images[i,paste0("folder_", image_type)], channel, "/")
}
extract_transect <- function (watershed, ref) {
    #' This function searches for all objects on an image and return the pixel intensity along the transect of each object
    #' Arguments:
    #'  watershed is the watershed image with object labels
    #'  ref is the original image with intensity data. This has to be single-channel
    oc <- ocontour(watershed) # all contour points
    transect <- NULL
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

        transect[[j]] <- tibble(x = lz$x, y = lz$y, Intensity = radial_gradient)
    }

    return(transect)
}
smooth_transect <- function (transect, span = .8) {
    #' This function smooth the transects using loess that fits a polynomial regression
    #' the parameter span controls the degree of smoothing. from span = 0 very rugged to span = 1 very smooth
    loess_custom <- function (formula, data) loess(formula, data, span = span)
    transect %>%
        # Smooth intensity
        nest(data = -ObjectID) %>%
        mutate(mod = map(data, loess_custom, formula = Intensity ~ DistanceToCenter),
               FittedIntensity = map(mod, `[[`, "fitted")) %>%
        select(-mod) %>%
        unnest(cols = c(data, FittedIntensity))
}
diff_transect <- function (transect_smooth) {
    #' This function scales the transect and calculate the derivative
    # Derivative. f'(x) = (f(x+h)-f(h)) / h
    transect_smooth %>%
        group_by(ObjectID) %>%
        # Scale the Distance
        mutate(ScaledDistanceToCenter = DistanceToCenter / max(DistanceToCenter)) %>%
        mutate(Derivative = c(NA, diff(FittedIntensity) / diff(ScaledDistanceToCenter))) %>%
        mutate(SecondDerivative = c(NA, diff(Derivative) / diff(ScaledDistanceToCenter))) %>%
        select(ObjectID, x, y, DistanceToCenter, ScaledDistanceToCenter, ends_with("Intensity"), ends_with("Derivative"))
}
calculate_transect <- function (transect_smooth) {
    #' This function calculates three types of transect features
    #' 1. the onset of the fist bumpt.
    #' 2. the number of bumps
    #' 3. the intensity of distance_to_center quantiles

    # Onset of the first bump. DIstance to the center is scaled to [0, 1]
    transect_onset_bump <- transect_smooth %>%
        group_by(ObjectID, .drop = F) %>%
        filter(SecondDerivative < 0) %>%
        slice(1) %>%
        select(ObjectID, OnsetBump = ScaledDistanceToCenter)


    # Number of bumps
    transect_n_bump <- transect_smooth %>%
        group_by(ObjectID, .drop = F) %>%
        filter(SecondDerivative < 0) %>%
        count(name = "Count")

    # Transect quantile by distance to center
    find_quantile <- function (x, q) abs(x - q) == min(abs(x - q))
    transect_q005 <- transect_smooth %>%
        select(ObjectID, ScaledDistanceToCenter, Intensity) %>%
        filter(find_quantile(ScaledDistanceToCenter, 0.05)) %>%
        slice(1) %>% # If there are two pixels equally close to 0.05, for example 0 and 0.1, then choose the first 1
        select(ObjectID, b.tran.q005 = Intensity)
    transect_q01 <- transect_smooth %>%
        select(ObjectID, ScaledDistanceToCenter, Intensity) %>%
        filter(find_quantile(ScaledDistanceToCenter, 0.1)) %>%
        slice(1) %>% # If there are two pixels equally close to 0.05, for example 0 and 0.1, then choose the first 1
        select(ObjectID, b.tran.q01 = Intensity)
    transect_q05 <- transect_smooth %>%
        select(ObjectID, ScaledDistanceToCenter, Intensity) %>%
        filter(find_quantile(ScaledDistanceToCenter, 0.5)) %>%
        slice(1) %>% # If there are two pixels equally close to 0.05, for example 0 and 0.1, then choose the first 1
        select(ObjectID, b.tran.q05 = Intensity)
    transect_q09 <- transect_smooth %>%
        select(ObjectID, ScaledDistanceToCenter, Intensity) %>%
        filter(find_quantile(ScaledDistanceToCenter, 0.9)) %>%
        slice(1) %>% # If there are two pixels equally close to 0.05, for example 0 and 0.1, then choose the first 1
        select(ObjectID, b.tran.q09 = Intensity)
    transect_q095 <- transect_smooth %>%
        select(ObjectID, ScaledDistanceToCenter, Intensity) %>%
        filter(find_quantile(ScaledDistanceToCenter, 0.95)) %>%
        slice(1) %>% # If there are two pixels equally close to 0.05, for example 0 and 0.1, then choose the first 1
        select(ObjectID, b.tran.q095 = Intensity)

    #
    transect_n_bump %>%
        left_join(transect_onset_bump, by = "ObjectID") %>%
        left_join(transect_q005, by = "ObjectID") %>%
        left_join(transect_q01, by = "ObjectID") %>%
        left_join(transect_q05, by = "ObjectID") %>%
        left_join(transect_q09, by = "ObjectID") %>%
        left_join(transect_q095, by = "ObjectID") %>%
        return()
}
plot_transect <- function (transect_smooth, title_name = NULL, scaled_distance = T) {
    p <- transect_smooth %>%
        pivot_longer(cols = c(Intensity, FittedIntensity, Derivative, SecondDerivative), names_to = "Measure", values_to = "Value") %>%
        mutate(Measure = ordered(Measure, c("Intensity", "FittedIntensity", "Derivative", "SecondDerivative"))) %>%
        ggplot() +
        facet_grid(Measure~., scales = "free_y") +
        theme_classic() +
        theme(panel.border = element_rect(color = 1, fill = NA)) +
        ggtitle(title_name)

    if (scaled_distance) {
        p <- p +
            geom_line(aes(x = ScaledDistanceToCenter, y = Value, group = ObjectID), lwd = .3) +
            geom_point(aes(x = ScaledDistanceToCenter, y = Value, group = ObjectID), size = .1)
    }

    if (scaled_distance == F) {
        p <- p +
            geom_line(aes(x = DistanceToCenter, y = Value, group = ObjectID), lwd = .3) +
            geom_point(aes(x = DistanceToCenter, y = Value, group = ObjectID), size = .1)
    }

    p <- p +
        geom_hline(data = tibble(Measure = factor(c("Intensity", "FittedIntensity", "Derivative", "SecondDerivative"), c("Intensity", "FittedIntensity", "Derivative", "SecondDerivative")),
                                 yintercept = c(NA, NA, 0, 0)),
                   aes(yintercept = yintercept), color = "red", linetype = 2)

    return(p)
}
draw_pixels <- function (img, pixel.x, pixel.y, color = "red") {
    #' This function takes an image and draw red on the assigned pixels
    stopifnot(length(pixel.x) == length(pixel.y))

    # If the image is gray-scale, render it into color mode
    if (length(dim(img)) == 2) img <- abind(img, img, img, along = 3) %>% Image(colormode = "Color")

    # Create a logical mask
    m <- Image(T, dim(img)[1:2])
    for (i in 1:length(pixel.x)) m[pixel.x[i], pixel.y[i]] <- F
    M <- abind(m, m, m, along = 3)

    # Combine with solid colored image, replace appropriate pixels in image
    mask <- Image(color, dim = dim(img)[1:2], colormode = "Color")
    ans <- img * M + mask * !M

    return(ans)
}
remove_outliers <- function (object_feature, multiplier = 2, features = c("b.sd", "b.mad", "b.mean", "b.q05", "b.q005", "b.tran.sd", "b.tran.mad")) {
    #' This function uses a interquantile rule to find outliers
    #' for a feature, if the data point falls outside the range of [Q1-2*IQR, Q3+2*IQR], it's a outlier

    stopifnot(features %in% colnames(object_feature))


    x <- 0
    for (feature in features) {
        drop_list <- object_feature %>%
            filter(!between(
                get(feature),
                quantile(get(feature), probs = 0.25, na.rm=TRUE) - (multiplier * IQR(get(feature), na.rm=TRUE)),
                quantile(get(feature), probs = 0.75, na.rm=TRUE) + (multiplier * IQR(get(feature), na.rm=TRUE))
            ))

        object_feature <- object_feature %>%
            filter(between(
                get(feature),
                quantile(get(feature), probs = 0.25, na.rm=TRUE) - (multiplier * IQR(get(feature), na.rm=TRUE)),
                quantile(get(feature), probs = 0.75, na.rm=TRUE) + (multiplier * IQR(get(feature), na.rm=TRUE))
        ))
        # Report the number of objects dropped
        x <- x + nrow(drop_list)
        if (nrow(drop_list) != 0) cat("\n", nrow(drop_list), " outlier object(s) dropped from feature", feature)
    }
    cat("\ntotal ", x, " outliers object(s) dropped\n")
    return(object_feature)
}

for (i in 1:nrow(list_images)) {
    image_name <- list_images$image_name[i]
    color_channel <- list_images$color_channel[i]

    # 6.0 gating the no-object images
    ## no watershed image
    if (image_name %in% plates_no_colony) {
        cat("\nno colony, no watershed image\t", image_name)
        next
    }

    ## Image with intensity
    image_rolled <- readImage(paste0(paste_folder_name("rolled", color_channel), image_name, ".tiff"))
    ## Read the watershed file only from the GREEN channel. This should contain one R object image_watershed2
    load(paste0(paste_folder_name("watershed", "green"), image_name, ".RData"))

    ## No object on the watershed image
    if (all(image_watershed2 == 0)) {
        cat("\tno object\t", i, "/", nrow(list_images), "\t", image_name)
        next
    }

    # 6.1 Extract transect data
    transect <- extract_transect(image_watershed2, image_rolled) %>%
        lapply(function(x) mutate(x, DistanceToCenter = 0:(length(x)-1))) %>%
        bind_rows(.id = "ObjectID")

    # 6.2 Calculate transect features
    ## Smooth the transect intensity
    transect_smooth <- transect %>%
        smooth_transect() %>%
        diff_transect()

    ## Transect features
    transect_feature <- transect_smooth %>%
        group_by(ObjectID) %>%
        # Summary statistics
        summarize(b.tran.mean = mean(Intensity),
                  b.tran.sd = sd(Intensity),
                  b.tran.mad = mad(Intensity),
                  b.center = Intensity[1], # center intensity
                  b.periphery = last(Intensity), # periphery intensity
        ) %>%
        mutate(b.diff.cp = b.periphery - b.center) %>%
        # Transect
        left_join(calculate_transect(transect_smooth), by = "ObjectID") %>%
        rename(t.bump.number = Count, # number of transect bumps
               t.bump.onset = OnsetBump) # onset of the first transect bum
    cat("\n\ntransect features")

    # 7.4 Compute feature. The output table is NULL if no object (no colony)
    object_feature <- computeFeatures(
        image_watershed2, image_rolled,
        methods.noref = c("computeFeatures.shape"),
        methods.ref = c("computeFeatures.basic", "computeFeatures.moment"),
        #basic.quantiles = c(0.01, 0.05, seq(0.1, 0.9, by = .1), 0.95, 0.99)
        basic.quantiles = c(0.01, 0.05, c(0.1, 0.2, 0.5, 0.8, 0.9), 0.95, 0.99)
    )

    # 7.5 Clean up object feature name
    object_feature <- object_feature %>%
        as_tibble(rownames = "ObjectID") %>%
        # Remove duplicatedly calculated properties
        select(ObjectID, starts_with("x.0"), starts_with("x.Ba")) %>%
        # Remove the redundant prefix
        rename_with(function(x) str_replace(x,"x.0.", ""), starts_with("x.0")) %>%
        rename_with(function(x) str_replace(x,"x.Ba.", ""), starts_with("x.Ba")) %>%
        # Join the transect  features
        left_join(transect_feature, by = "ObjectID") %>%
        select(ObjectID, starts_with("b."), starts_with("t."),  starts_with("s."), starts_with("m."))


    # 7.6 remove outliers, only when its an isolate images
    #' `DO NOT REMOVE THE OUTLIERS in the coculture iamges`
    if (length(str_split(image_name, "_")[[1]]) == 4) object_feature <- object_feature %>% remove_outliers()
    write_csv(object_feature, paste0(paste_folder_name("feature", color_channel), image_name, ".csv"))

    # 7.7 Mark the transect and contours
    ## Transects
    transect_smooth_outlier <- transect_smooth %>% filter(!(ObjectID %in% object_feature$ObjectID))
    transect_smooth <- transect_smooth %>% filter(ObjectID %in% object_feature$ObjectID)
    write_csv(transect_smooth, file = paste0(paste_folder_name("transect", color_channel), image_name, ".csv"))


    # 7.8 Mark the transects and contour
    ## Transect
    image_transect <- image_rolled %>%
        draw_pixels(transect_smooth$x, transect_smooth$y, color = "red") %>% # Draw transects of normal objects
        draw_pixels(transect_smooth_outlier$x, transect_smooth_outlier$y, color = "blue") # Draw transects of outliers
    #display(image_transect, method = "raster")
    ## Contour
    image_watershed3_outliers <- rmObjects(image_watershed2, unique(transect_smooth$ObjectID), reenumerate = F)
    image_watershed3 <- rmObjects(image_watershed2, unique(transect_smooth_outlier$ObjectID), reenumerate = F)
    image_transect <- paintObjects(image_watershed3, image_transect, col = "red") # Draw contours of normal objects
    image_transect <- paintObjects(image_watershed3_outliers, image_transect, col = "blue") # Draw contours of outliers
    #
    writeImage(image_transect, paste0(paste_folder_name("transect", color_channel), image_name, ".tiff"))
    cat("\tdraw contours and transects", i, "/", nrow(list_images), "\t", image_name)

}
