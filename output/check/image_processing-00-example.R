library(tidyverse)
library(EBImage)

# This main folder depends on your home directory and user name. Python somehow does not read ~/ instead I have to specify /Users/chang-yu/
#folder_main <- "/Users/chang-yu/Dropbox/lab/emergent-coexistence/data/raw/plate_scan/emergent_coexistence_plate_scan_check/"
folder_main <- "~/Dropbox/lab/emergent-coexistence/data/raw/plate_scan/emergent_coexistence_plate_scan_check/"
folder_script <- "~/Desktop/Lab/emergent-coexistence/output/check/"

#' Example pair
#'  D_T8_C1R2_1
#'  D_T8_C1R2_2
#'  D_T8_C1R2_5-95_1_2

# 00. generate two mapping files ----
## List of images
list_images <- tibble(
    #image_name = c("D_T8_C1R2_1", "D_T8_C1R2_2", "D_T8_C1R2_5-95_1_2"),
    #image_name = c("D_T8_C1R2_1", "D_T8_C1R2_3", "D_T8_C1R2_5-95_1_3"),
    image_name = c(c("D_T8_C1R2_1", "D_T8_C1R2_4", "D_T8_C1R2_5-95_1_4"),
                   c("D_T8_C1R6_2", "D_T8_C1R6_4", "D_T8_C1R6_5-95_2_4"),
                   c("D_T8_C1R6_1", "D_T0_C1R6_3", "D_T8_C1R6_5-95_1_3"))
) %>%
    mutate(
        folder_original = rep(paste0(folder_main, "check/D-00-original/"), n()),
        folder_green = rep(paste0(folder_main, "check/D-01-green_channel/"), n()),
        folder_green_rolled = rep(paste0(folder_main, "check/D-02-green_rolled/"), n()),
        # folder_green_threshold = rep(paste0(folder_main, "check/D-03-green_threshold/"), n()),
        # folder_green_round = rep(paste0(folder_main, "check/D-04-green_round/"), n()),
        # folder_green_watershed_file = rep(paste0(folder_main, "check/D-05-green_watershed_file/"), n()),
        folder_green_watershed_file = rep(paste0(folder_main, "check/D-05-green_watershed_file/"), n()),
        folder_green_watershed = rep(paste0(folder_main, "check/D-06-green_watershed/"), n()),
        folder_green_feature = rep(paste0(folder_main, "check/D-07-green_feature/"), n()),
        folder_green_transection = rep(paste0(folder_main, "check/D-08-green_transection/"), n()),
        folder_green_cluster = rep(paste0(folder_main, "check/D-09-green_cluster/"), n()),
    )


## List of image mapping between pair and isolate
list_image_mapping <- tibble(
    Batch = "D",
    Community = "C1R2",
    Isolate1 = "1",
    Isolate2 = "2",
    Freq1 = "95",
    Freq2 = "5",
    image_name_pair = "D_T8_C1R2_5-95_1_2",
    image_name_isolate1 = "D_T8_C1R2_1",
    image_name_isolate2 = "D_T8_C1R2_2"
)


# 01. green channel and grey scale ----

for (i in 1:nrow(list_images)) {
    image_name <- list_images$image_name[i]

    # 0. original image
    image_original <- readImage(paste0(list_images$folder_original[i], image_name, ".tiff"))

    # 1. Green channel
    temp <- image_original
    colorMode(temp) = Grayscale
    image_green <- temp[,,2]
    writeImage(image_green, paste0(list_images$folder_green[i], image_name, ".tiff"))
    cat("\ngreen channel\t", i, "/", nrow(list_images), "\t", list_images$image_name[i])
}


# 02. background subtraction  ----
# rolling ball
# This will take a while. Usually one image takes ~1min
library(reticulate)
# This rolling ball takes the image_processing-01-example.csv as input
#py_run_file("image_processing-02-rolling_ball.py")


# 03. segmentation ----
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
list_images <- read_csv("~/Desktop/Lab/emergent-coexistence/output/check/00-list_images-D.csv", show_col_types = F)
i = which(list_images$image_name == "D_T8_C1R2_5-95_2_4")
for (i in 1:nrow(list_images)) {
    image_name <- list_images$image_name[i]

    # Rolled image
    image_rolled <- readImage(paste0(list_images$folder_green_rolled[i], image_name, ".tiff"))

    # 3. Thresholding
    threshold <- otsu(image_rolled)
    #threshold <- 0.9
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
    #display(image_watershed2, method = "raster")
    save(image_watershed2, file = paste0(list_images$folder_green_watershed_file[i], image_name, ".RData")) # save watersed image object
    writeImage(colorLabels(image_watershed2), paste0(list_images$folder_green_watershed[i], image_name, ".tiff"))
    cat("\twatershed\t", i, "/", nrow(list_images), "\t", list_images$image_name[i])

}

# temp <- read_csv(paste0(folder_main, "check/D-07-green_feature/D_T8_C1R7_50-50_3_6.csv"), show_col_types = F)
# temp %>%
#     mutate(Circularity = 4 * pi * s.area / s.perimeter^2) %>%
#     filter(Circularity > 0.3) %>%
#     mutate(s.ratio = s.radius.sd/s.radius.mean) %>%
#     select(starts_with("s."), Circularity) %>%
#     view

# i for image_names
# j for objects
# k for pixels


# 04. calculate the features ----
library(EBImageExtra)
library(purrr)
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

for (i in 1:nrow(list_images)) {
    image_name <- list_images$image_name[i]
    image_rolled <- readImage(paste0(list_images$folder_green_rolled[i], image_name, ".tiff"))
    load(paste0(list_images$folder_green_watershed_file[i], image_name, ".RData")) # this should contain one R object image_watershed2

    # 7.1 Extract transection data
    #transections[[i]] <- extract_transection(image_watershed2, image_rolled) %>%
    transection <- extract_transection(image_watershed2, image_rolled) %>%
        lapply(function(x) mutate(x, DistanceToCenter = 1:length(x))) %>%
        bind_rows(.id = "ObjectID")
    write_csv(transection, file = paste0(list_images$folder_green_transection[i], image_name, ".csv"))
    cat("\ntransection")

    # Save the transection image
    image_transect <- draw_pixels(image_rolled, transection$x, transection$y)
    writeImage(image_transect, paste0(list_images$folder_green_transection[i], image_name, ".tiff"))
    cat("\ndraw transection")

    # 7.2 Smooth the curves using loess
    transection_smooth <- transections[[i]] %>%
        nest(data = -ObjectID) %>%
        mutate(mod = map(data, loess, formula = Intensity ~ DistanceToCenter),
               fitted = map(mod, `[[`, "fitted")) %>%
        select(-mod) %>%
        unnest(cols = c(data, fitted))
    cat("\nsmooth curve")

    # transection_smooth %>%
    #     #filter(ObjectID == 5) %>%
    #     ggplot() +
    #     geom_line(aes(x = DistanceToCenter, y = fitted, group = ObjectID), lwd = .3) +
    #     #geom_point(data = transections[[1]], aes(x = DistanceToCenter, y = Intensity, group = ObjectID)) +
    #     theme_classic() +
    #     guides(color = "none") +
    #     labs()

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
    #load(paste0(list_images$folder_green_watershed_file, image_name, ".RData")) # this should contain one R object image_watershed2
    object_feature <- computeFeatures(
        image_watershed2, image_rolled,
        methods.noref = c("computeFeatures.shape"),
        methods.ref = c("computeFeatures.basic", "computeFeatures.moment"),
        basic.quantiles = c(0.01, 0.05, seq(0.1, 0.9, by = .1), 0.95, 0.99)
    )

    ## Execute the name cleanup only if there is at least 1 object
    if (is_null(object_feature)) {
        #write_csv(object_feature, paste0(list_images$folder_green_feature[i], image_name, ".csv"))
        cat("\tno object\t", i, "/", nrow(list_images), "\t", list_images$image_name[i])
    }

    if (!is_null(object_feature)) {
        object_feature <- object_feature %>%
            as_tibble(rownames = "ObjectID") %>%
            # Remove duplicatedly calculated properties
            select(ObjectID, starts_with("x.0"), starts_with("x.Ba")) %>%
            # Remove the redundant prefix
            rename_with(function(x) str_replace(x,"x.0.", ""), starts_with("x.0")) %>%
            rename_with(function(x) str_replace(x,"x.Ba.", ""), starts_with("x.Ba")) %>%
            #
            select(ObjectID, starts_with("b."), starts_with("s."), starts_with("m.")) %>%
            # # Remove the round plate lid crack that has high intensity variability
            # filter(b.sd < 0.3)
            # Join the transection summary statistic features
            left_join(transection_feature, )

        write_csv(object_feature, paste0(list_images$folder_green_feature[i], image_name, ".csv"))
        cat("\tfeature\t", i, "/", nrow(list_images), "\t", list_images$image_name[i])
    }
}


"finish this
- smooth the transection curve -> done
also these values have to match the Object ID
1) fit the curve to have a value for concavity
2) normalize the size and calculate the roungh slope
3) the number of valleys
this can wait. Using the new feature, do the model selection for clustering

run a set of model selection algorithm
"


# 05. clustering ----

# 8. Clustering

#' 1. stepwise selection over the best subset selection: run regsubset to get the best models for p variables.
#' 3. compare models: divide the data into k-fold cross-validaiton. Find the one with lowest MSE
library(metafor) # a meta-analysis package
library(leaps) # for computing stepwise regression. But it only fits lm
library(glmulti) # an extension to include glm in feature selection
library(gridExtra) # for plotting tables on the graph


list_images <- read_csv("~/Desktop/Lab/emergent-coexistence/output/check/00-list_images-D.csv", show_col_types = F)
list_image_mapping <- read_csv(paste0(folder_script, "00-list_image_mapping-D.csv"), show_col_types = F)
list_image_mapping_folder <- list_image_mapping %>%
    left_join(select(list_images, image_name_pair = image_name, folder_feature_pair = folder_green_feature, folder_green_cluster), by = "image_name_pair") %>%
    left_join(select(list_images, image_name_isolate1 = image_name, folder_feature_isolate1 = folder_green_feature), by = "image_name_isolate1") %>%
    left_join(select(list_images, image_name_isolate2 = image_name, folder_feature_isolate2 = folder_green_feature), by = "image_name_isolate2")

## 8.1 read the feature files
i = 1
read_feature_combined <- function () {
    object_feature_pair <- paste0(
        list_image_mapping_folder$folder_feature_pair[i],
        list_image_mapping_folder$image_name_pair[i], ".csv"
    ) %>%
        read_csv(show_col_types = FALSE) %>%
        mutate(image_name = list_image_mapping_folder$image_name_pair[i],
               Group = "pair") # Label for supervised learning

    object_feature_isolate1 <- paste0(
        list_image_mapping_folder$folder_feature_isolate1[i],
        list_image_mapping_folder$image_name_isolate1[i], ".csv"
    ) %>%
        read_csv(show_col_types = FALSE) %>%
        mutate(image_name = list_image_mapping_folder$image_name_isolate1[i],
               Group = "isolate1") # Label for supervised learning

    object_feature_isolate2 <- paste0(
        list_image_mapping_folder$folder_feature_isolate2[i],
        list_image_mapping_folder$image_name_isolate2[i], ".csv"
    ) %>%
        read_csv(show_col_types = FALSE) %>%
        mutate(image_name = list_image_mapping_folder$image_name_isolate2[i],
               Group = "isolate2") # Label for supervised learning

    #
    object_feature_combined <- bind_rows(object_feature_pair, object_feature_isolate1, object_feature_isolate2) %>%
        dplyr::select(image_name, Group, everything())

    return(object_feature_combined)
}
list_image_mapping_folder$image_name_pair[i]

feature_candidates <- c("b.mean", "b.sd", "b.mad",
                        "s.area", "s.radius.mean", "s.radius.sd",
                        "b.tran.mean", "b.tran.sd", "b.tran.mad",
                        "b.center", "b.periphery", "b.diff.cp")

object_feature_combined <- read_feature_combined() %>%
    mutate(GroupBinary = case_when(
        Group == "isolate1" ~ 0,
        Group == "isolate2" ~ 1,
    )) %>%
    # Start with these parameters
    dplyr::select(image_name, ObjectID, Group, GroupBinary, all_of(feature_candidates))
object_feature_isolates <- object_feature_combined %>%
    filter(!is.na(GroupBinary))
object_feature_pair <- object_feature_combined %>%
    filter(is.na(GroupBinary))

# 8.2 exhaustive stepwise selection and cross-validation
best_subset <- glmulti(
    y = "GroupBinary",
    xr = feature_candidates,
    data = object_feature_isolates,
    level = 1,
    method = "h",
    crit = "aicc",
    minsize = 2,
    maxsize = 5,
    confsetsize = 100, # Keep the top n models
    plotty = F, report = F,
    fitfunction = glm
    #family = binomial # logistic
)

#plot(best_subset)
# Top models
weightable(best_subset) %>%  as_tibble()

# Multimodel inference, using the metafor package
eval(metafor:::.glmulti)
compute_importance <- function (best_subset) {
    #' This function is from https://www.metafor-project.org/doku.php/tips:model_selection_with_glmulti_and_mumin
    #' (see detailed interpretation in the webpage)
    #' where it is used to calculat the model-averaged values for all predictors
    #' across the 100 models from the best-subset selection.
    #' This should give a good inference about the predictors, not in the context
    #' of a single model that is declared to be the best, but across all possible models
    #' (taking their relative weights into consideration)
    mmi <- as.data.frame(coef(best_subset, varweighting="Johnson"))
    mmi <- data.frame(Estimate=mmi$Est, SE=sqrt(mmi$Uncond), Importance=mmi$Importance, row.names=row.names(mmi))
    mmi$z <- mmi$Estimate / mmi$SE
    mmi$p <- 2*pnorm(abs(mmi$z), lower.tail=FALSE)
    names(mmi) <- c("Estimate", "Std. Error", "Importance", "z value", "Pr(>|z|)")
    mmi$ci.lb <- mmi[[1]] - qnorm(.975) * mmi[[2]]
    mmi$ci.ub <- mmi[[1]] + qnorm(.975) * mmi[[2]]
    mmi <- mmi[order(mmi$Importance, decreasing=TRUE), c(1,2,4:7,3)]
    round(mmi, 4)
    return(as_tibble(mmi, rownames = "Feature"))
}
model_importance <- compute_importance(best_subset) %>% filter(Feature != "(Intercept)")
p1 <- model_importance %>%
    mutate(Feature = ordered(Feature, rev(Feature))) %>%
    ggplot() +
    geom_col(aes(x = Feature, y = Importance), color = 1, fill = grey(0.8)) +
    geom_hline(yintercept = 0.8, color = "red") +
    coord_flip() +
    theme_classic() +
    ggtitle("model-averaged importance of predictors")

# The final best model and table of coefficent
model <- best_subset@objects[[1]]
model_table <- broom::tidy(model) %>%
    mutate(estimate = round(estimate, 4),
           std.error = round(std.error, 4),
           statistic = round(statistic, 4))
p2 <- cowplot::ggdraw() +
    annotation_custom(tableGrob(model_table, theme = ttheme_default(base_size = 6)), xmin = 0, xmax = 1, ymin = 0, ymax = 1) +
    theme_void() +
    theme(legend.position = "none")


# Model predict
object_feature_predicted <- object_feature_pair %>%
    mutate(PredictedGroupProbability = predict(model, object_feature_pair)) %>%
    dplyr::select(image_name, ObjectID, PredictedGroupProbability) %>%
    # Criteria for catagorizing
    mutate(Group = case_when(
        PredictedGroupProbability < 0.4 ~ "predicted isolate1",
        PredictedGroupProbability > 0.6 ~ "predicted isolate2",
        PredictedGroupProbability > 0.4 & PredictedGroupProbability < 0.6 ~ "undecided"
    )) %>%
    mutate(Group = factor(Group, c("predicted isolate1", "predicted isolate2", "undecided"),))

object_feature_predicted_count <- object_feature_predicted %>%
    group_by(Group, .drop = F) %>%
    count(name = "Count") %>%
    ungroup() %>%
    mutate(Group = Group, Halign = c("right", "left", "center"), Valign = c(2, 2, 4))

p3 <- object_feature_predicted %>%
    ggplot() +
    geom_histogram(aes(x = PredictedGroupProbability), color = 1, fill = NA, bins = 30) +
    geom_vline(xintercept = c(0.4, 0.6), linetype = 2, color = "red") +
    geom_vline(xintercept = c(0,1), linetype = 2, color = "black") +
    geom_text(data = object_feature_predicted_count, x = c(.4,.6,.5), aes(hjust = Halign, vjust = Valign, label = paste0("  ", Group, ": ", Count, "  ")), y = Inf) +
    scale_x_continuous(expand = c(0,0.1), breaks = seq(-2,2, .2)) +
    scale_y_continuous(expand = c(0,0)) +
    theme_classic()

cat("\tprediction")

# Scatterplot for clustering intuition
set_color_names <- function () {
    c("#FF5A5F","#087E8B","#FF5A5F","#087E8B", "black") %>%
        setNames(c(
            paste0(list_image_mapping_folder$image_name_isolate1[i], " isolate1"),
            paste0(list_image_mapping_folder$image_name_isolate2[i], " isolate2"),
            paste0(list_image_mapping_folder$image_name_pair[i], " predicted isolate1"),
            paste0(list_image_mapping_folder$image_name_pair[i], " predicted isolate2"),
            paste0(list_image_mapping_folder$image_name_pair[i], " undecided")
        ))
}
set_fill_names <- function () {
    # c("white","white","#FF5A5F","#087E8B", "black") %>%
    # c("white","white","white","white", "white") %>%
    rep(NA, 5) %>%
        setNames(c(
            paste0(list_image_mapping_folder$image_name_isolate1[i], " isolate1"),
            paste0(list_image_mapping_folder$image_name_isolate2[i], " isolate2"),
            paste0(list_image_mapping_folder$image_name_pair[i], " predicted isolate1"),
            paste0(list_image_mapping_folder$image_name_pair[i], " predicted isolate2"),
            paste0(list_image_mapping_folder$image_name_pair[i], " undecided")
        ))
}
set_shape_names <- function () {
    # c(21,21,21,21,21) %>% # solid point with fill
    rep(1, 5) %>% # Hallow point
        setNames(c(
            paste0(list_image_mapping_folder$image_name_isolate1[i], " isolate1"),
            paste0(list_image_mapping_folder$image_name_isolate2[i], " isolate2"),
            paste0(list_image_mapping_folder$image_name_pair[i], " predicted isolate1"),
            paste0(list_image_mapping_folder$image_name_pair[i], " predicted isolate2"),
            paste0(list_image_mapping_folder$image_name_pair[i], " undecided")
        ))
}
set_alpha_names <- function () {
    c(.3,.3,1,1,1) %>% setNames(c(
        paste0(list_image_mapping_folder$image_name_isolate1[i], " isolate1"),
        paste0(list_image_mapping_folder$image_name_isolate2[i], " isolate2"),
        paste0(list_image_mapping_folder$image_name_pair[i], " predicted isolate1"),
        paste0(list_image_mapping_folder$image_name_pair[i], " predicted isolate2"),
        paste0(list_image_mapping_folder$image_name_pair[i], " undecided")
    ))
}
color_names <- set_color_names()
fill_names <- set_fill_names()
shape_names <- set_shape_names()
alpha_names <- set_alpha_names()

## combined the prediction and object label
object_feature_plot <- object_feature_predicted %>%
    select(image_name, ObjectID, Group) %>%
    left_join(select(object_feature_pair, -GroupBinary, -Group), by = c("image_name", "ObjectID")) %>%
    bind_rows(select(object_feature_isolates, -GroupBinary)) %>%
    mutate(ColorLabel = paste(image_name, Group),
           FillLabel = paste(image_name, Group),
           ShapeLabel = paste(image_name, Group),
           AlphaLabel = paste(image_name, Group))

p4 <-  object_feature_plot %>%
    arrange(ColorLabel) %>%
    ggplot() +
    geom_point(aes_string(x = model_importance$Feature[1], y = model_importance$Feature[2],
                          color = "ColorLabel", fill = "FillLabel", shape = "ShapeLabel", alpha = "AlphaLabel"),
               size = 2, stroke = .8) +
    scale_color_manual(values = color_names, name = "", label = names(color_names)) +
    scale_fill_manual(values = fill_names, name = "", label = names(color_names)) +
    scale_shape_manual(values = shape_names, name = "", label = names(color_names)) +
    scale_alpha_manual(values = alpha_names, name = "", label = names(color_names)) +
    theme_classic()
cat("\tclusterplot")

# PCA with the model selected variables
pcobj <- object_feature_plot %>%
    select(all_of(names(coef(model))[-1])) %>% # names(coef(model))[-1] # remove intercept term
    prcomp(center = TRUE, scale. = TRUE)
compute_pca_coord <- function (pcobj) {
    #' The function below comes from the source code of ggbiplot
    #' to replace the use of ggbiplot
    choices = 1:2 # which PCs to plot
    scale = 1
    obs.scale = 1 - scale
    var.scale = scale

    # Recover the SVD
    nobs.factor <- sqrt(nrow(pcobj$x) - 1)
    d <- pcobj$sdev
    u <- sweep(pcobj$x, 2, 1 / (d * nobs.factor), FUN = '*')
    v <- pcobj$rotation
    # Scores
    choices <- pmin(choices, ncol(u))
    df.u <- as.data.frame(sweep(u[,choices], 2, d[choices]^obs.scale, FUN='*'))
    # Directions
    v <- sweep(v, 2, d^var.scale, FUN='*')
    df.v <- as.data.frame(v[, choices])

    names(df.u) <- c('xvar', 'yvar')
    names(df.v) <- names(df.u)
    df.u <- df.u * nobs.factor

    # Variable Names
    df.v$varname <- rownames(v)
    # Variables for text label placement
    df.v$angle <- with(df.v, (180/pi) * atan(yvar / xvar))
    df.v$hjust = with(df.v, (1 - 1.5 * sign(xvar)) / 2)

    # Change the labels for the axes
    # Append the proportion of explained variance to the axis labels
    u.axis.labs <- paste('standardized PC', choices, sep='')
    u.axis.labs <- paste(u.axis.labs, sprintf('(%0.1f%% explained var.)', 100 * pcobj$sdev[choices]^2/sum(pcobj$sdev^2)))

    return(list(df.v = tibble(df.v), # Variables
                df.u = tibble(df.u), # Score
                u.axis.labs = u.axis.labs))
}

temp <- compute_pca_coord(pcobj)
df.u <- temp$df.u
df.v <- temp$df.v
u.axis.labs <- temp$u.axis.labs

p5 <- df.u %>%
    bind_cols(object_feature_plot) %>%
    arrange(ColorLabel) %>%
    ggplot() +
    # Draw directions of axis
    geom_segment(data = df.v, aes(x = 0, y = 0, xend = xvar, yend = yvar), arrow = arrow(length = unit(1/2, 'picas')), color = scales::muted('red')) +
    # Draw scores
    geom_point(aes(x = xvar, y = yvar,
                   color = ColorLabel, fill = FillLabel, shape = ShapeLabel, alpha = AlphaLabel),
               size = 2, stroke = .8) +
    # Label the variable axes
    geom_text(data = df.v, aes(label = varname, x = xvar, y = yvar, angle = angle, hjust = hjust), color = 'blue', size = 3) +
    scale_color_manual(values = color_names, name = "", label = names(color_names)) +
    scale_fill_manual(values = fill_names, name = "", label = names(color_names)) +
    scale_shape_manual(values = shape_names, name = "", label = names(color_names)) +
    scale_alpha_manual(values = alpha_names, name = "", label = names(color_names)) +
    theme_classic() +
    labs(x = u.axis.labs[1], y = u.axis.labs[2])
cat("\tpcaplot")


library(cowplot)
legend <- get_legend(p4 + theme(legend.box.margin = margin(0, 0, 0, 12)))
p <- plot_grid(
    p1, p2, p3,
    p4 + theme(legend.position = "none"),
    p5 + theme(legend.position = "none"),
    legend,
    nrow = 2, rel_widths = c(1,1,1),
    align = "h", axis = " tblr"
) + theme(plot.background = element_rect(fill = "white"))

ggsave(filename = paste0(list_images$folder_green_cluster[i], list_images$image_name[i], ".png"), plot = p, width = 10, height = 6)
cat("\nplot feature\t", i, "/", nrow(list_image_mapping_folder), "\t", list_images$image_name[i])



# Test plot for features selected
if (FALSE) {
    plot_one_feature <- function (object_feature_combined, feature1) {
        object_feature_combined %>%
            ggplot() +
            geom_histogram(aes_string(x = feature1, fill = "image_name"), color = 1, bins = 50, alpha = .5, position = "identity") +
            scale_fill_manual(values = color_names) +
            theme_classic() +
            theme(legend.position = "none")
    }
    plot_two_features <- function (object_feature_combined, feature1, feature2) {
        object_feature_combined %>%
            ggplot() +
            geom_point(aes_string(x = feature1, y = feature2, color = "image_name"), size = 1, shape = 1, stroke = .5) +
            scale_color_manual(values = color_names) +
            theme_classic()
    }
    subset_training_set <- function (x) {
        x %>%
            filter(Group != "pair") %>%
            # Dummy variable for logistic regression
            mutate(dummy = case_when(
                Group == "isolate1" ~ 0,
                Group == "isolate2" ~ 1
            ))
    }
    p1 <- plot_two_features(object_feature_combined, "b.mean", "b.sd")
    p2 <- plot_two_features(object_feature_combined, "b.mean", "b.mad")
    p3 <- plot_two_features(object_feature_combined, "b.sd", "b.mad")

    # Featured selected
    results %>%
        select(starts_with("b."), starts_with("s.")) %>%
        mutate(NumberOfPredictor = 1:n()) %>%
        pivot_longer(cols = -NumberOfPredictor, names_to = "Feature", values_to = "Value") %>%
        mutate(Feature = factor(Feature, names(object_feature_isolates)[-1])) %>%
        filter(Value) %>%
        ggplot() +
        geom_point(aes(x = NumberOfPredictor, y = Feature), size = 4) +
        scale_x_continuous(breaks = 1:10) +
        scale_y_discrete(drop = F) +
        #scale_y_discrete(breaks = names(object_feature_isolates)[-1]) +
        theme_classic() +
        theme(panel.grid.major = element_line(color = grey(0.9)))

    # fitting scores across the models
    results %>%
        mutate(NumberOfPredictors = 1:n()) %>%
        select(NumberOfPredictors, adj.r.squared, mallows_cp, BIC) %>%
        pivot_longer(cols = -NumberOfPredictors, names_to = "Statistic", values_to = "Value") %>%
        ggplot(aes(x = NumberOfPredictors, y = Value, color = Statistic)) +
        geom_line() +
        geom_point() +
        facet_wrap(~ Statistic, scale = "free") +
        theme_classic() +
        guides(color = "none")

}
# glmulti(
#     y = "GroupBinary",
#     #xr = c("b.mean"),
#     xr = names(object_feature_isolates)[-1],
#     #xr = str_subset(names(object_feature_isolates), "^b."),
#     # xr = names(object_feature_combined)[names(object_feature_isolates) != c("Group", "ObjectID")],
#     data = object_feature_isolates,
#     method = "h",
#     crit = "aic",
#     confsetsize = 5, # Keep 5 best models
#     plotty = F, report = F,
#     fitfunction = "glm"
# )
if (FALSE) {


    # 8.3 cross-validation
    # create training - testing data
    set.seed(5)
    sample <- sample(c(TRUE, FALSE), nrow(object_feature_isolates), replace = T, prob = c(0.6,0.4))
    train <- object_feature_isolates[sample, ]
    test <- object_feature_isolates[!sample, ]

    # perform best subset selection
    best_subset <- regsubsets(Salary ~ ., train, nvmax = 19)
    ## build an “X” matrix from data
    test_m <- model.matrix(Salary ~ ., data = test)
    ## create empty vector to fill with error values
    validation_errors <- vector("double", length = 19)

    for(j in 1:19) {
        j=1
        coef_x <- coef(best_subset, id = j)                     # extract coefficients for model size i
        pred_x <- test_m[ , names(coef_x)] %*% coef_x           # predict salary using matrix algebra
        validation_errors[j] <- mean((test$Salary - pred_x)^2)  # compute test error btwn actual & predicted salary
    }
    plot(validation_errors, type = "b")

}



# Quantile
if (FALSE) {

    plot_two_features(object_feature_combined, "b.q095", "b.q05")
    quantile_name_mapping <- tibble(Quantile = 1:999, QuantileColumn = str_replace(paste0("b.q", sprintf("%03d", 1:999)), "0+$", ""))

    object_feature_quantile <- object_feature_combined %>%
        select(image_name, ObjectID, Group, starts_with("b.q")) %>%
        pivot_longer(cols = starts_with("b."), names_to = "QuantileColumn", values_to = "Value") %>%
        left_join(quantile_name_mapping) %>% select(-QuantileColumn)


    object_feature_quantile %>%
        unite(col = "UniqueID", Group, ObjectID) %>%
        ggplot(aes(x = Quantile, y = Value, color = image_name, group = UniqueID)) +
        geom_point(alpha = 0.3, size = .5) +
        geom_line(alpha = 0.3, size = .5) +
        scale_color_manual(values = color_names) +
        theme_classic()

    object_feature_quantile %>%
        mutate(Quantile = factor(Quantile)) %>%
        ggplot() +
        geom_boxplot(aes(x = Quantile, y = Value, color = image_name)) +
        geom_point(aes(x = Quantile, y = Value, color = image_name), position = position_jitterdodge(jitter.width = .1), size = .2) +
        scale_color_manual(values = color_names) +
        theme_classic()
}
















