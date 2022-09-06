library(tidyverse)
library(EBImage)
library(EBImageExtra)
library(reticulate)
library(tidypredict)
library(gridExtra)

# This main folder depends on your home directory and user name. Python somehow does not read ~/ instead I have to specify /Users/chang-yu/
#folder_main <- "/Users/chang-yu/Dropbox/lab/emergent-coexistence/data/raw/plate_scan/emergent_coexistence_plate_scan_check/"
folder_main <- "~/Dropbox/lab/emergent-coexistence/data/raw/plate_scan/emergent_coexistence_plate_scan_check/"


#' Example pair
#'  D_T8_C1R2_1
#'  D_T8_C1R2_2
#'  D_T8_C1R2_5-95_1_2

# 0. Generate two mapping file
## List of images
list_images <- tibble(
    #image_name = c("D_T8_C1R2_1", "D_T8_C1R2_2", "D_T8_C1R2_5-95_1_2"),
    #image_name = c("D_T8_C1R2_1", "D_T8_C1R2_3", "D_T8_C1R2_5-95_1_3"),
    image_name = c("D_T8_C1R2_1", "D_T8_C1R2_4", "D_T8_C1R2_5-95_1_4"),
    folder_original = rep(paste0(folder_main, "check/D-00-original/"), 3),
    folder_green = rep(paste0(folder_main, "check/D-01-green_channel/"), 3),
    folder_green_rolled = rep(paste0(folder_main, "check/D-02-green_rolled/"), 3),
    #folder_green_threshold = rep(paste0(folder_main, "check/D-03-green_threshold/"), 3),
    #folder_green_round = rep(paste0(folder_main, "check/D-04-green_round/"), 3),
    folder_green_watershed_file = rep(paste0(folder_main, "check/D-05-green_watershed_file/"), 3),
    folder_green_watershed = rep(paste0(folder_main, "check/D-06-green_watershed/"), 3),
    folder_green_feature = rep(paste0(folder_main, "check/D-07-green_feature/"), 3),
    folder_green_cluster = rep(paste0(folder_main, "check/D-08-green_cluster/"), 3),
)



i = 1
image_name <- list_images$image_name[i]


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


# 2. Rolling ball. This will take a while. Usually one image takes ~1min
# This rolling ball takes the image_processing-01-example.csv as input
#py_run_file("image_processing-02-rolling_ball.py")


# 3-6. segmentation
detect_nonround_object <- function (image_object) {
        object_shape <- computeFeatures.shape(image_object) %>% as_tibble(rownames = "ObjectID")
        object_shape_round <- object_shape %>%
            # Area
            filter(s.area > 800 & s.area < 20000) %>%
            # Remove tape and label that has really large variation in radius
            filter(s.radius.sd < 10) %>%
            # Circularity = 1 means a perfect circle and goes down to 0 for non-cicular shapes
            mutate(Circularity = 4 * pi * s.area / s.perimeter^2) %>%
            filter(Circularity > 0.3) %>%
            arrange(desc(s.area))
        object_ID_nonround <- object_shape$ObjectID[!(object_shape$ObjectID %in% object_shape_round$ObjectID)]
        return(object_ID_nonround)
}

for (i in 1:nrow(list_images)) {

    image_name <- list_images$image_name[i]

    # Rolled image
    image_rolled <- readImage(paste0(list_images$folder_green_rolled[i], image_name, ".tiff"))

    # 3. Thresholding
    threshold <- otsu(image_rolled)
    image_thresholded <- image_rolled < threshold
    cat("\nthreshold")

    # 4. Detect round shaped object and remove super small size
    image_object <- bwlabel(image_thresholded)
    object_ID_nonround <- detect_nonround_object(image_object)
    image_round <- rmObjects(image_object, object_ID_nonround, reenumerate = T)
    cat("\tround object")

    # 5-6. Watershed
    image_distancemap <- distmap(image_round)
    cat("\tdistance map")
    image_watershed <- watershed(image_distancemap, tolerance = 1)
    ## After watershed, apply a second filter removing objects that are too small to be colonies
    object_ID_nonround2 <- detect_nonround_object(image_watershed)
    image_watershed2 <- rmObjects(image_watershed, object_ID_nonround2, reenumerate = T)
    save(image_watershed2, file = paste0(list_images$folder_green_watershed_file[i], image_name, ".RData")) # save watersed image object
    writeImage(colorLabels(image_watershed2), paste0(list_images$folder_green_watershed[i], image_name, ".tiff")) # save
    cat("\twatershed")

}



# i for image_names
# j for objects
# k for pixels
# 7. Calculate the features
## Transection

extract_transection <- function (watershed, ref) {
    #' Arguments:
    #'  watershed is the watershed image with object labels
    #'  ref is the original image with intensity data. This has to be single-channel
    oc <- ocontour(watershed)
    transection <- NULL # should be the number of objects
    #max_radius <- NULL # A test for object ID -> good. As long I use the oc output, they should be good
    for (j in 1:length(oc)) {
        z <- oc[[j]] # coordinate of contour pixels of a point
        cz <- apply(z, 2, mean) # object center
        radius = sqrt(rowSums((z - rep(cz, each=nrow(z)))^2)) # all radius from the center pixel to the contour pixels
        mz <- z[which.min(abs(radius - median(radius))),] # the coordinate of contour pixels with the median radius (or closest to the median)
        # max_radius[j] <- max(radius)
        # The coordinates of all along the line to the maximum line: this is arranged such that it's always from the center to the periphery
        lz <- bresenham(c(cz[1], mz[1]), c(cz[2], mz[2]))
        radial_gradient <- NULL
        for (k in 1:length(lz$x)) radial_gradient[k] <- ref@.Data[lz$x[k], lz$y[k]]
        transection[[j]] <- radial_gradient
    }
    return(transection)
}
transections <- NULL
for (i in 1:nrow(list_images)) {
    image_name <- list_images$image_name[i]
    image_rolled <- readImage(paste0(list_images$folder_green_rolled[i], image_name, ".tiff"))
    load(paste0(list_images$folder_green_watershed_file[i], image_name, ".RData")) # this should contain one R object image_watershed2

    transections[[i]] <- extract_transection(image_watershed2, image_rolled) %>%
        lapply(function(x) tibble(Intensity = x, DistanceToCenter = 1:length(x))) %>%
        bind_rows(.id = "ObjectID")
}

# This is 1-4
bind_rows(
    mutate(transections[[1]], Group = "isolate1"),
    mutate(transections[[2]], Group = "isolate2"),
    mutate(transections[[3]], Group = "pair")
) %>%
    group_by(Group, ObjectID) %>%
    mutate(DistanceToCenter = DistanceToCenter/max(DistanceToCenter)) %>%
    filter(Group != "pair") %>%
    ggplot() +
    geom_line(aes(x = DistanceToCenter, y = Intensity, color = Group, group = interaction(Group, ObjectID)), lwd = .3) +
    theme_classic() +
    #guides(color = "none") +
    labs()

transections[[1]] %>%
    group_by(ObjectID) %>%
    mutate(DistanceToCenter = DistanceToCenter/max(DistanceToCenter)) %>%
    ggplot() +
    geom_line(aes(x = DistanceToCenter, y = Intensity, group = ObjectID), lwd = .3) +
    theme_classic() +
    #guides(color = "none") +
    labs()



"
plot it on EBImage for each transection used
THis should be a line per colony
"

"finish this
- smooth the transection curve
also these values have to match the Object ID
1) fit the curve to have a value for concavity
2) normalize the size and calculate the slope
3) the number of valleys
"


# 7. Calculate feature
for (i in 1:nrow(list_images)) {
    i=1
    ## Compute feature. It is NULL if no object (no colony)
    load(paste0(list_images$folder_green_watershed_file, image_name, ".RData")) # this should contain one R object image_watershed2

    object_feature <- computeFeatures(
        image_watershed2, image_rolled,
        methods.noref = c("computeFeatures.shape"),
        methods.ref = c("computeFeatures.basic", "computeFeatures.moment")
    )

    ## Execute the name cleanup only if there is at least 1 object
    if (is_null(object_feature)) {
        #write_csv(object_feature, paste0(list_images$folder_green_feature, image_name, ".csv"))
        cat("\tno object\t", i, "/", nrow(list_images), "\t", list_images$image_name[i])
    }

    if (!is_null(object_feature)) {
        object_feature <- object_feature %>%
            as_tibble(rownames = "ObjectID") %>%
            # Remove duplicatedly calculated properties
            select(ObjectID, starts_with("x.0"), starts_with("x.Ba")) %>%
            # Remove the redundant prefix in the column names
            rename_with(function(x) str_replace(x,"x.0.", ""), starts_with("x.0")) %>%
            rename_with(function(x) str_replace(x,"x.Ba.", ""), starts_with("x.Ba")) %>%
            #
            select(ObjectID, starts_with("b."), starts_with("s."), starts_with("m."))
        write_csv(object_feature, paste0(list_images$folder_green_feature, image_name, ".csv"))
        cat("\tfeature\t", i, "/", nrow(list_images), "\t", list_images$image_name[i])
    }
}




# 8. Cluster
folder_script <- "/Users/chang-yu/Desktop/Lab/emergent-coexistence/output/check/"
batch <- "D"
list_images <- paste0(folder_script, "00-list_images-", batch, ".csv") %>% read_csv(show_col_types = F)
list_image_mapping <- paste0(folder_script, "00-list_image_mapping-", batch, ".csv") %>% read_csv(show_col_types = F)

list_image_mapping_folder <- list_image_mapping %>%
    left_join(select(list_images, image_name_pair = image_name, folder_feature_pair = folder_green_feature, folder_green_cluster), by = "image_name_pair") %>%
    left_join(select(list_images, image_name_isolate1 = image_name, folder_feature_isolate1 = folder_green_feature), by = "image_name_isolate1") %>%
    left_join(select(list_images, image_name_isolate2 = image_name, folder_feature_isolate2 = folder_green_feature), by = "image_name_isolate2")

## Read the feature files
#which(list_image_mapping_folder$image_name_pair == "D_T8_C1R2_5-95_1_3")
which(list_image_mapping_folder$image_name_pair == "D_T8_C1R4_5-95_2_3")
i = 54
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
        select(image_name, Group, everything())

    return(object_feature_combined)
}
set_color_names <- function() {
    color_names <- c(2, 3, 1) %>% setNames(
        c(list_image_mapping_folder$image_name_isolate1[i],
          list_image_mapping_folder$image_name_isolate2[i],
          list_image_mapping_folder$image_name_pair[i])
    )

}
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

object_feature_combined <- read_feature_combined()
color_names <- set_color_names()
#object_feature_combined$image_name %>% unique
# p0 <- plot_one_feature(object_feature_combined, "b.mean")
p1 <- plot_two_features(object_feature_combined, "b.mean", "b.sd")
p2 <- plot_two_features(object_feature_combined, "b.mean", "b.mad")
p3 <- plot_two_features(object_feature_combined, "b.sd", "b.mad")

# Logistic regression using the three main features
training <- subset_training_set(object_feature_combined)
model <- glm(dummy ~ b.mean + b.sd + b.mad, data = training, family = "binomial")

object_feature_predicted <- object_feature_combined %>%
    filter(Group == "pair") %>%
    tidypredict_to_column(model, vars = "PredictedGroupProbability") %>%
    select(image_name, ObjectID, PredictedGroupProbability)

object_feature_predicted_count <- object_feature_predicted %>%
    mutate(PredictedGroup = ifelse(PredictedGroupProbability < 0.5, 0, 1) %>% factor(c(0,1))) %>%
    group_by(PredictedGroup) %>%
    count(name = "Count") %>%
    ungroup() %>%
    mutate(Group = c("isolate1", "isolate2"), Align = c("right", "left"))

# Model predict
p4 <- object_feature_predicted %>%
    ggplot() +
    geom_histogram(aes(x = PredictedGroupProbability), color = 1, fill = NA, bins = 30) +
    geom_vline(xintercept = 0.5, linetype = 2, color = "red") +
    geom_text(data = object_feature_predicted_count, aes(x = 0.5, hjust = Align, label = paste0("  ", Group, ": ", Count, "  ")), y = Inf, vjust = 2) +
    scale_x_continuous(limits = c(0,1)) +
    theme_classic()

# Model fit
p5 <- ggplot() +
    geom_point() +
    annotation_custom(tableGrob(broom::tidy(model), theme = ttheme_default(base_size = 8)), xmin = .2, xmax = .8, ymin = .2, ymax = .8) +
    scale_x_continuous(limits = c(0,1)) +
    scale_y_continuous(limits = c(0,1)) +
    theme_void()


legend <- get_legend(p1 + theme(legend.box.margin = margin(0, 0, 0, 12)))
p <- plot_grid(p1 + theme(legend.position = "none"),
               p2 + theme(legend.position = "none"),
               p3 + theme(legend.position = "none"),
               p4 + theme(legend.position = "none"),
               p5 + theme(legend.position = "none"),
               legend, nrow = 2, align = "hv", axis = "lrtb") +
    theme(plot.background = element_rect(fill = "white"))
#p_output <- plot_grid(p, legend, rel_widths = c(2,1), nrow = 1) + theme(plot.background = element_rect(fill = "white"))
ggsave(filename = paste0(list_images$folder_green_cluster[i], image_name, ".png"), plot = p, width = 10, height = 6)
cat("\nplot feature\t", i, "/", nrow(list_image_mapping_folder), "\t", image_name)



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
















