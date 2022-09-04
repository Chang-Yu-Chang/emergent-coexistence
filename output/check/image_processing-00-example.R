library(tidyverse)
library(EBImage)
library(reticulate)
library(tidypredict)
library(gridExtra)

# This main folder depends on your home directory and user name. Python somehow does not read ~/ instead I have to specify /Users/chang-yu/
folder_main <- "/Users/chang-yu/Dropbox/lab/emergent-coexistence/data/raw/plate_scan/emergent_coexistence_plate_scan_check/"


#' Example pair
#'  D_T8_C1R2_1
#'  D_T8_C1R2_2
#'  D_T8_C1R2_5-95_1_2

# 0. Generate two mapping file
## List of images
list_images <- tibble(
    image_name = c("D_T8_C1R2_1", "D_T8_C1R2_2", "D_T8_C1R2_5-95_1_2"),
    folder_original = rep(paste0(folder_main, "check/D-00-original/"), 3),
    folder_green = rep(paste0(folder_main, "check/D-01-green_channel/"), 3),
    folder_rolled = rep(paste0(folder_main, "check/D-02-green_rolled/"), 3),
    folder_watershed = rep(paste0(folder_main, "check/D-06-watershed/"), 3),
    folder_feature = rep(paste0(folder_main, "check/D-07-feature/"), 3),
    folder_cluster = rep(paste0(folder_main, "check/D-08-cluster/"), 3),
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


for (i in 1:nrow(list_images)) {
    image_name <- list_images$image_name[i]
    folder_original <- list_images$folder_original[i]
    folder_green <- list_images$folder_green[i]
    folder_rolled <- list_images$folder_rolled[i]

    # 0. original image
    image_original <- readImage(paste0(folder_original, image_name, ".tiff"))

    # 1. Green channel
    temp <- image_original
    colorMode(temp) = Grayscale
    image_green <- temp[,,2]
    writeImage(image_green, paste0(folder_green, image_name, ".tiff"))
    cat("\ngreen channel\t", i, "/", nrow(list_images), "\t", list_images$image_name[i])
}


# 2. Rolling ball. This will take a while
# This rolling ball takes the image_processing-01-example.csv as input
py_run_file("image_processing-02-rolling_ball.py")


#
for (i in 1:nrow(list_images)) {
    image_name <- list_images$image_name[i]
    folder_rolled <- list_images$folder_rolled[i]
    folder_watershed <- list_images$folder_watershed[i]
    folder_feature <- list_images$folder_feature[i]

    # Rolled image
    image_rolled <- readImage(paste0(folder_rolled, image_name, ".tiff"))

    # 3. Thresholding
    threshold <- otsu(image_rolled)
    image_thresholded <- image_rolled < threshold
    cat("\nthreshold")

    # 4. Detect round shaped object and remove super small size
    image_object <- bwlabel(image_thresholded)
    detect_nonround_object <- function (image_object) {
        object_shape <- computeFeatures.shape(image_object) %>% as_tibble(rownames = "ObjectID")
        object_shape_round <- object_shape %>%
            # Area
            filter(s.area > 500 & s.area < 500000) %>%
            # Roundness = 1 means a perfect circle
            mutate(Roundness = (s.radius.max - s.radius.min)/2) %>%
            filter(Roundness > 0.1 & Roundness < 50) %>%
            # Circularity = 1 means a perfect circle and goes down to 0 for non-cicular shapes
            mutate(Circularity = 4 * pi * s.area / s.perimeter^2) %>%
            filter(Circularity > 0.3) %>%
            arrange(desc(s.area))
        object_ID_nonround <- object_shape$ObjectID[!(object_shape$ObjectID %in% object_shape_round$ObjectID)]
        return(object_ID_nonround)
    }
    object_ID_nonround <- detect_nonround_object(image_object)
    image_round <- rmObjects(image_object, object_ID_nonround, reenumerate = T)
    cat("\tround object")

    # 5. Distance map
    image_distancemap <- distmap(image_round)
    cat("\tdistance map")

    # 6. Watershed
    image_watershed <- watershed(image_distancemap, tolerance = 1)

    ## Second filter removing objects that are too small to be colonies
    object_ID_nonround2 <- detect_nonround_object(image_watershed)
    image_watershed2 <- rmObjects(image_watershed, object_ID_nonround2, reenumerate = T)
    writeImage(colorLabels(image_watershed2), paste0(folder_watershed, image_name, ".tiff"))
    cat("\twatershed")

    # 7. Calculate feature
    ## Compute feature. It is NULL if no object (no colony)
    object_feature <- computeFeatures(
        image_watershed2, image_rolled,
        methods.noref = c("computeFeatures.shape"),
        methods.ref = c("computeFeatures.basic", "computeFeatures.moment")
    )

    ## Execute the name cleanup only if there is at least 1 object
    if (is_null(object_feature)) {
        #write_csv(object_feature, paste0(folder_feature, image_name, ".csv"))
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
            select(ObjectID, starts_with("b."), starts_with("s."), starts_with("m."))
        write_csv(object_feature, paste0(folder_feature, image_name, ".csv"))
        cat("\tfeature\t", i, "/", nrow(list_images), "\t", list_images$image_name[i])
    }
}


# 8. Cluster
folder_script <- "/Users/chang-yu/Desktop/Lab/emergent-coexistence/output/check/"
batch <- "D"
list_images <- paste0(folder_script, "00-list_images-", batch, ".csv") %>% read_csv(show_col_types = F)
list_image_mapping <- paste0(folder_script, "00-list_image_mapping-", batch, ".csv") %>% read_csv(show_col_types = F)

list_image_mapping_folder <- list_image_mapping %>%
    left_join(select(list_images, image_name_pair = image_name, folder_feature_pair = folder_feature, folder_cluster), by = "image_name_pair") %>%
    left_join(select(list_images, image_name_isolate1 = image_name, folder_feature_isolate1 = folder_feature), by = "image_name_isolate1") %>%
    left_join(select(list_images, image_name_isolate2 = image_name, folder_feature_isolate2 = folder_feature), by = "image_name_isolate2")

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
ggsave(filename = paste0(folder_cluster, image_name, ".png"), plot = p, width = 10, height = 6)
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
















