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
    save(image_watershed2, file = paste0(list_images$folder_green_watershed_file[i], image_name, ".RData")) # save watersed image object
    writeImage(colorLabels(image_watershed2), paste0(list_images$folder_green_watershed[i], image_name, ".tiff"))
    cat("\twatershed\t", i, "/", nrow(list_images), "\t", list_images$image_name[i])

}



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
# library(caret) # for easy machine learning workflow
# library(bestglm) # extension to include glm in leaps
# library(glmnet) # for fitting glm via penalized maximum likelihood
library(leaps) # for computing stepwise regression. But it only fits lm
library(glmulti) # extension to include glm in leaps
#library(tidyverse)
library(tidypredict)
library(gridExtra)


list_images <- read_csv("~/Desktop/Lab/emergent-coexistence/output/check/00-list_images-D.csv", show_col_types = F)
list_image_mapping <- read_csv(paste0(folder_script, "00-list_image_mapping-D.csv"), show_col_types = F)
list_image_mapping_folder <- list_image_mapping %>%
    left_join(select(list_images, image_name_pair = image_name, folder_feature_pair = folder_green_feature, folder_green_cluster), by = "image_name_pair") %>%
    left_join(select(list_images, image_name_isolate1 = image_name, folder_feature_isolate1 = folder_green_feature), by = "image_name_isolate1") %>%
    left_join(select(list_images, image_name_isolate2 = image_name, folder_feature_isolate2 = folder_green_feature), by = "image_name_isolate2")

## 8.1 read the feature files
i = 2
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
#feature_candidates_test <- c("b.mean", "b.sd", "b.mad", "s.area", "s.radius.mean")
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
names(object_feature_isolates)

# 8.2 exhaustive stepwise selection and cross-validation
best_subset <- glmulti(
    y = "GroupBinary",
    #xr = c("b.mean"),
    xr = feature_candidates,
    #xr = str_subset(names(object_feature_isolates), "^b."),
    # xr = names(object_feature_combined)[names(object_feature_isolates) != c("Group", "ObjectID")],
    data = object_feature_isolates,
    level = 1,
    method = "h",
    crit = "aicc",
    minsize = 2,
    maxsize = 5,
    confsetsize = 5, # Keep the top n models
    plotty = F, report = T,
    fitfunction = "glm"
    #family = binomial # logistic
)
# plot(best_subset)
# top <- weightable(best_subset)
# top <- top[top$aicc <= min(top$aicc) + 2,]
# broom::tidy(best_subset@objects[[1]])
"plot this. or should be able to turn this into ggplot here"
plot(best_subset, type = "s")
model <- best_subset@objects[[1]]
broom::tidy(model)
# Linear Models: the absolute value of the t-statistic for each model parameter is used.
model_importance <- broom::tidy(model) %>%
    arrange(desc(abs(statistic))) %>%
    filter(term != "(Intercept)") %>%
    pull(term)

# object_feature_isolates %>%
#     ggplot() +
#     geom_point(aes(x = b.mean, y = b.sd, color = Group)) +
#     theme_classic()
object_feature_isolates %>%
    ggplot() +
    geom_point(aes_string(x = model_importance[1], y = model_importance[2], color = "Group")) +
    theme_classic()









if (FALSE) {

    predict_regsubsets <- function(model_object, newdata, id ,...) {
        form <- as.formula(model_object$call[[2]])
        mat <- model.matrix(form, newdata)
        coefi <- coef(model_object, id = id)
        xvars <- names(coefi)
        mat[, xvars] %*% coefi
    }
    n_feature_max = 10
    #best_subset <- regsubsets(GroupBinary ~ ., object_feature_isolates, nvmax = n_feature_max, method = "exhaustive")
    #results <- broom::tidy(best_subset)

    ## k-fold cross-validation
    number_k <- 10
    set.seed(1)
    folds <- sample(1:number_k, nrow(object_feature_isolates), replace = TRUE) # Not all training set will be of the size but roughly
    cv_errors <- matrix(NA, n_feature_max, number_k, dimnames = list(paste(1:n_feature_max), NULL)) # output cross-validation matrix
    for (k in 1:number_k) {
        # Perform best subset on row that are NOT in the test test k
        best_subset <- regsubsets(GroupBinary ~ . - image_name - ObjectID - Group, object_feature_isolates[folds != k,], nvmax = n_feature_max, method = "exhaustive")
        # Perform cross-validation
        for (j in 1:n_feature_max) {
            pred_x <- predict_regsubsets(best_subset, object_feature_isolates[folds == k,], id = j)
            cv_errors[j,k] <- mean((object_feature_isolates$GroupBinary[folds == k] - pred_x)^2)
        }
    }

    ## Find the best model
    best_model_index <- which.min(rowMeans(cv_errors))
    p1 <- tibble(NumberOfPredictors = 1:n_feature_max, CrossValidationError = rowMeans(cv_errors)) %>%
        ggplot(aes(x = NumberOfPredictors, y = CrossValidationError)) +
        geom_line() +
        geom_point(shape = 1, size = 3, stroke = 1) +
        geom_vline(xintercept = best_model_index, color = "red", linetype = 2) +
        scale_x_continuous(breaks = 1:n_feature_max) +
        theme_classic() +
        ggtitle("cross-validation for best model")

    ## Coefficients of the best model
    final_best_subset <- regsubsets(GroupBinary ~ . - image_name - ObjectID - Group, object_feature_isolates, nvmax = n_feature_max, method = "exhaustive")
    #plot(final_best_subset, scale = "adjr2") is the default alternative for the ggplot plot below
    p2 <- final_best_subset %>%
        broom::tidy() %>%
        mutate(NumberOfPredictor = 1:n()) %>%
        pivot_longer(cols = -c(NumberOfPredictor, r.squared, adj.r.squared, BIC, mallows_cp), names_to = "Feature", values_to = "Value") %>%
        mutate(adj.r.squared.discrete = factor(round(adj.r.squared,4))) %>%
        mutate(Feature = factor(Feature, c("(Intercept)", feature_candidates))) %>%
        ggplot() +
        geom_tile(aes(x = NumberOfPredictor, y = Feature, fill = Value, alpha = adj.r.squared)) +
        geom_vline(xintercept = c(1:n_feature_max)+0.5, color = grey(.8)) +
        geom_hline(yintercept = c(1:length(feature_candidates))+0.5, color = grey(.8)) +
        scale_fill_manual(values = c("TRUE" = "black", "FALSE" = "white")) +
        scale_x_continuous(expand = c(0,0), breaks = 1:n_feature_max) +
        scale_y_discrete(expand = c(0,0)) +
        theme_classic() +
        theme(panel.border = element_rect(fill = NA, color = 1)) +
        guides(fill = "none") +
        ggtitle("each column is a n-feature model")

    ## Coefficient importance by the abs(coef)
    best_model_coef <- coef(final_best_subset, best_model_index)
    best_model_coef_importance <- best_model_coef %>%
        enframe(name = "Feature", value = "Value") %>%
        mutate(Sign = ifelse(Value > 0, "positive", "negative"),
               Rank = rank(abs(Value))) %>%
        arrange(desc(Rank)) %>%
        mutate(Feature = ordered(Feature, rev(Feature)))

    p3 <- best_model_coef_importance %>%
        ggplot() +
        geom_col(aes(x = Feature, y = abs(Value), fill = Sign), color = 1) +
        #scale_x_reverse() +
        scale_fill_manual(values = c("positive" = "white", "negative" = "black")) +
        coord_flip() +
        theme_minimal() +
        theme(legend.position = c(.8, .3)) +
        labs(y = "coefficient")
    p3
}
if (FALSE) {
"there should a way to extract to extract the confidence interaval for each selected variables
now I can only extract coeff"

#best_model_coef <- coef(best_subset, best_model_index)
#names(best_model_coef)
predict_model_coef <- function (tibble_test, model_coef, ...) {
    #model_coef <- best_model_coef
    #    stopifnot(names(model_coef)[1] == "(Intercept)")
    #    tibble_test <- object_feature_pair
    coef_names <- names(model_coef)[-1]
    prediction <- rep(model_coef["(Intercept)"], nrow(tibble_test)) %>% unname
    for (i in 1:length(coef_names)) {
        prediction <- prediction + unname(unlist(tibble_test[, "b.mean"])) * model_coef[coef_names[i]]
    }

    return(prediction)
}

}

library(cowplot)
plot_grid(p1, p2, ncol = 1, align = "v", axis = "lr")
plot_grid(p1, p3, p2, NULL, ncol = 2, align = "v", axis = "lr")

# Model predict
object_feature_predicted <- object_feature_pair %>%
    mutate(PredictedGroupProbability = predict(model, object_feature_pair)) %>%
    dplyr::select(image_name, ObjectID, PredictedGroupProbability) %>%
    mutate(PredictedGroup = ifelse(PredictedGroupProbability < 0.5, 0, 1) %>% factor(c(0,1)))
#mutate(PredictedGroup = ifelse(PredictedGroupProbability < 0.5, 0, 1) %>% factor(c(0,1)))

object_feature_predicted_count <- object_feature_predicted %>%
    mutate(PredictedGroup = ifelse(PredictedGroupProbability < 0.5, 0, 1) %>% factor(c(0,1))) %>%
    group_by(PredictedGroup, .drop = F) %>%
    count(name = "Count") %>%
    ungroup() %>%
    mutate(Group = c("isolate1", "isolate2"), Align = c("right", "left"))

p4 <- object_feature_predicted %>%
    ggplot() +
    geom_histogram(aes(x = PredictedGroupProbability), color = 1, fill = NA, bins = 30) +
    geom_vline(xintercept = 0.5, linetype = 2, color = "red") +
    geom_text(data = object_feature_predicted_count, aes(x = 0.5, hjust = Align, label = paste0("  ", Group, ": ", Count, "  ")), y = Inf, vjust = 2) +
    scale_x_continuous(expand = c(0,1)) +
    scale_y_continuous(expand = c(0,0)) +
    theme_classic()

p4

# Scatterplot for clustering intuition
set_color_names <- function() {
    c(2,3,2,3) %>% setNames(c(
        paste0(list_image_mapping_folder$image_name_isolate1[i], " isolate1"),
        paste0(list_image_mapping_folder$image_name_isolate2[i], " isolate2"),
        paste0(list_image_mapping_folder$image_name_pair[i], " predicted isolate1"),
        paste0(list_image_mapping_folder$image_name_pair[i], " predicted isolate2")
    ))
}
set_shape_names <- function() c(16, 21) %>% setNames(c("pair", "isolate"))
color_names <- set_color_names()
shape_names <- set_shape_names()

p5 <- object_feature_predicted %>%
    select(image_name, ObjectID, PredictedGroup) %>%
    mutate(Group = case_when(
        PredictedGroup == 0 ~ paste0("predicted isolate1"),
        PredictedGroup == 1 ~ paste0("predicted isolate2")
    )) %>%
    left_join(select(object_feature_pair, -GroupBinary, -Group), by = c("image_name", "ObjectID")) %>%
    bind_rows(object_feature_isolates) %>%
    mutate(ColorLabel = paste(image_name, Group),
           ShapeLabel = ifelse(str_detect(Group, "predicted"), "pair", "isolate")) %>%
    # select(ColorLabel, ShapeLabel, image_name, Group) %>%
    # distinct()
    ggplot() +
    geom_point(aes_string(x = model_importance[1],
                          y = model_importance[2],
                          color = "ColorLabel", shape = "ShapeLabel"),
               size = 1, stroke = .8) +
    scale_color_manual(values = color_names) +
    scale_shape_manual(values = shape_names) +
    theme_classic()


# PCA with the model selected variables
library(ggfortify)
object_feature_isolates %>%
    select(all_of(model_importance)) %>%
    scale() %>%
    prcomp() %>%
    autoplot(data = object_feature_isolates %>% select(Group, all_of(feature_candidates)), colour = "Group") +
    theme_classic()



# Figures
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






# # Logistic regression using the three main features
# training <- subset_training_set(object_feature_combined)
# model <- glm(dummy ~ b.mean + b.sd + b.mad, data = training, family = "binomial")

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
















