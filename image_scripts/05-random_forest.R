#' This scripts run the random forest
#'
#' To use this Rscript, in bash environment:
#' Rscript 05-random_forest.R list_images.csv list_image_mapping.csv
#'
#' For example:
#' Rscript 05-random_forest.R mapping_files/00-list_images-B2-green.csv mapping_files/00-list_image_mapping-B2.csv

library(tidyverse)
library(cowplot)
library(caret) # for streamlining the model training process for complex classification and regression
library(randomForest) # For implementing random forest algorithm
source(here::here("processing_scripts/00-metadata.R"))

list_images <- read_csv(commandArgs(trailingOnly = T)[1], show_col_types = F)
list_image_mapping <- read_csv(commandArgs(trailingOnly = T)[2], show_col_types = F)

list_image_mapping_folder <- list_image_mapping %>%
    left_join(rename(list_images, image_name_pair = image_name), by = "image_name_pair") %>%
    left_join(select(list_images, image_name_isolate1 = image_name), by = "image_name_isolate1") %>%
    left_join(select(list_images, image_name_isolate2 = image_name), by = "image_name_isolate2")

#
read_feature <- function (type, color_channel) {
    stopifnot(type %in% c("pair", "isolate"))

    if (type == "pair") {
        object_feature_pair <- paste0(
            paste0(list_images[i,paste0("folder_feature")], color_channel, "/"),
            list_image_mapping_folder$image_name_pair[i], ".csv"
        ) %>%
            read_csv(show_col_types = FALSE) %>%
            mutate(image_name = list_image_mapping_folder$image_name_pair[i],
                   Group = NA) # Label for supervised learning
        object_feature <- object_feature_pair %>%
            select(image_name, ObjectID, Group, everything())
    }

    if (type == "isolate") {
        object_feature_isolate1 <- paste0(
            paste0(list_images[i,paste0("folder_feature")], color_channel, "/"),
            list_image_mapping_folder$image_name_isolate1[i], ".csv"
        ) %>%
            read_csv(show_col_types = FALSE) %>%
            mutate(image_name = list_image_mapping_folder$image_name_isolate1[i],
                   Group = "isolate1") # Label for supervised learning

        object_feature_isolate2 <- paste0(
            paste0(list_images[i,paste0("folder_feature")], color_channel, "/"),
            list_image_mapping_folder$image_name_isolate2[i], ".csv"
        ) %>%
            read_csv(show_col_types = FALSE) %>%
            mutate(image_name = list_image_mapping_folder$image_name_isolate2[i],
                   Group = "isolate2") # Label for supervised learning

        object_feature <- bind_rows(object_feature_isolate1, object_feature_isolate2) %>%
            select(image_name, ObjectID, Group, everything())
    }

    return(object_feature)
}
create_customRF <- function (){
    #' caret package does not take ntree as default for method="rf"
    #' This function is to customize the rf approach to include argument ntree
    customRF <- list(type = "Classification", library = "randomForest", loop = NULL)
    customRF$parameters <- data.frame(parameter = c("mtry", "ntree"), class = rep("numeric", 2), label = c("mtry", "ntree"))
    customRF$grid <- function(x, y, len = NULL, search = "grid") {}
    customRF$fit <- function(x, y, wts, param, lev, last, weights, classProbs, ...) {
        randomForest(x, y, mtry = param$mtry, ntree=param$ntree, ...)
    }
    customRF$predict <- function(modelFit, newdata, preProc = NULL, submodels = NULL)
        predict(modelFit, newdata)
    customRF$prob <- function(modelFit, newdata, preProc = NULL, submodels = NULL)
        predict(modelFit, newdata, type = "prob")
    customRF$sort <- function(x) x[order(x[,1]),]
    customRF$levels <- function(x) x$classes
    return(customRF)
}
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
    c(NA,NA,NA,NA, grey(.8)) %>%
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
    # rep(1, 5) %>% # Hallow point
    c(1,1,1,1,21) %>%
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

for (i in 1:nrow(list_image_mapping_folder)) {
    cat("\t", i)
    image_name <- list_image_mapping_folder$image_name_pair[i]

    ## Skip images with no colony
    if (list_image_mapping_folder$image_name_pair[i] %in% plates_no_colony) {cat("\nno colony, no watershed image\t", list_image_mapping_folder$image_name_pair[i]); next}


    # 9.1 Read object features from isolate and co-coculture images ----
    object_feature_isolates <- read_feature(type = "isolate", color_channel = "merged") %>%
        select(image_name, ObjectID, Group, all_of(feature_candidates)) %>%
        mutate(Group = factor(Group))
    object_feature_pair <- read_feature(type = "pair", color_channel = "merged") %>%
        select(image_name, ObjectID, Group, all_of(feature_candidates)) %>%
        filter(is.na(Group))
    cat("\nfeature")

    # 9.2 Train random forest model using repeated CV ----
    ## Repeated 10-fold cross-validation
    ctrl <- trainControl(
        method = "repeatedcv",
        repeats = 3,
        number = 10
    )
    ## Test an array of arguments for random forest
    mtry_range <- seq(2, (length(object_feature_isolates)-3), 5)
    rfGrid <- expand.grid(.mtry = mtry_range, .ntree = c(500))
    ## Customize RF to include argument ntree
    customRF <- create_customRF()

    ## Training
    set.seed(1)
    model <- train(
        Group ~ .,
        data = select(object_feature_isolates, -image_name, -ObjectID),
        method = customRF,
        metric = "Accuracy",
        tuneGrid = rfGrid,
        trControl = ctrl
    )
    cat("\tvalidation")
    # 9.3 Build the final model using the set of arguments (mtry) with highest validation accuracy ----
    ## Final model
    set.seed(1)
    model_final <- train(
        Group ~. ,
        data = select(object_feature_isolates, -image_name, -ObjectID),
        method = customRF,
        metric = "Accuracy",
        tuneGrid = expand.grid(.mtry = model$bestTune$mtry, .ntree = model$bestTune$ntree),
        trControl = ctrl
    )


    ## Plot the accuracy in validation and training sets
    p1 <- model$results %>%
        as_tibble() %>%
        mutate(ntree = factor(ntree)) %>%
        ggplot() +
        # Validation model accuracy
        geom_line(aes(x = mtry, y = Accuracy, group = ntree, color = "validation set"), size = 1) +
        geom_point(aes(x = mtry, y = Accuracy, group = ntree), size = 3, shape = 21, fill = "white", stroke = 2) +
        # Train model accuracy
        geom_point(aes(x = model_final$results$mtry, y = model_final$results$Accuracy, color = "training set"), size = 3, shape = 21, fill = "white", stroke = 2) +
        geom_hline(aes(yintercept = model_final$results$Accuracy, color = "training set"), size = 1) +
        scale_x_continuous(breaks = mtry_range) +
        scale_color_manual(values = c("validation set" = "black", "training set" = "blue")) +
        theme_cowplot() +
        theme(panel.grid.major = element_line(color = gray(0.8))) +
        labs(color = "")
    cat("\ttraining")

    # 9.4 Feature importance ----
    model_importance <- model_final$finalModel$importance %>%
        as_tibble(rownames = "Feature") %>%
        arrange(desc(MeanDecreaseGini))
    p2 <- model_importance %>%
        mutate(Feature = factor(Feature, rev(Feature))) %>%
        ggplot() +
        geom_col(aes(x = Feature, y = MeanDecreaseGini), color = 1, fill = NA) +
        coord_flip() +
        theme_classic() +
        theme(axis.text.y = element_text(size = 8))
    cat("\timportance")

    # 9.5 Prediction on the coculture images ----
    ## Build the model using all isolate data, using the selected arguments mtry, ntree
    object_feature_predicted <- object_feature_pair %>%
        # Predicted group
        mutate(Group = paste0("predicted ", predict(model_final, newdata = object_feature_pair))) %>%
        select(image_name, ObjectID, Group) %>%
        # Predicted probability of each class
        bind_cols(as_tibble(predict(model_final, newdata = object_feature_pair, type = "prob")) %>% rename(PredictedProbabilityIsolate1 = isolate1, PredictedProbabilityIsolate2 = isolate2)) %>%
        mutate(Group = factor(Group, c("predicted isolate1", "predicted isolate2")))

    object_feature_predicted_count <- object_feature_predicted %>%
        group_by(Group, .drop = F) %>%
        count(name = "Count") %>%
        ungroup() %>%
        mutate(Halign = c("right", "left"), Valign = c(2, 2)) %>%
        mutate(Group = as.character(Group) %>% str_replace("predicted isolate", "pd_iso"))


    ## Compute class probabilities of each object
    p3a <- predict(model, newdata = object_feature_pair, type = "prob") %>%
        as_tibble() %>%
        arrange(isolate2) %>%
        mutate(ObjectID = 1:n()) %>%
        pivot_longer(-ObjectID, names_to = "Class", values_to = "Probability") %>%
        mutate(Class = factor(Class, c("isolate1", "isolate2"))) %>%
        ggplot() +
        geom_col(aes(x = ObjectID, y = Probability, fill = Class), color = 1, size = .05) +
        geom_hline(yintercept = 0.5, color = "red", linetype = 2) +
        coord_flip() +
        scale_fill_manual(values = c("isolate1" = gray(0.3), "isolate2" = gray(0.9))) +
        scale_y_continuous(breaks = seq(0,1,.2), position = "right") +
        theme_classic() +
        theme(legend.title = element_blank())


    p3b <- predict(model, newdata = object_feature_pair, type = "prob") %>%
        as_tibble() %>%
        ggplot() +
        geom_histogram(aes(x = isolate2), color = 1, fill = NA, binwidth = 0.02, breaks = seq(0,1,.02)) +
        geom_vline(xintercept = c(.5), linetype = 2, color = "red") +
        geom_vline(xintercept = c(0,1), linetype = 2, color = "black") +
        geom_text(data = object_feature_predicted_count, x = c(.45,.55), aes(hjust = Halign, vjust = Valign, label = paste0("  ", Group, ": ", Count, "  ")), y = Inf) +
        scale_x_continuous(breaks = seq(0,1,.2), limits = c(0,1)) +
        scale_y_continuous(expand = c(0,0)) +
        theme_classic() +
        labs(x = "Probability of being isolate2")
    cat("\tprediction")

    # 9.6 Scatterplot for clustering intuition ----
    color_names <- set_color_names()
    fill_names <- set_fill_names()
    shape_names <- set_shape_names()
    alpha_names <- set_alpha_names()

    ## Combined the prediction and object label
    object_feature_plot <- object_feature_predicted %>%
        select(image_name, ObjectID, Group) %>%
        left_join(select(object_feature_pair, -Group), by = c("image_name", "ObjectID")) %>%
        bind_rows(object_feature_isolates) %>%
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
    cat("\tscatterplot")

    # 9.7 PCA with the model selected variables ----
    pcobj <- object_feature_plot %>%
        select(-image_name, -ObjectID, -Group, -ends_with("Label")) %>%
        prcomp(center = TRUE, scale. = TRUE)
    pca_coord <- compute_pca_coord(pcobj)

    p5 <- pca_coord$df.u %>%
        bind_cols(object_feature_plot) %>%
        arrange(ColorLabel) %>%
        ggplot() +
        # Draw scores
        geom_point(aes(x = xvar, y = yvar,
                       color = ColorLabel, fill = FillLabel, shape = ShapeLabel, alpha = AlphaLabel),
                   size = 2, stroke = .8) +
        scale_color_manual(values = color_names, name = "", label = names(color_names)) +
        scale_fill_manual(values = fill_names, name = "", label = names(color_names)) +
        scale_shape_manual(values = shape_names, name = "", label = names(color_names)) +
        scale_alpha_manual(values = alpha_names, name = "", label = names(color_names)) +
        theme_classic() +
        labs(x = pca_coord$u.axis.labs[1], y = pca_coord$u.axis.labs[2])
    cat("\tpcaplot")

    # 9.8 Combine all plots ----
    legend <- get_legend(p4 + theme(legend.box.margin = ggplot2::margin(0, 0, 0, 12), legend.direction = "horizontal"))
    p_left <- plot_grid(
        p1 + theme(legend.position = "top"),
        p2,
        align = "v", axis = "lr",
        nrow = 2, rel_heights = c(1, 1.5), scale = 0.9, labels = LETTERS[1:2]
    )

    p_middle <- plot_grid(
        p3a + theme(legend.position = "top"),
        p3b,
        align = "v", axis = "lr",
        nrow = 2, rel_heights = c(1.5, 1), scale = 0.9, labels = LETTERS[3:4]
    )

    p_right <- plot_grid(
        p4 + theme(legend.position = "none"),
        p5 + theme(legend.position = "none"),
        align = "v", axis = "lr",
        nrow = 2, rel_heights = c(1, 1), scale = 0.9, labels = LETTERS[5:6]
    )


    p <- plot_grid(
        p_left,
        p_middle,
        p_right,
        NULL,
        legend,
        NULL,
        nrow = 2, scale = 1, rel_heights = c(1,.1)
    ) + theme(plot.background = element_rect(fill = "white", color = NA))
    ggsave(filename = paste0(list_image_mapping_folder$folder_random_forest[i], list_image_mapping_folder$image_name_pair[i], ".png"), plot = p, width = 15, height = 8)

    # 9.9 Output csv ----
    ## Model selection. The data for making panel A
    bind_rows(
        mutate(model$results, FinalModel = F),
        mutate(model_final$results, FinalModel = T)
    ) %>%
        write_csv(paste0(list_image_mapping_folder$folder_random_forest[i], "validation-", list_image_mapping_folder$image_name_pair[i], ".csv"))

    ## Object prediction. The data for making panel C and D
    object_feature_predicted %>%
        write_csv(paste0(list_image_mapping_folder$folder_random_forest[i], list_image_mapping_folder$image_name_pair[i], ".csv"))

    cat("\t", i, "/", nrow(list_image_mapping_folder), "\t", list_image_mapping_folder$image_name_pair[i])

}








