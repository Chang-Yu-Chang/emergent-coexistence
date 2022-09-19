library(tidyverse)
library(cowplot)
library(leaps) # for computing stepwise regression. But it only fits lm
library(glmulti) # extension to include glm in leaps
library(gridExtra) # for making the table a grob


folder_script <- "~/Desktop/Lab/emergent-coexistence/output/check/"
list_images <- read_csv(commandArgs(trailingOnly = T)[1], show_col_types = F)
list_image_mapping <- read_csv(commandArgs(trailingOnly = T)[2], show_col_types = F)

# list_images <- read_csv(paste0(folder_script, "00-list_images-D.csv"), show_col_types = F)
# list_image_mapping <- read_csv(paste0(folder_script, "00-list_image_mapping-D.csv") , show_col_types = F)
# list_images <- read_csv(paste0(folder_script, "00-list_images-C2.csv"), show_col_types = F)
# list_image_mapping <- read_csv(paste0(folder_script, "00-list_image_mapping-C2.csv") , show_col_types = F)
# list_images <- read_csv(paste0(folder_script, "00-list_images-B2.csv"), show_col_types = F)
# list_image_mapping <- read_csv(paste0(folder_script, "00-list_image_mapping-B2.csv") , show_col_types = F)


list_image_mapping_folder <- list_image_mapping %>%
    left_join(select(list_images, image_name_pair = image_name, folder_feature_pair = folder_green_feature, folder_green_cluster), by = "image_name_pair") %>%
    left_join(select(list_images, image_name_isolate1 = image_name, folder_feature_isolate1 = folder_green_feature), by = "image_name_isolate1") %>%
    left_join(select(list_images, image_name_isolate2 = image_name, folder_feature_isolate2 = folder_green_feature), by = "image_name_isolate2")


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

    # Remove dark objects that are too dark to be colonies
    object_feature_combined <- object_feature_combined %>% filter(b.mean > 0)

    return(object_feature_combined)
}


plates_no_colony <- c(
    "B2_T8_C11R1_5-95_2_8",
    "B2_T8_C11R1_5-95_2_9",
    "B2_T8_C11R1_5-95_8_2",
    "B2_T8_C11R1_5-95_9_8",
    "B2_T8_C11R1_50-50_2_8",
    "B2_T8_C11R1_50-50_2_9",
    "C2_T8_C11R2_50-50_2_10",
    "C2_T8_C11R2_50-50_9_13"
)

i = which(list_image_mapping_folder$image_name_pair == "B2_T8_C11R1_5-95_1_3")
#for (i in temp_indices) {

for (i in 1:nrow(list_image_mapping_folder)) {
    #if (i < 6)  next
    # 8.1 read the feature files ----

    list_image_mapping_folder$image_name_pair[i]

    ## Skip images with no colony
    if (list_image_mapping_folder$image_name_pair[i] %in% plates_no_colony) {cat("\nno colony, no watershed image\t", list_image_mapping_folder$image_name_pair[i]); next}

    ## Read object features from isolate and co-coculture images
    feature_candidates <- c("b.mean", "b.sd", "b.mad",
                            "b.q005", "b.q05", "b.q095",
                            "s.area", "s.radius.mean", "s.radius.sd",
                            "b.tran.mean", "b.tran.sd", "b.tran.mad",
                            "b.center", "b.periphery", "b.diff.cp",
                            "b.tran.q005", "b.tran.q01", "b.tran.q05", "b.tran.q09", "b.tran.q095"
                            # "t.bump.number"
    )

    object_feature_combined <- read_feature_combined() %>%
        mutate(GroupBinary = case_when(
            Group == "isolate1" ~ 0,
            Group == "isolate2" ~ 1,
        )) %>%
        # Start with these parameters
        select(image_name, ObjectID, Group, GroupBinary, all_of(feature_candidates)) %>%
        {.}

    ## Remove outlier objects according to the criterion for each pair/isolate iamage. See this R script for more detail
    source("05a-remove_outlier_objects.R")

    object_feature_isolates <- object_feature_combined %>%
        filter(!is.na(GroupBinary))
    object_feature_pair <- object_feature_combined %>%
        filter(is.na(GroupBinary))
    cat("\nread feature")


    # 8.2 exhaustive feature selection  ----
    n_models <- glmulti(
        y = "GroupBinary",
        xr = feature_candidates,
        data = object_feature_isolates,
        level = 1,
        method = "d",
        crit = "aicc",
        minsize = 2,
        maxsize = 4,
        confsetsize = 100, # Keep the top n models
        plotty = F, report = F,
        fitfunction = glm
        #family = binomial # logistic
    ) %>% capture.output()

    cat("\tnumber of models:", str_replace(n_models[13], "\\[1\\] ", ""))

    best_subset <- glmulti(
        y = "GroupBinary",
        xr = feature_candidates,
        data = object_feature_isolates,
        level = 1,
        method = "h",
        crit = "aicc",
        minsize = 2,
        maxsize = 4,
        confsetsize = 100, # Keep the top n models
        plotty = F, report = F,
        fitfunction = glm# logistic
    )

    # Multimodel inference, using the metafor package
    #' The importance value for a particular predictor is equal to the sum of the weights/probabilities for the models in which the variable appears.
    #' these values can be regarded as the overall support for each variable across all models in the candidate set

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
        geom_hline(yintercept = 0.5, color = "red") +
        coord_flip() +
        theme_classic() +
        ggtitle("model-averaged importance of predictors")
    cat("\timportance")

    # 8.3 Use the selected important features to build a logit model ----
    selected_features <- model_importance %>% filter(Importance > 0.5) %>% pull(Feature)
    ## If only 0, 1, or 2 features pass the importance threshold 0.5, choose the top three features regardless of the importance
    if (length(selected_features) <= 2) selected_features <- model_importance %>% slice(1:3) %>% pull(Feature)

    ## Build model. Include all interaction terms between the selected features
    model_logit <- glm(GroupBinary ~ .^2,
                       data = select(object_feature_isolates, GroupBinary, all_of(selected_features)),
                       family = binomial, control = list(maxit = 50))

    ## Remove the non-significant interaction terms
    drop_bad_interactions <- function (model_logit) {
        bad_interactions <- broom::tidy(model_logit) %>%
            filter(p.value > 0.05 & str_detect(term, ":")) %>%
            pull(term)
        if (length(bad_interactions) > 0) for (k in 1:length(bad_interactions)) model_logit <- update(model_logit, as.formula(paste(".~.-", bad_interactions[k])))
        return(model_logit)
    }
    model_logit <- drop_bad_interactions(model_logit)
    model_table <- model_logit %>%
        broom::tidy() %>%
        mutate(estimate = round(estimate, 4),
               std.error = round(std.error, 4),
               statistic = round(statistic, 4))
    ## Plot the model in a table
    p2 <- cowplot::ggdraw() +
        annotation_custom(tableGrob(model_table, theme = ttheme_default(base_size = 8)), xmin = 0, xmax = 1, ymin = 0, ymax = 1) +
        theme_void() +
        theme(legend.position = "none")
    cat("\tcoefficient")


    # 8.4 Model prediction ----
    object_feature_predicted <- object_feature_pair %>%
        mutate(PredictedGroupProbability = predict(model_logit, object_feature_pair, type = "response")) %>%
        dplyr::select(image_name, ObjectID, PredictedGroupProbability) %>%
        # Criteria for catagorizing
        mutate(Group = case_when(
            PredictedGroupProbability < 0.2 ~ "predicted isolate1",
            PredictedGroupProbability > 0.8 ~ "predicted isolate2",
            PredictedGroupProbability > 0.2 & PredictedGroupProbability < 0.8 ~ "undecided"
        )) %>%
        mutate(Group = factor(Group, c("predicted isolate1", "predicted isolate2", "undecided"),))

    object_feature_predicted_count <- object_feature_predicted %>%
        group_by(Group, .drop = F) %>%
        count(name = "Count") %>%
        ungroup() %>%
        mutate(Halign = c("right", "left", "center"), Valign = c(2, 2, 4)) %>%
        mutate(Group = as.character(Group) %>% str_replace("predicted isolate", "pd_iso"))

    p3 <- object_feature_predicted %>%
        ggplot() +
        geom_histogram(aes(x = PredictedGroupProbability), binwidth = .02, color = 1, fill = "white", breaks = seq(0,1,.02)) +
        # geom_histogram(aes(x = PredictedGroupProbability), binwidth = 0.05, color = 1, fill = NA, bins = 30) +
        geom_vline(xintercept = c(0.2, 0.8), linetype = 2, color = "red") +
        geom_vline(xintercept = c(0,1), linetype = 2, color = "black") +
        geom_text(data = object_feature_predicted_count, x = c(.45,.55,.5), aes(hjust = Halign, vjust = Valign, label = paste0("  ", Group, ": ", Count, "  ")), y = Inf) +
        scale_x_continuous(expand = c(0,0.1), breaks = seq(0,1,.2)) +
        scale_y_continuous(expand = c(0,0)) +
        theme_classic()
    write_csv(object_feature_predicted, paste0(list_image_mapping_folder$folder_green_cluster[i], list_image_mapping_folder$image_name_pair[i], ".csv"))
    cat("\tprediction")

    # 8.5 Scatterplot for clustering intuition ----
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
    color_names <- set_color_names()
    fill_names <- set_fill_names()
    shape_names <- set_shape_names()
    alpha_names <- set_alpha_names()

    ## Combined the prediction and object label
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

    # 8.6 PCA with the model selected variables ----
    if (length(coef(model_logit)) >= 3) {
        # When more than two features picked
        features_to_plot <- names(coef(model_logit))[-1] %>%
            str_subset("^((?!:).)*$") # Remove interaction terms
        pcobj <- object_feature_plot %>%
            select(all_of(features_to_plot)) %>%
            prcomp(center = TRUE, scale. = TRUE)
    } else if (length(coef(model_logit)) == 2) {
        # When only one feature picked, use the top two important features
        pcobj <- object_feature_plot %>%
            select(all_of(model_importance$Feature[1:2])) %>%
            prcomp(center = TRUE, scale. = TRUE)
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

    # 8.7 Combine all plots
    legend <- get_legend(p4 + theme(legend.box.margin = margin(0, 0, 0, 12)))
    p <- plot_grid(
        p1, p2, p3,
        p4 + theme(legend.position = "none"),
        p5 + theme(legend.position = "none"),
        legend,
        nrow = 2, rel_widths = c(1,1,1), scale = .9,
        align = "hv", axis = " tblr"
    ) + theme(plot.background = element_rect(fill = "white", color = NA))

    ggsave(filename = paste0(list_image_mapping_folder$folder_green_cluster[i], list_image_mapping_folder$image_name_pair[i], ".png"), plot = p, width = 15, height = 8)
    cat("\t", i, "/", nrow(list_image_mapping_folder), "\t", list_image_mapping_folder$image_name_pair[i])

}
