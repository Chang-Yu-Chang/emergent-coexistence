library(tidyverse)
library(cowplot)
library(tidypredict)
library(gridExtra)

folder_script <- "/Users/chang-yu/Desktop/Lab/emergent-coexistence/output/check/"
list_images <- read_csv(commandArgs(trailingOnly = T)[1])
list_image_mapping <- read_csv(commandArgs(trailingOnly = T)[2])
# batch <- "D"
# list_images <- paste0(folder_script, "00-list_images-", batch, ".csv") %>% read_csv(show_col_types = F)
# list_image_mapping <- paste0(folder_script, "00-list_image_mapping-", batch, ".csv") %>% read_csv(show_col_types = F)

#
list_image_mapping_folder <- list_image_mapping %>%
    left_join(select(list_images, image_name_pair = image_name, folder_feature_pair = folder_blue_feature, folder_blue_cluster), by = "image_name_pair") %>%
    left_join(select(list_images, image_name_isolate1 = image_name, folder_feature_isolate1 = folder_blue_feature), by = "image_name_isolate1") %>%
    left_join(select(list_images, image_name_isolate2 = image_name, folder_feature_isolate2 = folder_blue_feature), by = "image_name_isolate2")


# which(list_image_mapping_folder$image_name_pair == "D_T8_C1R2_50-50_1_3")
i = 1
## These functions have the loop index i inside
read_feature_combined <- function () {
    object_feature_pair <- paste0(
        list_image_mapping_folder$folder_feature_pair[i],
        list_image_mapping_folder$image_name_pair[i], ".csv"
    ) %>%
        read_csv(show_col_types = FALSE) %>%
        mutate(image_name = list_image_mapping_folder$image_name_pair[i],
               Group = "pair")

    object_feature_isolate1 <- paste0(
        list_image_mapping_folder$folder_feature_isolate1[i],
        list_image_mapping_folder$image_name_isolate1[i], ".csv"
    ) %>%
        read_csv(show_col_types = FALSE) %>%
        mutate(image_name = list_image_mapping_folder$image_name_isolate1[i],
               Group = "isolate1")

    object_feature_isolate2 <- paste0(
        list_image_mapping_folder$folder_feature_isolate2[i],
        list_image_mapping_folder$image_name_isolate2[i], ".csv"
    ) %>%
        read_csv(show_col_types = FALSE) %>%
        mutate(image_name = list_image_mapping_folder$image_name_isolate2[i],
               Group = "isolate2")

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
        geom_point(aes_string(x = feature1, y = feature2, color = "image_name"), size = 2, shape = 1, stroke = .8, alpha = .7) +
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

for (i in 1:nrow(list_image_mapping_folder)) {
    image_name <- list_image_mapping_folder$image_name_pair[i]
    folder_blue_cluster <- list_image_mapping_folder$folder_blue_cluster[i]

    object_feature_combined <- read_feature_combined()
    color_names <- set_color_names()

    # Plot the three main features
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
        group_by(PredictedGroup, .drop = F) %>%
        count(name = "Count") %>%
        ungroup() %>%
        arrange(PredictedGroup) %>%
        mutate(Group = c("isolate1", "isolate2"), Align = c("right", "left"))

    # Model predict
    p4 <- object_feature_predicted %>%
        ggplot() +
        geom_histogram(aes(x = PredictedGroupProbability), color = 1, fill = NA, bins = 30) +
        geom_vline(xintercept = 0.5, linetype = 2, color = "red") +
        geom_text(data = object_feature_predicted_count, aes(x = 0.5, hjust = Align, label = paste0("  ", Group, ": ", Count, "  ")), y = Inf, vjust = 2) +
        scale_x_continuous(limits = c(-.1,1.1), breaks = seq(0, 1, by = .1)) +
        theme_classic()

    # Model parameter table
    p5 <- ggplot() +
        geom_point() +
        annotation_custom(tableGrob(broom::tidy(model), theme = ttheme_default(base_size = 8)), xmin = .2, xmax = .8, ymin = .2, ymax = .8) +
        scale_x_continuous(limits = c(0,1)) +
        scale_y_continuous(limits = c(0,1)) +
        theme_void()

    # Plot
    legend <- get_legend(p1 + theme(legend.box.margin = margin(0, 0, 0, 12)))
    p <- plot_grid(p1 + theme(legend.position = "none"),
                   p2 + theme(legend.position = "none"),
                   p3 + theme(legend.position = "none"),
                   p4 + theme(legend.position = "none"),
                   p5 + theme(legend.position = "none"),
                   legend, nrow = 2, align = "hv", axis = "lrtb") +
        theme(plot.background = element_rect(fill = "white", color = NA))
    ggsave(filename = paste0(folder_blue_cluster, image_name, ".png"), plot = p, width = 10, height = 6)
    cat("\nplot feature\t", i, "/", nrow(list_image_mapping_folder), "\t", image_name)


}


