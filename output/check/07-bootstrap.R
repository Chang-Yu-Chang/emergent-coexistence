#' This script reads the object probability from the random forest and bootstraps
#' the colony frequencies
#' Rscript 07-bootstrap.R 00-list_images-D-green.csv 00-list_image_mapping-D.csv
library(tidyverse)
library(cowplot)

folder_script <- "~/Desktop/Lab/emergent-coexistence/output/check/"
list_images <- read_csv(commandArgs(trailingOnly = T)[1], show_col_types = F)
list_image_mapping <- read_csv(commandArgs(trailingOnly = T)[2], show_col_types = F)

# list_images <- read_csv(paste0(folder_script, "00-list_images-D-green.csv"), show_col_types = F)
# list_image_mapping <- read_csv(paste0(folder_script, "00-list_image_mapping-D.csv") , show_col_types = F)

list_image_mapping_folder <- list_image_mapping %>%
    left_join(rename(list_images, image_name_pair = image_name), by = "image_name_pair") %>%
    left_join(select(list_images, image_name_isolate1 = image_name), by = "image_name_isolate1") %>%
    left_join(select(list_images, image_name_isolate2 = image_name), by = "image_name_isolate2")

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


set.seed(1)
n_bootstraps = 1000 # Number of bootstraps

for (i in 1:nrow(list_image_mapping_folder)) {
    ## Skip images with no colony
    if (list_image_mapping_folder$image_name_pair[i] %in% plates_no_colony) {cat("\nno colony, no watershed image\t", list_image_mapping_folder$image_name_pair[i]); next}
    # 10.0 Read the random forest object probabilities ----
    image_name <- list_image_mapping_folder$image_name_pair[i]
    object_feature_predicted <- read_csv(paste0(list_image_mapping_folder$folder_random_forest[i], image_name, ".csv"), show_col_types = F)
    cat("\n", nrow(object_feature_predicted), "objects")


    # 10.1 bootstrapping ----
    cat("\t bootstrap", n_bootstraps, " times")
    object_bootstrapped <- rep(list(NA), n_bootstraps)
    for (j in 1:n_bootstraps) {
        object_bootstrapped[[j]] <- object_feature_predicted %>%
            rowwise() %>%
            mutate(BootstrapID = j) %>%
            mutate(Group = sample(c("isolate1", "isolate2"), size = 1, replace = F, prob = c(PredictedProbabilityIsolate1, PredictedProbabilityIsolate2))) %>%
            select(image_name, BootstrapID, ObjectID, Group)
        if (j %% 100 == 0) cat(" ", j)
    }

    # 10.2 output the result ----
    object_bootstrapped <- bind_rows(object_bootstrapped) %>%
        pivot_wider(id_cols = ObjectID, names_from = BootstrapID, values_from = Group, names_prefix = "bootstrap_")
    write_csv(object_bootstrapped, paste0(list_image_mapping_folder$folder_bootstrap[i], image_name, ".csv"))

    cat("\t", i, "/", nrow(list_image_mapping_folder), "\t", list_image_mapping_folder$image_name_pair[i])
}


if (FALSE) {
    object_bootstrapped_count <- bind_rows(object_bootstrapped) %>%
        mutate(Group = factor(Group, c("isolate1", "isolate2"))) %>%
        group_by(BootstrapID, Group, .drop = F) %>%
        count(name = "Count") %>%
        group_by(BootstrapID) %>%
        mutate(Fraction = Count / sum(Count))


    object_bootstrapped_count %>%
        pivot_wider(id_cols = BootstrapID, names_from = Group, values_from = Fraction) %>%
        ggplot() +
        geom_histogram(aes(x = isolate1), color = 1, fill = NA) +
        sca
    theme_classic()

}
