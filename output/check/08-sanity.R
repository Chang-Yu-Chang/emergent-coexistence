#' This script reads the output from 07-bootstrap and compare the result to human
#' count result
library(tidyverse)
library(cowplot)
folder_script <- "~/Desktop/Lab/emergent-coexistence/output/check/"
folder_main <- "~/Dropbox/lab/emergent-coexistence/data/raw/plate_scan/emergent_coexistence_plate_scan_check/"
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

# 1. Read accuracy data ----
#batch_names <- c("D", "C2", "B2")
batch_names <- c("D", "C2", "B2", "C")

## Pairs mapping file
list_image_mapping_master <- rep(list(NA), length(batch_names))
for (j in 1:length(batch_names)) list_image_mapping_master[[j]] <- read_csv(paste0(folder_script, "00-list_image_mapping-", batch_names[j], ".csv") , show_col_types = F)
list_image_mapping_master <- bind_rows(list_image_mapping_master)
list_images_master <- rep(list(NA), length(batch_names))
for (j in 1:length(batch_names)) list_images_master[[j]] <- read_csv(paste0(folder_script, "00-list_images-", batch_names[j], "-green.csv") , show_col_types = F)
list_images_master <- bind_rows(list_images_master)
list_image_mapping_folder_master <- list_image_mapping_master %>%
    left_join(tibble(image_name_pair = plates_no_colony, Undecided = "no colony"), by = "image_name_pair") %>%
    left_join(rename(list_images_master, image_name_pair = image_name), by = "image_name_pair") %>%
    left_join(select(list_images_master, image_name_isolate1 = image_name), by = "image_name_isolate1") %>%
    left_join(select(list_images_master, image_name_isolate2 = image_name), by = "image_name_isolate2")


# 10.1 Count and combine bootstrapping result ----
boots <- rep(list(NA), nrow(list_image_mapping_folder_master))
for (i in 1:nrow(list_image_mapping_folder_master)) {
    i=1
#for (i in 1:3) {
    ## Skip images with no colony
    if (list_image_mapping_folder_master$image_name_pair[i] %in% plates_no_colony) {cat("\nno colony, no watershed image\t", list_image_mapping_folder_master$image_name_pair[i]); next}

    boots[[i]] <- paste0(list_image_mapping_folder_master$folder_bootstrap[i], list_image_mapping_folder_master$image_name_pair[i], ".csv") %>%
        read_csv(show_col_types = F) %>%
        pivot_longer(cols = -ObjectID, names_to = "BootstrapID", values_to = "Isolate", names_prefix = "bootstrap_") %>%
        mutate(Isolate = factor(Isolate, c("isolate1", "isolate2"))) %>%
        mutate(BootstrapID = factor(BootstrapID, 1:1000)) %>%
        group_by(BootstrapID, Isolate, .drop = F) %>%
        count(name = "Count") %>%
        group_by(BootstrapID) %>%
        mutate(TotalCount = sum(Count)) %>%
        # Only measure isolate1
        filter(Isolate == "isolate1") %>%
        ungroup() %>%
        mutate(image_name_pair = list_image_mapping_folder_master$image_name_pair[i]) %>%
        select(image_name_pair, BootstrapID, Isolate1Count = Count, TotalCount)
    cat("\n", i, "", nrow(list_image_mapping_folder_master))
}

bind_rows(boots) %>%
    write_csv(paste0(folder_main, "meta/bootstraps.csv"))















