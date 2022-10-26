#' This script creates the folder structure associated with the image processing pipeline
#' This script only need to run once

library(tidyverse)
source(here::here("analysis/00-metadata.R"))

temp_level1 <- rep(list(NA), length(batch_names))
temp_level2 <- rep(list(NA), length(batch_names))

# Create the list of folder directory
for (j in 1:length(batch_names)) {
    temp_level1[[j]] <- c(
        paste0(folder_pipeline, "examples/"),
        paste0(folder_pipeline, "meta/"),
        paste0(folder_pipeline, "images/", batch_names[[j]], "-", list_folders, "/")
    )
    temp_level2[[j]] <- c(
        paste0(folder_pipeline, "images/", batch_names[[j]], "-", list_folders[1], "/", list_channels, "/"),
        paste0(folder_pipeline, "images/", batch_names[[j]], "-", list_folders[2], "/", list_channels, "/"),
        paste0(folder_pipeline, "images/", batch_names[[j]], "-", list_folders[3], "/", "green/"),
        paste0(folder_pipeline, "images/", batch_names[[j]], "-", list_folders[4], "/", "green/"),
        paste0(folder_pipeline, "images/", batch_names[[j]], "-", list_folders[5], "/", "green/"),
        paste0(folder_pipeline, "images/", batch_names[[j]], "-", list_folders[6], "/", list_channels, "/"),
        paste0(folder_pipeline, "images/", batch_names[[j]], "-", list_folders[7], "/", list_channels, "/"),
        paste0(folder_pipeline, "images/", batch_names[[j]], "-", list_folders[7], "/", "merged")
    )
}

list_folder_level1 <- unlist(temp_level1) %>% unique
list_folder_level2 <- unlist(temp_level2) %>% unique

# Create folders
for (k in 1:length(list_folder_level1)) {
    if (!dir.exists(list_folder_level1[k])) {
        dir.create(list_folder_level1[k])
        cat("\n", list_folder_level1[k], "created")
    }
}

for (k in 1:length(list_folder_level2)) {
    if (!dir.exists(list_folder_level2[k])) {
        dir.create(list_folder_level2[k])
        cat("\n", list_folder_level2[k], "created")
    }
}


# Remove folders. DO NOT DO THIS UNLESS CERTAIN
# if (FALSE) {
#     for (k in 1:length(list_folder_level1)) {
#         unlink(list_folder_level1[k], recursive = TRUE)
#     }
# }










