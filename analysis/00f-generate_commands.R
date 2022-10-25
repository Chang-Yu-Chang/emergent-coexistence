#' This Rscript generates a bash command

library(tidyverse)
source(here::here("analysis/00-metadata.R"))

list_mapping_files <- list.files(paste0(folder_script, "mapping_files")) %>%
    str_subset("green|blue|red") %>%
    paste0("mapping_files/", .)

temp_commands <- rep(list(NA), length(list_pipeline_scipts))
for (i in 1:length(list_pipeline_scipts)) {
    if (list_pipeline_scipts[i] %in% c("01-channel.R", "04-feature.R")) {
        temp_commands[[i]] <- paste0("Rscript ", list_pipeline_scipts[i], " ", list_mapping_files)
    } else if (list_pipeline_scipts[i] == "02-rolling_ball.py") {
        temp_commands[[i]] <- paste0("python ", list_pipeline_scipts[i], " ", list_mapping_files)
    } else if (list_pipeline_scipts[i] %in% c("03-segmentation.R", "04a-merge_features.R", "05-random_forest.R")){
        temp_commands[[i]] <- paste0("Rscript ", list_pipeline_scipts[i], " ", str_subset(list_mapping_files, "green"))
    }
}


# Write the commands into a bash script
unlist(temp_commands) %>%
    writeLines(here::here("analysis/00g-commands.sh"))
