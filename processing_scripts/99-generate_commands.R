#' This Rscript generates a shell script containing bash commands

library(tidyverse)
source(here::here("analysis/00-metadata.R"))

create_commands <- function (batch_name) {
    c( paste0("# Batch ", batch_name),
        # Color channel
        paste0("Rscript ", list_pipeline_scripts[1], " ", folder_mapping_files, "00-list_images-", batch_names[j], "-", list_channels, ".csv"),
        # Rolling ball
        paste0("python ", list_pipeline_scripts[2], " ", folder_mapping_files, "00-list_images-", batch_names[j], "-", list_channels, ".csv"),
        # Segmentation
        paste0("Rscript ", list_pipeline_scripts[3], " ", folder_mapping_files, "00-list_images-", batch_names[j], "-green.csv"),
        # Features
        paste0("Rscript ", list_pipeline_scripts[4], " ", folder_mapping_files, "00-list_images-", batch_names[j], "-", list_channels, ".csv"),
        paste0("Rscript ", list_pipeline_scripts[5], " ", folder_mapping_files, "00-list_images-", batch_names[j], "-green.csv"),
        # Random forest
        paste0("Rscript ", list_pipeline_scripts[6], " ", folder_mapping_files, "00-list_images-", batch_names[j], "-green.csv ", folder_mapping_files, "00-list_image_mapping-", batch_names[j], ".csv")
    )
}
temp_commands <- rep(list(NA), length(batch_names))
for (j in 1:length(batch_names)) temp_commands[[j]] <- create_commands(batch_names[j])
commands <- unlist(temp_commands)


commands <- c(paste0("cd ", folder_script), commands)
# Write the commands into a bash script
writeLines(commands, here::here("analysis/99a-commands.sh"))


"
python 14-pairwise_16s_mismatch.py /Users/cychang/Dropbox/lab/emergent-coexistence/data/
"
