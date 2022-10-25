# This script stores the metadata shared by all scripts

library(tidyverse)

# This main folder depends on your home directory and user name. Python somehow does not read ~/ instead I have to specify /Users/chang-yu/

folder_script <- "~/Desktop/lab/emergent-coexistence/analysis/" # Enter the directory of analysis scripts
folder_main <- "~/Dropbox/lab/emergent-coexistence/plate_scan_pipeline/" # Enter the directory of data
# folder_main <- "/Users/chang-yu/Dropbox/lab/emergent-coexistence/plate_scan_pipeline/"
# folder_script <- "/Users/chang-yu//Desktop/lab/emergent-coexistence/analysis/"
folder_main <- "/Users/cychang/Dropbox/lab/emergent-coexistence/plate_scan_pipeline/"
folder_script <- "/Users/cychang/Desktop/lab/emergent-coexistence/analysis/"

# Metadata for shared by the analysis scripts
list_folders <- c("01-channel", "02-rolled", "03-threshold", "04-round", "05-watershed", "06-transet", "07-feature", "08-random_forest", "09-bootstrap", "10-images_and_random_forest")
list_channels <- c("red", "green", "blue")
list_pipeline_scipts <- c("01-channel.R", "02-rolling_ball.py", "03-segmentation.R", "04-feature.R", "04a-merge_features.R", "05-random_forest.R")
batch_names <- c("B2", "C", "C2", "D")

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
