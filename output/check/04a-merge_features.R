#' This script reads the RGB features and combines them together before feeding into models
library(tidyverse)


batch_names <- c("D", "C2", "B2", "C")
j=1

folder_main <- "~/Dropbox/lab/emergent-coexistence/data/raw/plate_scan/emergent_coexistence_plate_scan_check/"
folder_script <- "~/Desktop/Lab/emergent-coexistence/output/check/"



list_images <- read_csv("~/Desktop/Lab/emergent-coexistence/output/check/00-list_images-D-red.csv", show_col_types = F)



for (j in 1:length(batch_names)) {

}
