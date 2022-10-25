#' This script examines the number of file in each folder


library(tidyverse)
source(here::here("analysis/00-metadata.R"))
source(here::here("analysis/00a-folder_structure.R"))

count_tiff <- function (folder_directory) {
    folder_directory %>%
        list.files %>%
        str_subset(".tiff") %>%
        length
}

count_tiff(list_folder_level2[4])
