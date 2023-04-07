#' This script reads the RGB features and combines them together before feeding into models
#'
#' To use this Rscript, in bash environment
#' Rscript 04a-feature.R `mapping_files.csv`


library(tidyverse)
source(here::here("processing_scripts/00-metadata.R"))

list_images <- read_csv(commandArgs(trailingOnly = T)[1], show_col_types = F)

for (i in 1:nrow(list_images)) {
    if (list_images$image_name[i] %in% plates_no_colony) next
    temp <- NULL
    for (k in 1:3) {
        temp[[k]] <- paste0(list_images$folder_feature[i], list_channels[k], "/", list_images$image_name[i], ".csv") %>%
            read_csv(show_col_types = F) %>%
            mutate(ColorChannel = list_channels[k]) %>%
            select(ColorChannel, ObjectID, everything())
    }

    #' There are some missing NA, because outlier objectID is different across RGB
    #' Remove rows containing NA because random forest does not handle NA
    bind_rows(temp) %>%
        pivot_wider(id_cols = ObjectID, names_from = ColorChannel,
                    values_from = c(starts_with("s."), starts_with("m."), starts_with("b."), starts_with("t.")),
        ) %>%
        select(ObjectID, ends_with("green"), starts_with("b.mean"), starts_with("b.sd"), starts_with("b.mad")) %>%
        select(-contains("t.bump")) %>%
        # Drop rows that contains NA because of outliers
        drop_na() %>%
        write_csv(paste0(list_images$folder_feature[i], "merged/", list_images$image_name[i], ".csv"))

    cat("\t", i)

}















