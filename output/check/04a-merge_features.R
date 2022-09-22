#' This script reads the RGB features and combines them together before feeding into models
library(tidyverse)


folder_main <- "~/Dropbox/lab/emergent-coexistence/data/raw/plate_scan/emergent_coexistence_plate_scan_check/"
folder_script <- "~/Desktop/Lab/emergent-coexistence/output/check/"

#list_images <- read_csv("~/Desktop/Lab/emergent-coexistence/output/check/00-list_images-D-red.csv", show_col_types = F)
#batch_names <- c("D", "C2", "B2", "C")
batch_names <- c("C2")
color_channels <- c("red", "green", "blue")

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


j=1
k=1
# paste0(folder_script, "00-list_images-", batch_names[j], "-", color_channels[k], ".csv") %>%
#     read_csv(show_col_types = F)


for (j in 1:length(batch_names)) {
    cat("\nbatch ", batch_names[j])
    list_images <- list(
        read_csv(paste0(folder_script, "00-list_images-", batch_names[j], "-red.csv"), show_col_types = F),
        read_csv(paste0(folder_script, "00-list_images-", batch_names[j], "-green.csv"), show_col_types = F),
        read_csv(paste0(folder_script, "00-list_images-", batch_names[j], "-blue.csv"), show_col_types = F)
    )

    for (i in 1:nrow(list_images[[1]])) {
        if (list_images[[1]]$image_name[i] %in% plates_no_colony) next
        temp <- NULL
        for (k in 1:3) {
            temp[[k]] <- paste0(list_images[[k]]$folder_feature[i], color_channels[k], "/", list_images[[k]]$image_name[i], ".csv") %>%
                read_csv(show_col_types = F) %>%
                mutate(ColorChannel = color_channels[k]) %>%
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
            write_csv(paste0(list_images[[k]]$folder_feature[i], "merged/", list_images[[k]]$image_name[i], ".csv"))

        cat("\t", i)

    }
}














