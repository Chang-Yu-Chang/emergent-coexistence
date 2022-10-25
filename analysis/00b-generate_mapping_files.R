#' This scripts generates the master csv for image file name and directory

library(tidyverse)
source(here::here("analysis/00-metadata.R"))

#
for (j in 1:length(batch_names)) {
    # 0.1 list of existing original images ----
    folder_original <- paste0(folder_main, "images/", batch_names[j], "-00-original/")
    image_names <-
        list.files(folder_original) %>%
        # Remove all mixing pairs, for example with suffix 1_1, 2_2, etc
        #str_subset(paste0("^((?!", paste(paste0("_", 1:13, "_", 1:13), collapse = "|"), ").)*$")) %>%
        str_subset(paste0("^((?!", paste(paste0("_", 3:13, "_", 3:13), collapse = "|"), ").)*$")) %>%
        str_subset(paste0("^((?!_1_1\\.).)*$")) %>%
        str_subset(paste0("^((?!_2_2\\.).)*$")) %>%
        # Remove folders
        str_subset(".tiff") %>%
        #str_subset("T8") %>%
        str_replace(".tiff", "") %>%
        # Remove all that contain _-, which was a naming convention for the different dilution factors on the same plate
        str_subset("^((?!_-).)*$") %>%
        sort()

    # Manual key in plates using different naming convention
    if (batch_names[j] == "D") {
        image_names <- c(image_names, "D_T8_C4R1_50-50_1_3_-4") %>% sort
    } else if (batch_names[j] == "C2") {

    }

    # 0.2 list of image files and the folders to store them ----
    n_images <- length(image_names)

    # List of image file names and folders
    list_images <- tibble(
        image_name = image_names,
        folder_original = rep(paste0(folder_main, "images/", batch_names[j], "-00-original/"), n_images),
        folder_channel = rep(paste0(folder_main, "images/", batch_names[j], "-", list_folders[1], "/"), n_images),
        folder_rolled = rep(paste0(folder_main, "images/", batch_names[j], "-", list_folders[2],"/"), n_images),
        folder_threshold = rep(paste0(folder_main, "images/", batch_names[j], "-", list_folders[3], "/"), n_images),
        folder_round = rep(paste0(folder_main, "images/", batch_names[j], "-", list_folders[4], "/"), n_images),
        folder_watershed = rep(paste0(folder_main, "images/", batch_names[j], "-", list_folders[5], "/"), n_images),
        folder_transect = rep(paste0(folder_main, "images/", batch_names[j], "-", list_folders[6], "/"), n_images),
        folder_feature = rep(paste0(folder_main, "images/", batch_names[j], "-", list_folders[7], "/"), n_images),
        folder_random_forest = rep(paste0(folder_main, "images/", batch_names[j], "-", list_folders[8], "/"), n_images),
        folder_bootstrap = rep(paste0(folder_main, "images/", batch_names[j], "-", list_folders[9], "/"), n_images),
        folder_combined = rep(paste0(folder_main, "images/", batch_names[j], "-", list_folders[10], "/"), n_images)
    )

    # Repeat the rows 3 times for rgb channels
    for (color in c("red", "green", "blue")) {
        list_images %>%
            mutate(color_channel = color) %>%
            select(image_name, color_channel, everything()) %>%
            write_csv(paste0(folder_script, "mapping_files/00-list_images-", batch_names[j], "-", color, ".csv"))
        cat("\n", paste0(folder_script, "mapping_files/00-list_images-", batch_names[j], "-", color, ".csv"), "\tcreated")
    }


    # 0.3 Mapping file between isolate and pairs ----
    #' This section does not need rgb
    #' For example, D_T8_C1R7_3 has a length of 4 and is an isolate image
    #' whereas D_T8_C1R7_5-95_1_3 has a length of 6 and is a pair image
    name_length <- image_names %>% str_split("_") %>% sapply(length)
    # Include the plate image with addition dilution factor specified
    # ' This will spit out a warning message for additional pieces in 1 row. Ignore it.
    index_pair <- name_length == 6 | name_length == 7

    ## Parse file name
    list_image_pairs <- tibble(image_name_pair = image_names[index_pair]) %>%
        separate(col = image_name_pair, sep = "_", into = c("Batch", "temp", "Community", "Freqs", "Isolate1", "Isolate2"), remove = F) %>%
        # Reverse the frequency order here because the original frequencies on image file names and plate labels wer in a wrong order
        separate(col = Freqs, sep = "-", into = c("Freq2", "Freq1"), remove = T) %>%
        select(Batch, Community, Isolate1, Isolate2, Freq1, Freq2, image_name_pair)

    list_image_isolates <- tibble(image_name_isolate = image_names[!index_pair]) %>%
        separate(col = image_name_isolate, sep = "_", into = c("Batch", "temp", "Community", "Isolate"), remove = F) %>%
        select(Batch, Community, Isolate, image_name_isolate)

    ## Combine pair and isolate image lists
    list_image_mapping <- list_image_pairs %>%
        left_join(rename(list_image_isolates, Isolate1 = Isolate, image_name_isolate1 = image_name_isolate), by = c("Batch", "Community", "Isolate1")) %>%
        left_join(rename(list_image_isolates, Isolate2 = Isolate, image_name_isolate2 = image_name_isolate), by = c("Batch", "Community", "Isolate2"))

    write_csv(list_image_mapping, paste0(folder_script, "mapping_files/00-list_image_mapping-", batch_names[j], ".csv"))
    cat("\n", paste0(folder_script, "mapping_files/00-list_image_mapping-", batch_names[j], ".csv"), "\tcreated")

}

# A master mapping csv ----
list_image_mapping_master <- rep(list(NA), length(batch_names))
for (j in 1:length(batch_names)) list_image_mapping_master[[j]] <- read_csv(paste0(folder_script, "mapping_files/00-list_image_mapping-", batch_names[j], ".csv") , show_col_types = F)
list_image_mapping_master <- bind_rows(list_image_mapping_master)

pairs_freq_ID <- list_image_mapping_master %>%
    rename(Isolate1InitialODFreq = Freq1, Isolate2InitialODFreq = Freq2) %>%
    # Correct the isolate order
    rowwise() %>%
    mutate(Isolate1InitialODFreq = ifelse(Isolate1 > Isolate2, 5, Isolate1InitialODFreq),
           Isolate2InitialODFreq = ifelse(Isolate1 > Isolate2, 95, Isolate2InitialODFreq),
           FlipOrder = ifelse(Isolate1 > Isolate2, T, F)
    ) %>%
    mutate(temp = min(Isolate1,Isolate2), Isolate2 = max(Isolate1, Isolate2), Isolate1 = temp) %>%
    select(-temp) %>%
    mutate(Community = factor(Community, paste0("C", rep(1:12, each = 8), "R", rep(1:8, 12)))) %>%
    arrange(Community, Isolate1, Isolate2, Isolate1InitialODFreq) %>%
    ungroup() %>%
    select(-FlipOrder) %>%
    # Remove staph contamination
    filter(!(Community == "C11R2" & Isolate1 == 13)) %>%
    filter(!(Community == "C11R2" & Isolate2 == 13)) %>%
    filter(!(Batch == "B2" & Community == "C11R1" & Isolate1 == 1)) %>%
    filter(!(Batch == "C" & Community == "C11R1" & Isolate1 == 5))


## Two species pairs miss the 50-50 data
pairs_freq_ID %>%
    group_by(Batch, Community, Isolate1, Isolate2) %>%
    count() %>%
    filter(n != 3)
#   Batch Community Isolate1 Isolate2     n
#   <chr> <fct>        <dbl>    <dbl> <int>
# 1 C     C11R1            1        2     2
# 2 C     C11R1            1        3     2

## Append these two frequencies back
pairs_freq_ID <- tibble(Batch = c("C", "C"), Community = c("C11R1", "C11R1"),
                        Isolate1 = c(1,1), Isolate2 = c(2,3),
                        Isolate1InitialODFreq = c(50, 50), Isolate2InitialODFreq = c(50, 50)) %>%
    bind_rows(pairs_freq_ID) %>%
    mutate(Community = factor(Community, paste0("C", rep(1:12, each = 8), "R", rep(1:8, 12)))) %>%
    mutate(Isolate1 = factor(Isolate1, 1:13), Isolate2 = factor(Isolate2, 1:13)) %>%
    arrange(Community, Isolate1, Isolate2, Isolate1InitialODFreq)
write_csv(pairs_freq_ID, paste0(folder_main, "meta/00b-pairs_freq_ID.csv"))
cat("\n", paste0(folder_main, "meta/00b-pairs_freq_ID.csv"), "\tcreated")








