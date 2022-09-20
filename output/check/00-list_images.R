# Generate the master csv for image file name and directory

library(tidyverse)

# This main folder depends on your home directory and user name. Python somehow does not read ~/ instead I have to specify /Users/chang-yu/
# folder_main <- "~/Dropbox/lab/emergent-coexistence/data/raw/plate_scan/emergent_coexistence_plate_scan_check/"
# folder_script <- "~/Desktop/Lab/emergent-coexistence/output/check/"

folder_main <- "/Users/chang-yu/Dropbox/lab/emergent-coexistence/data/raw/plate_scan/emergent_coexistence_plate_scan_check/"
folder_script <- "/Users/chang-yu//Desktop/Lab/emergent-coexistence/output/check/"
# folder_main <- "/Users/cychang/Dropbox/lab/emergent-coexistence/data/raw/plate_scan/emergent_coexistence_plate_scan_check/"
# folder_script <- "/Users/cychang/Desktop/Lab/emergent-coexistence/output/check/"


#batch_names <- c("D", "C", "C2", "B2", "chromo")
batch_names <- c("D", "C2", "B2", "C")
j=1
for (j in 1:length(batch_names)) {

    # 0.1 list of exisitng original images ----
    folder_original <- paste0(folder_main, "check/", batch_names[j], "-00-original/")
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
        # 20220902 The chromogenic plates all have it. Come check later
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
        folder_original = rep(paste0(folder_main, "check/", batch_names[j], "-00-original/"), n_images),
        folder_channel = rep(paste0(folder_main, "check/", batch_names[j], "-01-channel/"), n_images),
        folder_rolled = rep(paste0(folder_main, "check/", batch_names[j], "-02-rolled/"), n_images),
        folder_threshold = rep(paste0(folder_main, "check/", batch_names[j], "-03-threshold/"), n_images),
        folder_round = rep(paste0(folder_main, "check/", batch_names[j], "-04-round/"), n_images),
        folder_watershed = rep(paste0(folder_main, "check/", batch_names[j], "-05-watershed/"), n_images),
        folder_transection = rep(paste0(folder_main, "check/", batch_names[j], "-06-transection/"), n_images),
        folder_feature = rep(paste0(folder_main, "check/", batch_names[j], "-07-feature/"), n_images),
        folder_logit = rep(paste0(folder_main, "check/", batch_names[j], "-08-logit/"), n_images),
        folder_random_forest = rep(paste0(folder_main, "check/", batch_names[j], "-09-random_forest/"), n_images)
    )

    # Repeat the rows 3 times for rgb channels
    for (color in c("red", "green", "blue")) {
        list_images %>%
            mutate(color_channel = color) %>%
            select(image_name, color_channel, everything()) %>%
            write_csv(paste0(folder_script, "00-list_images-", batch_names[j], "-", color, ".csv"))
    }


    ## if T8 isolate image has no-growth/containmination, use T0 image instead.
    #' Keep populating this list if found new ones
    #' D T0 C1R4 3
    #' D T0 C1R6 3
    #' D T1 C1R7 7
    # list_images %>%
    #     filter(str_detect(image_names, "T0") | str_detect(image_names, "T1")) %>%
    #     write_csv(paste0(folder_script, "00-list_images-no_growth_isolates.csv"))

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

    write_csv(list_image_mapping, paste0(folder_script, "00-list_image_mapping-", batch_names[j], ".csv"))

}




























