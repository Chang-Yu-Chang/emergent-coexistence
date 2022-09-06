# Generate the master csv for image file name and directory

library(tidyverse)

# This main folder depends on your home directory and user name. Python somehow does not read ~/ instead I have to specify /Users/chang-yu/
folder_main <- "/Users/chang-yu/Dropbox/lab/emergent-coexistence/data/raw/plate_scan/emergent_coexistence_plate_scan_check/"
folder_script <- "/Users/chang-yu/Desktop/Lab/emergent-coexistence/output/check/"

batch_names <- c("D", "C", "C2", "B2", "chromo")


for (j in 1:length(batch_names)) {
    folder_original <- paste0(folder_main, "check/", batch_names[j], "-00-original/")
    image_names <- list.files(folder_original) %>%
        # Remove folders
        str_subset(".tiff") %>%
        #str_subset("T8") %>%
        str_replace(".tiff", "") %>%
        # Remove all that contain _-, which was a naming convention for the different dilution factors on the same plate
        # 20220902 The chromogenic plates all have it. Come check later
        str_subset("^((?!_-).)*$") %>%
        # Remove all mixing pairs, for example with suffix 1_1, 2_2, etc
        str_subset(paste0("^((?!", paste(paste0("_", 1:13, "_", 1:13), collapse = "|"), ").)*$")) %>%
        sort()

    n_images <- length(image_names)

    # List of image file names and folders
    list_images <- tibble(
        image_name = image_names,
        folder_original = rep(paste0(folder_main, "check/", batch_names[j], "-00-original/"), n_images),

        folder_green = rep(paste0(folder_main, "check/", batch_names[j], "-01-green_channel/"), n_images),
        folder_green_rolled = rep(paste0(folder_main, "check/", batch_names[j], "-02-green_rolled/"), n_images),
        folder_green_watershed_file = rep(paste0(folder_main, "check/", batch_names[j], "-05-green_watershed_file/"), n_images),
        folder_green_watershed = rep(paste0(folder_main, "check/", batch_names[j], "-06-green_watershed/"), n_images),
        folder_green_feature = rep(paste0(folder_main, "check/", batch_names[j], "-07-green_feature/"), n_images),
        folder_green_cluster = rep(paste0(folder_main, "check/", batch_names[j], "-08-green_cluster/"), n_images),

        folder_red = rep(paste0(folder_main, "check/", batch_names[j], "-11-red_channel/"), n_images),
        folder_red_rolled = rep(paste0(folder_main, "check/", batch_names[j], "-12-red_rolled/"), n_images),
        folder_red_watershed_file = rep(paste0(folder_main, "check/", batch_names[j], "-15-red_watershed_file/"), n_images),
        folder_red_watershed = rep(paste0(folder_main, "check/", batch_names[j], "-16-red_watershed/"), n_images),
        folder_red_feature = rep(paste0(folder_main, "check/", batch_names[j], "-17-red_feature/"), n_images),
        folder_red_cluster = rep(paste0(folder_main, "check/", batch_names[j], "-18-red_cluster/"), n_images),

        folder_blue = rep(paste0(folder_main, "check/", batch_names[j], "-21-blue_channel/"), n_images),
        folder_blue_rolled = rep(paste0(folder_main, "check/", batch_names[j], "-22-blue_rolled/"), n_images),
        folder_blue_watershed_file = rep(paste0(folder_main, "check/", batch_names[j], "-25-blue_watershed_file/"), n_images),
        folder_blue_watershed = rep(paste0(folder_main, "check/", batch_names[j], "-26-blue_watershed/"), n_images),
        folder_blue_feature = rep(paste0(folder_main, "check/", batch_names[j], "-27-blue_feature/"), n_images),
        folder_blue_cluster = rep(paste0(folder_main, "check/", batch_names[j], "-28-blue_cluster/"), n_images)
    )

    write_csv(list_images, paste0(folder_script, "00-list_images-", batch_names[j], ".csv"))

    ## if T8 isolate image has no-growth/containmination, use T0 image instead.
    #' Keep populating this list if found new ones
    #' D T0 C1R4 3
    #' D T0 C1R6 3
    list_images %>%
        filter(str_detect(image_names, "T0")) %>%
        write_csv(paste0(folder_script, "00-list_images-no_growth_isolates.csv"))

    # Mapping file between isolate and pairs
    #' For example, D_T8_C1R7_3 has a length of 4 and is an isolate image
    #' whereas D_T8_C1R7_5-95_1_3 has a length of 6 and is a pair image
    name_length <- image_names %>% str_split("_") %>% sapply(length)
    index_pair <- name_length == 6

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




























