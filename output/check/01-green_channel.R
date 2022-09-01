#' This is the step 1 in making the colony morphologies distinct
#' TSA plate image, green channel, grey scale, subtract background
#'  The processed image files are stored in folders "Batch"_green/
#'  - check/B2_green/
#'  - check/D_green/
#'  - etc..
#' In this step, the majority of the pairs can be distinguished (XXX number out of 126 pairs)
#'
#' Use R EBImage to extract only the green channel from the original RGB image, and turn it into grey scale. This step is done in the R script (channel.R)
#' In ImageJ, remove the background using Process > Subtract Background > Light Background, with rolling ball radius 80 pixels -> click OK


"
Do not repeatedly run this script if the pro if the image has been processed in imageJ!
"


library(tidyverse)
library(EBImage)

# This is done in
list_batches <- c("C2")

for (j in 1:length(list_batches)) {
    folder_source <- paste0("transtivity-", list_batches[j],"/pairs/")
    folder_output <- paste0(list_batches[j],"_green/")
    list_img <- list.files(paste0("~/Dropbox/lab/emergent-coexistence/data/raw/plate_scan/emergent_coexistence_plate_scan_check/", folder_source)) %>%
        str_subset("T8") %>%
        str_replace(".tiff", "") %>%
        # Remove all that contain _-, which was a naming convention for the different dilution factors on the same plate
        str_subset("^((?!_-).)*$") %>%
        # Remove all mixing pairs, for example with suffix 1_1, 2_2, etc
        str_subset(paste0("^((?!", paste(paste0("_", 1:13, "_", 1:13), collapse = "|"), ").)*$"))


    for (i in 1:length(list_img)) {
        # File name
        img_name <- list_img[i]
        output_file <- paste0("~/Dropbox/lab/emergent-coexistence/data/raw/plate_scan/emergent_coexistence_plate_scan_check/check/", folder_output, img_name, ".tiff")
        cat("\n", i, "/", length(list_img), list_img[i])

        # If the output files already exist, dont overwrite it
        if (file.exists(output_file)) {
            cat("\tThe green channel image file already exists")

        } else if (!file.exists(output_file)) {
            # Read the original image
            p <- readImage(paste0("~/Dropbox/lab/emergent-coexistence/data/raw/plate_scan/emergent_coexistence_plate_scan_check/", folder_source, img_name, ".tiff"))
            colorMode(p) = Grayscale

            # Extract only the green channel
            p <- p[,,2]
            writeImage(p, output_file, quality = 85)
        }
    }
}




# Individual case where the T8 plate was contaminated. Instead, use the T0 or T1 plates with clear morphologies
# The order are: source folder, isolate plate image file name, and output folder
list_individuals <- list(
    c("transtivity-B2/single_isolates/", "B2_T1_C2R6_1_20180929", "check/B2_green/"),
    c("transtivity-D/single_isolates/", "D_T0_C1R4_3", "check/D_green/"),
    c("transtivity-D/single_isolates/", "D_T0_C1R6_3", "check/D_green/"),
)

for (i in 1:length(list_individuals)) {
    # File names
    folder_source <- list_individuals[[i]][1]
    img_name <- list_individuals[[i]][2]
    folder_output <- list_individuals[[i]][3]
    output_file <- paste0("~/Dropbox/lab/emergent-coexistence/data/raw/plate_scan/emergent_coexistence_plate_scan_check/", folder_output, img_name, ".tiff")
    cat("\n", i, "/", length(list_individuals), list_individuals[[i]][2])

    # If the output files already exist, dont overwrite it
    if (file.exists(output_file)) {
        cat("\tThe green channel image file already exists")

    } else if (!file.exists(output_file)) {
        # Read the original image
        p <- readImage(paste0("~/Dropbox/lab/emergent-coexistence/data/raw/plate_scan/emergent_coexistence_plate_scan_check/", folder_source, img_name, ".tiff"))
        colorMode(p) = Grayscale
        # Extract only the green channel
        p <- p[,,2]
        writeImage(p, output_file, quality = 85)
    }


}













