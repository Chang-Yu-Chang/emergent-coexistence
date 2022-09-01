library(EBImage)
library(tidyverse)

# Call the list of image files
folder_source <- "transtivity-B2/pairs/"
folder_output <- "B2/"
list_img <- list.files(paste0("~/Dropbox/lab/emergent-coexistence/data/raw/plan_scan/emergent_coexistence_plate_scan_check/", folder_source)) %>%
    str_subset("T8") %>%
    str_replace(".tiff", "") %>%
    # Remove all that contain _-, which was a naming convention for the different dilution factors on the same plate
    str_subset("^((?!_-).)*$") %>%
    # Remove all 1_1, 2_2, etc mixing pairs
    str_subset(paste0("^((?!", paste(paste0(1:10, "_", 1:10), collapse = "|"), ").)*$"))

#
for (i in 1:length(list_img)) {
    img_name <- list_img[i]
    p <- readImage(paste0("~/Dropibox/lab/emergent-coexistence/data/raw/plan_scan/emergent_coexistence_plate_scan_check/", folder_source, img_name, ".tiff"))
    colorMode(p) = Grayscale
    # Extract only the green channel
    p <- p[,,2]
    output_file <- paste0("~/Dropbox/lab/emergent-coexistence/data/raw/plan_scan/emergent_coexistence_plate_scan_check/check/", folder_output, img_name, ".tiff")
    if (!file.exists(output_file)) {
        writeImage(p, output_file, quality = 85)
    }
    cat("\n", i, "/", length(list_img), list_img[i])
}




# Individual case where the T8 plate was contaminated
list(
    c("transtivity-D/single_isolates/", "D_T0_C1R4_3"),
    c("transtivity-D/single_isolates/", "D_T0_C1R6_3"),
    c("transtivity-B2/single_isolates/", "D_T0_C1R6_3"),
)

img_name <- "D_T0_C1R6_3"
image_name <- "B2_T1_C2R6_1_20180929"
folder_name <- "transtivity-D/single_isolates/"
p <- readImage(paste0("~/Dropbox/lab/emergent-coexistence/data/raw/plan_scan/emergent_coexistence_plate_scan_check/", folder_name, img_name, ".tiff"))
colorMode(p) = Grayscale
# Extract only the green channel
p <- p[,,2]
writeImage(p, paste0("~/Dropbox/lab/emergent-coexistence/data/raw/plan_scan/emergent_coexistence_plate_scan_check/check/", img_name, ".tiff"), quality = 85)
