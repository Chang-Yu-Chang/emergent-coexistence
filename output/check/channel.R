library("EBImage")
library(tidyverse)



list_img <- list.files("~/Dropbox/lab/emergent-coexistence/data/raw/plan_scan/emergent_coexistence_plate_scan_check/transtivity-D/pairs") %>%
    str_subset("T8") %>%
    str_replace(".tiff", "")

for (i in 1:length(list_img)) {
    img_name <- list_img[i]
    p <- readImage(paste0("~/Dropbox/lab/emergent-coexistence/data/raw/plan_scan/emergent_coexistence_plate_scan_check/transtivity-D/pairs/", img_name, ".tiff"))
    colorMode(p) = Grayscale
    # Extract only the green channel
    p <- p[,,2]
    writeImage(p, paste0("~/Dropbox/lab/emergent-coexistence/data/raw/plan_scan/emergent_coexistence_plate_scan_check/check/", img_name, ".tiff"), quality = 85)
    cat("\n", i, "//", length(list_img), list_img[i])
}


# Individual case where the T8 plate was contaminated
#img_name <- "D_T0_C1R4_3"
img_name <- "D_T0_C1R6_3"
folder_name <- "transtivity-D/single_isolates/"
p <- readImage(paste0("~/Dropbox/lab/emergent-coexistence/data/raw/plan_scan/emergent_coexistence_plate_scan_check/", folder_name, img_name, ".tiff"))
colorMode(p) = Grayscale
# Extract only the green channel
p <- p[,,2]
writeImage(p, paste0("~/Dropbox/lab/emergent-coexistence/data/raw/plan_scan/emergent_coexistence_plate_scan_check/check/", img_name, ".tiff"), quality = 85)
