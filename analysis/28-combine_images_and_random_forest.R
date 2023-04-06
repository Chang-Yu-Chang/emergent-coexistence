#' This script combines the transect image with the random forest results
#' 1. Append the master mapping files with the file directory
#' 2. combine images and random forest results
#' 3. To merge png from a folder into a single pdf. Use imagemagick. Go to the folder and convert *.jpg XXX.pdf

library(tidyverse)
library(EBImage)
library(cowplot)
library(magick)
source(here::here("analysis/00-metadata.R"))

folder_mapping_files <- "mapping_files/"
# 1. Append the master mapping file with the file directory ----
list_image_mapping_master <- read_csv(paste0(folder_script, folder_mapping_files, "00-list_image_mapping_folder_master.csv"), show_col_types = F) %>%
    # Correct the isolate order
    rename(Isolate1InitialODFreq = Freq1, Isolate2InitialODFreq = Freq2) %>%
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
    select(-FlipOrder)

pairs_freq <- read_csv(paste0(folder_data, "temp/25-pairs_freq.csv"), show_col_types = F) %>%
    filter(Time == "T0") %>%
    select(PairFreqID, Batch, Community, Isolate1, Isolate2, Isolate1InitialODFreq)

list_image_mapping_master <- list_image_mapping_master %>%
    right_join(pairs_freq)

# 2. Combine image pngs with random forest results -----
for (i in 1:nrow(list_image_mapping_master)) {
    folder_original <- list_image_mapping_master$folder_original[i]
    folder_rolled <- paste0(list_image_mapping_master$folder_rolled[i], list_image_mapping_master$color_channel[i], "/")
    image_name <- list_image_mapping_master$image_name_pair[i]

    coculture_image_name <- paste0(folder_rolled, list_image_mapping_master$image_name_pair[i], ".tiff")
    random_forest_image_name <- paste0(list_image_mapping_master$folder_random_forest[i], list_image_mapping_master$image_name_pair[i], ".png")
    if (!file.exists(coculture_image_name)) {
        cat(coculture_image_name, "does not exist")
        next
    }

    if (!file.exists(random_forest_image_name)) {
        cat(random_forest_image_name, "doest not exist")
        next
    }

    p <- ggdraw() +
        theme(plot.background = element_rect(fill = "white")) +
        # Isolate1 original
        draw_image(paste0(folder_rolled, list_image_mapping_master$image_name_isolate1[i], ".tiff"),
                   x = -.33, y = .25, hjust = 0, vjust = 0, scale = .3) +
        # Isolate2 original
        draw_image(paste0(folder_rolled, list_image_mapping_master$image_name_isolate2[i], ".tiff"),
                   x = 0, y = .25, hjust = 0, vjust = 0, scale = .3) +
        # Coculture original
        draw_image(paste0(folder_rolled, list_image_mapping_master$image_name_pair[i], ".tiff"),
                   x = .33, y = .25, hjust = 0, vjust = 0, scale = .3) +
        # Random forest outcome
        draw_image(paste0(list_image_mapping_master$folder_random_forest[i], list_image_mapping_master$image_name_pair[i], ".png"),
                   x = 0, y = -.22, hjust = 0, vjust = 0, scale = .9) +
        # Annotation
        annotate("text", x = .02, y = .98, label = list_image_mapping_master$image_name_pair[i], size = 5, hjust = 0, fontface = "bold") +
        annotate("text", x = .17, y = .93, label = "Isolate1", size = 5, hjust = .5, fontface = "bold") +
        annotate("text", x = .5, y = .93, label = "Isolate2", size = 5, hjust = .5, fontface = "bold") +
        annotate("text", x = .83, y = .93, label = "Coculture", size = 5, hjust = .5, fontface = "bold") +
        theme(plot.background = element_rect(color = NA, fill = "white"))
    #annotate("text", x = .5, y = .45, label = "Random Forest Result", size = 5, hjust = .5, fontface = "bold")

    # ggsave(paste0(folder_pipeline, "images/", list_image_mapping_master$Batch[i], "-10-images_and_random_forest/", image_name, ".png"),
    #        plot = p, width = 10, height = 9, dpi = 300)
    ggsave(paste0(folder_pipeline, "random_forest/", image_name, ".png"), plot = p, width = 10, height = 9, dpi = 300)

    cat("\n", image_name, "\t", i, "/", nrow(list_image_mapping_master))
}

