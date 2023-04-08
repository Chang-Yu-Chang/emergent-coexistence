#' This script generates the list of pairs and frequency and their ID
#'
#' For the image processing pipeline, I processed all coculture images I have, including
#' those that may be contamination

library(tidyverse)
source(here::here("processing_scripts/00-metadata.R"))

list_image_mapping_folder_master <- read_csv(paste0("image_scripts/mapping_files/", "00-list_image_mapping_folder_master.csv"), show_col_types = F)

# 1. Isolates ----
# This one reads the old, hand-curated csv to generate ID
isolates_ID <- read_csv(paste0(folder_data, "raw/isolates1.csv"), col_types = cols())
isolates_ID$ID[isolates_ID$ExpID == "2.6.A.5"] <- "2.6.A.5"

nrow(isolates_ID) # 65 isolates
write_csv(isolates_ID, paste0(folder_data, "temp/00c-isolates_ID.csv"))
cat("\n", paste0(folder_data, "temp/00c-isolates_ID.csv"), "\tcreated")

# 2. Communities ----
communities_name <- c("C1R2", "C1R4", "C1R6", "C1R7", "C2R6", "C2R8", "C4R1", "C7R1", "C8R4", "C11R1", "C11R2", "C11R5")
communities_size <- c(4,5,5,7,4,4,3,4,3,9,12,5)
pp <- function(x) choose(x, 2)

communities <- data.frame(
    Community = communities_name,
    CommunitySize = communities_size,
    CommunityPairSize = pp(communities_size)
) %>%
    mutate(Community = factor(Community, communities_name))  %>%
    arrange(CommunitySize) %>%
    mutate(CommunityLabel = 1:12) %>%
    select(Community, CommunityLabel, everything())
nrow(communities)
write_csv(communities, paste0(folder_data, "temp/00c-communities.csv"))
cat("\n", paste0(folder_data, "temp/00c-communities.csv"), "\tcreated")






# 3. A mapping file for pairs and frequencies ----
pairs_freq_ID <- list_image_mapping_folder_master %>%
    rename(Isolate1InitialODFreq = Freq1, Isolate2InitialODFreq = Freq2) %>%
    # Correct the isolate order
    rowwise() %>%
    mutate(Isolate1InitialODFreq = ifelse(Isolate1 > Isolate2, 5, Isolate1InitialODFreq),
           Isolate2InitialODFreq = ifelse(Isolate1 > Isolate2, 95, Isolate2InitialODFreq),
           FlipOrder = ifelse(Isolate1 > Isolate2, T, F)
    ) %>%
    mutate(temp = min(Isolate1,Isolate2), Isolate2 = max(Isolate1, Isolate2), Isolate1 = temp) %>%
    select(-temp) %>%
    ungroup() %>%
    select(-FlipOrder) %>%
    bind_rows(tibble(
        ## Append two missing images to ensure the mapping files run ok
        Batch = c("C", "C"), Community = c("C11R1", "C11R1"),
        Isolate1 = c(1,1), Isolate2 = c(2,3),
        Isolate1InitialODFreq = c(50, 50), Isolate2InitialODFreq = c(50, 50)
    )) %>%
    mutate(Community = factor(Community, paste0("C", rep(1:12, each = 8), "R", rep(1:8, 12)))) %>%
    mutate(Isolate1 = factor(Isolate1, 1:12), Isolate2 = factor(Isolate2, 1:12)) %>%
    arrange(Community, Isolate1, Isolate2, Isolate1InitialODFreq)
nrow(pairs_freq_ID) # 159*3 = 477

pairs_ID <- pairs_freq_ID %>%
    distinct(Batch, Community, Isolate1, Isolate2) %>%
    mutate(PairID = 1:n()) %>%
    select(PairID, everything())
nrow(pairs_ID) # 159 pairs in the image pipeline and pairwise competiton

write_csv(pairs_freq_ID, paste0(folder_data, "temp/00c-pairs_freq_ID.csv"))
cat("\n", paste0(folder_data, "temp/00c-pairs_freq_ID.csv"), "\tcreated")

write_csv(pairs_ID, paste0(folder_data, "temp/00c-pairs_ID.csv"))
cat("\n", paste0(folder_data, "temp/00c-pairs_ID.csv"), "\tcreated")


