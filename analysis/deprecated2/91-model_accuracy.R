#' This scripts aggregates the random forest model prediction value

library(tidyverse)
library(cowplot)
source(here::here("analysis/00-metadata.R"))

# 1. Read accuracy data ----
## Pairs mapping file
list_image_mapping_folder_master <- read_csv(paste0(folder_script, folder_mapping_files, "00-list_image_mapping_folder_master.csv"), show_col_types = F)

## pairs of mismatch
pairs_mismatch <- read_csv(paste0(folder_data, "temp/15-pairs_mismatch.csv"), show_col_types = F)

## Same-bug pairs
pairs_samebug <- pairs_mismatch %>%
    filter(Mismatch == 0) %>%
    mutate(PairSameBug = T)
isolates_duplicate <- tibble(ID = c(462, 355, 356, 461, 452, 446, 305, 435, 444, 348, 460, 454), Duplicated = T)

## Random forest accuracy
temp <- rep(list(NA), length(batch_names))
for (j in 1:length(batch_names)) {
    list_images <- read_csv(paste0(folder_script, folder_mapping_files, "00-list_images-", batch_names[j], "-green.csv"), show_col_types = F)
    list_image_mapping <- read_csv(paste0(folder_script, folder_mapping_files, "00-list_image_mapping-", batch_names[j], ".csv") , show_col_types = F)

    list_image_mapping_folder <- list_image_mapping %>%
        left_join(rename(list_images, image_name_pair = image_name), by = "image_name_pair") %>%
        left_join(select(list_images, image_name_isolate1 = image_name), by = "image_name_isolate1") %>%
        left_join(select(list_images, image_name_isolate2 = image_name), by = "image_name_isolate2")

    accuracy_validation <- NULL
    # Read the prediction
    for (i in 1:nrow(list_image_mapping_folder)) {
        # Skip no colony results
        if (list_image_mapping_folder$image_name_pair[i] %in% plates_no_colony) next
        accuracy_validation[[i]] <- read_csv(paste0(list_image_mapping_folder$folder_random_forest[i], "validation-", list_image_mapping_folder$image_name_pair[i], ".csv"), show_col_types = F) %>%
            mutate(image_name_pair = list_image_mapping_folder$image_name_pair[i])
        cat(" ", i)
    }

    temp[[j]] <- bind_rows(accuracy_validation) %>%
        #rename(image_name_pair = image_name) %>%
        left_join(list_image_mapping_folder, by = "image_name_pair") %>%
        select(-starts_with("folder_"))
}

accuracy_validation_batch <- bind_rows(temp)


# 2. Clean up pairs of contamination/ duplication ----
## Append ID of Duplicated isolates to my internal ID
isolates_ID <- read_csv(paste0(folder_data, "temp/00c-isolates_ID.csv"), show_col_types = F) %>%
    select(ID, Community, Isolate) %>%
    left_join(isolates_duplicate, by = "ID") %>%
    replace_na(list(Duplicated = F))

## Append the column and remove duplicated isolates in the accuracy table
accuracy <- accuracy_validation_batch %>%
    filter(FinalModel) %>%
    # Remove C11R2 isolate 13, which is a Staph. It's also not included in isolates_ID_match
    filter(!((Community == "C11R2" & Isolate1 == 13) | ((Community == "C11R2" & Isolate2 == 13)))) %>%
    # Remove batch B2 C11R1 pairs that contains isolate 1
    filter(!(Batch == "B2" & Community == "C11R1" & (Isolate1 == 1 | Isolate2 == 1))) %>%
    # Keep batch C C11R1 pairs that contains isolate 1
    filter(!(Batch == "C" & Community == "C11R1" & (Isolate1 != 1 & Isolate2 != 1))) %>%
    # Label pairs containing duplicates
    left_join(rename(isolates_ID, ID1 = ID, Isolate1 = Isolate, Duplicated1 = Duplicated), by = c("Community", "Isolate1")) %>%
    left_join(rename(isolates_ID, ID2 = ID, Isolate2 = Isolate, Duplicated2 = Duplicated), by = c("Community", "Isolate2")) %>%
    mutate(PairType = case_when(
        Duplicated1 == F & Duplicated2 == F ~ "clean",
        Duplicated1 == T & Duplicated2 == F ~ "one duplicate",
        Duplicated1 == F & Duplicated2 == T ~ "one duplicate",
        Duplicated1 == T & Duplicated2 == T ~ "both duplicate"
    )) %>%
    mutate(PairType = factor(PairType, c("clean", "one duplicate", "both duplicate"))) %>%
    mutate(AccuracyPassThreshold = case_when(
        Accuracy >= 0.9 ~ T,
        Accuracy < 0.9 ~ F
    )) %>%
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
    select(-FlipOrder) %>%
    # Remove staph contamination
    filter(!(Community == "C11R2" & Isolate1 == 13)) %>%
    filter(!(Community == "C11R2" & Isolate2 == 13)) %>%
    filter(!(Batch == "B2" & Community == "C11R1" & Isolate1 == 1)) %>%
    filter(!(Batch == "C" & Community == "C11R1" & Isolate1 == 5)) %>%
    select(image_name_pair, Community, Isolate1, Isolate2, Isolate1InitialODFreq,
           Accuracy, AccuracyPassThreshold, PairType)




# 3. Model accuracy ----
accuracy_count <- accuracy %>%
    group_by(PairType) %>%
    count(name = "Count")

p1a <- accuracy %>%
    ggplot() +
    geom_histogram(aes(x = Accuracy), color = 1, binwidth = 0.01, breaks = seq(0.6,1,.01)) +
    geom_text(data = accuracy_count, x = -Inf, y = Inf, aes(label = paste0("N=",Count)), vjust = 2, hjust = -1) +
    geom_vline(xintercept = 0.9, color = "red", linetype = 2) +
    facet_grid(PairType~.) +
    theme_classic() +
    theme(panel.border = element_rect(color = 1, fill = NA)) +
    labs(x = "Accuracy", y = "Count")

p1b <- accuracy %>%
    ggplot() +
    geom_histogram(aes(x = Accuracy), color = 1, binwidth = 0.01, breaks = seq(0.6,1,.01)) +
    geom_text(data = accuracy_count, x = -Inf, y = Inf, aes(label = paste0("N=",Count)), vjust = 2, hjust = -1) +
    geom_vline(xintercept = 0.9, color = "red", linetype = 2) +
    scale_y_log10() +
    facet_grid(PairType~.) +
    theme_classic() +
    theme(panel.border = element_rect(color = 1, fill = NA)) +
    labs(x = "Accuracy", y = "Count")

p1 <- plot_grid(p1a, p1b, nrow = 1, axis = "tb", align = "h", scale = 0.9, labels = c("count", "log-scale")) + theme(plot.background = element_rect(fill = "white", color = NA))
ggsave(paste0(folder_data, "temp/91-random_forest-accuracy.png"), p1, width = 6, height = 6)


## How many of the clean pairs have accuracy lower than 0.9?
accuracy %>%
    group_by(PairType, AccuracyPassThreshold) %>%
    count(name = "Count") %>%
    group_by(PairType) %>%
    mutate(Fraction = Count / sum(Count))

#   PairType       AccuracyPassThreshold Count Fraction
#   <fct>          <lgl>                 <int>    <dbl>
# 1 clean          TRUE                    357    1
# 2 one duplicate  FALSE                    24    0.140
# 3 one duplicate  TRUE                    147    0.860
# 4 both duplicate FALSE                     3    0.143
# 5 both duplicate TRUE                     18    0.857


pairs_accuracy <- accuracy %>%
    # Remove pairs that have cocultures with no colony
    unite(col = "Pair", Community, Isolate1, Isolate2, sep = "_", remove = F) %>%
    filter(!(Pair %in% pairs_no_colony)) %>%
    group_by(Community, Isolate1, Isolate2) %>%
    summarize(AccuracyMean = mean(Accuracy), AccuracySd = sd(Accuracy))

write_csv(accuracy, paste0(folder_data, "temp/91-accuracy.csv"))
write_csv(pairs_accuracy, paste0(folder_data, "temp/91-pairs_accuracy.csv"))













