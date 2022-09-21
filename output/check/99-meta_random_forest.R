#' Meta analysis for the output generated in 09-random_forest
library(tidyverse)
library(cowplot)
folder_script <- "~/Desktop/Lab/emergent-coexistence/output/check/"
folder_main <- "~/Dropbox/lab/emergent-coexistence/data/raw/plate_scan/emergent_coexistence_plate_scan_check/"
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

# 1. Read data ----
batch_names <- c("D", "C2", "B2")
#batch_names <- c("D", "C2", "B2", "C")

## Pairs mapping file
list_image_mapping_master <- rep(list(NA), length(batch_names))
for (j in 1:length(batch_names)) list_image_mapping_master[[j]] <- read_csv(paste0(folder_script, "00-list_image_mapping-", batch_names[j], ".csv") , show_col_types = F)
list_image_mapping_master <- bind_rows(list_image_mapping_master)
list_image_mapping_master %>% left_join(tibble(image_name_pair = plates_no_colony, Undecided = "no colony"))


## pairs of mismatch
pairs_mismatch <- read_csv("~/Desktop/Lab/emergent-coexistence/output/check/pairs_mismatch.csv", show_col_types = F)

## Same-bug pairs
pairs_samebug <- pairs_mismatch %>%
    filter(Mismatch == 0) %>%
    mutate(PairSameBug = T)
isolates_duplicate <- tibble(ID = c(462, 355, 356, 461, 452, 446, 305, 435, 444, 348, 460, 454), Remove = T)
pairs_samebug_full <- pairs_samebug %>%
    mutate(revIsolate1 = Isolate2, revIsolate2 = Isolate1) %>%
    pivot_longer(cols = ends_with("Isolate1"), names_to = "temp", values_to = "Isolate1") %>%
    select(-temp) %>%
    pivot_longer(cols = ends_with("Isolate2"), names_to = "temp", values_to = "Isolate2") %>%
    select(-temp) %>%
    filter(Isolate1 != Isolate2)


## Pairs of challenging morphology
# C1R4	2,4	2,5	4,5
# C2R6 	1,3
# C2R8	3,4
# C4R1	1,3
# C7R1	1,3	1,4	3,4
# C8R4	1,3
# C11R1	1,3 	1,4	2,8	2,9	8,9
# C11R2	2,8	3,5	3,7	3,12	5,6	7,12	8,9
# C11R5	1,3	1,5	2,4	3,5
pairs_similar <- tibble(
    Community = c(rep("C1R4", 3),
                  "C2R6", "C2R8", "C4R1",
                  rep("C7R1", 3),
                  "C8R4",
                  rep("C11R1", 5),
                  rep("C11R2", 7),
                  rep("C11R5", 4)),
    Isolate1 = c(2, 2, 4,
                 1, 3, 1,
                 1, 1, 3,
                 1,
                 1,1,2,2,8,
                 2,3,3,3,5,7,8,
                 1,1,2,3),
    Isolate2 = c(4, 5, 5,
                 3, 4, 3,
                 3, 4, 4,
                 3,
                 3,4,8,9,9,
                 8,5,7,12,6,12,9,
                 3,5,4,5),
    SimilarByEye = T
)



## Object prediction
temp <- rep(list(NA), length(batch_names))

for (j in 1:length(batch_names)) {
    "
    # Have to edit this to include rgb
    "
    list_images <- read_csv(paste0(folder_script, "00-list_images-", batch_names[j], "-green.csv"), show_col_types = F)
    list_image_mapping <- read_csv(paste0(folder_script, "00-list_image_mapping-", batch_names[j], ".csv") , show_col_types = F)

    list_image_mapping_folder <- list_image_mapping %>%
        left_join(rename(list_images, image_name_pair = image_name), by = "image_name_pair") %>%
        left_join(select(list_images, image_name_isolate1 = image_name), by = "image_name_isolate1") %>%
        left_join(select(list_images, image_name_isolate2 = image_name), by = "image_name_isolate2")

    object_prediction <- NULL
    # Read the prediction
    for (i in 1:nrow(list_image_mapping_folder)) {
        # Skip no colony results
        if (list_image_mapping_folder$image_name_pair[i] %in% plates_no_colony) next
        object_prediction[[i]] <- read_csv(paste0(list_image_mapping_folder$folder_random_forest[i], list_image_mapping_folder$image_name_pair[i], ".csv"), show_col_types = F)
        cat(" ", i)
    }

    temp[[j]] <- bind_rows(object_prediction) %>%
        rename(image_name_pair = image_name) %>%
        left_join(list_image_mapping_folder, by = "image_name_pair") %>%
        select(-starts_with("folder_"))
}

object_prediction_batch <- bind_rows(temp)


## Random forest accuracy
for (j in 1:length(batch_names)) {
    list_images <- read_csv(paste0(folder_script, "00-list_images-", batch_names[j], "-green.csv"), show_col_types = F)
    list_image_mapping <- read_csv(paste0(folder_script, "00-list_image_mapping-", batch_names[j], ".csv") , show_col_types = F)

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

# 2. Remove duplicated isolates that have 0 mismatch and remove Staph ----

## Append ID of Duplicated isolates to my internal ID
isolates_ID_match <- read_csv("~/Dropbox/lab/emergent-coexistence/data/raw/pairwise_competition/isolates1.csv", col_types = cols()) %>%
    mutate(Assembly = "self_assembly") %>%
    select(ID, Community, Isolate) %>%
    left_join(isolates_duplicate, by = "ID") %>%
    replace_na(list(Remove = F))


## Append the column and remove duplicated isolates in the prediction tibble
object <- object_prediction_batch %>%
    # Remove C11R2 isolate 13, which is a Staph. It's also not included in isolates_ID_match
    filter(!((Community == "C11R2" & Isolate1 == 13) | ((Community == "C11R2" & Isolate2 == 13))))

object_prediction_batch %>%
    filter(((Community == "C11R2" & Isolate1 == 13) | ((Community == "C11R2" & Isolate2 == 13)))) %>%
     distinct(Community, Isolate1, Isolate2, Freq1) %>% nrow() # Number of coculture images containing C11R2 isolate 13

object_clean <- object %>%
    left_join(rename(isolates_ID_match, ID1 = ID, Isolate1 = Isolate, Remove1 = Remove), by = c("Community", "Isolate1")) %>%
    left_join(rename(isolates_ID_match, ID2 = ID, Isolate2 = Isolate, Remove2 = Remove), by = c("Community", "Isolate2")) %>%
    filter(Remove1 == F & Remove2 == F)
object_clean %>% distinct(Community, Isolate1, Isolate2, Freq1) %>% nrow() # Number of clean images containing no duplicate isolate

object_duplicate <- object %>%
    left_join(rename(isolates_ID_match, ID1 = ID, Isolate1 = Isolate, Remove1 = Remove), by = c("Community", "Isolate1")) %>%
    left_join(rename(isolates_ID_match, ID2 = ID, Isolate2 = Isolate, Remove2 = Remove), by = c("Community", "Isolate2")) %>%
    filter(Remove1 == T | Remove2 == T)
object_duplicate %>% distinct(Community, Isolate1, Isolate2, Freq1) %>% nrow() # Number of coculture images containting one or both duplicated isolates

## Append the column and remove duplicated isolates in the accuracy table
accuracy <- accuracy_validation_batch %>%
    # Remove C11R2 isolate 13, which is a Staph. It's also not included in isolates_ID_match
    filter(!((Community == "C11R2" & Isolate1 == 13) | ((Community == "C11R2" & Isolate2 == 13))))

accuracy_clean <- accuracy %>%
    left_join(rename(isolates_ID_match, ID1 = ID, Isolate1 = Isolate, Remove1 = Remove), by = c("Community", "Isolate1")) %>%
    left_join(rename(isolates_ID_match, ID2 = ID, Isolate2 = Isolate, Remove2 = Remove), by = c("Community", "Isolate2")) %>%
    filter(Remove1 == F & Remove2 == F)
accuracy_clean %>% distinct(Community, Isolate1, Isolate2, Freq1) %>% nrow() # Number of clean images containing no duplicate isolate

accuracy_duplicate <- accuracy %>%
    left_join(rename(isolates_ID_match, ID1 = ID, Isolate1 = Isolate, Remove1 = Remove), by = c("Community", "Isolate1")) %>%
    left_join(rename(isolates_ID_match, ID2 = ID, Isolate2 = Isolate, Remove2 = Remove), by = c("Community", "Isolate2")) %>%
    filter(Remove1 == T | Remove2 == T)
accuracy_duplicate %>% distinct(Community, Isolate1, Isolate2, Freq1) %>% nrow() # Number of coculture images containting one or both duplicated isolates

accuracy_train <- bind_rows(
    accuracy_clean %>% filter(FinalModel) %>% mutate(CleanPair = T),
    accuracy_duplicate %>% filter(FinalModel) %>% mutate(CleanPair = F)
)


# 3. Model accuracy
accuracy_train_count <- accuracy_train %>%
    mutate(CleanPair = factor(CleanPair, c(T, F))) %>%
    group_by(CleanPair) %>%
    count(name = "Count")

p1a <- accuracy_train %>%
    mutate(CleanPair = factor(CleanPair, c(T, F))) %>%
    ggplot() +
    geom_histogram(aes(x = Accuracy), color = 1, binwidth = 0.01, breaks = seq(0.6,1,.01)) +
    geom_text(data = accuracy_train_count, x = -Inf, y = Inf, aes(label = paste0("N=",Count)), vjust = 2, hjust = -1) +
    geom_vline(xintercept = 0.9, color = "red", linetype = 2) +
    facet_grid(CleanPair~., labeller = labeller(CleanPair = c("TRUE" = "clean pair", "FALSE" = "duplicated isolates"))) +
    theme_classic() +
    theme(panel.border = element_rect(color = 1, fill = NA)) +
    labs(x = "Accuracy", y = "Count")

p1b <- accuracy_train %>%
    mutate(CleanPair = factor(CleanPair, c(T, F))) %>%
    ggplot() +
    geom_histogram(aes(x = Accuracy), color = 1, binwidth = 0.01, breaks = seq(0.6,1,.01)) +
    geom_text(data = accuracy_train_count, x = -Inf, y = Inf, aes(label = paste0("N=",Count)), vjust = 2, hjust = -1) +
    geom_vline(xintercept = 0.9, color = "red", linetype = 2) +
    scale_y_log10() +
    facet_grid(CleanPair~., labeller = labeller(CleanPair = c("TRUE" = "clean pair", "FALSE" = "duplicated isolates"))) +
    theme_classic() +
    theme(panel.border = element_rect(color = 1, fill = NA)) +
    labs(x = "Accuracy", y = "Count")

p1 <- plot_grid(p1a, p1b, nrow = 1, axis = "tb", align = "h", scale = 0.9, labels = c("count", "log-scale")) + theme(plot.background = element_rect(fill = "white", color = NA))
ggsave(paste0(folder_main, "examination/random_forest-accuracy.png"), p1, width = 10, height = 6)


## How many of the clean pairs have accurracy lower than 0.9?
accuracy_train %>%
    #filter(CleanPair) %>%
    mutate(AccuracyPassThreshold = case_when(
        Accuracy >= 0.9 ~ T,
        Accuracy < 0.9 ~ F
    )) %>%
    group_by(CleanPair, AccuracyPassThreshold) %>%
    count(name = "Count") %>%
    group_by(CleanPair) %>%
    mutate(Fraction = Count / sum(Count))

# Groups:   CleanPair [2]
#   CleanPair AccuracyPassThreshold Count Fraction
#   <lgl>     <lgl>                 <int>    <dbl>
# 1 FALSE     FALSE                    28   0.145
# 2 FALSE     TRUE                    165   0.855
# 3 TRUE      FALSE                     6   0.0168
# 4 TRUE      TRUE                    352   0.983


## Pairs that have accuracy < 0.9
temp <- accuracy_train %>%
    filter(CleanPair) %>%
    filter(Accuracy < 0.9) %>%
    select(image_name_pair, Batch, Community, Isolate1, Isolate2, Freq1, Freq2, Accuracy)


# 4. Prediction overview ----


# Object prediction counts
object_count <- object_clean %>%
    mutate(Group = factor(Group)) %>%
    group_by(image_name_pair, Community, Isolate1, Isolate2, Group, .drop = F) %>%
    count(name = "Count") %>%
    ungroup()


# Object prediction counts, pairs ordered by the fraction of undecided colonies
object_count_ordered <- object_count %>%
    group_by(image_name_pair) %>%
    mutate(Fraction = Count / sum(Count)) %>%
    pivot_wider(id_cols = c(image_name_pair, Community, Isolate1, Isolate2), names_from = Group, values_from = Fraction) %>%
    arrange(`predicted isolate1`, `predicted isolate2`)  %>%
    ungroup() %>%
    # Assign sorting ID to the pairs
    mutate(PairOrder = 1:n()) %>%
    mutate(image_name_pair = ordered(image_name_pair, image_name_pair)) %>%
    pivot_longer(cols = c(`predicted isolate1`, `predicted isolate2`), names_to = "Group", values_to = "Fraction")


# Plot the result
fill_names <- c("#FF5A5F","#087E8B", "black") %>% setNames(c("predicted isolate1", "predicted isolate2", "undecided"))
p1 <- object_count_ordered %>%
    mutate(AccuracyPassThreshold = ifelse(image_name_pair %in% temp$image_name_pair, "low accuracy", "high accuracy")) %>%
    ggplot() +
    geom_col(aes(x = image_name_pair, y = Fraction, fill = Group, color = AccuracyPassThreshold), alpha = .5) +
    geom_hline(yintercept = 0.05, linetype = 3, color = 1) +
    scale_fill_manual(values = fill_names) +
    scale_color_manual(values = c("low accuracy" = "black", "high accuracy" = NA)) +
    scale_y_continuous(expand = c(0,0)) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
ggsave(paste0(folder_main, "examination/random_forest-prediction.png"), p1, width = 40, height = 10)


# number of sorted out / low undecided / high undecided
object_count_undecided <- object_count_ordered %>%
    pivot_wider(names_from = Group, values_from = Fraction) %>%
    mutate(Undecided = case_when(
        undecided == 0 ~ "sorted out",
        undecided < 0.05 ~ "low undecided",
        undecided >= 0.05 ~ "high undecided"
    ))

object_count_undecided %>%
    group_by(Undecided) %>%
    count(name = "Count") %>%
    ungroup() %>%
    mutate(Fraction = Count / sum(Count))































