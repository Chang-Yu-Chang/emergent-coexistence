#' Meta analysis for the output generated in 08-logit
library(tidyverse)
library(cowplot)
folder_script <- "~/Desktop/lab/emergent-coexistence/output/check/"
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
    list_images <- read_csv(paste0(folder_script, "00-list_images-", batch_names[j], ".csv"), show_col_types = F)
    list_image_mapping <- read_csv(paste0(folder_script, "00-list_image_mapping-", batch_names[j], ".csv") , show_col_types = F)

    list_image_mapping_folder <- list_image_mapping %>%
        left_join(select(list_images, image_name_pair = image_name, folder_feature_pair = folder_green_feature, folder_transection_pair = folder_green_transection, folder_green_cluster), by = "image_name_pair") %>%
        left_join(select(list_images, image_name_isolate1 = image_name, folder_feature_isolate1 = folder_green_feature, folder_transection_isolate1 = folder_green_transection), by = "image_name_isolate1") %>%
        left_join(select(list_images, image_name_isolate2 = image_name, folder_feature_isolate2 = folder_green_feature, folder_transection_isolate2 = folder_green_transection), by = "image_name_isolate2")

    object_prediction <- NULL
    # Read the prediction
    for (i in 1:nrow(list_image_mapping_folder)) {
        # Skip no colony results
        if (list_image_mapping_folder$image_name_pair[i] %in% plates_no_colony) next
        object_prediction[[i]] <- read_csv(paste0(list_image_mapping_folder$folder_green_cluster[i], list_image_mapping_folder$image_name_pair[i], ".csv"), show_col_types = F)
        cat(" ", i)
    }

    temp[[j]] <- bind_rows(object_prediction) %>%
        rename(image_name_pair = image_name) %>%
        left_join(list_image_mapping_folder, by = "image_name_pair") %>%
        select(-starts_with("folder_"))
}

object_prediction_batch <- bind_rows(temp)

# 2. Remove duplicated isolates that have 0 mismatch and remove Staph ----

## Append ID of Duplicated isolates to my internal ID
isolates_ID_match <- read_csv("~/Dropbox/lab/emergent-coexistence/data/raw/pairwise_competition/isolates1.csv", col_types = cols()) %>%
    mutate(Assembly = "self_assembly") %>%
    select(ID, Community, Isolate) %>%
    left_join(isolates_duplicate, by = "ID") %>%
    replace_na(list(Remove = F))

## Append the column and remove duplicated isolates
object <- object_prediction_batch %>%
    # Remove C11R2 isolate 13, which is a Staph. It's also not included in isolates_ID_match
    filter(!((Community == "C11R2" & Isolate1 == 13) | ((Community == "C11R2" & Isolate2 == 13))))

object_prediction_batch %>%
    filter(((Community == "C11R2" & Isolate1 == 13) | ((Community == "C11R2" & Isolate2 == 13)))) %>%
     distinct(Community, Isolate1, Isolate2, Freq1) %>% nrow() # Number of coculture images containing C11R2 isolate 13


## Remove duplicated isolates
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

# 3. Prediction overview ----

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
    arrange(undecided, `predicted isolate1`, `predicted isolate2`)  %>%
    ungroup() %>%
    # Assign sorting ID to the pairs
    mutate(PairOrder = 1:n()) %>%
    mutate(image_name_pair = ordered(image_name_pair, image_name_pair)) %>%
    pivot_longer(cols = c(undecided, `predicted isolate1`, `predicted isolate2`), names_to = "Group", values_to = "Fraction")


# Plot the result
fill_names <- c("#FF5A5F","#087E8B", "black") %>% setNames(c("predicted isolate1", "predicted isolate2", "undecided"))
p1 <- object_count_ordered %>%
    ggplot() +
    geom_col(aes(x = image_name_pair, y = Fraction, fill = Group), alpha = .5) +
    geom_hline(yintercept = 0.05, linetype = 3, color = 1) +
    scale_fill_manual(values = fill_names) +
    scale_y_continuous(expand = c(0,0)) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
ggsave(paste0(folder_main, "examination/prediction.png"), p1, width = 40, height = 10)


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

# 20220909
#   Undecided Count Fraction
#   <chr>           <int>    <dbl>
# 1 high undecided     43    0.239
# 2 low undecided      32    0.178
# 3 sorted out        105    0.583


# 20220912 including threee batches, only linear models
#   Undecided Count Fraction
#   <chr>           <int>    <dbl>
# 1 high undecided    166    0.270
# 2 low undecided      91    0.148
# 3 sorted out        357    0.581

# 20220912 including threee batches, including interaction terms
#   Undecided Count Fraction
#   Undecided      Count Fraction
#   <chr>          <int>    <dbl>
# 1 high undecided   129    0.234
# 2 low undecided     84    0.152
# 3 sorted out       338    0.613

# 20220914 including threee batches, including interaction terms. Remove the duplicated isolates
# So the total number of pairs are 358
#   Undecided      Count Fraction
#   <chr>          <int>    <dbl>
# 1 high undecided    72    0.201
# 2 low undecided     54    0.151
# 3 sorted out       232    0.648


"
among 358 clean images, how many of them are sorted / low undecided /high undecided?

among these images, how many unique pairs are there?
how many of the unique pairs are all sorted? how many of them have 3 pairs highly undecided?
"

# Undecided group by unique species pairs

## Clean up the column name and freq order
object_unique_pairs <- object_count_undecided %>%
    # Append the frequency data
    left_join(list_image_mapping_master, by = c("image_name_pair", "Community", "Isolate1", "Isolate2")) %>%
    # Unique species pair
    rowwise() %>%
    mutate(temp_Isolate1 = min(Isolate1, Isolate2), temp_Isolate2 = max(Isolate1, Isolate2)) %>%
    mutate(temp_Freq1 = ifelse(Isolate1 < Isolate2, Freq1, Freq2), temp_Freq2 = ifelse(Isolate1 < Isolate2, Freq2, Freq1)) %>%
    unite(col = "pair_name", Batch, Community, temp_Isolate1, temp_Isolate2, remove = F) %>%
    # Clean up columns
    select(pair_name, Community, starts_with("temp"), undecided, `predicted isolate1`, `predicted isolate2`, Undecided) %>%
    rename_with(~str_replace(., "temp_", ""), starts_with("temp")) %>%
    mutate(Community = factor(Community, paste0("C", rep(1:12, each = 8), "R", 1:8))) %>%
    arrange(Community, pair_name, Freq1)


object_unique_pairs_undecied <- object_unique_pairs %>%
    group_by(pair_name, Undecided) %>%
    count(name = "Count") %>%
    # Sort the pairs
    pivot_wider(names_from = Undecided, values_from = Count, values_fill = 0) %>%
    arrange(`high undecided`, `low undecided`, `sorted out`) %>%
    ungroup() %>%
    mutate(PairOrder = 1:n(), pair_name = factor(pair_name, pair_name)) %>%
    pivot_longer(cols = c("sorted out", "low undecided", "high undecided"), names_to = "Undecided", values_to = "Count") %>%
    mutate(Undecided = factor(Undecided, c("sorted out", "low undecided", "high undecided")))


fill_names <- c("#EFF1ED", "#717744", "#373D20") %>% setNames(c("sorted out", "low undecided", "high undecided"))
p3 <- object_unique_pairs_undecied %>%
    ggplot() +
    geom_col(aes(x = pair_name, y = Count, fill = Undecided)) +
    scale_fill_manual(values = fill_names) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5))

ggsave(paste0(folder_main, "examination/pairs_undecided.png"), p3, width = 20, height = 5)


# Groups of
object_unique_pairs_undecied_group <- object_unique_pairs_undecied %>%
    pivot_wider(names_from = Undecided, values_from = Count) %>%
    mutate(UniquePairGroup = case_when(
        `sorted out` == 3 | (`sorted out` == 2 & `low undecided` == 0 & `high undecided` == 0) ~ "g1",
        (`sorted out` == 1 | `sorted out` == 2) & `low undecided` == 1 & `high undecided` == 0 ~ "g2",
         `low undecided` >= 2 | `high undecided` <= 2 ~ "g3",
        `sorted out` == 0 & `low undecided` == 0 & `high undecided` == 3 ~ "g4",
    ))

object_unique_pairs_undecied_group %>%
    group_by(UniquePairGroup) %>%
    count(name = "Count") %>%
    ungroup() %>%
    mutate(Fraction = Count / sum(Count))

# 20220914
#   UniquePairGroup Count Fraction
#   <chr>           <int>    <dbl>
# 1 g1                 55    0.458
# 2 g2                 26    0.217
# 3 g3                 22    0.183
# 4 g4                 17    0.142

object_unique_pairs_undecied_group_by_eye <-  object_unique_pairs_undecied_group %>%
    separate(pair_name, into = c("Batch", "Community", "Isolate1", "Isolate2"), remove = F, convert = T) %>%
    left_join(pairs_similar, by = c("Community", "Isolate1", "Isolate2")) %>%
    replace_na(list(SimilarByEye = F)) %>%
    # Sort the pairs
    arrange(desc(SimilarByEye), `high undecided`, `low undecided`, `sorted out`) %>%
    ungroup() %>%
    mutate(PairOrder = 1:n(), pair_name = factor(pair_name, pair_name)) %>%
    pivot_longer(cols = c("sorted out", "low undecided", "high undecided"), names_to = "Undecided", values_to = "Count")

p4 <- object_unique_pairs_undecied_group_by_eye %>%
    filter(UniquePairGroup %in% c("g3", "g4")) %>%
    ggplot() +
    geom_col(aes(x = pair_name, y = Count, fill = Undecided, color = SimilarByEye)) +
    scale_fill_manual(values = fill_names) +
    scale_color_manual(values = c("TRUE" = 1, "FALSE" = "white")) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
ggsave(paste0(folder_main, "examination/pairs_undecided_byeye.png"), p4, width = 15, height = 5)


#
object_unique_pairs_undecied_group_by_eye %>%
    filter(UniquePairGroup %in% c("g3", "g4")) %>%
    pivot_wider(names_from = Undecided, values_from = Count) %>%
    group_by(SimilarByEye) %>%
    count(name = "Count") %>%
    ungroup() %>%
    mutate(Fraction = Count / sum(Count))

#   SimilarByEye Count Fraction
#   <lgl>        <int>    <dbl>
# 1 FALSE           34    0.872
# 2 TRUE             5    0.128


object_unique_pairs_undecied_group_by_eye %>%
    filter(UniquePairGroup %in% c("g4")) %>%
    pivot_wider(names_from = Undecided, values_from = Count) %>%
    filter(SimilarByEye == F) %>%
    pull(pair_name) %>%
    as.character() %>%
    paste(collapse = "\n") %>%
    cat()




# # Remove the same-bug pairs ----
# ## Number of non-same-bug pairs
# object_prediction_count %>%
#     filter(PairSameBug == F) %>%
#     group_by(Undecided) %>%
#     count()
# #filter(Undecided == "high undecided") %>%
# left_join(pairs_samebug_full, by = c("Community", "Isolate1", "Isolate2")) %>%
#     replace_na(list(PairSameBug = F)) %>%
#     arrange(desc(PairSameBug)) %>%
#     mutate(image_name_pair = factor(image_name_pair, image_name_pair)) %>%
#     #filter(PairSameBug == F) %>%
#     pivot_longer(cols = c(undecided, `predicted isolate1`, `predicted isolate2`), names_to = "Group", values_to = "Fraction")
# p2 <- object_prediction_count_undecided_nonsamebug  %>%
#     ggplot() +
#     geom_col(aes(x = image_name_pair, y = Fraction, fill = Group, color = PairSameBug), alpha = .5) +
#     scale_color_manual(values = c("TRUE" = "black", "FALSE" = NA)) +
#     scale_fill_manual(values = fill_names) +
#     scale_y_continuous(expand = c(0,0)) +
#     theme_classic() +
#     theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
#
#
#
# object_prediction_count_undecided_nonsamebug %>%
#     filter(PairSameBug == F)
# pivot_wider(names_from = Group, values_from = Fraction) %>%
#     rowwise() %>%
#     mutate(newIsolate1 = min(Isolate1, Isolate2), newIsolate2 = max(Isolate1, Isolate2)) %>%
#     select(Community, Isolate1 = newIsolate1, Isolate2 = newIsolate2) %>%
#     group_by(Community, Isolate1, Isolate2) %>%
#     arrange(Community, Isolate1, Isolate2) %>%
#     count()
#







# Plot the transects of isolates and pair ----
image_tocheck <- "D_T8_C1R2_5-95_2_3"
#image_tocheck_index = which(list_image_mapping_folder$image_name_pair %in% image_tocheck)
i = image_tocheck_index[1]

for (i in image_tocheck_index) {
    list_image_mapping_folder$image_name_pair[i]
    transection_isolate1 <- read_csv(paste0(list_image_mapping_folder$folder_transection_isolate1[i], list_image_mapping_folder$image_name_isolate1[i], ".csv"), show_col_types = F)
    transection_isolate2 <- read_csv(paste0(list_image_mapping_folder$folder_transection_isolate2[i], list_image_mapping_folder$image_name_isolate2[i], ".csv"), show_col_types = F)
    transection_pair <- read_csv(paste0(list_image_mapping_folder$folder_transection_pair[i], list_image_mapping_folder$image_name_pair[i], ".csv"), show_col_types = F)

    color_names <- c("#FF5A5F","#087E8B", "black") %>% setNames(c("isolate1", "isolate2", "pair"))
    p3 <- bind_rows(
        mutate(transection_isolate1, Group = "isolate1"),
        mutate(transection_isolate2, Group = "isolate2"),
        mutate(transection_pair, Group = "pair")
    ) %>%
        ggplot() +
        geom_line(aes(x = ScaledDistanceToCenter, y = Intensity, color = Group, group = interaction(Group, ObjectID)), lwd = .1) +
        scale_color_manual(values = color_names) +
        facet_wrap(~Group) +
        theme_classic() +
        theme(panel.border = element_rect(color = 1, fill = NA)) +
        ggtitle(list_image_mapping_folder$image_name_pair[i])
    ggsave(paste0(folder_main, "examination/transection-", list_image_mapping_folder$image_name_pair[i],".png"), p3, width = 15, height = 6)
}

