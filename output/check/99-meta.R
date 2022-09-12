#' Meta analysis for the output generated in 09-cluster
library(tidyverse)
library(cowplot)
folder_script <- "~/Desktop/Lab/emergent-coexistence/output/check/"
folder_main <- "~/Dropbox/lab/emergent-coexistence/data/raw/plate_scan/emergent_coexistence_plate_scan_check/"
list_images <- read_csv(paste0(folder_script, "00-list_images-D.csv"), show_col_types = F)
list_image_mapping <- read_csv(paste0(folder_script, "00-list_image_mapping-D.csv") , show_col_types = F)

list_image_mapping_folder <- list_image_mapping %>%
    left_join(select(list_images, image_name_pair = image_name, folder_feature_pair = folder_green_feature, folder_transection_pair = folder_green_transection, folder_green_cluster), by = "image_name_pair") %>%
    left_join(select(list_images, image_name_isolate1 = image_name, folder_feature_isolate1 = folder_green_feature, folder_transection_isolate1 = folder_green_transection), by = "image_name_isolate1") %>%
    left_join(select(list_images, image_name_isolate2 = image_name, folder_feature_isolate2 = folder_green_feature, folder_transection_isolate2 = folder_green_transection), by = "image_name_isolate2")


# Prediction overview ----
i = 1
object_prediction <- NULL
# Read the prediction
for (i in 1:nrow(list_image_mapping_folder[i])) {
    object_prediction[[i]] <- read_csv(paste0(list_image_mapping_folder$folder_green_cluster[i], list_image_mapping_folder$image_name_pair[i], ".csv"), show_col_types = F)
    cat(" ", i)
}


# Object prediction results
object_prediction <- bind_rows(object_prediction) %>%
    rename(image_name_pair = image_name) %>%
    left_join(list_image_mapping_folder, by = "image_name_pair") %>%
    select(-starts_with("folder_"))


# Object prediction counts
object_prediction_count <- object_prediction %>%
    mutate(Group = factor(Group)) %>%
    group_by(image_name_pair, Community, Isolate1, Isolate2, Group, .drop = F) %>%
    count(name = "Count")

# Object prediction counts, pairs ordered by the fraction of undecided colonies
object_prediction_count_ordered <- object_prediction_count %>%
    group_by(image_name_pair) %>%
    mutate(Fraction = Count / sum(Count)) %>%
    pivot_wider(id_cols = c(image_name_pair, Community, Isolate1, Isolate2), names_from = Group, values_from = Fraction) %>%
    arrange(undecided, `predicted isolate1`, `predicted isolate2`)  %>%
    ungroup() %>%
    # Assign sorting ID to the pairs
    mutate(PairOrder = 1:n()) %>%
    mutate(image_name_pair = ordered(image_name_pair, image_name_pair)) %>%
    pivot_longer(cols = c(undecided, `predicted isolate1`, `predicted isolate2`), names_to = "Group", values_to = "Fraction")


#
fill_names <- c("#FF5A5F","#087E8B", "black") %>% setNames(c("predicted isolate1", "predicted isolate2", "undecided"))
if (FALSE) {
    p1 <- object_prediction_count %>%
        ggplot() +
        geom_col(aes(x = image_name_pair, y = Count, fill = Group), alpha = .5) +
        scale_fill_manual(values = fill_names) +
        scale_y_continuous(expand = c(0,0)) +
        theme_classic() +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5))

    p2 <- object_prediction_count %>%
        ggplot() +
        geom_col(aes(x = image_name_pair, y = Count, fill = Group), alpha = .5, position = position_fill()) +
        scale_fill_manual(values = fill_names) +
        scale_y_continuous(expand = c(0,0)) +
        theme_classic() +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5))

}

p1 <- object_prediction_count_ordered %>%
    ggplot() +
    geom_col(aes(x = image_name_pair, y = Fraction, fill = Group), alpha = .5) +
    geom_hline(yintercept = 0.05, linetype = 3, color = 1) +
    scale_fill_manual(values = fill_names) +
    scale_y_continuous(expand = c(0,0)) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
ggsave(paste0(folder_main, "examination/prediction.png"), p1, width = 40, height = 10)


# number of sorted out / low undecided / high undecided
object_prediction_count_undecided <- object_prediction_count_ordered %>%
    pivot_wider(names_from = Group, values_from = Fraction) %>%
    mutate(undecided_group = case_when(
        undecided == 0 ~ "sorted out",
        undecided < 0.05 ~ "low undecided",
        undecided >= 0.05 ~ "high undecided"
    ))

##
object_prediction_count_undecided %>%
    group_by(undecided_group) %>%
    count(name = "Count") %>%
    ungroup() %>%
    mutate(Fraction = Count / sum(Count))

# 20220909
#   undecided_group Count Fraction
#   <chr>           <int>    <dbl>
# 1 high undecided     43    0.239
# 2 low undecided      32    0.178
# 3 sorted out        105    0.583


# 20220911 after adding one more feature
#   undecided_group Count Fraction
#   <chr>           <int>    <dbl>
# 1 high undecided     45    0.25
# 2 low undecided      35    0.194
# 3 sorted out        100    0.556

# Check it the pair is composed of the same bug
pairs_samebug <- tibble(
    Community = c("C1R4", rep("C11R5", 4)),
    Isolate1 = c(4, 1, 1, 3, 2),
    Isolate2 = c(5, 3, 5, 5, 4),
    PairSameBug = T
)
pairs_samebug_full <- pairs_samebug %>%
    mutate(revIsolate1 = Isolate2, revIsolate2 = Isolate1) %>%
    pivot_longer(cols = ends_with("Isolate1"), names_to = "temp", values_to = "Isolate1") %>%
    select(-temp) %>%
    pivot_longer(cols = ends_with("Isolate2"), names_to = "temp", values_to = "Isolate2") %>%
    select(-temp) %>%
    filter(Isolate1 != Isolate2)

object_prediction_count_undecided %>%
    filter(undecided_group == "high undecided") %>%
    rowwise() %>%
    mutate(newIsolate1 = min(Isolate1, Isolate2), newIsolate2 = max(Isolate1, Isolate2)) %>%
    select(Community, Isolate1 = newIsolate1, Isolate2 = newIsolate2) %>%
    left_join(pairs_samebug, by = c("Community", "Isolate1", "Isolate2")) %>%
    replace_na(list(PairSameBug = F)) %>%
    group_by(PairSameBug, .drop = F) %>%
    count(name = "Count") %>%
    ungroup() %>%
    mutate(Fraction = Count / sum(Count))

# 202220909 among the highly undecided pairs
#   PairSameBug Count Fraction
#   <lgl>       <int>    <dbl>
# 1 FALSE          28    0.651
# 2 TRUE           15    0.349


# 20222091 among the highly undecided pairs
#   PairSameBug Count Fraction
#   <lgl>       <int>    <dbl>
# 1 FALSE          32    0.711
# 2 TRUE           13    0.289


# coculture images that are not the same-bug pairs
object_prediction_count_undecided_nonsamebug <- object_prediction_count_undecided %>%
    #filter(undecided_group == "high undecided") %>%
    left_join(pairs_samebug_full, by = c("Community", "Isolate1", "Isolate2")) %>%
    replace_na(list(PairSameBug = F)) %>%
    arrange(desc(PairSameBug)) %>%
    mutate(image_name_pair = factor(image_name_pair, image_name_pair)) %>%
    #filter(PairSameBug == F) %>%
    pivot_longer(cols = c(undecided, `predicted isolate1`, `predicted isolate2`), names_to = "Group", values_to = "Fraction")
p2 <- object_prediction_count_undecided_nonsamebug  %>%
    ggplot() +
    geom_col(aes(x = image_name_pair, y = Fraction, fill = Group, color = PairSameBug), alpha = .5) +
    geom_hline(yintercept = 0.05, linetype = 3, color = 1) +
    scale_color_manual(values = c("TRUE" = "black", "FALSE" = NA)) +
    scale_fill_manual(values = fill_names) +
    scale_y_continuous(expand = c(0,0)) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5))

ggsave(paste0(folder_main, "examination/prediction_undecided.png"), p2, width = 20, height = 8)
#p <- plot_grid(p1, p2, p3, ncol = 1, align = "v", axis = "lr")


# Undecided gorup by unique species pairs
pairs_freq <- object_prediction_count_undecided_nonsamebug %>%
    pivot_wider(names_from = Group, values_from = Fraction) %>%
    select(-undecided, -`predicted isolate1`, -`predicted isolate2`, -PairOrder) %>%
    separate(image_name_pair, into = c("Batch", "temp1", "temp2", "Freq1", "Freq2", "temp3", "temp4")) %>%
    select(-starts_with("temp"))

temp_index <- pairs_freq$Isolate1 > pairs_freq$Isolate2
temp <- pairs_freq$Isolate1[temp_index]
pairs_freq$Isolate1[temp_index] <- pairs_freq$Isolate2[temp_index]
pairs_freq$Isolate2[temp_index] <- temp
temp <- pairs_freq$Freq1[temp_index]
pairs_freq$Freq1[temp_index] <- pairs_freq$Freq2[temp_index]
pairs_freq$Freq2[temp_index] <- temp

fill_names <- c("#EFF1ED", "#717744", "#373D20") %>% setNames(c("sorted out", "low undecided", "high undecided"))
p3 <- pairs_freq %>%
    filter(PairSameBug == F) %>%
    # Clean up names
    select(-Freq2) %>%
    unite(col = "pair_name", Batch, Community, Isolate1, Isolate2) %>%
    group_by(pair_name, undecided_group) %>%
    count(name = "Count") %>%
    # Sort the pairs
    pivot_wider(names_from = undecided_group, values_from = Count, values_fill = 0) %>%
    arrange(`high undecided`, `low undecided`, `sorted out`) %>%
    ungroup() %>%
    mutate(PairOrder = 1:n(), pair_name = factor(pair_name, pair_name)) %>%
    pivot_longer(cols = c("sorted out", "low undecided", "high undecided"), names_to = "undecided_group", values_to = "Count") %>%
    #
    mutate(undecided_group = factor(undecided_group, c("sorted out", "low undecided", "high undecided"))) %>%
    ggplot() +
    geom_col(aes(x = pair_name, y = Count, fill = undecided_group)) +
    scale_fill_manual(values = fill_names) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5))

ggsave(paste0(folder_main, "examination/pairs_undecided.png"), p3, width = 10, height = 3)

# Check the community pairs with high undecided number of coculture images
object_prediction_count_undecided_nonsamebug %>%
    # Toggle this T/F to check the pairs
    filter(PairSameBug ==  F) %>%
    pivot_wider(names_from = Group, values_from = Fraction) %>%
    rowwise() %>%
    mutate(newIsolate1 = min(Isolate1, Isolate2), newIsolate2 = max(Isolate1, Isolate2)) %>%
    select(Community, Isolate1 = newIsolate1, Isolate2 = newIsolate2) %>%
    group_by(Community, Isolate1, Isolate2) %>%
    arrange(Community, Isolate1, Isolate2) %>%
    count()

# Same-bug pairs
#   Community Isolate1 Isolate2     n
#   <chr>        <dbl>    <dbl> <int>
# 1 C11R5            1        3     3
# 2 C11R5            1        5     3
# 3 C11R5            2        4     3
# 4 C11R5            3        5     3
# 5 C1R4             4        5     3
#
# Different-bug pairs
#    Community Isolate1 Isolate2     n
#    <chr>        <dbl>    <dbl> <int>
#  1 C1R2             2        4     2
#  2 C1R4             1        2     1
#  3 C1R4             2        3     1
#  4 C1R6             1        2     2
#  5 C1R6             1        3     2
#  6 C1R6             2        3     3
#  7 C1R6             3        5     3
#  8 C1R7             1        7     3
#  9 C1R7             3        5     3
# 10 C1R7             3        6     2
# 11 C1R7             3        7     1
# 12 C1R7             5        6     2
# 13 C4R1             1        3     3


# Remove the same-bug pairs ----
## Number of non-same-bug pairs
object_prediction_count %>%

    filter(PairSameBug == F) %>%
    group_by(undecided_group) %>%
    count()
#filter(undecided_group == "high undecided") %>%
left_join(pairs_samebug_full, by = c("Community", "Isolate1", "Isolate2")) %>%
    replace_na(list(PairSameBug = F)) %>%
    arrange(desc(PairSameBug)) %>%
    mutate(image_name_pair = factor(image_name_pair, image_name_pair)) %>%
    #filter(PairSameBug == F) %>%
    pivot_longer(cols = c(undecided, `predicted isolate1`, `predicted isolate2`), names_to = "Group", values_to = "Fraction")
p2 <- object_prediction_count_undecided_nonsamebug  %>%
    ggplot() +
    geom_col(aes(x = image_name_pair, y = Fraction, fill = Group, color = PairSameBug), alpha = .5) +
    scale_color_manual(values = c("TRUE" = "black", "FALSE" = NA)) +
    scale_fill_manual(values = fill_names) +
    scale_y_continuous(expand = c(0,0)) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5))



object_prediction_count_undecided_nonsamebug %>%
    filter(PairSameBug == F)
pivot_wider(names_from = Group, values_from = Fraction) %>%
    rowwise() %>%
    mutate(newIsolate1 = min(Isolate1, Isolate2), newIsolate2 = max(Isolate1, Isolate2)) %>%
    select(Community, Isolate1 = newIsolate1, Isolate2 = newIsolate2) %>%
    group_by(Community, Isolate1, Isolate2) %>%
    arrange(Community, Isolate1, Isolate2) %>%
    count()








# Plot the transection of isoaltes and pair ----

draw_transection <- function (transection, smooth = F, span = .75) {
    color_names <- c("#FF5A5F","#087E8B", "black") %>% setNames(c("isolate1", "isolate2", "pair"))

    # Smooth the curve
    if (smooth) {
        loess_custom <- function (formula, data) loess(formula, data, span = span)
        transection <- transection %>%
            nest(data = -ObjectID) %>%
            mutate(mod = map(data, loess_custom, formula = Intensity ~ DistanceToCenter),
                   FittedIntensity = map(mod, `[[`, "fitted")) %>%
            select(-mod) %>%
            unnest(cols = c(data, FittedIntensity))
    }

    transection %>%
        ggplot() +
        geom_line(aes(x = DistanceToCenter, y = FittedIntensity, group = ObjectID)) +
        scale_color_manual(values = color_names) +
        theme_classic() +
        theme(panel.border = element_rect(color = 1, fill = NA))
}
image_tocheck <- "D_T8_C1R6_5-95_3_5"
image_tocheck_index = which(list_image_mapping_folder$image_name_pair %in% image_tocheck)
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
        draw_transection()

        # #filter(Group != "pair") %>%
        # ggplot() +
        # geom_line(aes(x = DistanceToCenter, y = fitted, color = Group, group = interaction(Group, ObjectID))) +
        # scale_color_manual(values = color_names) +
        # facet_wrap(~Group) +
        # theme_classic() +
        # theme(panel.border = element_rect(color = 1, fill = NA)) +
        # ggtitle(list_image_mapping_folder$image_name_pair[i])

    ggsave(paste0(folder_main, "examination/transection-", list_image_mapping_folder$image_name_pair[i],".png"), p3, width = 15, height = 6)
}






