#' This scripts reads the T0 vs. T8 frequencies and determine the competition outcomes
#' for each unique species pair by computing the fitness functions
#' 1. Read the cleaned-up frequencies data from 094
#' 2. Bootstrap frequency changes
#' 3. Map the frequency changes to the fitness function
#' 4. Determine the competition outcomes

library(tidyverse)
library(cowplot)

folder_script <- "~/Desktop/Lab/emergent-coexistence/output/check/"
folder_main <- "~/Dropbox/lab/emergent-coexistence/data/raw/plate_scan/emergent_coexistence_plate_scan_check/"
pairs_freq_ID <- read_csv(paste0(folder_main, "meta/pairs_freq_ID.csv"), show_col_types = F)
pairs_T0_boots <- read_csv(paste0(folder_main, "meta/pairs_T0_boots.csv"), show_col_types = F)
pairs_T8_boots <- read_csv(paste0(folder_main, "meta/pairs_T8_boots.csv"), show_col_types = F)
pairs_no_colony <- c(
    "C11R1_2_8",
    "C11R1_2_9",
    "C11R1_8_9",
    "C11R2_2_10"
)
plates_no_colony <- c(
    "B2_T8_C11R1_5-95_2_8",
    "B2_T8_C11R1_5-95_2_9",
    "B2_T8_C11R1_5-95_8_2",
    "B2_T8_C11R1_5-95_9_8",
    "B2_T8_C11R1_50-50_2_8",
    "B2_T8_C11R1_50-50_2_9",
    "C2_T8_C11R2_50-50_2_10",
    "C2_T8_C11R2_50-50_9_13",
    "C_T8_C11R1_50-50_1_2", # no plate
    "C_T8_C11R1_50-50_1_3" # no plate
)


# 1. Plot the frequency change for each bootstrap for each coculture----
pairs_ID <- distinct(pairs_freq_ID, Batch, Community, Isolate1, Isolate2)
#pairs_ID <- pairs_ID %>% slice(1:2)
#temp <- rep(list(NA), nrow(pairs_ID))

for (i in 1:nrow(pairs_ID)) {
    community <- pairs_ID$Community[i]
    isolate1 <- pairs_ID$Isolate1[i]
    isolate2 <- pairs_ID$Isolate2[i]
    pair_name <- paste0(community, "_", isolate1, "_", isolate2)
    if (pair_name %in% pairs_no_colony) {cat("\nT8 has no colony, skip pair\t", pair_name); next}

    # Increase or decrease sign
    pair_boots <- bind_rows(
        filter(pairs_T0_boots, Community == community, Isolate1 == isolate1, Isolate2 == isolate2),
        filter(pairs_T8_boots, Community == community, Isolate1 == isolate1, Isolate2 == isolate2),
    ) %>%
        select(Community, Isolate1, Isolate2, Isolate1InitialODFreq, Time, BootstrapID, Isolate1CFUFreq) %>%
        pivot_wider(names_from = Time, values_from = Isolate1CFUFreq) %>%
        mutate(FreqChange = T8-T0) %>%
        mutate(FreqChangeSign = case_when(
            FreqChange > 0 ~ "increase",
            FreqChange == 0 ~ "same",
            FreqChange < 0 ~ "decrease",
        )) %>%
        pivot_longer(cols = c(T0, T8), names_to = "Time", values_to = "Isolate1CFUFreq")

    # Table of increase and decrease
    pair_boots_table <- pair_boots %>%
        filter(Time == "T0") %>%
        mutate(FreqChangeSign = factor(FreqChangeSign, c("increase", "same", "decrease"))) %>%
        group_by(Community, Isolate1, Isolate2, Isolate1InitialODFreq, FreqChangeSign, .drop = F) %>%
        count(name = "Count") %>%
        group_by(Community, Isolate1, Isolate2, Isolate1InitialODFreq) %>%
        mutate(Fraction = Count / sum(Count))

    #
    color_names <- c("increase" = "#EF798A", "same" = grey(0.8), "decrease" = "#7D82B8")
    fill_names <- c("increase" = "#EF798A", "same" = grey(0.8), "decrease" = "#7D82B8")

    # T0 vs. T8 interaction plots
    p1 <- pair_boots %>%
        ggplot() +
        #geom_hline(yintercept = c(0, 1), linetype = 1) +
        geom_point(aes(x = Time, y = Isolate1CFUFreq, color = FreqChangeSign), shape = 21) +
        ## Add color on increase or decrease
        geom_line(aes(x = Time, y = Isolate1CFUFreq, color = FreqChangeSign, group = BootstrapID), size = .1) +
        scale_color_manual(values = color_names) +
        scale_y_continuous(limits = c(0,1), expand = c(0, 0)) +
        facet_grid(.~Isolate1InitialODFreq) +
        theme_classic() +
        theme(panel.border = element_rect(color = 1, fill = NA)) +
        ggtitle(pair_name)

    # Barplot
    p2 <- pair_boots_table %>%
        ggplot() +
        geom_hline(yintercept = c(0, 1), linetype = 1) +
        geom_col(aes(x = FreqChangeSign, y = Fraction, fill = FreqChangeSign), color = 1) +
        geom_hline(yintercept = c(0.05, 0.95), linetype = 2, color = "red") +
        scale_y_continuous(expand = c(0, 0)) +
        scale_fill_manual(values = color_names) +
        facet_grid(.~Isolate1InitialODFreq) +
        theme_classic() +
        theme(panel.border = element_rect(color = 1, fill = NA))


    p <- plot_grid(p1, p2, nrow = 2, axis = "lr", align = "v")
    ggsave(paste0(folder_main, "check/D-11-T0_T8_frequencies/", pair_name, ".png"), p, width = 10, height = 5)
    write_csv(pair_boots_table, paste0(folder_main, "check/D-11-T0_T8_frequencies/", pair_name, ".csv"))
    cat("\n", pair_name, "\t", i, "/", nrow(pairs_ID))
}


#pairs_boots_table <- bind_rows(temp)
#write_csv(pairs_boots_table, paste0(folder_main, "meta/pairs_boots_table.csv"))


# 2. Calculate the significance of frequency change ----
temp <- rep(list(NA), nrow(pairs_ID))
for (i in 1:nrow(pairs_ID)) {
    community <- pairs_ID$Community[i]
    isolate1 <- pairs_ID$Isolate1[i]
    isolate2 <- pairs_ID$Isolate2[i]
    pair_name <- paste0(community, "_", isolate1, "_", isolate2)
    if (pair_name %in% pairs_no_colony) {cat("\nT8 has no colony, skip pair\t", pair_name); next}

    #ggsave(paste0(folder_main, "check/D-11-T0_T8_frequencies/", pair_name, ".png"), p, width = 10, height = 5)
    temp[[i]] <- read_csv(paste0(folder_main, "check/D-11-T0_T8_frequencies/", pair_name, ".csv"), show_col_types = F)
}
pairs_boots_table <- bind_rows(temp[which(!is.na(temp))])

pairs_fitness <- pairs_boots_table %>%
    group_by(Community) %>%
    # mutate(FreqChangeSign = factor(FreqChangeSign, c("increase", "decrease", "same"))) %>%
    # group_by(Community, Isolate1, Isolate2, Isolate1InitialODFreq, FreqChangeSign, .drop = F) %>%
    # count(name = "Count") %>%
    mutate(Community = factor(Community, paste0("C", rep(1:12, each = 8), "R", rep(1:8, 12)))) %>%
    arrange(Community, Isolate1, Isolate2, Isolate1InitialODFreq) %>%
    group_by(Community, Isolate1, Isolate2, Isolate1InitialODFreq, .drop = F) %>%
    mutate(Fraction = Count / sum(Count)) %>%
    # C11R1 isolate 6
    filter(!is.na(FreqChangeSign)) %>%
    pivot_wider(id_cols = c(Community, Isolate1, Isolate2, Isolate1InitialODFreq), names_from = FreqChangeSign, values_from = Fraction)

make_interaction_type <- function () {

    interaction_type_three <- tibble(
        FromRare = rep(c(1, -1, 0), each = 9),
        FromMedium = rep(rep(c(1, -1, 0), each = 3), 3),
        FromAbundant = rep(c(1, -1, 0), 9),
        InteractionType = NA,
        InteractionTypeFiner = NA
    )
    interaction_type_two <- tibble(
        FromRare = rep(c(1, -1, 0), each = 3),
        FromMedium = rep(NA, 9),
        FromAbundant = rep(c(1, -1, 0), 3),
        InteractionType = NA,
        InteractionTypeFiner = NA
    )

    interaction_type <- bind_rows(interaction_type_three, interaction_type_two)

    ## Assign interaction types to combinations of frequency changes signs
    interaction_type$InteractionType[c(1, 14, 28, 32, 10, 13, 31)] <- "exclusion"
    interaction_type$InteractionType[c(2, 3, 5, 8, 9, 23, 26, 29, 30, 33, 4, 11, 12, 15, 17, 20, 34, 35)] <- "coexistence"
    interaction_type$InteractionType[c(27, 36)] <- "coexistence"

    ## Assign finer interaction types to combinations of frequency changes signs
    interaction_type$InteractionTypeFiner[c(1, 14, 28, 32)] <- "competitive exclusion"
    interaction_type$InteractionTypeFiner[c(10, 13, 31)] <- "mutual exclusion"
    interaction_type$InteractionTypeFiner[c(2, 3, 5, 8, 9, 23, 26, 29, 30, 33)] <- "stable coexistence"
    interaction_type$InteractionTypeFiner[c(4, 11, 12, 15, 17, 20, 34, 35)] <- "frequency-dependent coexistence"
    interaction_type$InteractionTypeFiner[c(27,36)] <- "neutrality"
    interaction_type <- interaction_type %>%  mutate(FitnessFunction = paste(FromRare, FromMedium, FromAbundant, sep = "_"))
}
interaction_type <- make_interaction_type()

pairs_outcome <- pairs_fitness %>%
    mutate(FitnessChange = case_when(
        increase > 0.95 ~ 1,
        decrease > 0.95 ~ -1,
        increase < 0.95 & decrease < 0.95 ~ 0
    )) %>%
    pivot_wider(id_cols = c(Community, Isolate1, Isolate2), names_from = Isolate1InitialODFreq, values_from = FitnessChange) %>%
    rename(FromRare = `5`, FromMedium = `50`, FromAbundant = `95`) %>%
    left_join(interaction_type) %>%
    select(Community, Isolate1, Isolate2, InteractionType, InteractionTypeFiner, FitnessFunction)

# 3. Plot the competition outcomes ----

assign_interaction_color <- function (level = "simple") {
    if (level == "simple") {
        interaction_type <- c("exclusion", "coexistence")
        interaction_color <- c("#DB7469", "#557BAA")
        names(interaction_color) <- interaction_type
        return(interaction_color)
    }
    if (level == "matrix") {
        interaction_type <- c("exclusion", "coexistence", "exclusion violating rank", "bistability", "neutrality", "self", "undefined")
        interaction_color <- c("#DB7469", "#557BAA", "#8CB369", "#EECF6D", "#8650C4", "black", "grey80")
        names(interaction_color) <- interaction_type
        return(interaction_color)
    }
    if (level == "finer") {
        interaction_type <- c("competitive exclusion", "stable coexistence", "mutual exclusion", "frequency-dependent coexistence", "neutrality", "exclusion violating rank")
        interaction_color <- c("#DB7469", "#557BAA", "#FFBC42", "#B9FAF8", "#8650C4", "#8CB369")
        names(interaction_color) <- interaction_type
        return(interaction_color)
    }
}
fill_names <- assign_interaction_color()

pairs_outcome %>%
    ungroup() %>%
    #unite(col = "FitnessFunction", FromRare, FromMedium, FromAbundant, sep = "_") %>%
    group_by(FitnessFunction, InteractionType) %>%
    count(name = "Count")


# Exclusion vs. coexistence
pairs_outcome %>%
    ggplot() +
    geom_bar(aes(x = Community, fill = InteractionType), position = "fill") +
    scale_fill_manual(values = fill_names) +
    theme_classic()

# Finer scale
pairs_outcome %>%
    ggplot() +
    geom_bar(aes(x = Community, fill = InteractionTypeFiner), position = "fill") +
    scale_fill_manual(values = assign_interaction_color(level = "finer")) +
    theme_classic()


# 4. Compare the pairs_T0_boots result to random forest classified result

# 4 aggregate the predicted result ----
temp <- rep(list(NA), nrow(pairs_freq_ID))
#for (i in 1:nrow(pairs_freq_ID)) {
for (i in 220:nrow(pairs_freq_ID)) {
    # Skip images with no colony
    if (pairs_freq_ID$image_name_pair[i] %in% plates_no_colony) {cat("\nno colony\t", pairs_freq_ID$image_name_pair[i]); next}
    # Skip no-plate coculture
    if (is.na(pairs_freq_ID$image_name_pair[i])) {cat("\nno plate\t", pairs_freq_ID$image_name_pair[i]); next}
    # Read the random forest predicted object
    image_name <- pairs_freq_ID$image_name_pair[i]
    object_feature_predicted <- read_csv(paste0(folder_main, "check/", pairs_freq_ID$Batch[i], "-09-random_forest/", pairs_freq_ID$image_name_pair[i], ".csv"), show_col_types = F)
    cat("\n", nrow(object_feature_predicted), "objects")
    #
    temp[[i]] <- object_feature_predicted %>%
        select(image_name_pair = image_name, ObjectID, Group)
    cat("\t", i, "/", nrow(pairs_freq_ID), "\t", pairs_freq_ID$image_name_pair[i])
}

object_predicted <- bind_rows(temp[which(!is.na(temp))])


pairs_T8 <- object_predicted %>%
    mutate(Group = factor(Group, c("predicted isolate1", "predicted isolate2"))) %>%
    group_by(image_name_pair, Group, .drop = F) %>%
    count(name = "Count") %>%
    group_by(image_name_pair) %>%
    mutate(TotalCount = sum(Count), Frequency = Count/TotalCount) %>%
    filter(Group == "predicted isolate1") %>%
    pivot_wider(names_from = Group, values_from = Count) %>%
    rename(Isolate1Count = `predicted isolate1`, Isolate1CFUFreq = Frequency) %>%
    # Append jrow information
    left_join(pairs_freq_ID, by = "image_name_pair") %>%
    select(image_name = image_name_pair, Batch, Community, Isolate1, Isolate2, Isolate1InitialODFreq, Isolate1Count, TotalCount, Isolate1CFUFreq) %>%
    ungroup()

"
append T0 boots to T8
"




#     ggplot() +
#     geom_histogram(aes(x = Isolate1CFUFreq, fill = Time), color = 1) +
#     facet_grid(Isolate1InitialODFreq ~ Time, scales = "free_y") +
#     coord_flip() +
#     theme_bw() +
#     theme(panel.border = element_rect(color = 1, fill = NA))

# pairs_T0_boots %>%
#     filter(Community == "C1R2", Isolate1 == 2, Isolate2 == 3) %>%
#     #group_by(Isolate1InitialODFreq) %>%
#     ggplot() +
#     geom_histogram(aes(x = Isolate1CFUFreq)) +
#     facet_grid(Isolate1InitialODFreq ~ .) +
#     theme_classic

if (FALSE) {
    # 2. Append T0 bootstraps to T8 bootstraps ----
pairs_fitness_ID <- pairs_ID %>% select(-Batch)
temp <- rep(list(NA), nrow(pairs_fitness_ID))

i = which(pairs_fitness_ID$Community == "C2R8" & pairs_fitness_ID$Isolate1 == 3 & pairs_fitness_ID$Isolate2 == 4 &
          pairs_fitness_ID$Isolate1InitialODFreq == 50)

for (i in 1:nrow(pairs_fitness_ID)) {
    community <- pairs_fitness_ID$Community[i]
    isolate1 <- pairs_fitness_ID$Isolate1[i]
    isolate2 <- pairs_fitness_ID$Isolate2[i]
    isolate1_freq <- pairs_fitness_ID$Isolate1InitialODFreq[i]
    isolate2_freq <- pairs_fitness_ID$Isolate2InitialODFreq[i]

    pairs_T0_boot <- pairs_T0_boots %>%
        filter(Community == community, Isolate1 == isolate1, Isolate2 == isolate2,
               Isolate1InitialODFreq == isolate1_freq, Isolate2InitialODFreq == isolate2_freq)

    pairs_T8_boot <- pairs_T8_boots %>%
        filter(Community == community, Isolate1 == isolate1, Isolate2 == isolate2,
               Isolate1InitialODFreq == isolate1_freq, Isolate2InitialODFreq == isolate2_freq)

    if (nrow(pairs_T8_boot) == 0) {
        cat("\nNo colony on T8 for", paste(community, isolate1, isolate2, isolate1_freq, isolate2_freq, sep = " "))
        next
    }

    temp[[i]] <- bind_rows(pairs_T0_boot, pairs_T8_boot) %>%
        select(Community, Isolate1, Isolate2, Isolate1InitialODFreq, Time, BootstrapID, Isolate1CFUFreq) %>%
        pivot_wider(names_from = Time, values_from = Isolate1CFUFreq) %>%
        mutate(FreqChange = T8-T0) %>%
        mutate(FreqChangeSign = case_when(
            FreqChange > 0 ~ "increase",
            FreqChange == 0 ~ "same",
            FreqChange < 0 ~ "decrease",
        ))
    cat(" ", i)
}

pairs_fitness <- bind_rows(temp[which(!is.na(temp))])




"
check the below T8
C11R1 1 2 50-50
C11R1 1 3 50-50
No colony on T8 for C11R1 2 8 5 95
No colony on T8 for C11R1 2 8 50 50
No colony on T8 for C11R1 2 8 95 5
No colony on T8 for C11R1 2 9 50 50
No colony on T8 for C11R1 2 9 95 5
No colony on T8 for C11R1 8 9 5 95
No colony on T8 for C11R2 2 10 50 50
that should have colonies
"



}
