# This script stores the metadata shared by all scripts

library(tidyverse)

# 1. Directories
# This main folder depends on your home directory and user name
folder_script <- "/Users/cychang/Desktop/lab/emergent-coexistence/analysis/" # Enter the directory of analysis scripts
folder_pipeline <- "/Users/cychang/Dropbox/lab/emergent-coexistence/pipeline/" # Enter the directory of image processing pipeline
folder_data <- "/Users/cychang/Dropbox/lab/emergent-coexistence/data/" # Enter the directory of data
folder_simulation <-  "/Users/cychang/Dropbox/lab/emergent-coexistence/simulation/" # Enter the directory of simulation data

# 2. For image processing pipeline ----
list_folders <- c("01-channel", "02-rolled", "03-threshold", "04-round", "05-watershed", "06-transect", "07-feature", "08-random_forest", "09-bootstrap", "10-images_and_random_forest")
list_channels <- c("red", "green", "blue")
list_pipeline_scripts <- c("01-channel.R", "02-rolling_ball.py", "03-segmentation.R", "04-feature.R", "04a-merge_features.R", "05-random_forest.R")
batch_names <- c("B2", "C", "C2", "D")

# List of no colony cocultures. They are excluded from the loop in image processing pipeline to avoid error
pairs_no_colony <- c(
    "C11R1_2_8",
    "C11R1_2_9",
    "C11R1_8_9",
    "C11R2_2_10",
    "C11R1_1_2",
    "C11R1_1_3"
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

# Random forest
feature_candidates <- c(
    paste0(c(
        "s.area", "s.radius.mean", "s.radius.sd",
        "s.perimeter", "s.radius.mean", "s.radius.sd", "s.radius.min", "s.radius.max",
        "m.cx", "m.cy", "m.majoraxis", "m.eccentricity", "m.theta",
        "b.mean", "b.sd", "b.mad",
        "b.q001", "b.q005", "b.q01", "b.q02", "b.q05", "b.q08", "b.q09", "b.q095", "b.q099",
        "b.tran.mean", "b.tran.sd", "b.tran.mad",
        "b.center", "b.periphery", "b.diff.cp",
        "b.tran.q005", "b.tran.q01", "b.tran.q05", "b.tran.q09", "b.tran.q095"
        #"t.bump.number"
    ), "_green"),
    paste0(c("b.mean", "b.sd", "b.mad"), rep(c("_red", "_blue"), each = 3))
)


# For determining the competition outcomes
interaction_type_finer <- c(
    "competitive exclusion", "stable coexistence",
    "mutual exclusion", "frequency-dependent coexistence",
    "coexistence at 5%", "coexistence at 95%",
    "2-freq neutrality", "3-freq neutrality",
    "unknown"
)


# 3. For plotting ----
assign_interaction_color <- function (level = "simple") {
    if (level == "simple") {
        interaction_type <- c("exclusion", "coexistence", "unknown")
        interaction_color <- c("#DB7469", "#557BAA", grey(0.5))
        names(interaction_color) <- interaction_type
        return(interaction_color)
    }

    if (level == "NA data") {
        interaction_type <- c("exclusion", "coexistence", "unknown", "no colony or low accuracy")
        interaction_color <- c("#DB7469", "#557BAA", grey(0.5), grey(0.8))
        names(interaction_color) <- interaction_type
        return(interaction_color)

    }

    if (level == "hierarchy") {
        interaction_type <- c("exclusion", "coexistence", "exclusion violating rank", "bistability", "neutrality", "self", "unknown")
        interaction_color <- c("#DB7469", "#557BAA", "#8CB369", "#EECF6D", "#8650C4", "black", "grey50")
        names(interaction_color) <- interaction_type
        return(interaction_color)
    }
    if (level == "finer") {
        interaction_type <- interaction_type_finer
        #interaction_type <- c("competitive exclusion", "stable coexistence", "mutual exclusion", "frequency-dependent coexistence", "neutrality", "exclusion violating rank")
        interaction_color <- c("#DB7469", "#557BAA",
                               "#FFBC42", "#B9FAF8",
                               "lightblue", "cyan",
                               "#8650C4", "purple",
                               grey(0.5))
        names(interaction_color) <- interaction_type
        return(interaction_color)
    }
}
interaction_color <- assign_interaction_color()
frequency_color <- c("95"="#292F36", "50"="#9F87AF", "5"="#7D7C7C")
category_color <- c(sugar = "#ED6A5A", acid = "#03CEA4", fermenter = "#8A89C0", respirator = "#FFCB77")
paint_white_background <- function () theme(plot.background = element_rect(fill = "white", color = NA))

#
# family_colors <- c(
#     Others = grey(0.5),
#     Enterobacteriaceae = "#397eb8",
#     Pseudomonadaceae = "#e21e26",
#     Aeromonadaceae = "#4fb148",
#     Sphingobacteriaceae = "#984e9e",
#     Moraxellaceae = "firebrick",
#     Comamonadaceae = "yellow",
#     Alcaligenaceae = "darkorchid2"
# )
#
# genus_colors <- c(
#     Others = grey(0.5),
#     Enterobactor1 = "#225ea8",
#     Klebsiella1 = "#3eb6c5",
#     Raoultella1 = "#a3d6b2",
#     Citrobacter1 = "#fcf8cf",
#     Pseudomonas1 = "#7d1517",
#     Pseudomonas2 = "#b31e24",
#     Pseudomonas3 = "#d63226",
#     Pseudomonas4 = "#e44b34",
#     Pseudomonas5 = "#ec6448",
#     Pseudomonas6 = "#f68d5c",
#     Aeromonas1 = "#8fd1c6"
# )

# outcome_colors <- c("1-exclusion" = "#CC79A7",
#                     "2-exclusion" = "#D55E00",
#                     "3-coexistence" = "#0072B2",
#                     "4-coexistence" = "#F0E442",
#                     "5-inconclusive" = "#999999")

outcome_colors <- c("1-exclusion" = "firebrick",
                    "2-exclusion" = "pink",
                    "3-coexistence" = "royalblue4",
                    "4-coexistence" = "lightblue",
                    "5-inconclusive" = "#999999")

outcome_labels <- c("exclusion with extinction",
                    "exclusion without extinction",
                    "coexistence with MIC",
                    "coexistence without MIC",
                    "inconclusive")

# Process

remove_ineligible_pairs <- function(pairs) {
    pairs %>%
        arrange(outcome, PairID) %>%
        mutate(PairID = factor(PairID, unique(PairID))) %>%
        # Remove no-colony pairs. six pairs
        drop_na(outcome) %>%
        # Remove low-accuracy model pairs. nine pairs
        filter(AccuracyMean > 0.9)
}

# flip_winner_species_freq <- function (pairs_freq) {
#     temp_index <- which(pairs_freq$Isolate1IsLoser)
#     if (length(temp_index) !=0) {
#     pairs_freq_flipped <- pairs_freq[temp_index, ] %>%
#         rename(temp = Isolate1, Isolate1 = Isolate2) %>%
#         rename(Isolate2 = temp) %>%
#         mutate(
#             Isolate1CFUFreqMean = 1-Isolate1CFUFreqMean,
#             Isolate1InitialODFreq = 100-Isolate1InitialODFreq,
#             Isolate1CFUFreqMedian = 1-Isolate1CFUFreqMedian,
#             Isolate1CFUFreqPercentile5 = 1-Isolate1CFUFreqPercentile5,
#             Isolate1CFUFreqPercentile95 = 1-Isolate1CFUFreqPercentile95
#         )
#
#     bind_rows(pairs_freq[-temp_index, ], pairs_freq_flipped) %>%
#         arrange(Time, PairID) %>%
#         return()
#     } else if (length(temp_index) == 0) {
#         return(pairs_freq)
#     }
# }






