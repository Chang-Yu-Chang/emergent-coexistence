# This script stores the metadata shared by all scripts

# 1. Directories
# This main folder depends on your home directory
folder_script <- "processing_scripts/" # Enter the directory of processing scripts
folder_pipeline <- "pipeline/" # Enter the directory of image processing pipeline
folder_data <- "data/" # Enter the directory of data


# 2. For image processing pipeline ----
list_folders <- c("01-channel", "02-rolled", "03-threshold", "04-round", "05-watershed", "06-transect", "07-feature", "08-random_forest", "09-bootstrap")
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
    "C_T8_C11R1_50-50_1_2",
    "C_T8_C11R1_50-50_1_3"
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
    ), "_green"),
    paste0(c("b.mean", "b.sd", "b.mad"), rep(c("_red", "_blue"), each = 3))
)

# 3. For plotting ----
paint_white_background <- function () theme(plot.background = element_rect(fill = "white", color = NA))
frequency_color <- c("95"="#292F36", "50"="#9F87AF", "5"="#7D7C7C")

outcome_colors <- c("1-exclusion" = "firebrick",
                    "2-exclusion" = "hotpink2",
                    "3-coexistence" = "dodgerblue3",
                    "4-coexistence" = "lightblue",
                    "5-inconclusive" = "#999999")

outcome_labels <- c("competitive exclusion",
                    "on the path to\ncompetitive exclusion",
                    "stable coexistence\n(mutual invasibility)",
                    "coexistence without\nevidence of mutual invasibility",
                    "inconclusive")
family_colors <- c(
    Others = grey(0.5),
    Enterobacteriaceae = "#397eb8",
    Pseudomonadaceae = "#e21e26",
    Aeromonadaceae = "#4fb148",
    Sphingobacteriaceae = "#984e9e",
    Moraxellaceae = "firebrick",
    Comamonadaceae = "yellow",
    Alcaligenaceae = "darkorchid2"
)


remove_ineligible_pairs <- function(pairs) {
    pairs %>%
        arrange(outcome, PairID) %>%
        mutate(PairID = factor(PairID, unique(PairID))) %>%
        # Remove no-colony pairs. six pairs
        drop_na(outcome) %>%
        # Remove low-accuracy model pairs. nine pairs
        filter(AccuracyMean > 0.9)
}



