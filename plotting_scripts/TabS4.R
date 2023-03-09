library(tidyverse)
library(cowplot)
library(officer)
library(flextable)
source(here::here("analysis/00-metadata.R"))

isolates <- read_csv(paste0(folder_data, "output/isolates.csv"), show_col_types = F)
pairs <- read_csv(paste0(folder_data, "output/pairs.csv"), show_col_types = F)
communities <- read_csv(paste0(folder_data, "temp/00c-communities.csv"), show_col_types = F)


# Table S4 fitness function ----
make_interaction_type <- function () {
    #' This function generates the fitness function table.
    #' There are a total of 27 possibilities
    interaction_type <- tibble(
        FromRare = rep(c(1, -1, 0), each = 9),
        FromMedium = rep(rep(c(1, -1, 0), each = 3), 3),
        FromAbundant = rep(c(1, -1, 0), 9),
        InteractionType = NA,
        InteractionTypeFiner = NA
    )
    # Assign interaction types to combinations of frequency changes signs
    interaction_type$InteractionType[c(1,10,13,14)] <- "exclusion"
    interaction_type$InteractionType[c(2:6,8,11,20,23, 9,18,21,24,25,26,27)] <- "coexistence"

    # Assign finer interaction types to combinations of frequency changes signs
    interaction_type$InteractionTypeFiner[c(1,14)] <- "competitive exclusion"
    interaction_type$InteractionTypeFiner[c(10,13)] <- "mutual exclusion"
    interaction_type$InteractionTypeFiner[c(2,5,8)] <- "stable coexistence"
    interaction_type$InteractionTypeFiner[c(4,6,11,20)] <- "frequency-dependent coexistence"
    interaction_type$InteractionTypeFiner[c(3)] <- "coexistence at 95%"
    interaction_type$InteractionTypeFiner[c(23)] <- "coexistence at 5%"
    interaction_type$InteractionTypeFiner[c(9,18,21,24:26)] <- "neutrality"
    interaction_type$InteractionTypeFiner[c(9,18,21,24:26)] <- "2-freq neutrality"
    interaction_type$InteractionTypeFiner[c(27)] <- "3-freq neutrality"
    interaction_type <- interaction_type %>%  mutate(FitnessFunction = paste(FromRare, FromMedium, FromAbundant, sep = "_"))

    # Fill in unknown
    interaction_type %>%
        replace_na(list(InteractionType = "unknown", InteractionTypeFiner = "unknown"))
}
reformat <- function(x) {
    x <- ifelse(x == "1", "+", x)
    x <- ifelse(x == "-1", "-", x)
    x <- factor(x, c("+", "-", "0"))
    return(x)
}
interaction_type <- make_interaction_type()

pairs_interaction <- pairs %>%
    filter(!is.na(InteractionType)) %>%
    filter(AccuracyMean > 0.9) %>%
    group_by(InteractionType, InteractionTypeFiner, FitnessFunction) %>%
    count(name = "Count")

ft <- interaction_type %>%
    left_join(pairs_interaction, by = c("InteractionType", "InteractionTypeFiner", "FitnessFunction")) %>%
    replace_na(list(Count = 0)) %>%
    mutate(across(starts_with("From"), reformat)) %>%
    select(-FitnessFunction) %>%
    rename(`Outcome` = InteractionType, `Finer outcome` = InteractionTypeFiner, `From rare` = FromRare, `From medium` = FromMedium, `From abundant` = FromAbundant) %>%
    mutate(` ` = 1:n()) %>%
    select(` `, everything()) %>%
    janitor::adorn_totals() %>%
    flextable() %>%
    width(j = 1, width = .5) %>%
    width(j = 2:5, width = 1.3) %>%
    width(j = 6, width = 3) %>%
    vline(j = 1, border = NULL, part = "all") %>%
    hline(i = 27, border = fp_border(color = "black", style = "solid", width = 2))

save_as_image(ft, here::here("plots/TableS4-fitness_function.png"), webshot = "webshot2")




























