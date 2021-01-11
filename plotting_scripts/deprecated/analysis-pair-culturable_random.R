# Random culturable isolates

suppressWarnings(suppressMessages(library(tidyverse)))
suppressWarnings(suppressMessages(library(data.table)))
#options(dplyr.summarise.inform = FALSE)
suppressWarnings(suppressMessages(library(cowplot)))

#args = commandArgs(trailingOnly = T)
#input_independent <- fread(args[[1]])
#input_independent <- fread("../data/raw/simulation/mapping_files/input_independent.csv")
#input_independent_pairs <- input_independent %>% filter(grepl("pair-culturable_isolates", exp_id))


#
read_pair_list <- function (pair_culturable_list) {
    # Mapping list of pairs, well, inital frequency, and isolate IDs
    pair_culturable_list %>%
        filter(Type == "consumer") %>%
        select(Well, Pair, InitialFrequency, ID) %>%
        group_by(Well) %>%
        mutate(Isolate = 1:2) %>%
        pivot_wider(names_from = Isolate, names_prefix = "Isolate", values_from = ID)
}
read_pair_competition <- function (pair_culturable_data, pair_list) {
    pair_culturable_data %>%
        filter(Type == "consumer") %>%
        left_join(select(pair_list, Well, Pair, InitialFrequency), by = "Well") %>%
        group_by(Pair, InitialFrequency, Transfer) %>%
        mutate(ID = factor(ID)) %>%
        mutate(TotalAbundance = sum(Abundance), RelativeAbundance = Abundance/TotalAbundance) %>%
        select(Pair, InitialFrequency, Transfer, ID, RelativeAbundance)
}
determine_pair_outcome <- function(pair_competition, pair_list) {
    frequency_changes <- pair_competition %>%
        left_join(pair_list, by = c("Pair", "InitialFrequency")) %>%
        filter(ID == Isolate1) %>% select(-ID) %>%
        group_by(Pair) %>%
        pivot_wider(names_from = c(Transfer), names_prefix = "T", values_from = RelativeAbundance) %>%
        mutate(FrequencyChange = ifelse((T5 - T0)>0, "T", "F")) %>%  # TRUE = Isolate1 increases, FALSE = Isolate1 decreases
        select(Pair, Isolate1, Isolate2, FrequencyChange) %>%
        group_by(Pair, Isolate1, Isolate2) %>%
        summarise(FrequencyChangePattern = paste0(FrequencyChange, collapse = "-"))

    map_frequency_pattern <- tibble(FrequencyChangePattern = c("T-T-T", "F-F-F", "T-F-F", "T-T-F", "T-F-T", "F-T-F", "F-F-T", "F-T-T"),
        InteractionType = c("exclusion", "exclusion", "coexistence", "coexistence", "coexistence", "coexistence", "mutual exclusion", "mutual exclusion"))

    frequency_changes$From = NA
    frequency_changes$To = NA
    frequency_changes$InteractionType = NA

    for (i in 1:nrow(frequency_changes)) {
        if (frequency_changes$FrequencyChangePattern[i] == "T-T-T") {
            from <- frequency_changes$Isolate1[i]
            to <- frequency_changes$Isolate2[i]
            interaction_type <- "exclusion"
        } else if (frequency_changes$FrequencyChangePattern[i] == "F-F-F") {
            from <- frequency_changes$Isolate2[i]
            to <- frequency_changes$Isolate1[i]
            interaction_type <- "exclusion"
        } else if (frequency_changes$FrequencyChangePattern[i] %in% c("T-F-F", "T-T-F", "T-F-T", "F-T-F")) {
            from <- frequency_changes$Isolate1[i]
            to <- frequency_changes$Isolate2[i]
            interaction_type <- "coexistence"
        } else {
            from <- NA
            to <- NA
            interaction_type <- NA
        }

        frequency_changes$From[i] <- from
        frequency_changes$To[i] <- to
        frequency_changes$InteractionType[i] <- interaction_type
    }

    frequency_changes %>%
        select(Pair, Isolate1, Isolate2, InteractionType, From, To) %>%
        return()
}

# temp_list <- rep(list(NA), nrow(input_independent_pairs))
# for (i in 1:nrow(input_independent_pairs)) {
#     cat("\nexp_id = ", input_independent_pairs$exp_id[i])
#     cat(",\tseed = ", input_independent_pairs$seed[i])
#     df_pair_list <- fread(paste0("../data/raw/simulation/pair-culturable-", i, ".txt")) %>%
#         read_pair_list()
#
#     df_pair_competition <-
#         fread(paste0("../data/raw/simulation/pair-culturable_isolates-", i, "_composition.txt")) %>%
#         read_pair_competition(df_pair_list)
#
#     df_pair_outcome <- determine_pair_outcome(df_pair_competition, df_pair_list)
#
#     temp_list[[i]] <- df_pair_outcome
# }
#
# df_pair_outcomes <- bind_rows(temp_list, .id = "Seed")
#
# #mapping_seed <- tibble(Seed = factor(1:4), l = c(0.2, 0, 0, 0.2), q = c(0.8, 0, 0.8, 0))
#
# p1 <- df_pair_outcomes %>%
# #    left_join(mapping_seed) %>%
#     group_by(Seed, InteractionType) %>%
#     summarise(Count = n()) %>%
#     ggplot() +
#     geom_bar(aes(x = Seed, y = Count, fill = InteractionType), stat = "identity") +
#     scale_x_discrete(expand = c(0,0)) +
#     scale_y_continuous(expand = c(0,0)) +
# #    facet_wrap(l~q, labeller = label_both, nrow = 1, scales = "free_x") +
#     theme_cowplot() +
#     theme(legend.position = "top", legend.title = element_blank()) +
#     panel_border(color = 1) +
#     ggtitle("Random pairs of culturable isolates")
# #cat("\n Made the pairwise outcome across seeds")
# p1
# ggsave("../plots/Fig2A.png", plot = p1, width = 5, height = 5)
#
#
#
#
#
#
#
#
#
