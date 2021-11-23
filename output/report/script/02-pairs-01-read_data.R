#' Read and combine data for pairwise competition
library(tidyverse)

communities <- read_csv(here::here("data/output/communities.csv"))

# Pairwise only data ----
pairs_ID <- read_csv(here::here("data/temp/pairs_ID.csv"))

## From 02D
pairs_freq <- read_csv(here::here("data/output/pairs_freq.csv"))
pairs_interaction <- read_csv(here::here("data/temp/pairs_interaction.csv"))

## From 02E
pairs_taxonomy <- read_csv(here::here("data/temp/pairs_taxonomy.csv"))
pairs_16S <- read_csv(here::here("data/temp/pairs_16S.csv"))
#pairs_tree_distance <- read_csv(here::here("data/temp/pairs_tree_distance.csv"))

## Combine data by assigning pairs information.186 pairs
pairs <- pairs_ID %>%
    left_join(pairs_interaction, by = c("Community", "Isolate1", "Isolate2")) %>%
    left_join(pairs_taxonomy, by = c("Community", "Isolate1", "Isolate2")) %>%
    left_join(pairs_16S, by = c("Community", "Isolate1", "Isolate2"))


# Flip the sign so that Isolate1 is the dominant and Isolate2 is subdominant
#' 50-50 frequency changes. To determine the dominance in coexistence pairs.
#' The isolate that increases frequency in the 50-50 competition is the dominant strain
pairs_coexist_dominant <- pairs %>%
    filter(InteractionType == "coexistence") %>%
    select(Community, Isolate1, Isolate2) %>%
    left_join(pairs_freq) %>%
    filter(Isolate1InitialODFreq == 50) %>%
    select(Community, Isolate1, Isolate2, Time, Isolate1MeasuredFreq) %>%
    pivot_wider(names_from = "Time", values_from = "Isolate1MeasuredFreq") %>%
    mutate(Isolate1Dominant = T8 > T0) %>%
    select(Community, Isolate1, Isolate2, Isolate1Dominant)

pairs <- left_join(pairs, pairs_coexist_dominant) %>%
        mutate(Assembly = ifelse(str_detect(Community, "C\\d+R\\d+"), "self_assembly",
                                 ifelse(str_detect(Community, "AcrAss"), "across_community",
                                        ifelse(str_detect(Community, "RanAss"), "random_assembly", NA)))) %>%
    select(Assembly, everything())


#
"
remove this line if the random network data is completed
"
pairs <- pairs %>% filter(!(is.na(InteractionType) & Assembly != "self_assembly"))

#
pairs_meta <- pairs


# Swap so that fermenter is always isolate1
for (i in 1:nrow(pairs_meta)) {
    if (is.na(pairs_meta$Fermenter1[i]) | is.na(pairs_meta$Fermenter2[i])) next
    if (!pairs_meta$Fermenter1[i] & pairs_meta$Fermenter2[i]) {
        print(i)
        temp <- pairs_meta$Isolate1[i]
        pairs_meta$Isolate1[i] <- pairs_meta$Isolate2[i]
        pairs_meta$Isolate2[i] <- temp
        temp <- pairs_meta$Fermenter1[i]
        pairs_meta$Fermenter1[i] <- pairs_meta$Fermenter2[i]
        pairs_meta$Fermenter2[i] <- temp
        temp <- pairs_meta$Family1[i]
        pairs_meta$Family1[i] <- pairs_meta$Family2[i]
        pairs_meta$Family2[i] <- temp
        temp <- pairs_meta$Genus1[i]
        pairs_meta$Genus1[i] <- pairs_meta$Genus2[i]
        pairs_meta$Genus2[i] <- temp
    }
}
# Swap so that winner is also isolate1
for (i in 1:nrow(pairs_meta)) {
    #if (is.na(pairs_meta$Fermenter1[i]) | is.na(pairs_meta$Fermenter2[i])) next
    if (pairs_meta$InteractionType[i] == "exclusion" & pairs_meta$From[i] != pairs_meta$Isolate1[i]) {
        print(i)
        temp <- pairs_meta$Isolate1[i]
        pairs_meta$Isolate1[i] <- pairs_meta$Isolate2[i]
        pairs_meta$Isolate2[i] <- temp
        temp <- pairs_meta$Fermenter1[i]
        pairs_meta$Fermenter1[i] <- pairs_meta$Fermenter2[i]
        pairs_meta$Fermenter2[i] <- temp
        temp <- pairs_meta$Family1[i]
        pairs_meta$Family1[i] <- pairs_meta$Family2[i]
        pairs_meta$Family2[i] <- temp
        temp <- pairs_meta$Genus1[i]
        pairs_meta$Genus1[i] <- pairs_meta$Genus2[i]
        pairs_meta$Genus2[i] <- temp
    }
}
# Swap so that in the coexistence pairs_meta, isolate1 is the dominant strain
for (i in 1:nrow(pairs_meta)) {
    if (is.na(pairs_meta$Isolate1Dominant[i])) next
    if (pairs_meta$Isolate1Dominant[i] == FALSE & pairs_meta$InteractionType[i] == "coexistence") {
        print(i)
        temp <- pairs_meta$Isolate1[i]
        pairs_meta$Isolate1[i] <- pairs_meta$Isolate2[i]
        pairs_meta$Isolate2[i] <- temp
        temp <- pairs_meta$Fermenter1[i]
        pairs_meta$Fermenter1[i] <- pairs_meta$Fermenter2[i]
        pairs_meta$Fermenter2[i] <- temp
        temp <- pairs_meta$Family1[i]
        pairs_meta$Family1[i] <- pairs_meta$Family2[i]
        pairs_meta$Family2[i] <- temp
        temp <- pairs_meta$Genus1[i]
        pairs_meta$Genus1[i] <- pairs_meta$Genus2[i]
        pairs_meta$Genus2[i] <- temp
        temp <- pairs_meta$From[i]
        pairs_meta$From[i] <- pairs_meta$To[i]
        pairs_meta$To[i] <- temp
    }
}


# Isolates data from 01-isolates ----
isolates_to_be_joint <- read_csv(here::here("data/output/isolates.csv")) %>%
    select(Community, Isolate, ID, starts_with("rmid_"), starts_with("rmax_"), starts_with("X_"), starts_with("pH_"), starts_with("OD620_"), ends_with("CS"))

pairs_meta <- pairs_meta %>%
    left_join(rename_with(isolates_to_be_joint, ~ paste0(., "1"), !contains("Community"))) %>%
    left_join(rename_with(isolates_to_be_joint, ~ paste0(., "2"), !contains("Community"))) %>%
    mutate(
        # Difference in glucose, acetate, lactate, succinate growth rate
        rmid_glu_d = rmid_glucose1 - rmid_glucose2, rmid_ace_d = rmid_acetate1 - rmid_acetate2,
        rmid_lac_d = rmid_lactate1 - rmid_lactate2, rmid_suc_d = rmid_succinate1 - rmid_succinate2,
        rmax_glu_d = rmax_glucose1 - rmax_glucose2, rmax_ace_d = rmax_acetate1 - rmax_acetate2,
        rmax_lac_d = rmax_lactate1 - rmax_lactate2, rmax_suc_d = rmax_succinate1 - rmax_succinate2,
        # Matched CS
        MatchedCS = PreferredCS1 == PreferredCS2,
        # Difference in pH
        pH_0hr_d = pH_0hr1 - pH_0hr2, pH_16hr_d = pH_16hr1 - pH_16hr2,
        pH_28hr_d = pH_28hr1 - pH_28hr2, pH48hr_d = pH_48hr1 - pH_48hr2,
        # Difference total acid production
        X_sum_0hr_d = X_sum_0hr1 - X_sum_0hr2, X_sum_16hr_d = X_sum_16hr1 - X_sum_16hr2,
        X_sum_28hr_d = X_sum_28hr1 - X_sum_28hr2, X_sum_48hr_d = X_sum_48hr1 - X_sum_48hr2
    )

#
write_csv(pairs, here::here("data/output/pairs.csv"))
write_csv(pairs_meta, file = here::here("data/output/pairs_meta.csv"))


