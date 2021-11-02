# Plot the difference in growth rates on glucose and acetate
library(tidyverse)
library(cowplot)

# Data ----
isolates <- read_csv(here::here("data/output/isolates.csv"))
isolates_preference <- read_csv(here::here("data/temp/isolates_preference.csv"))
pairs <- read_csv(here::here("data/output/pairs.csv")); pairs$InteractionType[pairs$InteractionType == "neutrality"] <- "coexistence"
pairs_freq <- read_csv(here::here("data/output/pairs_freq.csv"))
# Growth rate data from Jean
isolates_growth <- read_csv(here::here("data/raw/growth_rate/Growthcurver.csv"))
# Growth rate data from Sylvie
isolates_growth_syl <- read_csv(here::here("data/raw/growth_rate/Estrela_2021_isolates_grmax.csv"))
# Byproduct measurement on glucose. Data from Sylvie
isolates_byproduct <- read_csv(here::here("data/raw/growth_rate/By_Products_Glucose.csv")) %>%
    select(OD620_16h = OD620, ID = SangerID, Glucose_perc, acetate_mM, succinate_mM, lactate_mM, gluconate_mM, ketogluconate_mM)
isolates_byproduct_time <- read_csv(here::here("data/raw/growth_rate/Estrela_2021_isolates_ph_OAs.csv")) %>%
    select(ID = SangerID, Time = time_hours, OD620, pH, Glucose_perc, acetate_mM, succinate_mM, lactate_mM)
# OD data from Jean. Filter the 16hr data only
isolates_OD_DW <- read_csv(here::here("data/raw/growth_rate/OD_Data_MH2.csv")) %>%
    filter(Time == 16) %>% select(ID = Sequence, CS = Carbon_Source, OD) %>%
    mutate(CS = sub("[LD]-", "", CS) %>% tolower() %>% paste0("OD620_16h_", .))


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
    select(Community, Isolate1, Isolate2, Isolate1Dominant, T0, T8)
pairs <- left_join(pairs, pairs_coexist_dominant)


# Swap so that fermenter is always isolate1
for (i in 1:nrow(pairs)) {
    if (is.na(pairs$Fermenter1[i]) | is.na(pairs$Fermenter2[i])) next
    if (!pairs$Fermenter1[i] & pairs$Fermenter2[i]) {
        temp <- pairs$Isolate1[i]
        pairs$Isolate1[i] <- pairs$Isolate2[i]
        pairs$Isolate2[i] <- temp
        temp <- pairs$Fermenter1[i]
        pairs$Fermenter1[i] <- pairs$Fermenter2[i]
        pairs$Fermenter2[i] <- temp
        temp <- pairs$Family1[i]
        pairs$Family1[i] <- pairs$Family2[i]
        pairs$Family2[i] <- temp
        temp <- pairs$Genus1[i]
        pairs$Genus1[i] <- pairs$Genus2[i]
        pairs$Genus2[i] <- temp
    }
}
# Swap so that winner is also isolate1
for (i in 1:nrow(pairs)) {
    #if (is.na(pairs$Fermenter1[i]) | is.na(pairs$Fermenter2[i])) next
    if (pairs$InteractionType[i] == "exclusion" & pairs$From[i] != pairs$Isolate1[i]) {
        temp <- pairs$Isolate1[i]
        pairs$Isolate1[i] <- pairs$Isolate2[i]
        pairs$Isolate2[i] <- temp
        temp <- pairs$Fermenter1[i]
        pairs$Fermenter1[i] <- pairs$Fermenter2[i]
        pairs$Fermenter2[i] <- temp
        temp <- pairs$Family1[i]
        pairs$Family1[i] <- pairs$Family2[i]
        pairs$Family2[i] <- temp
        temp <- pairs$Genus1[i]
        pairs$Genus1[i] <- pairs$Genus2[i]
        pairs$Genus2[i] <- temp
    }
}
# Swap so that in the coexistence pairs, isolate1 is the dominant strain
for (i in 1:nrow(pairs)) {
    if (is.na(pairs$Isolate1Dominant[i])) next
    if (pairs$Isolate1Dominant[i] == FALSE) {
        temp <- pairs$Isolate1[i]
        pairs$Isolate1[i] <- pairs$Isolate2[i]
        pairs$Isolate2[i] <- temp
        temp <- pairs$Fermenter1[i]
        pairs$Fermenter1[i] <- pairs$Fermenter2[i]
        pairs$Fermenter2[i] <- temp
        temp <- pairs$Family1[i]
        pairs$Family1[i] <- pairs$Family2[i]
        pairs$Family2[i] <- temp
        temp <- pairs$Genus1[i]
        pairs$Genus1[i] <- pairs$Genus2[i]
        pairs$Genus2[i] <- temp
        temp <- pairs$From[i]
        pairs$From[i] <- pairs$To[i]
        pairs$To[i] <- temp
    }
}

# Jean's growth curve data. Use the fitted Rmid
isolates_growth_w <- isolates_growth %>%
    separate(col = SID, sep = "_", into  = c("ID", "CS"), convert = T) %>%
    select(-SangerID, -Family) %>%
    select(ID, CS, RMid) %>%
    filter(CS %in% c("D-Glucose", "Acetate", "D-Lactate", "Succinate", "Gluconate", "2-Ketogluconate")) %>%
    pivot_wider(names_from = CS, values_from = RMid) %>%
    right_join(select(isolates, Community, Isolate, ID))

pairs_meta <- pairs %>%
    left_join(rename_with(isolates_growth_w, ~ paste0(., "1"), !contains("Community"))) %>%
    left_join(rename_with(isolates_growth_w, ~ paste0(., "2"), !contains("Community"))) %>%
    # Difference in glucose, acetate, lactate, succinate growth rate
    mutate(r_glu_d = `D-Glucose1` - `D-Glucose2`, r_ace_d = Acetate1 - Acetate2,
           r_lac_d = `D-Lactate1` - `D-Lactate2`, r_suc_d = Succinate1 - Succinate2)

# Isolate OD and glucose byproduct at 16hr. Sylvie's data
pairs_meta <- pairs_meta %>%
    left_join(rename_with(isolates_byproduct, ~ paste0(., "1"), !contains("Community"))) %>%
    left_join(rename_with(isolates_byproduct, ~ paste0(., "2"), !contains("Community")))

# Isolate's OD at 16 hr in DW96 plate. Jean's data
isolates_OD_DW_w <- isolates_OD_DW %>%
    filter(!is.na(ID), ID != "N/A") %>%
    mutate(ID = as.numeric(ID)) %>%
    pivot_wider(names_from = CS, values_from = OD)
pairs_meta <- pairs_meta %>%
    left_join(rename_with(isolates_OD_DW_w, ~ paste0(., "1"), !contains("Community"))) %>%
    left_join(rename_with(isolates_OD_DW_w, ~ paste0(., "2"), !contains("Community")))
# Match the preferred CS
pairs_meta <- pairs_meta %>%
    left_join(rename_with(isolates_preference, ~ paste0(., "1"))) %>%
    left_join(rename_with(isolates_preference, ~ paste0(., "2"))) %>%
    mutate(MatchedCS = PreferredCS1 == PreferredCS2)
# isolate pH
isolates_pH <- isolates_byproduct_time %>%
    select(ID, Time, pH) %>%
    pivot_wider(names_from = Time, values_from = pH, names_glue = "{Time}hr")
pairs_meta <- pairs_meta %>%
    left_join(rename_with(isolates_pH, ~ paste0(., "1"))) %>%
    left_join(rename_with(isolates_pH, ~ paste0(., "2"))) %>%
    mutate(pH0hr_d = `0hr1` - `0hr2`, pH16hr_d = `16hr1` - `16hr2`,
           pH28hr_d = `28hr1` - `28hr2`, pH48hr_d = `48hr1` - `48hr2`)
# Total amount of secretion
isolates_byproduct_time_sum <- isolates_byproduct_time %>%
    group_by(ID, Time) %>%
    mutate(ByproductSum = sum(2 * acetate_mM, 4 * succinate_mM, 3* lactate_mM)) %>%
    select(ID, Time, ByproductSum) %>%
    pivot_wider(names_from = Time, values_from = ByproductSum, names_glue = "ByproductSum_{Time}hr")
pairs_meta <- pairs_meta %>%
    left_join(rename_with(isolates_byproduct_time_sum, ~ paste0(., "1"), !contains("Time"))) %>%
    left_join(rename_with(isolates_byproduct_time_sum, ~ paste0(., "2"), !contains("Time"))) %>%
    mutate(byproduct_sum_16hr_d = ByproductSum_16hr1 - ByproductSum_16hr2,
           byproduct_sum_28hr_d = ByproductSum_28hr1 - ByproductSum_28hr2,
           byproduct_sum_48hr_d = ByproductSum_48hr1 - ByproductSum_48hr2)

write_csv(pairs_meta, file = here::here("data/output/pairs_meta.csv"))
