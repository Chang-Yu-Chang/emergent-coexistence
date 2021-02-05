#' This script reads and merge the colony counts of pairs and their diluton factor while plating
library(tidyverse)
library(data.table)
switch_pairwise_column <- function (df, bypair = T) {
  if (any(is.factor(df$Isolate1))) df$Isolate1 <- as.numeric(df$Isolate1); df$Isolate2 <- as.numeric(df$Isolate2)
  if ("Isolate1FreqPredicted" %in% colnames(df)) {
    if (bypair == T) {
      temp_index <- df$Isolate1 > df$Isolate2
      df[temp_index, c("Isolate1", "Isolate2", "Isolate1Freq", "Isolate2Freq", "Isolate1FreqPredicted", "Isolate2FreqPredicted")] <-
        df[temp_index, c("Isolate2", "Isolate1", "Isolate2Freq", "Isolate1Freq", "Isolate2FreqPredicted", "Isolate1FreqPredicted")]

      df %>% arrange(Isolate1, Isolate2, Isolate1Freq) %>% return()
    } else if (bypair == F) {
      temp_index <- df$Isolate1Freq == 5
      df[temp_index, c("Isolate1", "Isolate2", "Isolate1Freq", "Isolate2Freq", "Isolate1FreqPredicted", "Isolate2FreqPredicted")] <-
        df[temp_index, c("Isolate2", "Isolate1", "Isolate2Freq", "Isolate1Freq", "Isolate2FreqPredicted", "Isolate1FreqPredicted")]

      df %>% arrange(Isolate1Freq, Isolate1, Isolate2) %>% return()
    }
  } else {

    if (bypair == T) {
      temp_index <- df$Isolate1 > df$Isolate2
      df[temp_index, c("Isolate1", "Isolate2", "Isolate1Freq", "Isolate2Freq")] <-
        df[temp_index, c("Isolate2", "Isolate1", "Isolate2Freq", "Isolate1Freq")]

      df %>% arrange(Isolate1, Isolate2, Isolate1Freq) %>% return()
    } else if (bypair == F) {
      temp_index <- df$Isolate1Freq == 5
      df[temp_index, c("Isolate1", "Isolate2", "Isolate1Freq", "Isolate2Freq")] <-
        df[temp_index, c("Isolate2", "Isolate1", "Isolate2Freq", "Isolate1Freq")]

      df %>% arrange(Isolate1Freq, Isolate1, Isolate2) %>% return()
    }
  }
}

# Read colony counts ----
## TSA plates
pairs_CFU <-
  fread(here::here("data/raw/pairwise_competition/result_pairwise_competition_arranged.csv"), stringsAsFactors = T) %>%
  # Remove contamination C11R2 isolate 13, which is a Staphylococcus
  filter(!(Community == "C11R2" & Isolate1 == 13), !(Community == "C11R2" & Isolate2 == 13)) %>%
#  filter(!(Community == "C11R2" & Isolate1 == 8 &Isolate2 == 9 & ColonyCount <= 5)) %>%
  # Remove contamination batch B2 community C11R1 isolate 1, which I repeate the experiment at batch C
  filter(!(Experiment == "Transitivity_B2" & Community == "C11R1" & Isolate1 == 1)) %>%
  # Configurate the attributes
  mutate(
    Community = ordered(Community, levels = communities_name),
    Isolate1 = ordered(Isolate1, 1:12), Isolate2 = ordered(Isolate2, 1:12),
    Isolate1Freq = as.character(Isolate1Freq), Isolate2Freq = as.character(Isolate2Freq),
    # Colony counts of isolate 2
    ColonyCount2 = ColonyCount - ColonyCount1
  ) %>%
  select(Experiment, Community, Isolate1, Isolate2, Isolate1Freq, Isolate2Freq, ColonyCount, ColonyCount1, ColonyCount2, Contamination) %>%
  as_tibble() %>%
  # Switch isolate1 and isolate2 so that the order of isolate1 is always smaller than isolate2
  switch_pairwise_column(bypair = T) %>%
  arrange(Community, Isolate1, Isolate2, Isolate1Freq, Isolate2Freq)

## Chromogentic plates; they already have dilution factor in the df
pairs_competition_chromo <-
  fread(here::here("data/raw/pairwise_competition/result_pairwise_competition_chromogenic.csv")) %>%
  # Remove the pairs that I still cannot tell on chromogenic agar plates
  filter(!is.na(ColonyCount)) %>%
  # Plating from glyerol stock, add 0.5 to dilution factor so to process later
  mutate(DilutionFactor = DilutionFactor + 0.5) %>%
  # Configurate the attributes
  mutate(
    Experiment = rep("Chromogenic", nrow(.)),
    Community = ordered(Community, levels = communities_name),
    Isolate1 = ordered(Isolate1, 1:12), Isolate2 = ordered(Isolate2, 1:12),
    Isolate1Freq = as.character(Isolate1Freq), Isolate2Freq = as.character(Isolate2Freq),
    # Colony counts of isolate 2
    ColonyCount2 = ColonyCount - ColonyCount1
  ) %>%
  select(Experiment, Community, Isolate1, Isolate2, Isolate1Freq, Isolate2Freq, ColonyCount, ColonyCount1, ColonyCount2, DilutionFactor) %>%
  as_tibble() %>%
  # Switch isolate1 and isolate2 so that the order of isolate1 is always smaller than isolate2
  switch_pairwise_column(bypair = T) %>%
  arrange(Community, Isolate1, Isolate2, Isolate1Freq, Isolate2Freq)

# Read dilution factors from the csv for colony counts at T8 ----
pairs_DF <-
  fread(here::here("data/raw/pairwise_competition/colony_count_dilution_factor.csv")) %>%
  # Remove duplicate plating
  filter(!(Isolate2 %in% c("1a", "1b"))) %>%
  # Remove monoculuture
  filter(Isolate1Freq != "monoculture", Isolate1 != Isolate2) %>%
  # Remove contamination C11R2 isolate 13, which is a Staphylococcus
  filter(Isolate1 != 13, Isolate2 != 13) %>%
  # Remove contamination batch B2 community C11R1 isolate 1, which I repeate the experiment at batch C
  filter(!(Community == "C11R1" & Batch == "B2" & (Isolate1 == 1 | Isolate2 == 1))) %>%
  # Remove duplicate that I used as controls in  batch C
  filter(
    !(Community == "C11R1" & Batch == "C" & ((Isolate1 == 5 & Isolate2 == 6) | (Isolate1 == 6 & Isolate2 == 5))),
    !(Community == "C11R1" & Batch == "C" & ((Isolate1 == 5 & Isolate2 == 7) | (Isolate1 == 7 & Isolate2 == 5)))) %>%
  # Configurate the attributes
  mutate(
    Experiment = paste0("Transtivity_", Batch),
    Community = ordered(Community, levels = communities_name),
    Isolate1 = ordered(Isolate1, 1:12), Isolate2 = ordered(Isolate2, 1:12),
    Isolate1Freq = as.character(Isolate1Freq), Isolate2Freq = as.character(Isolate2Freq)
  ) %>%
  select(Experiment, Community, Isolate1, Isolate2, Isolate1Freq, Isolate2Freq, DilutionFactor) %>%
  as_tibble() %>%
  switch_pairwise_column(bypair = T) %>%
  arrange(Community, Isolate1, Isolate2, Isolate1Freq, Isolate2Freq)

# Merge df of colony counts and df of dilution factors ----
pairs_competition <- pairs_CFU %>%
  left_join(pairs_DF, by = c("Experiment", "Community", "Isolate1", "Isolate2", "Isolate1Freq", "Isolate2Freq")) %>%
  switch_pairwise_column(bypair = T) %>%
  mutate(Community = ordered(Community, levels = communities_name)) %>%
  arrange(Community, Isolate1, Isolate2, Isolate1Freq, Isolate2Freq)


# Fill in the pairs specified by chromogenic plates ----
for (i in 1:nrow(pairs_competition_chromo)) {
  index_temp <- which(
    pairs_competition$Community == pairs_competition_chromo$Community[i] &
      pairs_competition$Isolate1 == pairs_competition_chromo$Isolate1[i] &
      pairs_competition$Isolate2 == pairs_competition_chromo$Isolate2[i] &
      pairs_competition$Isolate1Freq == pairs_competition_chromo$Isolate1Freq[i] &
      pairs_competition$Isolate2Freq == pairs_competition_chromo$Isolate2Freq[i]
  )
  pairs_competition[index_temp, c("ColonyCount", "ColonyCount1", "ColonyCount2", "DilutionFactor")] <-
    pairs_competition_chromo[i, c("ColonyCount", "ColonyCount1", "ColonyCount2", "DilutionFactor")] %>%
    mutate(DilutionFactor = as.integer(DilutionFactor))
}

#
pairs_ID <- pairs_competition %>%
  select(Community, Isolate1, Isolate2) %>%
  distinct()
#
fwrite(pairs_competition, file = here::here("data/temp/pairs_competition.csv"))
fwrite(pairs_ID, file = here::here("data/temp/pairs_ID.csv"))







