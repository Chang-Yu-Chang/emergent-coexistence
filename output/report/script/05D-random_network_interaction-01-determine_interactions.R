#' Determine interaction types in CASEU random netwoork results
#'
#' 1. Read and plot the data of frequency changes
#' 2. Compute the sign of frequency changes
#' 3. Determine the interaction types by the fitness function
#' 4. Combine with the isolates information

library(tidyverse)
library(data.table)
library(cowplot)

# CASEU result
reshape_df <- function(x, com, tim = NULL) {
  x %>%
    filter(Community == com) %>%
    {if (tim %in% c("T0", "T3")) mutate(., Time = tim) else  .} %>%
    mutate(RawDataType = "Sanger", Contamination = F, ErrorIsolate1MeasuredFreq = NA) %>%
    select(Community, Isolate1, Isolate2,
           Isolate1InitialODFreq = Isolate1Freq, Isolate2InitialODFreq = Isolate2Freq,
           Time,
           Isolate1MeasuredFreq = Isolate1FreqPredicted, ErrorIsolate1MeasuredFreq,
           RawDataType, Contamination)
}

## RanAss1 in T0 C P2
RanAss1_T0 <- fread(here::here("data/temp/CASEU_RN3.csv")) %>% as_tibble() %>%
  reshape_df(com = "RanAss1", tim = "T0")

## RanAss1 in T3 C P2
RanAss1_T3 <- fread(here::here("data/temp/CASEU_RN2.csv")) %>% as_tibble() %>%
  reshape_df(com = "RanAss1", tim = "T3")
RanAss1 <- bind_rows(RanAss1_T0, RanAss1_T3)

## AcrAss1 in T0 AD P2 and T3 AD P2
AcrAss1 <- fread(here::here("data/temp/CASEU_RN3.csv")) %>% as_tibble() %>% # T3 AD P2
  reshape_df(com = "AcrAss1", tim = NA)

## RanAss2 in T0 AD P2 and T3 AD P2
RanAss2 <- fread(here::here("data/temp/CASEU_RN3.csv")) %>% as_tibble() %>% # T3 AD P2
  reshape_df(com = "RanAss2", tim = NA)

## Isolates
isolates_random <- fread(here::here("data/output/isolates_random.csv")) %>%
  select(Community = AssemblyCommunity, Isolate = AssemblyIsolate, ExpID, Family, Genus)

"
Combine the result I have for AcrAss1 T0 T3, RanAss T0 T3, half of RanAss2 T0 T3
"


# Save the relative abundance at T3; Reshape to the pairwise competition data format ----
random_network_pairs_melted <- bind_rows(AcrAss1, RanAss1, RanAss2)

# Sign of frequency changes ----
random_network_pairs_frequency <-
  random_network_pairs_melted %>%
  pivot_wider(names_from = Time, values_from = Isolate1MeasuredFreq) %>%
  select(Community, Isolate1, Isolate2, InitialFrequency = Isolate1InitialODFreq, T0, T3) %>%
  group_by(Community) %>%
  arrange(Community, Isolate1, Isolate2, InitialFrequency) %>%
  mutate(FrequencyChange = (T3 - T0) > 0) %>%
  select(-T0, -T3)

# Determine interactions ----
interaction_type <- tibble(
  `0.05` = rep(c(T,F), each = 2),
  `0.95` = rep(c(T,F), 2),
  InteractionType = c("exclusion", "coexistence", "bistability", "exclusion")
)

random_network_pairs_fitness <-
  random_network_pairs_frequency %>%
  pivot_wider(names_from = InitialFrequency, values_from = FrequencyChange) %>%
  filter(!is.na(`0.05`), !is.na(`0.95`)) %>%
  left_join(interaction_type, by = c("0.05", "0.95")) %>%
  mutate(From = Isolate1, To = Isolate2)

# Switch arrow sign when isolate2 outcompetes isolate1
temp_index <- random_network_pairs_fitness$`0.05` == F & random_network_pairs_fitness$`0.95` == F
random_network_pairs_fitness[temp_index, c("From", "To")] <- random_network_pairs_fitness[temp_index, c("To", "From")]
random_network_pairs_interaction <- random_network_pairs_fitness %>%
  select(Community, Isolate1, Isolate2, InteractionType, From, To) %>%
  mutate(Isolate1 = factor(Isolate1), Isolate2 = factor(Isolate2))


# Join isolate information ----
match_isolates_pairs <- function (pairs_df, isolates_df, variable_name = "Epsilon") {
  if (!("Isolate" %in% names(isolates_df))) {
    stop("Isolate df does not have a column of Isolate")
  } else if (!("Isolate1" %in% names(pairs_df)) | !("Isolate1" %in% names(pairs_df))) {
    stop("Pair df does not have columns of Isolate1 and Isolate2")
  }

  # Congifurate isolate df
  temp <- isolates_df %>%
    mutate(Isolate1 = ordered(Isolate, levels = 1:12), Isolate2 = ordered(Isolate, levels = 1:12)) %>%
    data.table::as.data.table()
  temp[,paste0("Isolate", 1:2)] <- temp[,c("Isolate", "Isolate")]
  temp[,paste0(variable_name, 1:2)] <- temp[,rep(variable_name, 2), with = FALSE] # temp is a data.table object
  temp <- as_tibble(temp) %>% mutate(Isolate1 = ordered(Isolate1, levels = 1:12), Isolate2 = ordered(Isolate2, levels = 1:12))

  # Match isolate 1 and isolate 2 by join isolates df to pairs df
  pairs_df %>%
    mutate(Isolate1 = ordered(Isolate1, levels = 1:12), Isolate2 = ordered(Isolate2, levels = 1:12)) %>%
    left_join(temp[,c("Community", "Isolate1", paste0(variable_name, 1))], by = c("Community", "Isolate1")) %>%
    left_join(temp[,c("Community", "Isolate2", paste0(variable_name, 2))], by = c("Community", "Isolate2")) %>%
    invisible() %>%
    return()
}
random_network_pairs <- random_network_pairs_interaction %>%
  match_isolates_pairs(isolates_df = isolates_random, variable_name = "Family") %>%
  match_isolates_pairs(isolates_df = isolates_random, variable_name = "Genus")

fwrite(random_network_pairs, file = here::here("data/output/random_network_pairs.csv"))
fwrite(random_network_pairs_frequency, file = here::here("data/temp/random_network_pairs_frequency.csv"))


















