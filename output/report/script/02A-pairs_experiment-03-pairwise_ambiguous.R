#' Subset the ambigious pairs that are hard to distinguish on TSA agar plates
#' and subset the isolates that are in these pairs.
library(tidyverse)
library(data.table)

# Ambiguous pairs
pairs_ambiguous <- fread(here::here("data/temp/pairs_competition.csv")) %>%
  filter(is.na(ColonyCount1)) %>%
  select(-starts_with("ColonyCount"), -DilutionFactor)

# Ambiguous isolates
isolates_ambiguous <-
  pairs_ambiguous %>%
  select(Community, Isolate1, Isolate2) %>% # Select columns
  gather(key = "XX", value = "Isolate", -Community) %>% dplyr::select(-XX) %>% # Melt the df
  distinct() %>% arrange(Community, Isolate)  %>% # Remove duplicate
  left_join(isolates_ID_match,  by = c("Community", "Isolate")) %>% # Match isolate ID
#  select(ExpID, ID, Community, Isolate, Experiment) %>% # Rearrange column order
  arrange(Community, Isolate) %>%
#  filter(!(Community == "C11R2" & Isolate == 13) & !(Experiment == "Transitivity_B2" & Community == "C11R1" & Isolate == 1)) %>% # Remove contamination
  distinct(ExpID, .keep_all = T) # Remove duplicate across batch


fwrite(pairs_ambiguous, file = here::here("data/temp/pairs_ambiguous.csv"))
fwrite(isolates_ambiguous, file = here::here("data/temp/isolates_ambiguous.csv"))
