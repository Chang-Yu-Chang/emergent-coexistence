#' Read and match isolates' taxonomy (family, genus, and fermenter) to pairs
library(tidyverse)

# Read data ----
## Read isolates' RDP and ID match
isolates_RDP <- read_csv(here::here("data/temp/isolates_RDP.csv"))
isolates_ID_match <- read_csv(here::here("data/temp/isolates_ID_match.csv")) %>% select(ExpID, ID, Community, Isolate)

## Read pairs
pairs_ID <- read_csv(here::here("data/temp/pairs_ID.csv"))

## Match isolates' information
isolates_RDP_ID <- isolates_ID_match %>% left_join(isolates_RDP, by = c("ID"))

# Match the isolates taxonomic information to pairs
pairs_taxonomy <- pairs_ID %>%
    left_join(rename_with(isolates_RDP_ID, ~ paste0(., "1"), !contains("Community"))) %>%
    left_join(rename_with(isolates_RDP_ID, ~ paste0(., "2"), !contains("Community")))

if (FALSE) {
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

pairs_taxonomy <- pairs_ID %>%
  match_isolates_pairs(isolates_df = isolates_RDP_ID, variable_name = "Fermenter") %>%
  match_isolates_pairs(isolates_df = isolates_RDP_ID, variable_name = "Family") %>%
  match_isolates_pairs(isolates_df = isolates_RDP_ID, variable_name = "Genus") %>%
  as_tibble()

}


# Assign taxonomy
pairs_taxonomy$PairFermenter <- NA
pairs_taxonomy$PairFamily <- NA
## Fermenter
pairs_taxonomy$PairFermenter[pairs_taxonomy$Fermenter1 == T & pairs_taxonomy$Fermenter2 == T] <- "FF"
pairs_taxonomy$PairFermenter[pairs_taxonomy$Fermenter1 == F & pairs_taxonomy$Fermenter2 == F] <- "NN"
pairs_taxonomy$PairFermenter[(pairs_taxonomy$Fermenter1 == T & pairs_taxonomy$Fermenter2 == F) | pairs_taxonomy$Fermenter1 == F & pairs_taxonomy$Fermenter2 == T] <- "FN"

## Family
pairs_taxonomy$PairFamily[pairs_taxonomy$Family1 == "Enterobacteriaceae" & pairs_taxonomy$Family2 == "Enterobacteriaceae"] <- "EE"
pairs_taxonomy$PairFamily[pairs_taxonomy$Family1 == "Pseudomonadaceae" & pairs_taxonomy$Family2 == "Pseudomonadaceae"] <- "PP"
pairs_taxonomy$PairFamily[(pairs_taxonomy$Family1 == "Enterobacteriaceae" & pairs_taxonomy$Family2 == "Pseudomonadaceae") | (pairs_taxonomy$Family1 == "Pseudomonadaceae" & pairs_taxonomy$Family2 == "Enterobacteriaceae")] <- "EP"

#
write_csv(pairs_taxonomy, file = here::here("data/temp/pairs_taxonomy.csv"))

