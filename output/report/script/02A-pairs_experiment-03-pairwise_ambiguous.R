#' Subset the ambiguous pairs that are hard to distinguish on TSA agar plates
#' and subset the isolates that are in these pairs.
library(tidyverse)

isolates_ID_match <- read_csv("~/Dropbox/lab/emergent-coexistence/data/temp/isolates_ID_match.csv", col_types = cols())

# Ambiguous pairs
pairs_ambiguous <- read_csv("~/Dropbox/lab/emergent-coexistence/data/temp/pairs_competition.csv", col_types = cols()) %>%
    filter(is.na(ColonyCount1)) %>%
    select(-starts_with("ColonyCount"), -DilutionFactor)

# Ambiguous isolates
isolates_ambiguous <- pairs_ambiguous %>%
    select(Community, Isolate1, Isolate2) %>%
    pivot_longer(-Community, names_to = "temp", values_to = "Isolate") %>%
    select(-temp) %>%
    distinct(Community, Isolate) %>%
    arrange(Community, Isolate) %>%
    left_join(isolates_ID_match,  by = c("Community", "Isolate")) %>%
    # Remove duplicate across batch
    distinct(ExpID, .keep_all = T) %>%
    select(Assembly, Community, Isolate, ExpID, ID)


write_csv(pairs_ambiguous, "~/Dropbox/lab/emergent-coexistence/data/temp/pairs_ambiguous.csv")
write_csv(isolates_ambiguous, "~/Dropbox/lab/emergent-coexistence/data/temp/isolates_ambiguous.csv")
