library(tidyverse)
isolates_ID_match <- read_csv("output/check/temp/isolates_ID_match.csv")
pairs_ID <- read_csv("output/check/temp/pairs_ID.csv")

# Read the pair ID and isolate ID
pairs_ID_isolate <- pairs_ID %>%
    mutate(Community = factor(Community, paste0("C", rep(1:12, each = 8), "R", rep(1:8, 12)))) %>%
    rename(ID1 = s1, ID2 = s2) %>%
    left_join(select(isolates_ID_match, ID1 = ID, Isolate1 = CYC_isolate)) %>%
    left_join(select(isolates_ID_match, ID2 = ID, Isolate2 = CYC_isolate)) %>%
    arrange(Community, Isolate1, Isolate2)

# Switch the isolate's ID columns such that isolate1 always has smaller number than isolate2
temp_index <- pairs_ID_isolate$Isolate1 > pairs_ID_isolate$Isolate2

temp <- pairs_ID_isolate$ID1[temp_index]
pairs_ID_isolate$ID1[temp_index] <- pairs_ID_isolate$ID2[temp_index]
pairs_ID_isolate$ID2[temp_index] <- temp

temp <- pairs_ID_isolate$Isolate1[temp_index]
pairs_ID_isolate$Isolate1[temp_index] <- pairs_ID_isolate$Isolate2[temp_index]
pairs_ID_isolate$Isolate2[temp_index] <- temp

pairs_ID_isolate <- pairs_ID_isolate %>% arrange(Community, Isolate1, Isolate2)

## Dont Uncomment this trunk as the csv is later manually edited
# pairs_ID_isolate %>%
#     mutate(Pair = 1:n()) %>%
#     select(Pair, everything()) %>%
#     write_csv("output/check/temp/pairs_ID_isolate.csv")



#pairs_ID_isolate %>% filter(Community == "C11R5")


isolates_ID_match %>%
    filter(Community == "C1R2")

library(EBImage)
img = readImage(f)
