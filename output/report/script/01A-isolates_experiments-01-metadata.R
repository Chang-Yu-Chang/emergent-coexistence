#' Generate isolate and community metadata files
#' Two sets of isolates used: self-assembled community isolates and Jean's natural isolates
library(tidyverse)

# Self-assembled communities
isolates_ID_match <- read_csv(here::here("data/raw/pairwise_competition/isolates1.csv")) %>%
    mutate(Assembly = "self_assembly") %>%
    select(Assembly, ExpID, ID, Community, Isolate, Family, Genus)

# Random communities. Pick random isolates for experiment
jean_isolates <- read_csv(here::here("data/raw/pairwise_competition/jean_isolates.csv")) %>%
    select(ExpID, ID, Family, Genus)
## Randomly pick isolates. 16 isolates from self-assembled community isolates. 8 isolates in each of the new assembled communities
## Randomly pick 8 from 13 communities, and pick one isolate in each picked community. Two communities
temp_df <- data.frame(ExpID = c("1.4.A.1", "1.6.A.5", "2.6.A.3", "2.8.A.1", "8.4.A.2", "10.2.A.1", "11.1.B.1", "11.2.B.2"), stringsAsFactors = F)
temp_df2 <- data.frame(ExpID = c("1.6.A.6", "1.7.A.1.2", "2.8.A.3", "4.1.A.2", "8.4.A.1", "10.2.A.2", "11.1.B.3", "1.2.A.5"), stringsAsFactors = F)
## Randomly pick another 16 isolates from Jean's isolates. Two communities
jean_isolates <- jean_isolates %>%
    filter(Genus %in% c("Enterobacter", "Pseudomonas", "Pantoea", "Citrobacter", "Raoultella", "Kluyvera")) %>%
    mutate(ExpID = ordered(ExpID, paste0("JVN", 1:100)))
temp_df3 <-  jean_isolates %>%
    filter(ExpID %in% paste0("JVN", c(1,2,7,8,21,25,30, 46))) %>% # 30 and 46 are Pseudomonas
    arrange(ExpID) %>% as_tibble()
temp_df4 <- jean_isolates %>%
    filter(ExpID %in% paste0("JVN", c(3, 11, 12, 14, 23, 28, 35, 41))) %>% # 14, 23, 35 are Pseudomonas
    arrange(ExpID) %>% as_tibble()
## bind_row the new isolates list for each community assembly
isolates_assembly_across <-
    temp_df %>%
    bind_rows(temp_df2) %>%
    left_join(isolates_ID_match, by = c("ExpID")) %>%
    mutate(Assembly = "across_community",
           Community = rep(paste0("AcrAss", 1:2), each = 8)) %>%
    mutate(Isolate = rep(1:8, 2)) %>%
    select(Assembly, ExpID, ID, Community, Isolate, Family, Genus)
isolates_assembly_random <-
    temp_df3 %>%
    bind_rows(temp_df4) %>%
    mutate(Assembly = rep("random_assembly", 16),
           Community = rep(paste0("RanAss", 1:2), each = 8),
           Isolate = rep(1:8, 2)) %>%
    select(Assembly, ExpID, ID, Community, Isolate, Family, Genus)

isolates_random <- isolates_assembly_across %>%
    bind_rows(isolates_assembly_random) %>%
    as_tibble()

isolates_ID_match <- isolates_ID_match %>% bind_rows(isolates_random)
write_csv(isolates_ID_match, here::here("data/temp/isolates_ID_match.csv"))

# Communities
communities_name <- c("C1R2", "C1R4", "C1R6", "C1R7", "C2R6", "C2R8", "C4R1", "C7R1", "C8R4", "C10R2", "C11R1", "C11R2", "C11R5", "AcrAss1", "AcrAss2", "RanAss1", "RanAss2")
communities_size <- c(4,5,5,7,4,4,3,4,3,3,9,12,5, 8, 8, 8, 8)
pp <- function(x) choose(x, 2)
tt <- function(x) choose(x, 3)

communities <- data.frame(
    Community = communities_name,
    CommunitySize = communities_size,
    CommunityPairSize = pp(communities_size),
    CommunitiyMotifSize = tt(communities_size)
)

# Pairs
pairs_ID <- communities %>%
    select(Community, CommunitySize) %>%
    split.data.frame(f = .$Community) %>%
    lapply(function(x) {
        combn(x$CommunitySize, 2) %>% t() %>%
            as_tibble() %>%
            setNames(c("Isolate1", "Isolate2"))
    }) %>%
    bind_rows(.id = "Community") %>%


#
write_csv(pairs_ID, here::here("data/temp/pairs_ID.csv"))
write_csv(communities, here::here("data/output/communities.csv"))











