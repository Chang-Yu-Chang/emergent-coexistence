#' Generate isolate and community metadata files
#' Two sets of isolates used: self-assembled community isolates and Jean's natural isolates
library(tidyverse)

# Isolates from self-assembled communities
isolates_ID_match <- read_csv("~/Dropbox/lab/emergent-coexistence/data/raw/pairwise_competition/isolates1.csv", col_types = cols()) %>%
    mutate(Assembly = "self_assembly") %>%
    select(Assembly, ExpID, ID, Community, Isolate)

# Communities
communities_name <- c("C1R2", "C1R4", "C1R6", "C1R7", "C2R6", "C2R8", "C4R1", "C7R1", "C8R4", "C10R2", "C11R1", "C11R2", "C11R5")
communities_size <- c(4,5,5,7,4,4,3,4,3,3,9,12,5)
pp <- function(x) choose(x, 2)

communities <- data.frame(
    Community = communities_name,
    CommunitySize = communities_size,
    CommunityPairSize = pp(communities_size)
) %>%
    mutate(Community = factor(Community, communities_name))  %>%
    arrange(CommunitySize) %>%
    mutate(CommunityLabel = 1:13) %>%
    select(Community, CommunityLabel, everything())

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
    left_join(communities, by = "Community") %>%
    arrange(CommunityLabel, Isolate1, Isolate2) %>%
    mutate(PairID = 1:n()) %>%
    select(PairID, Community, Isolate1, Isolate2)



#
write_csv(isolates_ID_match, "~/Dropbox/lab/emergent-coexistence/data/temp/isolates_ID_match.csv")
write_csv(pairs_ID, "~/Dropbox/lab/emergent-coexistence/data/temp/pairs_ID.csv")
write_csv(communities, "~/Dropbox/lab/emergent-coexistence/data/output/communities.csv")











