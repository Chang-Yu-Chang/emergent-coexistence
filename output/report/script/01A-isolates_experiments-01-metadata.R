#' Generate isolate and community metadata files

# Isolates
isolates_ID_match <- fread(here::here("data/raw/pairwise_competition", "isolates1.csv")) %>%
    select(ExpID, ID, Community, Isolate)
fwrite(isolates_ID_match, here::here("data/temp/isolates_ID_match.csv"))

# Communities
communities_name <- c("C1R2", "C1R4", "C1R6", "C1R7", "C2R6", "C2R8", "C4R1", "C7R1", "C8R4", "C10R2", "C11R1", "C11R2", "C11R5")
communities_size <- c(4,5,5,7,4,4,3,4,3,3,9,12,5)
pp <- function(x) choose(x, 2)
tt <- function(x) choose(x, 3)

communities <- data.frame(
    Community = communities_name,
    CommunitySize = communities_size,
    CommunityPairSize = pp(communities_size),
    CommunitiyMotifSize = tt(communities_size)
)

fwrite(communities, here::here("data/output/communities.csv"))
