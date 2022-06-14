#' Tournament ranks
library(tidyverse)
source(here::here("plotting_scripts/misc.R"))

isolates_ID_match <- read_csv("~/Dropbox/lab/emergent-coexistence/data/temp/isolates_ID_match.csv", col_types = cols())
pairs_interaction <- read_csv("~/Dropbox/lab/emergent-coexistence/data/output/pairs_interaction.csv", col_types = cols()) %>%
    filter(Set == "CFUandCASEU")
communities <- read_csv("~/Dropbox/lab/emergent-coexistence/data/output/communities.csv", col_types = cols())

isolates_tournament <- communities %>%
    select(comm = Community, everything()) %>%
    rowwise() %>%
    mutate(pairs_comm = pairs_interaction %>% filter(Community == comm) %>% list()) %>%
    mutate(tournaments_comm = pairs_comm %>% tournament_rank() %>% list()) %>%
    select(Community = comm, tournaments_comm) %>%
    unnest(cols = tournaments_comm)

write_csv(isolates_tournament, "~/Dropbox/lab/emergent-coexistence/data/temp/isolates_tournament.csv")
