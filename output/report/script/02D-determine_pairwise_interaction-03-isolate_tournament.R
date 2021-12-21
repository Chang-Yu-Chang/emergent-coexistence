#' Tournament ranks
library(tidyverse)
source(here::here("plotting_scripts/network_functions.R"))

isolates <- read_csv(here::here("data/output/isolates.csv"))
pairs <- read_csv(here::here("data/output/pairs.csv"))
communities <- read_csv(here::here("data/output/communities.csv")) %>% filter(str_detect(Community, "C\\d"))

isolates_tournament <- communities %>%
    select(comm = Community, everything()) %>%
    rowwise() %>%
    mutate(pairs_comm = pairs %>% filter(Community == comm) %>% list()) %>%
    mutate(tournaments_comm = pairs_comm %>% tournament_rank() %>% list()) %>%
    select(Community = comm, tournaments_comm) %>%
    unnest(cols = tournaments_comm)

write_csv(isolates_tournament, file = here::here("data/temp/isolates_tournament.csv"))
