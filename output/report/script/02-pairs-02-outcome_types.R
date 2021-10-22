#' Types of pairwise competition outcomes
library(tidyverse)
library(data.table)
communities <- read_csv(here::here("data/output/communities.csv"))

# Read data
pairs_freq <- read_csv(here::here("data/output/pairs_freq.csv"))
pairs <- read_csv(here::here("data/output/pairs.csv"))
pairs_interaction_fitness <- read_csv(here::here("data/temp/pairs_interaction_fitness.csv")) %>% mutate(Community = ordered(Community, communities$Community))

# Example for plotting interaction types ----
## Four main examples; mainly for figure 1
pairs_example_outcomes <- pairs %>%
  group_by(InteractionType) %>%
  arrange(desc(Community), Isolate1) %>%
  slice(1) %>%
  select(Community, Isolate1, Isolate2, InteractionType, InteractionTypeFiner, From, To)


## Examples of outcomes in a finer scale; plot for figure 2C
pairs_example_outcomes_finer <-
  pairs_interaction_fitness %>%
  group_by(InteractionTypeFiner) %>%
  slice(1) %>%
  ungroup() %>%
  select(Community, Isolate1, Isolate2, InteractionType, InteractionTypeFiner, From, To)

#
write_csv(pairs_example_outcomes, here::here("data/output/pairs_example_outcomes.csv"))
write_csv(pairs_example_outcomes_finer, here::here("data/output/pairs_example_outcomes_finer.csv"))
