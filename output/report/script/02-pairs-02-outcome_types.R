#' Types of pairwise competition outcomes
library(tidyverse)
communities <- read_csv(here::here("data/output/communities.csv"))

# Read data
pairs_freq <- read_csv(here::here("data/output/pairs_freq.csv"))
pairs <- read_csv(here::here("data/output/pairs.csv"))
#pairs_interaction_fitness <- read_csv(here::here("data/temp/pairs_interaction_fitness.csv")) %>% mutate(Community = ordered(Community, communities$Community))

# Example for plotting interaction types ----
## Stable coexistence and coexistence ain examples
if (FALSE) {
pairs_example_outcomes <- pairs %>%
    filter(Assembly == "self_assembly") %>%
    group_by(InteractionType) %>%
    arrange(desc(Community), Isolate1) %>%
    slice(1) %>%
    select(Community, Isolate1, Isolate2, InteractionType, InteractionTypeFiner, From, To)

}
pairs_example_outcomes <- pairs %>%
    filter((Community == "C1R2" & Isolate1 == 1 & Isolate2 == 2) |
               (Community == "C1R7" & Isolate1 == 3 & Isolate2 == 6) |
               (Community == "C2R8" & Isolate1 == 1 & Isolate2 == 2))


## Examples of outcomes in a finer scale
pairs_example_outcomes_finer <- pairs %>%
    filter((Community == "C1R2" & Isolate1 == 1 & Isolate2 == 2) |
               (Community == "C1R7" & Isolate1 == 3 & Isolate2 == 6) |
               (Community == "C2R8" & Isolate1 == 1 & Isolate2 == 2) |
               (Community == "C11R1" & Isolate1 == 1 & Isolate2 == 5) |
               (Community == "C2R6" & Isolate1 == 1 & Isolate2 == 4))

#
write_csv(pairs_example_outcomes, here::here("data/output/pairs_example_outcomes.csv"))
write_csv(pairs_example_outcomes_finer, here::here("data/output/pairs_example_outcomes_finer.csv"))
