#' Subset the example pairs

library(tidyverse)
communities <- read_csv("~/Dropbox/lab/emergent-coexistence/data/output/communities.csv", col_types = cols())
pairs_freq <- read_csv("~/Dropbox/lab/emergent-coexistence/data/output/pairs_freq.csv", col_types = cols())
pairs <- read_csv("~/Dropbox/lab/emergent-coexistence/data/output/pairs.csv", col_types = cols())

# Example for plotting interaction types ----
## Stable coexistence and coexistence ain examples
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
write_csv(pairs_example_outcomes, "~/Dropbox/lab/emergent-coexistence/data/output/pairs_example_outcomes.csv")
write_csv(pairs_example_outcomes_finer, "~/Dropbox/lab/emergent-coexistence/data/output/pairs_example_outcomes_finer.csv")
