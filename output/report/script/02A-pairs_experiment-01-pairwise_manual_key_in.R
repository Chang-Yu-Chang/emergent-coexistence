#' Generate the empty csv table for manual key-in
library(tidyverse)
library(data.table)

# This function creates an empty csv file for later manual key-in
interaction_table_make <- function (experiment, community, community_size) {
  stopifnot(length(community) == length(community_size))

  # Make an empty list
  table_all <- rep(list(NA), length(community))

  # Loop through communities
  for (i in 1:length(community_size)) {
    table_all[[i]] <-
      cbind(rep(community[i], choose(community_size[i], 2) * 3),
        rbind(t(combn(community_size[i], 2)),
          t(combn(community_size[i], 2)),
          t(combn(community_size[i], 2)))) %>%
      as.data.frame() %>%
      setNames(c("Community", "Isolate1", "Isolate2")) %>%
      arrange(Community, Isolate1, Isolate2) %>%
      mutate(Isolate1Freq = rep(c(5, 50, 95), choose(community_size[i], 2)),
        Isolate2Freq = rep(c(95, 50, 5), choose(community_size[i], 2)))

  }

  #
  rbindlist(table_all) %>%
    arrange(Community, Isolate1, Isolate2, Isolate1Freq, Isolate2Freq) %>%
    mutate(Experiment = experiment) %>%
    return()
}


# Generate the empty table in csv file
temp <- rbind(
  interaction_table_make("Transitivity_B2",c("C11R1", "C2R6", "C2R8", "C7R1", "C8R4", "C10R2"), c(9, 4, 4, 4, 3, 3)),
  interaction_table_make("Transitivity_C", c("C11R1"), 9),
  interaction_table_make("Transitivity_C2", c("C11R2"), c(13)),
  interaction_table_make("Transitivity_D",c("C1R7", "C11R5", "C1R4", "C1R6", "C1R2", "C4R1"), c(7, 5, 5, 5, 4, 3))) %>%
  filter(!(Experiment == "Transitivity_B2" & Community == "C11R1" & (Isolate1 == 1 | Isolate2 == 1))) %>% # Remove C11R1 isolate 1 in batch B2
 # filter(!(Experiment == "Transitivity_C" & (Isolate1 != 1 & Isolate2 != 1))) %>% # Add C11R1 isolate 1 in batch C
  mutate(Isolate1 = ordered(Isolate1, 1:13), Isolate2 = ordered(Isolate2, 1:13)) %>%
  mutate(Isolate1 = ordered(Isolate1, 1:13), Isolate2 = ordered(Isolate2, 1:13)) %>%
  mutate(Community = factor(Community, c("C1R2", "C1R4", "C1R6", "C1R7", "C2R6", "C2R8", "C4R1", "C7R1", "C8R4", "C10R2", "C11R1", "C11R2", "C11R5"))) %>%
  arrange(Community, Isolate1Freq, Isolate2Freq, Isolate1, Isolate2) %>%
  mutate(ColonyCount=NA, ColonyCount1=NA, ColonyCount2=NA
  )

if (FALSE) fwrite(temp, root$find_file("data/raw/pairwise_competition/result_pairwise_competition_keyin.csv"), row.names = F)

#' Manual key-in
#' Use the code below for competition result.
"
Code | Result
-|-
  A positve integer | winner isolates
`0` | coexistence
`-1`| manual error
`-2`| similar morphology
`-3`| need to re-plate; low colony number; no growth
"
