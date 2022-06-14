# CFU and CASEU comparison
library(tidyverse)
library(cowplot)
library(gridExtra)
source(here::here("plotting_scripts/misc.R"))

isolates <- read_csv("~/Dropbox/lab/emergent-coexistence/data/output/isolates.csv", col_types = cols())
pairs <- read_csv("~/Dropbox/lab/emergent-coexistence/data/output/pairs.csv", col_types = cols()) %>% mutate(InteractionType = ifelse(InteractionType == "neutrality", "coexistence", InteractionType))
pairs_freq <- read_csv("~/Dropbox/lab/emergent-coexistence/data/output/pairs_freq.csv", col_types = cols()) %>% mutate(Time = factor(Time, c("Tini", "Tend")))
pairs_example_outcomes_finer <- read_csv("~/Dropbox/lab/emergent-coexistence/data/output/pairs_example_outcomes_finer.csv", col_types = cols())
communities <- read_csv("~/Dropbox/lab/emergent-coexistence/data/output/communities.csv", col_types = cols())
community_factor <- communities$Community
communities_size <- communities$CommunitySize
communities_hierarchy <- read_csv("~/Dropbox/lab/emergent-coexistence/data/output/communities_hierarchy.csv", col_types = cols())
load("~/Dropbox/lab/emergent-coexistence/data/output/communities_network.Rdata") # communities_network


pairs_interaction <- read_csv("~/Dropbox/lab/emergent-coexistence/data/output/pairs_interaction.csv", col_types = cols()) %>%
    mutate(InteractionType = ifelse(InteractionType == "neutrality", "coexistence", InteractionType))
pairs_ID <- read_csv("~/Dropbox/lab/emergent-coexistence/data/temp/pairs_ID.csv", col_types = cols())


pairs_interaction %>%
    mutate(Set = factor(Set, c("CFUonly", "CASEUonly", "CFUandCASEU"))) %>%
    select(Set, PairID, InteractionType) %>%
    pivot_wider(names_from = Set, values_from = InteractionType) %>%
    mutate(ResultMatch = CASEUonly == CFUandCASEU) %>%
    group_by(ResultMatch) %>%
    summarize(Count = n())

