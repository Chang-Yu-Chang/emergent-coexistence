# Supplementary figures

library(tidyverse)
library(tidymodels)
library(cowplot)
library(ggraph)
library(tidygraph)
library(ggsci)
library(ggprism)
source(here::here("plotting_scripts/network_functions.R"))

sequences_abundance <- read_csv(here::here("data/temp/sequences_abundance.csv"))
communities <- read_csv(here::here("data/output/communities.csv"))
pairs <- read_csv(here::here("data/output/pairs.csv")) %>% mutate(InteractionType = ifelse(InteractionType == "neutrality", "coexistence", InteractionType))
pairs_example_outcomes <- read_csv(here::here("data/output/pairs_example_outcomes.csv"))
pairs_example_outcomes_finer <- read_csv(here::here("data/output/pairs_example_outcomes_finer.csv"))
pairs_freq <- read_csv(here::here("data/output/pairs_freq.csv"))
networks_motif <- read_csv(here::here("data/output/networks_motif.csv"))
networks_motif_randomized <- read_csv(here::here("data/output/networks_motif_randomized.csv"))
networks_motif_randomized_percentile <- read_csv(here::here("data/output/networks_motif_randomized_percentile.csv")) %>% mutate(Community = factor(Community, communities$Community)) %>% arrange(Community)
load(here::here("data/output/network_community.Rdata"))






















