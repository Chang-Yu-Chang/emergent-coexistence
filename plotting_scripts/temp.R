# Scripts for figures
library(tidyverse)
library(cowplot)
library(officer)
library(tidygraph)
library(ggraph)
library(gridExtra)
library(officer)
library(flextable)
source(here::here("plotting_scripts/network_functions.R"))

sequences_abundance <- read_csv(here::here("data/temp/sequences_abundance.csv"), col_types = cols())
communities <- read_csv(here::here("data/output/communities.csv"), col_types = cols())
isolates <- read_csv(here::here("data/output/isolates.csv"), col_types = cols())
pairs <- read_csv(here::here("data/output/pairs.csv"), col_types = cols()) %>% mutate(InteractionType = ifelse(InteractionType == "neutrality", "coexistence", InteractionType))
pairs_freq <- read_csv(here::here("data/output/pairs_freq.csv"), col_types = cols())
load("~/Dropbox/lab/invasion-network/data/output/network.Rdata")
load("~/Dropbox/lab/invasion-network/data/output/network_randomized.Rdata")
community_factor <- c(communities %>% filter(str_detect(Community, "C\\d")) %>% arrange(CommunitySize) %>% pull(Community),
                      communities %>% filter(str_detect(Community, "Ass")) %>% pull(Community))
communities_size <- communities %>% mutate(Community = factor(Community, community_factor)) %>% arrange(Community) %>% pull(CommunitySize)


pairs_freq %>%
    filter(str_detect(Community, "C\\d+R\\d+")) %>%
    filter(Time == "T8") %>%
    group_by(Community, Isolate1, Isolate2) %>%
    summarize(RawDataType = ifelse(any(RawDataType == "Sanger"), "Sanger", "CFU")) %>%
    group_by(RawDataType) %>%
    count()
