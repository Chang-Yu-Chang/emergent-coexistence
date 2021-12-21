#' Relative abundances within pairs at T8
#' relative abundance

library(tidyverse)
library(cowplot)
isolates <- read_csv(here::here("data/output/isolates.csv")) %>% filter(Assembly == "self_assembly")
pairs <- read_csv(here::here("data/output/pairs.csv"))
pairs_freq <- read_csv(here::here("data/output/pairs_freq.csv")) %>% mutate(Community_pairs = Community) %>% select(-Community)
