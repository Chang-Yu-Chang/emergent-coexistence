#' This scripts generate the figures for community trios and pairs

library(tidyverse)
library(cowplot)
source(here::here("analysis/00-metadata.R"))
source(here::here("simulation/01-generate_input.R"))



# 0. parameters ----
input_parameters <- read_csv(here::here("simulation/01-input_parameters.csv"), col_types = cols())
input_poolPairs <- read_csv(here::here("simulation/03a-input_poolPairs.csv"), col_types = cols())
input_withinCommunityPairs <- read_csv(here::here("simulation/03b-input_withinCommunityPairs.csv"), col_types = cols())
