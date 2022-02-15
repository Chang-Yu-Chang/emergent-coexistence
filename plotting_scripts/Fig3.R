# Simulation. Random trios in species pool

library(cowplot)
library(tidygraph)
library(ggraph)
library(tidyverse)
library(data.table)
source("network_functions.R")
#source("analysis-pair-culturable_random.R")
#source("analysis-trio-culturable_random.R")
#source("analysis-pair-from_top_down_community.R")
whether_save_temp_file <- T

# input_independent <- fread("../data/raw/simulation/mapping_files/input_independent.csv")
# input_independent_trios <- input_independent %>% filter(grepl("trio-culturable_isolates", exp_id)) %>% filter(seed == 2)
# interaction_color <- assign_interaction_color()

