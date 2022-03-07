#' Count the node degree
library(tidyverse)
library(tidygraph)
library(igraph)
source(here::here("plotting_scripts/network_functions.R"))
isolates <- read_csv(here::here("data/output/isolates.csv"))
pairs <- read_csv(here::here("data/output/pairs.csv"))
communities <- read_csv(here::here("data/output/communities.csv"))
load("~/Dropbox/lab/invasion-network/data/output/network_community.Rdata")
load("~/Dropbox/lab/invasion-network/data/output/network_randomized.Rdata") # Randomized networks


networks_degree <- net_list[1:13] %>%
    lapply(function(x) {tibble(Degree = count_degree(x)) %>% mutate(Isolate = 1:n())}) %>%
    bind_rows(.id = "Community")

networks_degree_randomized <- rep(list(NA), 13)
for (i in 1:length(networks_degree_randomized)) {
    networks_degree_randomized[[i]] <- net_randomized_list[[i]] %>%
        lapply(function(x) {tibble(Degree = count_degree(x)) %>% mutate(Isolate = 1:n())}) %>%
        bind_rows(.id = "Replicate")
    cat(i)
}

networks_degree_randomized <- networks_degree_randomized %>%
    setNames(names(net_randomized_list)) %>%
    bind_rows(.id = "Community")


# Save the result
write_csv(networks_degree, file = here::here("data/output/networks_degree.csv"))
write_csv(networks_degree_randomized, file = here::here("data/output/networks_degree_randomized.csv"))















