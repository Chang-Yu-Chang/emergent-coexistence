#' Count the node degree
library(tidyverse)
library(tidygraph)
library(igraph)
source(here::here("plotting_scripts/misc.R"))
communities <- read_csv("~/Dropbox/lab/emergent-coexistence/data/output/communities.csv", col_types = cols())
load("~/Dropbox/lab/emergent-coexistence/data/output/communities_network.Rdata")
load("~/Dropbox/lab/emergent-coexistence/data/output/communities_network_randomized.Rdata")

communities_degree <- communities_network %>%
    rowwise() %>%
    mutate(Degree = count_degree(Network) %>% list) %>%
    unnest_longer(col = Degree, indices_to = "Isolate", values = "Count") %>%
    select(Community, Isolate, Count)


communities_degree_randomized <- communities_network_randomized %>%
    select(Community, NetworkRandomized) %>%
    unnest_longer(col = NetworkRandomized, indices_to = "Replicate", values_to = "Network") %>%
    rowwise() %>%
    mutate(Degree = count_degree(Network) %>% list) %>%
    unnest_longer(col = Degree, indices_to = "Isolate", values = "Count") %>%
    select(Community, Isolate, Count)

if (FALSE) {
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

}


# Save the result
write_csv(communities_degree, "~/Dropbox/lab/emergent-coexistence/data/output/communities_degree.csv")
write_csv(communities_degree_randomized, "~/Dropbox/lab/emergent-coexistence/data/output/communities_degree_randomized.csv")















