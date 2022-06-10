#' Count the network components
library(tidyverse)
library(tidygraph)
library(igraph)
source(here::here("plotting_scripts/misc.R"))
isolates <- read_csv("~/Dropbox/lab/emergent-coexistence/data/output/isolates.csv", col_types = cols())
pairs <- read_csv("~/Dropbox/lab/emergent-coexistence/data/output/pairs.csv", col_types = cols())
communities <- read_csv("~/Dropbox/lab/emergent-coexistence/data/output/communities.csv", col_types = cols())
load("~/Dropbox/lab/emergent-coexistence/data/output/communities_network.Rdata")
load("~/Dropbox/lab/emergent-coexistence/data/output/communities_network_randomized.Rdata")

# Observation
communities_component <- communities_network %>%
    select(Community, Network) %>%
    rowwise() %>%
    mutate(Component = count_component(Network)) %>%
    unnest(Component) %>%
    select(Community, Component)


# Randomized networks
communities_component_randomized <- communities_network_randomized %>%
    select(Community, NetworkRandomized) %>%
    unnest_longer(col = NetworkRandomized, indices_to = "Replicate", values_to = "Network") %>%
    rowwise() %>%
    mutate(Component = count_component(Network)) %>%
    unnest(Component) %>%
    select(Community, Replicate, Component)

if (FALSE) {
networks_component <- net_list[1:13] %>%
    lapply(count_component) %>%
    bind_rows(.id = "Community")

networks_component_randomized <- rep(list(NA), 13)
for (i in 1:length(networks_component_randomized)) {
    networks_component_randomized[[i]] <- net_randomized_list[[i]] %>%
        lapply(count_component) %>%
        bind_rows(.id = "Replicate")
    cat(i)
}

networks_component_randomized <- networks_component_randomized %>%
    setNames(names(net_randomized_list)) %>%
    bind_rows(.id = "Community")

}



# Save the result
write_csv(communities_component, "~/Dropbox/lab/emergent-coexistence/data/output/communities_component.csv")
write_csv(communities_component_randomized, "~/Dropbox/lab/emergent-coexistence/data/output/communities_component_randomized.csv")














