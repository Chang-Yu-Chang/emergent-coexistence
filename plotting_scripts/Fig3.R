# Plot simulation data


library(tidyverse)
library(data.table)

# Pairs with fixed initial frequencies

df_comp <- list.files("../data/raw/simulation/", full.names = T) %>%
    lapply(fread) %>%
    rbindlist()

#
plot_temporal_comp <- function(df_comp) {
    df_comp %>% 
    filter(Type == "consumer") %>% 
    mutate(ID = factor(ID)) %>% 
    filter(Well %in% paste0("W", 0:19)) %>%
    ggplot() +
    geom_bar(aes(x = Transfer, y = Abundance, fill = ID), position = "stack", stat = "identity") +
    facet_wrap(Well~.) +
    guides(fill = F) +
    theme_bw()
}
plot_temporal_resource <- function(df_comp) {
    df_comp %>% 
        filter(Type == "resource") %>% 
        mutate(ID = factor(ID)) %>% 
        filter(Well %in% paste0("W", 0:19)) %>%
        ggplot() +
        geom_bar(aes(x = Transfer, y = Abundance, fill = ID), position = "stack", stat = "identity") +
        facet_wrap(Well~.) +
        guides(fill = F) +
        theme_bw()
}
count_richness <- function(df_comp){
    df_comp %>%
        filter(Type == "consumer") %>% 
        mutate(ID = factor(ID)) %>% 
        filter(Transfer == max(Transfer)) %>% 
        filter(Abundance >= 1e-6) %>% 
        group_by(Well) %>%
        summarize(Richness = n()) %>%
        pull(Richness) %>% table()
}

df_comp <- fread("../data/raw/simulation/monoculture-random_isolates-1_composition.txt")
df_comp <- fread("../data/raw/simulation/pair-random_isolates-1_composition.txt")

plot_temporal(df_comp)
count_richness(df_comp)


plot_temporal_resource(df_comp)

df_comp %>%
    filter(Well == "W6") %>%
    filter(Type == "resource")
    
        