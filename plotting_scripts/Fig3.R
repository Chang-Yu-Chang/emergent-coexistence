# Plot simulation data


library(tidyverse)
library(data.table)

# Pairs with fixed initial frequencies

df_comp <- list.files("../data/raw/simulation/", full.names = T) %>%
    lapply(fread) %>%
    rbindlist()


df_comp %>%
    filter(Type == "consumer") %>% 
    mutate(ID = factor(ID)) %>% 
    filter(Transfer >= 10)
#    filter(Seed == 1, Well %in% paste0("W", 0:9)) %>% 
    ggplot() +
    geom_bar(aes(x = Transfer, y = Abundance, fill = ID), position = "stack", stat = "identity") +
    facet_wrap(Well~.) +
    guides(fill = F) +
    theme_bw()
