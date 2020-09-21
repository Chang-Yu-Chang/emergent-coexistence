# Generate initial plate condition for reconstituted pairs of top-down assembled communities

library(data.table)
library(tidyverse)

df_top_down <- list.files("../data/raw/simulation/", pattern = "top_down", full.names = T) %>% 
    lapply(fread) %>% 
    bind_rows()

df_top_down <- fread("../data/raw/simulation/community-top_down-1_composition.txt")
#count_richness(df_top_down) 

df_top_down_comp <- df_top_down %>% filter(Transfer == max(Transfer), Type == "consumer")
community_list <- df_top_down_comp %>%
    group_by(Well) %>% 
    mutate(SumAbundance = sum(Abundance), RelativeAbundance = Abundance/SumAbundance) %>% 
    filter(RelativeAbundance >= 0.01) %>%  # Isolate with more than 1% relative abundance
    summarize(Richness = n()) %>% 
    filter(Richness >= 3) # Communities with at least three

temp_list <- rep(list(NA), nrow(community_list))
names(temp_list) <- paste0(unique(df_top_down_comp$exp_id), "-", community_list$Well)

for (i in 1:nrow(community_list)) {
    temp_community_ID <- community_list$Well[i]
    temp_consumer <- df_top_down_comp %>% filter(Well == temp_community_ID)
    temp_resource <- df_top_down %>% filter(Well == temp_community_ID, Transfer == 0, Type %in% c("resource", "R0"))
    temp_exp_id <- unique(temp_consumer$exp_id)
    n_pairs <- choose(length(temp_consumer$ID), 2)
    
    df_consumer <- combn(temp_consumer$ID, 2) %>% 
        as_tibble() %>% 
        setnames(paste0("W", 1:n_pairs)) %>% 
        pivot_longer(cols = starts_with("W"), names_to = "Well", values_to = "ID") %>% 
        arrange(Well) %>%
        mutate(Transfer = 0, Type = "consumer", Abundance = 0.5) %>% 
        mutate(exp_id = paste0(temp_exp_id, "-", temp_community_ID, "-reconstituted_pairs")) %>% 
        select(exp_id, Transfer, Type, ID, Well, Abundance) 
        
    df_resource <- temp_resource[rep(1:2, each = n_pairs),] %>%
        group_by(Type) %>% 
        mutate(Well = paste0("W", 1:n_pairs)) %>% 
        ungroup() %>% 
        mutate(exp_id = paste0(temp_exp_id, "-", temp_community_ID, "-reconstituted_pairs"))
    
        
    temp_list[[i]] <- bind_rows(df_consumer, df_resource)
        
}

for (i in 1:length(temp_list)) fwrite(temp_list[[i]], file = paste0("../data/raw/simulation/reconstituted_pairs/", unique(temp_list[[i]]$exp_id), ".txt"))

 # Have to add tin the csv, the n_wells corrspond to the number of pairs tested
# Write the csv below
#community_list






df_mono <- fread("../../community-selection/Data/Raw/f1_additive-monoculture-26_function.txt")

df_mono %>% filter(Transfer == 40) %>%
    filter(Biomass >= 1) %>% 
    ggplot(aes(x = Biomass)) +
    geom_histogram()
    









