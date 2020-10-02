# Plot simulation data

library(data.table)
library(tidyverse)

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
make_edge_list <- function(df_comp_pairs) {
    df_pairs_T0 <- df_comp_pairs %>%
        filter(Transfer == 0) %>% 
        filter(Type == "consumer") %>%
        group_by(Well) %>% 
        mutate(Species = paste0("S", 1:2))
    
    list_pairs <- df_pairs_T0 %>% 
        select(exp_id, ID, Well, Species)
    
    df_pairs_T10 <- df_comp_pairs %>%
        filter(Transfer == max(Transfer)) %>% 
        filter(Type == "consumer") %>%
        group_by(Well) %>%
        right_join(mutate(list_pairs, Transfer = max(.$Transfer), Type = "consumer"),
            by = c("exp_id", "Transfer", "Type", "ID", "Well")) %>% 
        replace_na(list(Abundance = 0)) %>% 
        arrange(Well)
    
    df_edge_list <- df_pairs_T10 %>% 
        pivot_wider(names_from = Species, values_from = c(ID, Abundance)) 
    
    # Determine the competition outcomes
    df_edge_list$InteractionType <- NA
    df_edge_list$from <- NA
    df_edge_list$to <- NA
    df_edge_list$from_abundance <- NA
    df_edge_list$to_abundance <- NA
    
    for (i in 1:nrow(df_edge_list)) {
        ID_S1 <- df_edge_list$ID_S1[i]
        ID_S2 <- df_edge_list$ID_S2[i]
        abundance_S1 <- df_edge_list$Abundance_S1[i]
        abundance_S2 <- df_edge_list$Abundance_S2[i]
        
        if (df_edge_list$Abundance_S1[i] <= 0 & df_edge_list$Abundance_S2[i] <= 0) {
            df_edge_list$InteractionType[i] <- "no-growth"
        } else if (df_edge_list$Abundance_S1[i] > 0 & df_edge_list$Abundance_S2[i] > 0) {
            df_edge_list$InteractionType[i] <- "coexistence"
            df_edge_list$from[i] <- ID_S1
            df_edge_list$to[i] <- ID_S2
            df_edge_list$from_abundance[i] <- abundance_S1
            df_edge_list$to_abundance[i] <- abundance_S2
        } else if (df_edge_list$Abundance_S1[i] > 0 & df_edge_list$Abundance_S2[i] <= 0) {
            df_edge_list$InteractionType[i] <- "exclusion"
            df_edge_list$from[i] <- ID_S1
            df_edge_list$to[i] <- ID_S2
            df_edge_list$from_abundance[i] <- abundance_S1
            df_edge_list$to_abundance[i] <- abundance_S2
        } else if (df_edge_list$Abundance_S1[i] <= 0 & df_edge_list$Abundance_S2[i] > 0) {
            df_edge_list$InteractionType[i] <- "exclusion"
            df_edge_list$from[i] <- ID_S2
            df_edge_list$to[i] <- ID_S1
            df_edge_list$from_abundance[i] <- abundance_S2
            df_edge_list$to_abundance[i] <- abundance_S1
            
        }
    }
    
    df_edge_list %>%
        ungroup() %>% 
        select(from, to, InteractionType, from_abundance, to_abundance) %>% 
        #filter(InteractionType != "no-growth") %>% 
        arrange(from, to) %>%
        mutate(exp_id = unique(df_comp_pairs$exp_id)) %>% 
        select(exp_id, everything())
    
}
count_pairs <- function(edge_list) {
    edge_list %>%
        group_by(exp_id, InteractionType) %>%
        summarize(Count = n())
}
plot_pairs_leakage <- function(df) {
    df %>%
        ggplot() +
        geom_bar(aes(x = Leakage, y = Count, fill = InteractionType), position = "fill", stat = "identity") +
        scale_x_continuous(breaks = seq(0, 0.9, 0.1), expand = c(0,0)) +
        scale_y_continuous(expand = c(0,0)) +
        theme_bw() +
        theme(legend.position = "top", legend.title = element_blank()) +
        labs(x = "Leakage level", y = "Fraction of random pairs") 
}

# Panel A: random pairs
pair_counts_list <- list.files("../data/raw/simulation/", pattern = "leakage", full.names = T) %>% 
    lapply(function (x) {
        fread(x) %>% 
            make_edge_list() %>% 
            count_pairs()
    })

pair_counts <- pair_counts_list %>% 
    bind_rows() %>% 
    separate(exp_id, into = c("temp1", "temp2", "Leakage", "Specialist", "Seed"), sep = "-") %>% 
    unite("Treatment", temp1:temp2) %>% 
    mutate(Leakage = as.numeric(gsub("leakage", "", Leakage))/100) %>% 
    mutate(Specialist = as.numeric(gsub("specialist", "", Specialist))/100)

p1 <- pair_counts %>% 
    plot_pairs_leakage() +
    facet_grid(Specialist~.)

#ggsave("../plots/Fig3A_temp.png", plot = p1, width = 10, height = 10)


# Panel B: reconstituted pairs from the top-down assembled communities
pair_counts_list <- list.files("../data/raw/simulation/", pattern = "reconstituted_pairs_", full.names = T) %>% 
    lapply(function (x) {
        fread(x) %>% 
            make_edge_list() %>% 
            count_pairs()
    })


df_top_down <- fread("../data/raw/simulation/community-top_down-1_composition.txt") %>% 
    filter(Transfer == max(Transfer), Type == "consumer") %>% 
    mutate(Community = Well) %>%
    group_by(Community) %>% 
    summarize(Richness = n()) %>% 
    filter(Richness >= 3)

pair_counts <- pair_counts_list %>% 
    bind_rows() %>% 
    separate(exp_id, into = c("temp1", "temp2", "Seed", "Community", "temp3"), sep = "-") %>% 
    unite("Treatment", starts_with("temp"), sep = "-") %>% 
    left_join(df_top_down)
    
pair_counts <- pair_counts %>%
    group_by(Community) %>%
    mutate(SumCount = sum(Count), FracCount = Count / SumCount)
    
p2 <- pair_counts %>% 
    ggplot() +
    geom_jitter(aes(x = Richness, y = FracCount, color = InteractionType), shape = 21, size = 3, width = 0.05) +
    theme_bw() +
    theme(legend.position = "top", legend.title = element_blank()) +
    labs(x = "Richness", y = "Fraction of pairwise coexistence") +
    ggtitle("l=0.5, q=0.8")

p2
#ggsave("../plots/Fig3B_temp.png", plot = p2, width = 5, height = 5)


# Panel community
df_top_down <- fread("../data/raw/simulation/community-top_down-1_composition.txt") %>% 
    filter(Transfer == max(Transfer), Type == "consumer") %>% 
    mutate(Community = Well) %>%
    group_by(Community) %>% 
    mutate(SumAbundance = sum(Abundance), RelativeAbundance = Abundance/SumAbundance) %>% 
    filter(RelativeAbundance >= 0.01) %>% 
    summarize(Richness = n()) 

plot_temporal_comp(df_top_down)
#    filter(Richness >= 3)




if (FALSE) {
    
    #df_comp <- fread("../data/raw/simulation/monoculture-random_isolates-1_composition.txt")
    df_comp <- fread("../data/raw/simulation/pair-random_isolates-leakage0-1_composition.txt")
    plot_temporal_comp(df_comp)
    count_richness(df_comp)
    
    
    
    df_comp <- list.files("../data/raw/simulation/", pattern = "leakage", full.names = T) %>% 
        lapply(fread)
    
    df_comp %>%
        lapply(count_richness)
    
    #plot_temporal_resource(df_comp)
    
}