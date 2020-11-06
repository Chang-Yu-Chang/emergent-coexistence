#' Run randomization of the observed networks 


library(tidygraph)
library(tidyverse)
options(dplyr.summarise.inform = FALSE)
library(data.table)
source("network_functions.R")

# Read data
communities_abundance <- fread("../data/temp/communities_abundance.csv")
communities <- fread("../data/temp/communities.csv")
community_sizes <- communities$CommunitySize
community_names <- communities$Community
community_names_ordered_by_size <- communities %>% arrange(CommunitySize) %>% pull(Community)
isolates <- fread("../data/output/isolates.csv")
isolates_melted <- fread("../data/output/isolates_melted.csv")
pairs <- fread("../data/output/pairs.csv")
#pairs_melted <- fread("../data/output/pairs_melted.csv")


# Make networks 
graph_list <- rep(list(NA), length(community_names))
names(graph_list) <- community_names
for (i in 1:length(community_names)) {
    isolates_subset <- isolates %>% filter(Community == community_names[i])
    pairs_subset <- pairs %>% filter(Community == community_names[i])
    graph_list[[i]] <- make_network(isolates = isolates_subset, pairs = pairs_subset)
}

save(graph_list, file = "../data/temp/graph_list.Rdata")

# Make summary motif counts
observed_motif_count <- graph_list %>%
    lapply(summarize_network_motif) %>%
    bind_rows(.id = "Community")

fwrite(observed_motif_count, file = "../data/temp/observed_motif_counts.txt")

# Randomize networks and count
b = 1000

cat("\n Randomizing the empirical graphs")
temp_list <- rep(list(rep(list(NA), b)), length(graph_list))
names(temp_list) <- names(graph_list)
for (j in 1:length(graph_list)){
    cat("\ngraph:", names(graph_list)[j], "\n")
    for (i in 1:b) {
        temp_list[[j]][[i]] <- count_motif(randomize_network(graph_list[[j]]))
        if (i%%10 == 0) cat(i, " ")
    }
}

random_motif_counts <- temp_list %>%
    lapply(function(x) {
        lapply(x, function(y) {tibble(Motif = factor(1:7), Count = y)}) %>%
            rbindlist(idcol = "Seed")
    }) %>% 
    bind_rows(.id = "Community")

random_motif_counts_p95 <- random_motif_counts %>% 
    group_by(Community, Motif) %>% 
    filter(Count >= quantile(Count, 0.95)) %>% 
    distinct(Community, Motif, Count) %>% 
    arrange(Community, Motif, Count) %>% 
    slice_min(Count) %>% 
    mutate(Percentile = "p95")
random_motif_counts_p05 <- random_motif_counts %>% 
    group_by(Community, Motif) %>% 
    filter(Count <= quantile(Count, 0.05)) %>% 
    distinct(Community, Motif, Count) %>% 
    arrange(Community, Motif, Count) %>% 
    slice_max(Count) %>% 
    mutate(Percentile = "p05")
random_motif_counts_percentile <- bind_rows(random_motif_counts_p05, random_motif_counts_p95) %>% 
    mutate(Community = ordered(Community, levels = community_names_ordered_by_size))

fwrite(random_motif_counts, file = "../data/temp/random_motif_counts.txt")
fwrite(random_motif_counts_percentile, file = "../data/temp/random_motif_counts_percentile.txt")




if (FALSE) {
    # Join the final fraction of pairs
    
    pairs_freq <- pairs_melted %>%
        unite(col = "InitialFrequency", Isolate1InitialODFreq:Isolate2InitialODFreq, sep = "-") %>% 
        filter(Time == "T8") %>% 
        group_by(Community, Isolate1, Isolate2) %>% 
        summarize(Isolate1Freq = mean(Isolate1MeasuredFreq)) %>% 
        split.data.frame(f = .$Community) %>% 
        lapply(function(x) {
            reversed_edges <- tibble(Community = x$Community, Isolate1 = x$Isolate2, Isolate2 = x$Isolate1, Isolate1Freq = 1-x$Isolate1Freq)
            bind_rows(x, reversed_edges)
        }) %>%
        rbindlist()
    
    pairs_freq_comm <- filter(pairs_freq, Community == "C7R1")
    
    compute_hierarchy_score <- function(pairs_freq_comm) {
        community_size <- length(unique(c(pairs_freq_comm$Isolate1, pairs_freq_comm$Isolate2)))
        competitive_score <- pairs_freq_comm %>% 
            group_by(Isolate1) %>% 
            summarise(Isolate1CompetitiveScore = mean(Isolate1Freq))
        
        pairs_freq_comm %>% 
            left_join(competitive_score, by = "Isolate1") %>% 
            left_join(select(competitive_score, Isolate2 = Isolate1, Isolate2CompetitiveScore = Isolate1CompetitiveScore), by = "Isolate2") %>% 
            filter(Isolate1CompetitiveScore > Isolate2CompetitiveScore) %>%
            summarise(HierarchyScore = sum(Isolate1Freq)/choose(community_size, 2)) %>%
            return()
    }
    
    pairs_freq %>%
        split.data.frame(.$Community) %>% 
        lapply(compute_hierarchy_score) %>% 
        rbindlist(idcol = "Community")
    compute_hierarchy_score()
    
    competitive_score <- pairs_freq %>% 
        group_by(Community, Isolate1) %>% 
        summarise(Isolate1CompetitiveScore = mean(Isolate1Freq))
    
    pairs_freq %>% 
        left_join(competitive_score) %>% 
        left_join(select(competitive_score, Community, Isolate2 = Isolate1, Isolate2CompetitiveScore = Isolate1CompetitiveScore)) %>% 
        filter(Isolate1CompetitiveScore > Isolate2CompetitiveScore) %>%
        group_by(Community) %>% 
        left_join(communities) %>% 
        summarise(HierarchyScore = sum(Isolate1Freq)/CommunitySize^2) %>% view()
    
    pairs <- pairs %>% left_join(pairs_freq)
    
    # 
    graph <- graph_list[[1]]
    
    graph %>% 
        activate(edges) %>% 
        select(Community, Isolate1, Isolate2, InitialFrequency, Time, Isolate1MeasuredFreq, RawDataType) %>% 
        filter(Isolate1 == temp_isolates[i] | Isolate2 == temp_isolates[i]) %>%
        mutate(FocalIsolateFreq = ifelse(Isolate1 == temp_isolates[i], Isolate1MeasuredFreq, 1-Isolate1MeasuredFreq)) %>% 
        
        
        pairs_freq <- pairs_melted %>% 
        unite(col = "InitialFrequency", Isolate1InitialODFreq:Isolate2InitialODFreq, sep = "-") %>% 
        select(Community, Isolate1, Isolate2, InitialFrequency, Time, Isolate1MeasuredFreq, RawDataType) %>% 
        as_tibble()
    
    temp_comm <- pairs_freq %>%
        filter(Community == "C1R2") %>% 
        filter(Time == "T8")
    
    temp_isolates <- unique(c(temp_comm$Isolate1, temp_comm$Isolate2))
    temp_list <- rep(list(NA), length(temp_isolates))
    for (i in 1:length(temp_isolates)) { 
        temp_list[[i]] <- temp_comm %>% 
            filter(Isolate1 == temp_isolates[i] | Isolate2 == temp_isolates[i]) %>%
            mutate(FocalIsolateFreq = ifelse(Isolate1 == temp_isolates[i], Isolate1MeasuredFreq, 1-Isolate1MeasuredFreq)) %>% 
            mutate(Isolate = temp_isolates[i]) %>% 
            select(Isolate, FocalIsolateFreq)
    }
    
    temp_list %>% bind_rows %>%
        group_by(Isolate) %>% 
        summarize(CompetitiveScore = mean(FocalIsolateFreq))
    
}




