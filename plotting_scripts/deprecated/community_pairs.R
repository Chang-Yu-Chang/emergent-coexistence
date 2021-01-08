library(tidygraph)
library(ggraph)
library(tidyverse)
library(data.table)
library(cowplot)
source("network_functions.R")


# Read data
communities <- fread("../data/temp/communities.csv")
community_names_ordered_by_size <- communities %>% arrange(CommunitySize) %>% pull(Community)
simulated_motif_counts <- fread("../data/temp/simulated_motif_counts.txt")
observed_motif_counts <- fread("../data/temp/observed_motif_counts.txt")
random_motif_counts <- fread("../data/temp/random_motif_counts.txt")
random_motif_counts_percentile <- fread("../data/temp/random_motif_counts_percentile.txt")
load("../data/temp/graph_list.Rdata") # Load observed networks graph_list
load("../data/temp/example_motif_list.Rdata") # Load example motif graphs example_motifs

pairs <- fread("../data/output/pairs.csv") %>%
    as_tibble()


pairs %>%
    filter(PairFermenter %in% c("FF", "FN", "NN")) %>%
    ggplot() +
    geom_bar(aes(x = PairFermenter, fill = InteractionType), stat = "count") +
    #scale_x_discrete(expand = c(0,0)) +
    scale_y_continuous(expand = c(0,0)) +
    theme_cowplot() +
    theme(legend.title = element_blank()) +
    labs(x = "Pairs")


pairs %>%
    filter(PairFermenter %in% c("FF", "FN", "NN")) %>%
    group_by(PairFermenter, InteractionType) %>%
    summarize(Count = n()) %>%
    group_by(PairFermenter) %>%
    mutate(RelativeCount = Count / sum(Count)) %>%
    ggplot() +
    geom_bar(aes(x = PairFermenter, y = RelativeCount, fill = InteractionType), stat = "identity") +
    scale_y_continuous(expand = c(0,0)) +
    theme_cowplot() +
    theme(legend.title = element_blank()) +
    panel_border(color = 1) +
    labs(x = "Pairs")


