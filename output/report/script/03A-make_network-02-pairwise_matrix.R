#' Compute the tournament ranks and plot pairwise interaction matrices
#' There are R functions that I wrote for convenience. The source code is saved in `plotting_scripts/network_function.R`
#' 1. tournament_rank() computes the isolates' ranks
#' 2. adj_matrix_plot() plots the pairwise interaction matrix. It also uses other self-use
#' functions such as adj_from_net() and tournament_rank().
#' Make pairwise network from outcome of pairwise competitions
library(tidyverse)
library(data.table)
library(tidygraph)
library(ggraph)
source(here::here("plotting_scripts/network_functions.R"))

# Load Rdata which has a R list for community network `net_list` -----
load(file = here::here("data/output/network_community.Rdata"))
communities <- read_csv(here::here("data/output/communities.csv"))
pairs <- read_csv(here::here("data/output/pairs.csv")); pairs$InteractionType[pairs$InteractionType == "neutrality"] <- "coexistence"
isolates <- read_csv(here::here("data/output/isolates.csv"))
interaction_type <- c("exclusion", "coexistence", "lose", "bistability", "neutrality", "self", "undefined")
interaction_color = c("#DB7469", "#557BAA", "#73C966", "#EECF6D", "#8650C4", "black", "grey80")
names(interaction_color) <- interaction_type


for (i in 1:length(net_list)) {
    net_list[[i]] <- net_list[[i]] %>%
        activate(nodes) %>%
        left_join(isolates %>% filter(Community == communities$Community[i]) %>% select(Community, Isolate, PlotRank, ID))
}

## Plot pairwise matrix
p_net_matrix_plot_list <- lapply(net_list, plot_adjacent_matrix) %>% setNames(communities$Community)

plot_adjacent_matrix(net_list[[1]], show.axis = T)
plot_adjacent_matrix(net_list[[2]], show.axis = T)

pairs %>%
    left_join(rename_with(isolates, ~ paste0(., "1"), !contains("Community"))) %>%
    left_join(rename_with(isolates, ~ paste0(., "2"), !contains("Community"))) %>%
    filter(Community == "C1R2") %>%
    select(ID1, ID2, Isolate1, Isolate2, From, To, InteractionType)


label_ID <- as.character(isolates$ID[isolates$Community == "C1R2"]) %>% setNames(as.character(isolates$PlotRank[isolates$Community == "C1R2"]))

isolates %>%
    split.data.frame("Community") %>%
    lapply(function(x) {
    pull(x, ID) %>%
    combn(2) %>% t() %>%
    as_tibble() %>% setNames(c("ID1", "ID2"))
    }) %>%
    setNames(communities$Community)
    data.table::rbindlist()

pairs %>%
    left_join(rename_with(isolates, ~ paste0(., "1"), !contains("Community"))) %>%
    left_join(rename_with(isolates, ~ paste0(., "2"), !contains("Community"))) %>%
    filter(Community == "C1R2") %>%
    mutate(ID1 = factor(ID1), ID2 = factor(ID2)) %>%
    #filter(PlotRank1 >= PlotRank2) %>%
    ggplot() +
    #geom_tile(aes(x = ID1, y = ID2, fill = InteractionType), width = 0.9, height = 0.9) +
    geom_tile(aes(x = ID2, y = ID1, fill = InteractionType), width = 0.9, height = 0.9) +
    scale_fill_manual(values = interaction_color) +
    theme_classic() +
    NULL



if (FALSE) {
  # Trait ranks ----
  ## OD620
  p_matrix_plot_list_OD620 <- rep(list(NA), length(communities$Community)) %>% setNames(communities$Community)
  for (i in 1:length(communities$Community)) {
    p_matrix_plot_list_OD620[[i]] <-
      adj_matrix_plot(net_list[[i]],
                      axis_order_by = "OD620",
                      isolate_label = c("name", "ExpID", "Family", "Genus", "Fermenter"))
  }


  ## Plot pairwise matrix
  p_matrix_plot_list_growth_rate <- rep(list(NA), length(communities$Community)) %>% setNames(communities$Community)
  for (i in 1:length(communities$Community)) {
    p_matrix_plot_list_growth_rate[[i]] <-
      adj_matrix_plot(net_list[[i]],
                      axis_order_by = "glucose_MaxGrowthRate",
                      isolate_label = c("name", "ExpID", "Family", "Genus", "Fermenter"))
  }


  ## Plot pairwise matrix
  p_matrix_plot_list_lag <- rep(list(NA), length(communities$Community)) %>% setNames(communities$Community)
  for (i in 1:length(communities$Community)) {
    p_matrix_plot_list_lag[[i]] <-
      adj_matrix_plot(net_list[[i]],
                      axis_order_by = "glucose_LagTime",
                      isolate_label = c("name", "ExpID", "Family", "Genus", "Fermenter"))
  }

  #Plot the communities' pairwise matrix, whose axis are ordered by isolates' traits
  # Yield in glucose
  pdf(here::here("output/figure/communities_pairwise_matrix_OD620.pdf"), width = 12, height = 8); for (i in 1:length(net_list)) print(p_matrix_plot_list_OD620[[i]]); invisible(dev.off())

  # Maximal growth rate in glucose
  pdf(here::here("output/figure/communities_pairwise_matrix_growth_rate.pdf"), width = 12, height = 8); for (i in 1:length(net_list)) print(p_matrix_plot_list_growth_rate[[i]]); invisible(dev.off())

  # Lag time in glucose
  pdf(here::here("output/figure/communities_pairwise_matrix_lag.pdf"), width = 12, height = 8); for (i in 1:length(net_list)) print(p_matrix_plot_list_lag[[i]]); invisible(dev.off())
}
