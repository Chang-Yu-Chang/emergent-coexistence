#' Compute the tournament ranks and plot pairwise interaction matrices
#' There are R functions that I wrote for convenience. The source code is saved in `misc/network_function.R`
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
communities <- fread(here::here("data/output/communities.csv"))
community_name <- communities$Community
community_size <- communities$CommunitySize

## Plot pairwise matrix
p_net_matrix_plot_list <- rep(list(NA), length(community_name)) %>% setNames(community_name)
for (i in 1:length(community_name)) {
  p_net_matrix_plot_list[[i]] <- plot_adjacent_matrix(net_list[[i]])
  #      axis_order_by = "competitive_rank",
  #      isolate_label = c("name", "ExpID", "Family", "Genus", "Fermenter"))
}

if (FALSE) {
  # Trait ranks ----
  ## OD620
  p_matrix_plot_list_OD620 <- rep(list(NA), length(community_name)) %>% setNames(community_name)
  for (i in 1:length(community_name)) {
    p_matrix_plot_list_OD620[[i]] <-
      adj_matrix_plot(net_list[[i]],
                      axis_order_by = "OD620",
                      isolate_label = c("name", "ExpID", "Family", "Genus", "Fermenter"))
  }


  ## Plot pairwise matrix
  p_matrix_plot_list_growth_rate <- rep(list(NA), length(community_name)) %>% setNames(community_name)
  for (i in 1:length(community_name)) {
    p_matrix_plot_list_growth_rate[[i]] <-
      adj_matrix_plot(net_list[[i]],
                      axis_order_by = "glucose_MaxGrowthRate",
                      isolate_label = c("name", "ExpID", "Family", "Genus", "Fermenter"))
  }


  ## Plot pairwise matrix
  p_matrix_plot_list_lag <- rep(list(NA), length(community_name)) %>% setNames(community_name)
  for (i in 1:length(community_name)) {
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
