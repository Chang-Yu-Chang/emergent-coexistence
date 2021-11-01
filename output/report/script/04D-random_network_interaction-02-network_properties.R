#' Determine network propetris
#' 1. Make network
#'    - Plot network
#'    - Plot matrix
#' 2. Network motif
#'    - Make random networks
#' 3. Competitive hierarchy
#'    - Compute

library(tidyverse)
library(data.table)
library(tidygraph)
library(ggraph)
source(here::here("plotting_scripts/network_functions.R"))
tournament_rank <- function(pairs) {
  isolate_name <- pairs %>% select(Isolate1, Isolate2) %>% unlist %>% unique %>% sort()
  # Isolates' ranks in the tournament
  tour_rank <- data.frame(
    Isolate = isolate_name,
    # Win
    Win = filter(pairs, InteractionType == "exclusion") %>%
      select(From) %>% unlist() %>% factor(isolate_name) %>% table() %>% as.vector(),
    # Lose
    Lose = filter(pairs, InteractionType == "exclusion") %>%
      select(To) %>% unlist() %>% factor(isolate_name) %>% table() %>% as.vector(),
    # Draw; Note that I consider neturality and bistability as draw in the tournament
    Draw = filter(pairs, InteractionType %in% c("coexistence", "neutrality", "bistability")) %>%
      select(From, To) %>% unlist() %>% factor(isolate_name) %>% table() %>% as.vector())

  # Arrange the df by score
  tour_rank <- tour_rank %>%
    mutate(Score = Win - Lose + 0 * Draw, Game = Win + Lose + Draw) %>%
    arrange(desc(Score))

  # Calculate rank by score; same scores means the same ranks
  temp_score <- ordered(tour_rank$Score, levels = sort(unique(tour_rank$Score), decreasing = T))
  temp_score_table <- table(temp_score)
  temp <- NULL; temp_counter = 1
  for (i in 1:length(temp_score_table)) {
    temp <- c(temp, rep(temp_counter, temp_score_table[i]))
    temp_counter <- temp_counter + temp_score_table[i]
  }

  tour_rank$Rank <- temp
  tour_rank$PlotRank <- 1:nrow(tour_rank)
  return(tour_rank)
}

pairs_random <- fread(here::here("data/output/pairs_random.csv"))
communities <- fread(here::here("data/output/communities.csv"))
communities_name <- communities$Community

isolates_random_tournament <-
  pairs_random %>%
  split.data.frame(f = .$Community) %>%
  lapply(tournament_rank) %>%
  rbindlist(idcol = "Community")

isolates_random <- fread(here::here("data/output/isolates_random.csv")) %>%
  select(Community, Isolate, ExpID, Family, Genus) %>%
  left_join(isolates_random_tournament)

### Only use the isolates in the pairs
# temp <- pairs_random %>%
#   select(Community, Isolate1, Isolate2) %>%
#   group_by(Community) %>%
#   pivot_longer(-Community, names_to = "temp", values_to = "Isolate") %>%
#   distinct(Isolate) %>%
#   select(Community, Isolate)
#isolates_random <- isolates_random %>% right_join(temp)


# Make network ----
random_network_communities_name <- c("AcrAss1", "RanAss1") #random_network_communities_name <- c("AcrAss1", "AcrAss2", "RanAss1", "RanAss2")
net_random_list <- rep(list(NA), length(random_network_communities_name))
#net_random_list <- rep(list(NA), length(random_network_communities_name))
names(net_random_list) <- random_network_communities_name
for (i in 1:length(net_random_list)) {
  net_random_list[[i]] <-
    make_network(isolates = filter(isolates_random, Community == random_network_communities_name[i]),
                 pairs = filter(pairs_random, Community == random_network_communities_name[i]))
}


# Plot network
# Graph
p_net_random_list <- rep(list(NA), length(net_random_list))
names(p_net_random_list) <- random_network_communities_name
for (i in 1:length(net_random_list)) p_net_random_list[[i]] <- plot_competitive_network(net_random_list[[i]])

# Matrix
p_net_random_matrix_list <- rep(list(NA), length(net_random_list))
names(p_net_random_matrix_list) <- random_network_communities_name
for (i in 1:length(net_random_list)) p_net_random_matrix_list[[i]] <- plot_adjacent_matrix(net_random_list[[i]])

## Save network result
save(net_random_list, p_net_random_list, p_net_random_matrix_list, file = here::here("data/output/network_random.Rdata"))


# Make randomized networks

# Read networks
load(here::here("data/output/network_random.Rdata")) # net_random_list

# Observation
networks_random_motif <-
  net_random_list %>%
  lapply(function (net) {
    data.frame(Motif = paste0("Motif", 1:7), CountMotif = count_motif(net))
  }) %>%
  rbindlist(idcol = "Community")

fwrite(networks_random_motif, here::here("data/output/networks_random_motif.csv"))

# Make random networks and save them in a list ----
b <- 1000 # Number of boostrapping/randomization

# Make an empty two-layer R list
net_random_randomized_list <- rep(list(rep(list(NA), b)), length(net_random_list))
names(net_random_randomized_list) <- random_network_communities_name

tt <- proc.time()
for (i in 1:length(net_random_list)) {
  temp_tt <- proc.time()
  for (b_loop_index in 1:b) {
    set.seed(b_loop_index)
    net_random_randomized_list[[i]][[b_loop_index]] <- randomize_network(net_random_list[[i]])
    if (b_loop_index %% 1000 == 0) cat("\n boostrap =", b_loop_index)
  }
  # Print
  cat("\n\n", random_network_communities_name[i])
  cat("\n", (proc.time() - temp_tt)[3], "seconds")
  if (i == length(net_random_list)) cat("\n\n total time:", (proc.time() - tt)[3], "seconds\n\n")

}


# Save the data file
save(net_random_randomized_list, file = here::here("data/output/network_random_randomized.Rdata"))

# Motif count----
# Count the motifs in observed networks
networks_random_motif <- lapply(net_random_list, function (net) {
  data.frame(Motif = paste0("Motif", 1:7), CountMotif = count_motif(net))
}) %>%
  rbindlist(idcol = "Community")

# Count the motifs in randomized networks
net_random_motif_randomized_list <- rep(list(NA), length(net_random_list))
names(net_random_motif_randomized_list) <- random_network_communities_name

tt <- proc.time()
for (i in 1:length(net_random_list)) {
  temp_tt <- proc.time()
  net_random_motif_randomized_list[[i]] <- lapply(net_random_randomized_list[[i]], function (net) {
    data.frame(Motif = paste0("Motif", 1:7), CountMotif = count_motif(net))
  }) %>%
    rbindlist(idcol = "Randomization")
  # Print
  cat("\n\n", communities_name[i])
  cat("\n", (proc.time() - temp_tt)[3], "seconds")
  if (i == length(net_random_list)) cat("\n\n total time:", (proc.time() - tt)[3], "seconds")
}

net_random_motif_randomized <- rbindlist(net_random_motif_randomized_list, idcol = "Community")

## Find 5th and 95th percentiles for each motif within a community
b = 1000
net_random_motif_randomized_percentile <-
  net_random_motif_randomized %>%
  group_by(Community, Motif) %>%
  arrange(CountMotif) %>%
  select(-Randomization) %>%
  slice(c(b*0.05, b*0.95)) %>%
  mutate(Percentile = c("p5", "p95")) %>%
  #  spread(Percentile, MotifCount) %>%
  {.}

# Plot observations and 5% and 95% percentiles ----
p_net_random_motif_randomized <-
  ggplot() +
  # 5% and 95% percentiles in randomized networks
  geom_point(data = net_random_motif_randomized_percentile, aes(x = Motif, y = CountMotif, group = Motif), col = "black") +
  geom_segment(data = spread(net_random_motif_randomized_percentile, Percentile, CountMotif),
               aes(x = Motif, xend = Motif, y = p5, yend = p95), col = "black") +
  # Observations
  geom_point(data = networks_random_motif, aes(x = Motif, y = CountMotif), col = "red") +
  #  scale_colour_manual(name="Error Bars",values=c("black", "red")) +
  facet_wrap(Community ~., scale = "free_y") +
  theme_bw() +
  labs(x = "", y = "Count")






