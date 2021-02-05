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

random_network_pairs <- fread(here::here("data/temp/random_network_pairs.csv"))
random_network_pairs_freq <- fread(here::here("data/temp/random_network_pairs_freq.csv"))
communities <- fread(here::here("data/output/communities.csv"))
communities_name <- communities$Community

isolates_random_tournament <-
  random_network_pairs %>%
  split.data.frame(f = .$Community) %>%
  lapply(tournament_rank) %>%
  rbindlist(idcol = "Community")

isolates_random <- fread(here::here("data/output/isolates_random.csv")) %>%
  select(Community = AssemblyCommunity, Isolate = AssemblyIsolate, ExpID, Family, Genus) %>%
  left_join(isolates_random_tournament)

### Only use the isolates in the pairs
temp <- random_network_pairs %>%
  select(Community, Isolate1, Isolate2) %>%
  group_by(Community) %>%
  pivot_longer(-Community, names_to = "temp", values_to = "Isolate") %>%
  distinct(Isolate) %>%
  select(Community, Isolate)
isolates_random <- isolates_random %>% right_join(temp)


# Make network ----
random_network_communities_name <- c("AcrAss1", "RanAss1") #random_network_communities_name <- c("AcrAss1", "AcrAss2", "RanAss1", "RanAss2")
random_network_list <- rep(list(NA), length(random_network_communities_name))
#random_network_list <- rep(list(NA), length(random_network_communities_name))
names(random_network_list) <- random_network_communities_name
for (i in 1:length(random_network_list)) {
  random_network_list[[i]] <-
    make_network(isolates = filter(isolates_random, Community == random_network_communities_name[i]),
                 pairs = filter(random_network_pairs, Community == random_network_communities_name[i]))
}


# Plot network
# Graph
p_random_network_list <- rep(list(NA), length(random_network_list))
names(p_random_network_list) <- random_network_communities_name
for (i in 1:length(random_network_list)) p_random_network_list[[i]] <- plot_competitive_network(random_network_list[[i]])

# Matrix
p_random_network_matrix_list <- rep(list(NA), length(random_network_list))
names(p_random_network_matrix_list) <- random_network_communities_name
for (i in 1:length(random_network_list)) p_random_network_matrix_list[[i]] <- plot_adjacent_matrix(random_network_list[[i]])

## Save network result
save(random_network_list, p_random_network_list, p_random_network_matrix_list, file = here::here("data/temp/random_network_community.Rdata"))


# Make randomized networks

# Read networks
load(here::here("data/temp/random_network_community.Rdata")) # random_network_list

# Observation
random_networks_motif <-
  random_network_list %>%
  lapply(function (net) {
    data.frame(Motif = paste0("Motif", 1:7), CountMotif = count_motif(net))
  }) %>%
  rbindlist(idcol = "Community")

fwrite(random_networks_motif, here::here("data/output/random_networks_motif.csv"))

# Make random networks and save them in a list ----
b <- 1000 # Number of boostrapping/randomization

# Make an empty two-layer R list
random_network_randomized_list <- rep(list(rep(list(NA), b)), length(random_network_list))
names(random_network_randomized_list) <- random_network_communities_name

tt <- proc.time()
for (i in 1:length(random_network_list)) {
  temp_tt <- proc.time()
  for (b_loop_index in 1:b) {
    set.seed(b_loop_index)
    random_network_randomized_list[[i]][[b_loop_index]] <- randomize_network(random_network_list[[i]])
    if (b_loop_index %% 1000 == 0) cat("\n boostrap =", b_loop_index)
  }
  # Print
  cat("\n\n", random_network_communities_name[i])
  cat("\n", (proc.time() - temp_tt)[3], "seconds")
  if (i == length(random_network_list)) cat("\n\n total time:", (proc.time() - tt)[3], "seconds\n\n")

}


# Save the data file
save(random_network_randomized_list, file = here::here("data/output/random_network_randomized.Rdata"))

# Motif count----
# Count the motifs in observed networks
random_networks_motif <- lapply(random_network_list, function (net) {
  data.frame(Motif = paste0("Motif", 1:7), CountMotif = count_motif(net))
}) %>%
  rbindlist(idcol = "Community")

# Count the motifs in randomized networks
random_networks_motif_randomized_list <- rep(list(NA), length(random_network_list))
names(random_networks_motif_randomized_list) <- random_network_communities_name

tt <- proc.time()
for (i in 1:length(random_network_list)) {
  temp_tt <- proc.time()
  random_networks_motif_randomized_list[[i]] <- lapply(random_network_randomized_list[[i]], function (net) {
    data.frame(Motif = paste0("Motif", 1:7), CountMotif = count_motif(net))
  }) %>%
    rbindlist(idcol = "Randomization")
  # Print
  cat("\n\n", communities_name[i])
  cat("\n", (proc.time() - temp_tt)[3], "seconds")
  if (i == length(random_network_list)) cat("\n\n total time:", (proc.time() - tt)[3], "seconds")
}

random_networks_motif_randomized <- rbindlist(random_networks_motif_randomized_list, idcol = "Community")

## Find 5th and 95th percentiles for each motif within a community
b = 1000
random_networks_motif_randomized_percentile <-
  random_networks_motif_randomized %>%
  group_by(Community, Motif) %>%
  arrange(CountMotif) %>%
  select(-Randomization) %>%
  slice(c(b*0.05, b*0.95)) %>%
  mutate(Percentile = c("p5", "p95")) %>%
  #  spread(Percentile, MotifCount) %>%
  {.}

# Plot observations and 5% and 95% percentiles ----
p_random_networks_motif_randomized <-
  ggplot() +
  # 5% and 95% percentiles in randomized networks
  geom_point(data = random_networks_motif_randomized_percentile, aes(x = Motif, y = CountMotif, group = Motif), col = "black") +
  geom_segment(data = spread(random_networks_motif_randomized_percentile, Percentile, CountMotif),
               aes(x = Motif, xend = Motif, y = p5, yend = p95), col = "black") +
  # Observations
  geom_point(data = random_networks_motif, aes(x = Motif, y = CountMotif), col = "red") +
  #  scale_colour_manual(name="Error Bars",values=c("black", "red")) +
  facet_wrap(Community ~., scale = "free_y") +
  theme_bw() +
  labs(x = "", y = "Count")

p_random_networks_motif_randomized





