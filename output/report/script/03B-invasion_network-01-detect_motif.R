#' Detect the motifs in an invasion network
#' 1. Detect motif counts in observed network
#' 2. Detect motif counts in randomized networks
#' 3. Compute 5% and 95% percentiles of motif counts in randomized networks
#' 4. Plot the distribution of randomized network motifs
library(tidyverse)
library(data.table)
library(igraph)
source(here::here("plotting_scripts/network_functions.R"))

# Data
isolates <- fread(here::here("data/output/isolates.csv"))
pairs <- fread(here::here("data/output/pairs.csv"))
communities <- fread(here::here("data/output/communities.csv"))
communities_name <- communities$Community
load(here::here("data/output/network_community.Rdata")) # net_list
load(here::here("data/output/network_randomized.Rdata")) # net_randomized_list

# R function for detecting competitive network motifs
# Only extract 7 motifs that appears in a fully connected competitive network
#count_motif <- function (net) triad_census(net)[c(10, 9, 12, 14, 13, 15, 16)]


# Count motifs ----
# Count the motifs in observed networks
networks_motif <- lapply(net_list, function (net) {
  data.frame(Motif = paste0("Motif", 1:7), CountMotif = count_motif(net))
}) %>%
  rbindlist(idcol = "Community")

# Count the motifs in randomized networks
networks_motif_randomized_list <- rep(list(NA), length(net_list))
names(networks_motif_randomized_list) <- communities_name

tt <- proc.time()
for (i in 1:length(net_list)) {
  temp_tt <- proc.time()
  networks_motif_randomized_list[[i]] <- lapply(net_randomized_list[[i]], function (net) {
    data.frame(Motif = paste0("Motif", 1:7), CountMotif = count_motif(net))
  }) %>%
    rbindlist(idcol = "Randomization")
  # Print
  cat("\n\n", communities_name[i])
  cat("\n", (proc.time() - temp_tt)[3], "seconds")
  if (i == length(net_list)) cat("\n\n total time:", (proc.time() - tt)[3], "seconds")
}

networks_motif_randomized <- rbindlist(networks_motif_randomized_list, idcol = "Community")

## Find 5th and 95th percentiles for each motif within a community
b = 1000
networks_motif_randomized_percentile <-
  networks_motif_randomized %>%
  group_by(Community, Motif) %>%
  arrange(CountMotif) %>%
  select(-Randomization) %>%
  slice(c(b*0.05, b*0.95)) %>%
  mutate(Percentile = c("p5", "p95")) %>%
  {.}

fwrite(networks_motif, file = here::here("data/output/networks_motif.csv"))
fwrite(networks_motif_randomized, file = here::here("data/temp/networks_motif_randomized.csv"))
fwrite(networks_motif_randomized_percentile, here::here("data/temp/networks_motif_randomized_percentile.csv"))







# Old code
if (FALSE) {

  # Detect motifs in observations ----
  networks_motif <-
    data.frame(Community = rep(community_name, each = 7),
               Motif = rep(1:7, length(community_name)),
               MotifCount = NA) %>%
    mutate(Community = ordered(Community, levels = community_name))

  for (i in 1:length(community_name)) networks_motif[networks_motif$Community == names(net_list)[i], "MotifCount"] <- count_motif(net_list[[i]])


  # Detect motifs in randomized networks ----
  #' This may take ~5 minutes for 1000 boostrapping
  b <- 1000 # Number of boostrapping

  ## Empty df
  networks_motif_random <-
    data.frame(
      Community = rep(community_name, each = 7 * b), # 7 motifs
      Bootstrap = rep(1:b, each = 7),
      Motif = rep(1:7, length(community_name)),
      MotifCount = NA) %>%
    mutate(Community = ordered(Community, levels = community_name))

  ## Randomization
  temp_time <- proc.time()
  for (i in 1:length(community_name)) {
    for (b_loop_index in 1:b) {
      # Randomize network
      net_random <- network_randomize(net_list[[i]], method = "shuffle_interaction")

      # Update simulation result
      temp_index <- networks_motif_random$Community == community_name[i] & networks_motif_random$Bootstrap == b_loop_index
      networks_motif_random[temp_index,"MotifCount"] <- count_motif(net_random$net_rand)
    }
    cat(community_name[i], " ")
    cat((proc.time() - temp_time)[3], "seconds\n\n")
  }

  ## Find 5th and 95th percentiles for each motif within a community
  networks_motif_random_percentile <-
    networks_motif_random %>%
    group_by(Community, Motif) %>%
    arrange(MotifCount) %>%
    select(-Bootstrap) %>%
    slice(c(b*0.05, b*0.95)) %>%
    mutate(Percentile = c("p5", "p95")) %>%
    #  spread(Percentile, MotifCount) %>%
    {.}


  ## Save the randomized network motif count in Rdata ----
  if (TRUE) save(networks_motif, networks_motif_random_percentile, networks_motif_random, file = here::here("data/temp/networks_motif_random.Rdata"))



  # Plot the observed motif counts and randomized network motif counts ----
  load(here::here("data/temp/networks_motif_random.Rdata"))
  b <- max(networks_motif_random$Bootstrap) # Extract the number of bootstrapping

  ## Plot observations and 5% and 95% percentiles
  p_motif_random <-
    ggplot() +
    # 5% and 95% percentiles in randomized networks
    geom_point(data = networks_motif_random_percentile, aes(x = Motif, y = MotifCount, group = Motif), col = "black") +
    geom_segment(data = spread(networks_motif_random_percentile, Percentile, MotifCount),
                 aes(x = Motif, xend = Motif, y = p5, yend = p95), col = "black") +
    # Observations
    geom_point(data = networks_motif, aes(x = Motif, y = MotifCount), col = "red") +
    #  scale_colour_manual(name="Error Bars",values=c("black", "red")) +
    facet_wrap(Community ~., scale = "free_y") +
    theme_bw() +
    ggtitle("Motif counts in the invasion networks") +
    labs(x = "motif", y = "Count")

  ## Plot histogram
  p_motif_random_distribution_list <- rep(list(NA), length(communities$Community)) %>%
    setNames(communities$Community)

  for (i in 1:length(communities$Community)) {
    networks_motif_comm <- filter(networks_motif, Community == communities$Community[i])
    networks_motif_random_comm <- filter(networks_motif_random, Community == communities$Community[i])
    p_motif_random_distribution_list[[i]] <- motif_distribution_plot(networks_motif_comm, networks_motif_random_comm)
  }

}

