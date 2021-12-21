#' Count the number of pairwise coexistence as a function of distance to diagonal (|i-j|)
library(tidyverse)
library(data.table)
library(tidygraph)
library(igraph)
source(here::here("plotting_scripts/network_functions.R"))
isolates <- read_csv(here::here("data/output/isolates.csv"))
pairs <- read_csv(here::here("data/output/pairs.csv"))
communities <- read_csv(here::here("data/output/communities.csv"))
load(here::here("data/output/network_community.Rdata"))
load(here::here("data/output/network_randomized.Rdata"))
# 20211215 remove random assembly network as it's not completed yet
net_list <- net_list[names(net_list) %>% str_detect("C\\d") %>% which]
net_randomized_list <- net_randomized_list[names(net_randomized_list) %>% str_detect("C\\d") %>% which]

# R function for calculating number of coexistence as the function of distance to diagonal
adj_from_net <- function(net) {
  # Get adjacent matrix
  net_m <- get.adjacency(net, attr = "InteractionType", sparse = F)

  # Exclusion: win or lose
  temp_index <- which(net_m=="exclusion", arr.ind = T) %>% as.data.frame()
  for(i in 1:nrow(temp_index)) net_m[temp_index$col[i], temp_index$row[i]] <- "lose"
  net_m[net_m=="exclusion"] <- "win"

  # Bistability
  temp_index2 <- which(net_m=="bistability", arr.ind = T) %>% as.data.frame()
  for(i in 1:nrow(temp_index2)) net_m[temp_index2$col[i], temp_index2$row[i]] <- "bistability"
  net_m[net_m=="bistability"] <- "bistability"

  # Neutrality
  temp_index2 <- which(net_m=="neutrality", arr.ind = T) %>% as.data.frame()
  for(i in 1:nrow(temp_index2)) net_m[temp_index2$col[i], temp_index2$row[i]] <- "neutrality"
  net_m[net_m=="neutrality"] <- "neutrality"

  # Diagonal
  diag(net_m) <- "self"

  # Undefined
  temp_index <- which(net_m=="undefined", arr.ind = T) %>% as.data.frame()
  for(i in 1:nrow(temp_index)) net_m[temp_index$col[i], temp_index$row[i]] <- "undefined"

  # NA
  net_m[net_m == "" | is.na(net_m)] <- NA

  return(net_m)
}
## Find the fraction of coexistence as a function of distance to diagonal
diag_distance <- function (net) {
  # Convert network to matrix
  temp_matrix <- adj_from_net(net)

  # Re-order the matrix axis by the isolates' competitive rank
  temp_rank <- tournament_rank(net)$Isolate
  temp_matrix <- temp_matrix[temp_rank, temp_rank]

  # Order the matrix axis by competitive score
  which(temp_matrix == "coexistence", arr.ind = T) %>%
    as_tibble() %>%
    filter(col > row) %>%
    mutate(DistanceToDiagonal = abs(row - col)) %>%
    group_by(DistanceToDiagonal) %>%
    summarize(CountCoexistence = n())
}

# Count the distance to diagonal in observed networks
networks_diag <- lapply(net_list, diag_distance) %>% bind_rows(.id = "Community")

# Count the distance to diagonal in randomized networks
networks_diag_randomized_list <- rep(list(NA), length(net_list))
names(networks_diag_randomized_list) <- names(net_list)

tt <- proc.time()
for (i in 1:length(net_list)) {
  temp_tt <- proc.time()
  networks_diag_randomized_list[[i]] <- lapply(net_randomized_list[[i]], diag_distance) %>%
    bind_rows(.id = "Replicate")
  # Print
  cat("\n\n", communities$Community[i])
  cat("\n", (proc.time() - temp_tt)[3], "seconds")
  if (i == length(net_list)) cat("\n\n total time:", (proc.time() - tt)[3], "seconds")
}

networks_diag_randomized <- bind_rows(networks_diag_randomized_list, .id = "Community")

# Save the result
write_csv(networks_diag, file = here::here("data/output/networks_diag.csv"))
write_csv(networks_diag_randomized, file = here::here("data/output/networks_diag_randomized.csv"))














