#' Plot the result from dioagonal analysis
library(tidyverse)
library(data.table)

communities <- fread(here::here("data/output/communities.csv"))
communities_name <- communities$Community
networks_diag <- fread(here::here("data/output/networks_diag.csv")) %>%
  mutate(Community = ordered(Community, communities_name))
networks_diag_randomized <- fread(here::here("data/output/networks_diag_randomized.csv")) %>%
  mutate(Community = ordered(Community, communities_name))

# Specify the diagonal distances that have zero-entries of coexistence count ----
## Make temporary df that has the full diagonal distances
temp_list <- rep(list(NA), length(communities_name))
for (i in 1:length(communities_name)) temp_list[[i]] <- data.frame(Community = communities$Community[i], DistanceToDiagonal = 1:(communities$CommunitySize[i]-1))
temp_df <- rbindlist(temp_list) %>% as_tibble # For 13 observed networks
temp_df2 <- rep(list(temp_df), 10000) %>% rbindlist(idcol = "Randomization") %>% as_tibble()  # For 10000 randomization

## Joint the temp df to data
networks_diag <- temp_df %>%
  left_join(networks_diag, by = c("Community", "DistanceToDiagonal")) %>%
  replace_na(list(CountCoexistence = 0))
networks_diag_randomized <- temp_df2 %>%
  left_join(networks_diag_randomized, by = c("Randomization", "Community", "DistanceToDiagonal")) %>%
  replace_na(list(CountCoexistence = 0))



# Plot ----

# R function for plotting boxplot
plot_boxplot <- function(networks_diag, networks_diag_randomized, facet_by_community = TRUE) {
  ggplot() +
    # Random networks
    geom_boxplot(data = networks_diag_randomized, aes(x = DistanceToDiagonal, y = CountCoexistence, group = DistanceToDiagonal)) +
    # Observed networks
    geom_point(data = networks_diag, aes(x = DistanceToDiagonal, y = CountCoexistence, group = DistanceToDiagonal), col = "red", size = 2) +
    geom_line(data = networks_diag, aes(x = DistanceToDiagonal, y = CountCoexistence), col = "red") +
    {if (facet_by_community) facet_wrap(Community~., scale = "free") } +
    # Asethetic setup
    theme_bw() +
    theme(panel.grid.minor = element_blank()) +
    labs(x = "distance to diagional (|i-j|)", y = "count of coexistence")
}

# Individual community
p_networks_diag <- plot_boxplot(networks_diag, networks_diag_randomized, facet_by_community = T) +
  ggtitle("13 communities")

# Network pool
## Pool all communities
networks_diag_pool <- networks_diag %>% group_by(DistanceToDiagonal) %>% summarize(CountCoexistence = sum(CountCoexistence))
networks_diag_randomized_pool <- networks_diag_randomized %>% group_by(Randomization, DistanceToDiagonal) %>% summarize(CountCoexistence = sum(CountCoexistence))
p_networks_diag_pool1 <- plot_boxplot(networks_diag_pool, networks_diag_randomized_pool, facet_by_community = F) +
  ggtitle("13 communities pooled together")

## Pool small communities
networks_diag_pool <- networks_diag %>% filter(!Community %in% c("C11R1", "C11R2")) %>% group_by(DistanceToDiagonal) %>% summarize(CountCoexistence = sum(CountCoexistence))
networks_diag_randomized_pool <- networks_diag_randomized %>% filter(!Community %in% c("C11R1", "C11R2")) %>% group_by(Randomization, DistanceToDiagonal) %>% summarize(CountCoexistence = sum(CountCoexistence))
p_networks_diag_pool2 <- plot_boxplot(networks_diag_pool, networks_diag_randomized_pool, facet_by_community = F) +
  ggtitle("11 small communities pooled together")

## Pool small communities
networks_diag_pool <- networks_diag %>% filter(!Community %in% c("C1R7", "C11R1", "C11R2")) %>% group_by(DistanceToDiagonal) %>% summarize(CountCoexistence = sum(CountCoexistence))
networks_diag_randomized_pool <- networks_diag_randomized %>% filter(!Community %in% c("C1R7", "C11R1", "C11R2")) %>% group_by(Randomization, DistanceToDiagonal) %>% summarize(CountCoexistence = sum(CountCoexistence))
p_networks_diag_pool3 <- plot_boxplot(networks_diag_pool, networks_diag_randomized_pool, facet_by_community = F) +
  ggtitle("10 small communities pooled together")













