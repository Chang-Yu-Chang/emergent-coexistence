#' Read the motifs in three community pools
#' 1. All 13 communities pooled together
#' 2. Only large communities (C1R7, C11R1, C11R2)
#' 3. Only samll communities


# Load randomized network motif counts
networks_motif <- fread(root$find_file("data/temp/networks_motif.csv")) %>% mutate(Community = ordered(Community, communities_name))
networks_motif_rand <- fread(root$find_file("data/temp/networks_motif_rand.csv")) %>% mutate(Community = ordered(Community, communities_name))
networks_motif_rand_percentile <- fread(root$find_file("data/temp/networks_motif_rand_percentile.csv")) %>% mutate(Community = ordered(Community, communities_name))
b <- max(networks_motif_rand$Randomization) # Extract the number of Randomization

# R function for plotting the motif distribution in randomized networks  ----
motif_distribution_plot <- function (
  net_motif_comm,  # Observation
  net_motif_random_comm # Randomized networks
) {
  # Compute the 5% and 95% percentiles of randomized network motif counts
  net_motif_random_percentile_comm <-
    net_motif_random_comm %>%
    group_by(Community, Motif) %>%
    arrange(Count) %>%
    select(-Randomization) %>%
    slice(c(b*0.05, b*0.95)) %>%
    mutate(Percentile = c("p5", "p95")) %>%
    #  spread(Percentile, CountMotif) %>%
    {.}

  # Plot the histogram
  p_histogram <-
    net_motif_random_comm %>%
    ggplot() +
    # Distribution of randomized network motif counts
    geom_histogram(aes(x = CountMotif), col = "black", fill = "white") +
    # 5% and 95% percentiles of randomized network motif counts
    geom_vline(data = net_motif_random_percentile_comm, aes(xintercept = CountMotif, group = Motif), col = "black", alpha = 0.3) +
    geom_rect(data = spread(net_motif_random_percentile_comm, Percentile, CountMotif),
      aes(xmin = p5, xmax = p95), ymin = -Inf, ymax = Inf, fill = "black", alpha = 0.3) +
    # Observation
    geom_vline(data = net_motif_comm, aes(xintercept = CountMotif), col = "red") +
    facet_grid(Community ~ Motif, scale = "free_y") +
    coord_flip() +
    theme_bw() +
    theme(axis.text.y = element_blank(),
      strip.text = element_blank(),
      strip.background = element_blank(),
      axis.text.x = element_text(angle = 90)) +
    labs(x = "motif count") +
    ggtitle(unique(networks_motif_comm$Community))

  return(p_histogram)
}


motif_distribution_plot_pool <- function (
  net_motif_pool,  # Observation
  net_motif_random_pool # Randomized networks
) {
  # Compute the 5% and 95% percentiles of randomized network motif counts
  net_motif_random_pool_percentile <-
    net_motif_random_pool %>%
    group_by(Motif) %>%
    arrange(CountMotif) %>%
    select(-Randomization) %>%
    slice(c(b*0.05, b*0.95)) %>%
    mutate(Percentile = c("p5", "p95")) %>%
    #  spread(Percentile, CountMotif) %>%
    {.}


  # Plot the histogram
  p_histogram <-
    net_motif_random_pool %>%
    ggplot() +
    # Histogram
    geom_histogram(aes(x = CountMotif), col = "black", fill = "white") +
    # 5% and 95% percentiles
    geom_vline(data = net_motif_random_pool_percentile, aes(xintercept = CountMotif), color = "black", alpha = 0.3) +
    geom_rect(data = spread(net_motif_random_pool_percentile, Percentile, CountMotif),
      aes(xmin = p5, xmax = p95), ymin = -Inf, ymax = Inf, fill = "black", alpha = 0.3) +
    # Observation
    geom_vline(data = net_motif_pool, aes(xintercept = CountMotif), color = "red") +
    facet_grid(.~Motif, scale = "free_y") +
    coord_flip() +
    theme_bw() +
    theme(axis.text.y = element_blank(),
      strip.text = element_blank(),
      strip.background = element_blank(),
      axis.text.x = element_text(angle = 90)) +
    labs(x = "motif count")

  return(p_histogram)
}


# Pool 1: all communities ----
## Filter for pool 1 and pool motif counts
networks_motif_pool1 <-
  networks_motif %>%
  group_by(Motif) %>%
  summarize(CountMotif = sum(CountMotif))

networks_motif_rand_pool1 <-
  networks_motif_rand %>%
  group_by(Randomization, Motif) %>%
  summarize(CountMotif = sum(CountMotif))


# Pool 2: three large communities ----
networks_motif_pool2 <-
  networks_motif %>%
  filter(Community %in% c("C1R7", "C11R1", "C11R2")) %>%
  group_by(Motif) %>%
  summarize(CountMotif = sum(CountMotif))

networks_motif_rand_pool2 <-
  networks_motif_rand %>%
  filter(Community %in% c("C1R7", "C11R1", "C11R2")) %>%
  group_by(Randomization, Motif) %>%
  summarize(CountMotif = sum(CountMotif))


# Pool 3: 10 small communities ----
networks_motif_pool3 <-
  networks_motif %>%
  filter(!(Community %in% c("C1R7", "C11R1", "C11R2"))) %>%
  group_by(Motif) %>%
  summarize(CountMotif = sum(CountMotif))

networks_motif_rand_pool3 <-
  networks_motif_rand %>%
  filter(!(Community %in% c("C1R7", "C11R1", "C11R2"))) %>%
  group_by(Randomization, Motif) %>%
  summarize(CountMotif = sum(CountMotif))



# Plot histogram
p_motif_random_distribution_pool_list <- rep(list(NA), 3)
p_motif_random_distribution_pool_list[[1]] <- motif_distribution_plot_pool(networks_motif_pool1, networks_motif_rand_pool1)
p_motif_random_distribution_pool_list[[2]] <- motif_distribution_plot_pool(networks_motif_pool2, networks_motif_rand_pool2)
p_motif_random_distribution_pool_list[[3]] <- motif_distribution_plot_pool(networks_motif_pool3, networks_motif_rand_pool3)


# Combine the motif to histograms
p_motif_pool_list <- rep(list(NA), 3)
titles <- c("Pool 1: all 13 communities", "Pool 2: three large communities", "Pool 3: 10 small communities")

for (i in 1:3) {
  # Title
  title <- ggdraw() +
    draw_label(titles[i], fontface = 'bold', x = 0, hjust = 0) +
    theme(plot.margin = margin(0, 0, 0, 7))
  # Combine figures
   p_motif_pool_list[[i]] <-
plot_grid(title, p_motif_example, p_motif_random_distribution_pool_list[[i]],
  ncol = 1, rel_heights = c(0.1, 0.5, 1.5), align = "v", axis = "lr")
}






