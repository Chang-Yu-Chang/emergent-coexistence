#' Types of pairwise competition outcomes
library(tidyverse)
library(data.table)

# Read data
pairs_freq <- fread(here::here("data/temp/pairs_freq.csv"))
pairs <- fread(here::here("data/output/pairs.csv"))
pairs_interaction_fitness <- fread(here::here("data/temp/pairs_interaction_fitness.csv")) %>% mutate(Community = ordered(Community, communities_name)) %>% as_tibble()

# R function for plotting frequency changes ----
plot_pairs_freq <- function(pairs_list, pairs_freq_df, show_strip = TRUE) {
  # A list of pairs with Community, Isolate1, and Isolate2
  pairs_list %>%
    # Joint the df `pairs_freq` that saves frequency data
    left_join(pairs_freq_df, by = c("Community", "Isolate1", "Isolate2")) %>%
    ungroup() %>%
    mutate(Isolate1InitialODFreq = factor(Isolate1InitialODFreq),
      InteractionType = ordered(InteractionType, c("exclusion", "bistability", "coexistence", "neutrality")),
      InteractionTypeFiner = ordered(InteractionTypeFiner, c("competitive exclusion", "mutual exclusion",
        "stable coexistence", "frequency-dependent coexistence", "neutrality"))) %>%
    ggplot(aes(x = Time, y = Isolate1MeasuredFreq, color = Isolate1InitialODFreq, group = Isolate1InitialODFreq)) +
    geom_point(size = 2) +
    geom_line(size = 1) +
    scale_y_continuous(breaks = seq(0,1,0.5)) +
    theme_bw() +
    theme(
      axis.title.x = element_blank(),
      panel.spacing = unit(5, "mm"),
      panel.grid.minor = element_blank()
      ) +
      {if (show_strip == FALSE) theme(strip.background = element_blank(), strip.text = element_blank())} +
    guides(color = F) +
    labs(x = "transfer", y = "isolate relative abundance")
}

# R function for plotting bars
plot_bar <- function(pairs, bar_by = "InteractionType", fill_by = "InteractionType", show_fill_legend = F) {
  # Modify bars
  pairs$InteractionTypeFiner[pairs$InteractionTypeFiner == "frequency-dependent coexistence"] <- "frequency-dependent\ncoexistence"

  # Color palettes
  interaction_type <- c("exclusion", "coexistence", "neutrality", "mutual exclusion", "frequency-dependent\ncoexistence")
  myColor = c("#DB7469", "#557BAA", "#8650C4", "red", "blue")
  names(myColor) <- interaction_type

  # Plot
  pairs %>%
    ungroup() %>%
    mutate(
      InteractionType = ordered(InteractionType, c("exclusion", "coexistence", "neutrality")),
      InteractionTypeFiner = ordered(InteractionTypeFiner, c("competitive exclusion", "mutual exclusion",
        "stable coexistence", "frequency-dependent\ncoexistence", "neutrality"))) %>%
    group_by(.dots = list(bar_by, fill_by)) %>%
    summarize(Count = n(), Fraction = n()/nrow(.)) %>%
    ggplot() +
    geom_bar(aes_string(x = bar_by, y = "Fraction", fill = fill_by), stat = "identity", color = 1, position = "dodge") +
    scale_y_continuous(limits = c(0, 1)) +
    scale_fill_manual(values = myColor) +
    facet_grid(as.formula(paste0("~", bar_by)), scale = "free_x") +
    { if (show_fill_legend == FALSE) guides(fill = F) } +
    theme_cowplot() +
    theme(axis.title.x = element_blank(),
      panel.spacing = unit(0, "mm"),
      panel.grid.major.x = element_blank(),
      legend.position = "top",
      legend.title = element_blank(),
      strip.text = element_blank())
}

# Example for plotting interaction types ----
## Four main examples; mainly for figure 1
pairs_example_outcomes <- pairs %>%
  group_by(InteractionType) %>%
  arrange(desc(Community), Isolate1) %>%
  slice(1) %>%
  select(Community, Isolate1, Isolate2, InteractionType, InteractionTypeFiner, From, To)


## Examples of outcomes in a finer scale; plot for figure 2C
pairs_example_outcomes_finer <-
  pairs_interaction_fitness %>%
  group_by(InteractionTypeFiner) %>%
  slice(1) %>%
  ungroup() %>%
  select(Community, Isolate1, Isolate2, InteractionType, InteractionTypeFiner, From, To)









