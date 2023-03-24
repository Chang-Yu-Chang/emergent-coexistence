library(tidyverse)
library(cowplot)
source(here::here("analysis/00-metadata.R"))
source(here::here("simulation/01-generate_input.R"))

# 0. parameters ----
input_parameters <- read_csv(here::here("simulation/01-input_parameters.csv"), col_types = cols())
input_poolPairs <- read_csv(here::here("simulation/03a-input_poolPairs.csv"), col_types = cols())
input_withinCommunityPairs <- read_csv(here::here("simulation/03b-input_withinCommunityPairs.csv"), col_types = cols())

# 1. Pool pairs ----
poolPairs_N_freq <- read_csv(paste0(folder_simulation, "aggregated/13-poolPairs_N_freq.csv"), col_types = cols()) %>%
    mutate(Time = factor(Time, c("init", "end"))) %>%
    mutate(Pair = factor(Pair, paste0("P", 1:1000))) %>%
    mutate(InitialFrequency = factor(InitialFrequency, c(5, 50, 95)))

# Barplot
monocultureSets_richness <- read_csv(paste0(folder_simulation, "aggregated/03-monocultureSets_richness.csv"), col_types = cols()) %>%
    mutate(Community = factor(Community, paste0("W", 0:(nrow(input_poolPairs)-1)))) %>%
    mutate(PairSize = choose(Richness, 2)) %>%
    arrange(desc(Richness)) %>%
    slice(1:20) %>%
    mutate(CommunityLabel = 1:20)
poolPairs_N_outcome <- read_csv(paste0(folder_simulation, "aggregated/13-poolPairs_N_outcome.csv"), col_types = cols()) %>%
    mutate(Community = factor(Community, paste0("W", 0:(nrow(input_poolPairs)-1)))) %>%
    # Top 20 richness
    filter(Community %in% monocultureSets_richness$Community)

p1 <- poolPairs_N_outcome %>%
    group_by(Community, InteractionType, .drop = F) %>%
    summarize(Count = n()) %>%
    arrange(Community, InteractionType) %>%
    mutate(Fraction = Count / sum(Count)) %>%
    left_join(monocultureSets_richness) %>%
    ggplot() +
    geom_col(aes(x = CommunityLabel, y = Fraction, fill = InteractionType), color = 1) +
    annotate("text", x = 1:20, y = 1.15, label = monocultureSets_richness$Richness, size = 4) +
    annotate("text", x = 21, y = 1.15, label = c("n. of species"), size = 4, hjust = 0) +
    annotate("segment", x = .5, xend = 25, y = 1.1, yend = 1.1, color = "black") +
    annotate("text", x = 21, y = 1.05, label = c("n. of tested pairs"), size = 4, hjust = 0) +
    annotate("text", x = 1:20, y = 1.05, label = monocultureSets_richness$PairSize, size = 4) +
    scale_fill_manual(values = interaction_color) +
    scale_x_continuous(breaks = 1:20) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 1.2), breaks = seq(0, 1, .2)) +
    coord_cartesian(xlim = c(0.5, 20.5), ylim = c(0, 1), clip = "off") +
    theme_classic() +
    theme(plot.margin = unit(c(2,.5,.5,.5), "cm"))  +
    guides(color = "none", fill = guide_legend(title = "")) +
    labs(x = "simulated random species set") +
    ggtitle("Pool pairs")


# 2. Community pairs ----
withinCommunityPairs_N_freq <- read_csv(paste0(folder_simulation, "aggregated/13-withinCommunityPairs_N_freq.csv"), col_types = cols()) %>%
    mutate(Time = factor(Time, c("init", "end"))) %>%
    mutate(Pair = factor(Pair, paste0("P", 1:1000))) %>%
    mutate(InitialFrequency = factor(InitialFrequency, c(5, 50, 95)))

communities_richness <- read_csv(paste0(folder_simulation, "aggregated/12-communities_richness.csv"), col_types = cols()) %>%
    mutate(Community = factor(Community, paste0("W", 0:(nrow(input_withinCommunityPairs)-1)))) %>%
    mutate(PairSize = choose(Richness, 2)) %>%
    arrange(desc(Richness)) %>%
    slice(1:20) %>%
    mutate(CommunityLabel = 1:20)
withinCommunityPairs_N_outcome <- read_csv(paste0(folder_simulation, "aggregated/13-withinCommunityPairs_N_outcome.csv"), col_types = cols()) %>%
    mutate(Community = factor(Community, paste0("W", 0:(nrow(input_withinCommunityPairs)-1)))) %>%
    # Top 20 richness
    filter(Community %in% communities_richness$Community)


p2 <- withinCommunityPairs_N_outcome %>%
    group_by(Community, InteractionType, .drop = F) %>%
    summarize(Count = n()) %>%
    arrange(Community, InteractionType) %>%
    mutate(Fraction = Count / sum(Count)) %>%
    left_join(communities_richness) %>%
    ggplot() +
    geom_col(aes(x = CommunityLabel, y = Fraction, fill = InteractionType), color = 1) +
    #geom_text(aes(x = Community, label = Richness), y = 0.9) +
    annotate("text", x = 1:20, y = 1.15, label = communities_richness$Richness, size = 4) +
    annotate("text", x = 21, y = 1.15, label = c("n. of species"), size = 4, hjust = 0) +
    annotate("segment", x = .5, xend = 25, y = 1.1, yend = 1.1, color = "black") +
    annotate("text", x = 21, y = 1.05, label = c("n. of tested pairs"), size = 4, hjust = 0) +
    annotate("text", x = 1:20, y = 1.05, label = communities_richness$PairSize, size = 4) +
    scale_fill_manual(values = interaction_color) +
    scale_x_continuous(breaks = 1:20) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 1.2), breaks = seq(0, 1, .2)) +
    coord_cartesian(xlim = c(0.5, 20.5), ylim = c(0, 1), clip = "off") +
    theme_classic() +
    theme(plot.margin = unit(c(2,.5,.5,.5), "cm"))  +
    guides(color = "none", fill = guide_legend(title = "")) +
    labs(x = "simulated community") +
    ggtitle("Within-community pairs")

p <- plot_grid(p1, p2, nrow = 2, align = "v", axis = "rl", labels = LETTERS[c(1,2)], hjust = 0, label_x = 0.01) + paint_white_background()
#p <- p2
ggsave(here::here("plots/FigS19-pairwise_outcomes.png"), p, width = 8, height = 8)




#
poolPair_fraction <- poolPairs_N_outcome %>%
    mutate(Community = factor(Community, monocultureSets_richness$Community)) %>%
    group_by(Community, InteractionType, .drop = F) %>%
    summarize(Count = n()) %>%
    arrange(Community, InteractionType) %>%
    mutate(Fraction = Count / sum(Count)) %>%
    left_join(monocultureSets_richness) %>%
    select(Community, InteractionType, Fraction) %>%
    filter(InteractionType == "exclusion") %>%
    ungroup()

withinCommunityPairs_fraction <- withinCommunityPairs_N_outcome %>%
    mutate(Community = factor(Community, communities_richness$Community)) %>%
    mutate(InteractionType = factor(InteractionType, c("exclusion", "coexistence"))) %>%
    group_by(Community, InteractionType, .drop = F) %>%
    summarize(Count = n()) %>%
    drop_na() %>%
    arrange(Community, InteractionType) %>%
    mutate(Fraction = Count / sum(Count)) %>%
    left_join(communities_richness) %>%
    select(Community, InteractionType, Fraction) %>%
    filter(InteractionType == "exclusion") %>%
    ungroup()

t.test(poolPair_fraction$Fraction, withinCommunityPairs_fraction$Fraction) %>%
    broom::tidy()

wilcox.test(poolPair_fraction$Fraction, withinCommunityPairs_fraction$Fraction)

stderror <- function(x) sd(x)/sqrt(length(x))
mean(poolPair_fraction$Fraction)
stderror(poolPair_fraction$Fraction)

mean(withinCommunityPairs_fraction$Fraction)
stderror(withinCommunityPairs_fraction$Fraction)



