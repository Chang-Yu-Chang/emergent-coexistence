library(tidyverse)
library(cowplot)
source(here::here("analysis/00-metadata.R"))

pairs <- read_csv(paste0(folder_data, "output/pairs.csv"), show_col_types = F)
communities <- read_csv(paste0(folder_data, "temp/00c-communities.csv"), show_col_types = F)

# 0. Clean up column factors ----
# Arrange communities by size
communities <- communities %>%
    arrange(CommunitySize) %>%
    mutate(Community = factor(Community, Community))

# Clean up the pairs data ----
pairs <- pairs %>%
    # Remove no-colony pairs
    unite(col = "Pair", Community, Isolate1, Isolate2, sep = "_", remove = F) %>%
    filter(!(Pair %in% pairs_no_colony)) %>%
    # Remove low-accuracy model pairs
    filter(AccuracyMean > 0.9)


# Figure S6 pairwise competition between highly abundant species ----
# This csv keeps all Sangers that match to one ESV with up to 0-2 mismatch
isolates_abundance_loose <- read_csv(paste0(folder_data, "temp/14-isolates_abundance_loose.csv"), show_col_types = F)

# Highly abundant isolate species have >0.05 relative ESV abundance
isolates_abundant <- isolates_abundance_loose %>%
    select(Community, Isolate, RelativeAbundance) %>%
    filter(RelativeAbundance > 0.05)

isolates_abundant_richness <- isolates_abundant %>%
    group_by(Community) %>%
    summarize(Richness = n()) %>%
    left_join(communities) %>%
    arrange(CommunityLabel)

#
pairs_abundant <- pairs %>%
    select(PairID, Community, Isolate1, Isolate2, InteractionType) %>%
    left_join(rename(isolates_abundant, Isolate1 = Isolate, RelativeAbundance1 = RelativeAbundance)) %>%
    left_join(rename(isolates_abundant, Isolate2 = Isolate, RelativeAbundance2 = RelativeAbundance)) %>%
    filter(!is.na(RelativeAbundance1), !is.na(RelativeAbundance2)) %>%
    left_join(pairs)

pairs_abundant_richness <- pairs_abundant %>%
    group_by(Community) %>%
    summarize(Richness = n()) %>%
    left_join(communities) %>%
    arrange(CommunityLabel)

#
p <- pairs_abundant %>%
    filter(!is.na(FitnessFunction)) %>%
    group_by(Community, InteractionType) %>%
    count(name = "Count") %>%
    group_by(Community) %>% mutate(Fraction = Count / sum(Count), TotalCount = sum(Count)) %>%
    ungroup() %>%
    mutate(Community = factor(Community, communities$Community)) %>%
    arrange(Community) %>%
    left_join(communities, by = "Community") %>%
    mutate(CommunityLabel = factor(CommunityLabel)) %>%
    replace_na(list(InteractionType = "unknown")) %>%
    ggplot() +
    geom_col(aes(x = CommunityLabel, fill = InteractionType, y = Fraction), color = 1, width = .8, linewidth = .5) +
    annotate("text", x = 1:13, y = 1.15, label = isolates_abundant_richness$Richness, size = 4) +
    annotate("text", x = 14, y = 1.15, label = "n. of species", size = 4, hjust = 0) +
    annotate("segment", x = .5, xend = 18, y = 1.1, yend = 1.1, color = "black") +
    annotate("text", x= 1:13, y = 1.05, label = pairs_abundant_richness$Richness, size = 4) +
    annotate("text", x = 14, y = 1.05, label = "n. of tested pairs", size = 4, hjust = 0) +
    scale_fill_manual(values = assign_interaction_color(), breaks = c("coexistence", "exclusion", "unknown")) +
    scale_x_discrete(breaks = 1:13, expand = c(0.01, 0)) +
    scale_y_continuous(breaks = seq(0,1,0.2), limit = c(0, 1.3), expand = c(0,0)) +
    coord_cartesian(xlim = c(0.5, 13.5), ylim = c(0, 1), clip = "off") +
    theme_classic() +
    theme(legend.text = element_text(size = 12),
          axis.text = element_text(color = 1, size = 12),
          #axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          axis.title = element_text(color = 1, size = 12),
          legend.title = element_blank(),
          legend.position = "right",
          plot.margin = unit(c(2,.5,.5,.5), "cm")
    ) +
    labs(x = "Community", y = "Fraction")

ggsave(here::here("plots/FigS10-pairwise_competition_abundant.png"), p, width = 6, height = 3)


# Stats
pairs_abundant %>%
    group_by(InteractionType) %>%
    count(name = "Count") %>%
    ungroup() %>%
    mutate(Fraction = Count / sum(Count))
