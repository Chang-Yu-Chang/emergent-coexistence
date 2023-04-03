library(tidyverse)
library(cowplot)
source(here::here("analysis/00-metadata.R"))

communities <- read_csv(paste0(folder_data, "output/communities_remained.csv"), show_col_types = F)
isolates <- read_csv(paste0(folder_data, "output/isolates_remained.csv"), show_col_types = F)
pairs <- read_csv(paste0(folder_data, "output/pairs_remained.csv"), show_col_types = F)
communities_abundance_T12 <- read_csv(paste0(folder_data, "temp/32-communities_abundance_T12.csv"), show_col_types = F)

#
isolates_abundant <- isolates %>% filter(RelativeAbundance > 0.05)
nrow(isolates_abundant) # 51 isolates matched to abundant ESVs

isolates_abundant_richness <- isolates_abundant %>%
    group_by(Community) %>%
    summarize(CommunitySizeAbundant = n())

pairs_abundant <- pairs %>%
    select(PairID, Community, Isolate1, Isolate2, outcome) %>%
    left_join(select(isolates_abundant, Community, Isolate1 = Isolate, RelativeAbundance1 = RelativeAbundance)) %>%
    left_join(select(isolates_abundant, Community, Isolate2 = Isolate, RelativeAbundance2 = RelativeAbundance)) %>%
    filter(!is.na(RelativeAbundance1), !is.na(RelativeAbundance2))
nrow(pairs_abundant) # 85 pairs

pairs_abundant_richness <- pairs_abundant %>%
    group_by(Community) %>%
    summarize(CommunityPairSizeAbundant = n())

n_ESVs <- communities_abundance_T12 %>%
    group_by(Community) %>%
    count(name = "ESVRichness")

communities_abundant <- communities %>%
    left_join(isolates_abundant_richness) %>%
    left_join(pairs_abundant_richness) %>%
    left_join(n_ESVs) %>%
    select(CommunityLabel, CommunitySize, CommunitySizeAbundant, CommunityPairSize, CommunityPairSizeAbundant, ESVRichness)

sum(communities_abundant$CommunitySize) # 64 isolates
sum(communities_abundant$CommunitySizeAbundant) # 51 abundant isolates
sum(communities_abundant$CommunityPairSize) # 145 pairs
sum(communities_abundant$CommunityPairSizeAbundant) # 85 pairs of abundant isolates

#
p <- pairs_abundant %>%
    group_by(Community, outcome) %>%
    count(name = "Count") %>%
    # Total count
    group_by(Community) %>% mutate(Fraction = Count / sum(Count), TotalCount = sum(Count)) %>%
    left_join(communities, by = "Community") %>%
    #mutate(outcome = factor(outcome, c("5-inconclusive", "1-exclusion", "2-exclusion", "3-coexistence", "4-coexistence"))) %>%
    ungroup() %>%
    ggplot() +
    geom_col(aes(x = CommunityLabel, fill = outcome, y = Fraction), color = 1, width = .8, linewidth = .5, position = position_stack(reverse = T)) +
    # Number of distcint ESVs
    annotate("text", x = 14, y = 1.25, label = "n. of ESVs", size = 4, hjust = 0) +
    annotate("text", x = 1:13, y = 1.25, label = communities_abundant$ESVRichness, size = 4) +
    annotate("segment", x = .5, xend = 18, y = 1.2, yend = 1.2, color = "black") +
    # Number of isolates
    annotate("text", x = 14, y = 1.15, label = "n. of isolates", size = 4, hjust = 0) +
    annotate("text", x = 1:13, y = 1.15, label = communities_abundant$CommunitySizeAbundant, size = 4) +
    annotate("segment", x = .5, xend = 18, y = 1.1, yend = 1.1, color = "black") +
    # Number of tested pairs
    annotate("text", x = 14, y = 1.05, label = "n. of tested pairs", size = 4, hjust = 0) +
    annotate("text", x = 1:13, y = 1.05, label = communities_abundant$CommunityPairSizeAbundant, size = 4) +
    scale_fill_manual(values = outcome_colors, breaks = names(outcome_colors), labels = outcome_labels) +
    scale_x_continuous(breaks = 1:13, expand = c(0.01, 0)) +
    scale_y_continuous(breaks = seq(0,1,0.2), limit = c(0, 1.45), expand = c(0,0)) +
    coord_cartesian(xlim = c(0.5, 13.5), ylim = c(0, 1), clip = "off") +
    theme_classic() +
    theme(
        legend.text = element_text(size = 12),
        legend.title = element_blank(),
        #legend.key.size = unit(5, "mm"),
        legend.key.width = unit(8, "mm"),
        legend.key.height = unit(8, "mm"),
        legend.spacing.y = unit(8, "mm"),
        #legend.margin = margin(10,10,1,0, unit = "mm"),
        legend.position = "right",
        panel.border = element_rect(color = 1, fill = NA),
        axis.text = element_text(color = 1, size = 12),
        axis.title = element_text(color = 1, size = 12),
        plot.margin = unit(c(40, 20, 5, 5), "mm")
    ) +
    #guides(fill = guide_legend(byrow = T, ncol = 2)) +
    guides(fill = guide_legend(byrow = T, ncol = 1)) +
    labs(x = "community", y = "fraction")

ggsave(here::here("plots/FigS12-pairwise_competition_abundant.png"), p, width = 10, height = 6)


# Stats
nrow(pairs_abundant)
pairs_abundant %>%
    group_by(outcome) %>%
    count(name = "Count") %>%
    ungroup() %>%
    mutate(Fraction = Count / sum(Count), TotalCount = sum(Count))






