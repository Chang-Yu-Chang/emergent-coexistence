library(tidyverse)
library(cowplot)
source(here::here("analysis/00-metadata.R"))

communities <- read_csv(paste0(folder_data, "temp/00c-communities.csv"), show_col_types = F)
pairs <- read_csv(paste0(folder_data, "output/pairs.csv"), show_col_types = F)
pairs <- remove_ineligible_pairs(pairs)

# Pairwise competition between highly abundant species. Use the Islate sanger where each matches to one ESV
isolates_abundance_all_sanger <- read_csv(paste0(folder_data, "temp/32-isolates_abundance_all_sanger.csv"), show_col_types = F)

isolates_abundance_all_sanger <- isolates_abundance_all_sanger %>%
    # Remove the two isolates with bad alignment. Threshold = 8
    filter(BasePairMismatch <= 4) %>%
    select(Community, Isolate, RelativeAbundance) %>%
    # Highly abundant isolate species means they have  >0.05 relative abundance
    filter(RelativeAbundance > 0.05)

isolates_abundant_richness <- isolates_abundance_all_sanger %>%
    group_by(Community) %>%
    summarize(Richness = n()) %>%
    left_join(communities) %>%
    arrange(CommunityLabel)

pairs_abundant <- pairs %>%
    select(PairID, Community, Isolate1, Isolate2, outcome) %>%
    left_join(rename(isolates_abundance_all_sanger, Isolate1 = Isolate, RelativeAbundance1 = RelativeAbundance)) %>%
    left_join(rename(isolates_abundance_all_sanger, Isolate2 = Isolate, RelativeAbundance2 = RelativeAbundance)) %>%
    filter(!is.na(RelativeAbundance1), !is.na(RelativeAbundance2)) %>%
    left_join(pairs)

pairs_abundant_richness <- pairs_abundant %>%
    group_by(Community) %>%
    summarize(Richness = n()) %>%
    left_join(communities) %>%
    arrange(CommunityLabel)


p <- pairs_abundant %>%
    group_by(Community, outcome) %>%
    count(name = "Count") %>%
    # Total count
    group_by(Community) %>% mutate(Fraction = Count / sum(Count), TotalCount = sum(Count)) %>%
    left_join(communities, by = "Community") %>%
    ungroup() %>%
    ggplot() +
    geom_col(aes(x = CommunityLabel, fill = outcome, y = Fraction), color = 1, width = .8, linewidth = .5) +
    annotate("text", x = 1:13, y = 1.15, label = communities$CommunitySize, size = 4) +
    annotate("text", x = 14, y = 1.15, label = "n. of species", size = 4, hjust = 0) +
    annotate("segment", x = .5, xend = 18, y = 1.1, yend = 1.1, color = "black") +
    geom_text(aes(x = CommunityLabel, y = 1.05, label = TotalCount), size = 4) +
    annotate("text", x = 14, y = 1.05, label = "n. of tested pairs", size = 4, hjust = 0) +
    scale_fill_manual(values = outcome_colors, labels = outcome_labels) +
    scale_x_continuous(breaks = 1:13, expand = c(0.01, 0)) +
    scale_y_continuous(breaks = seq(0,1,0.2), limit = c(0, 1.3), expand = c(0,0)) +
    coord_cartesian(xlim = c(0.5, 13.5), ylim = c(0, 1), clip = "off") +
    theme_classic() +
    theme(
        legend.text = element_text(size = 10),
        legend.title = element_blank(),
        legend.key.size = unit(.5, "cm"),
        legend.spacing.y = unit(.3, "cm"),
        legend.position = "right",
        panel.border = element_rect(color = 1, fill = NA),
        axis.text = element_text(color = 1, size = 10),
        axis.title = element_text(color = 1, size = 10),
        plot.margin = unit(c(1,.5,.5,.5), "cm")
    ) +
    guides(fill = guide_legend(byrow = TRUE)) +
    labs(x = "community", y = "fraction")

ggsave(here::here("plots/FigS10-pairwise_competition_abundant.png"), p, width = 6, height = 3)


# Stats
nrow(pairs_abundant)
pairs_abundant %>%
    group_by(outcome) %>%
    count(name = "Count") %>%
    ungroup() %>%
    mutate(Fraction = Count / sum(Count), TotalCount = sum(Count))
