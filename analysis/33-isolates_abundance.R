library(tidyverse)
library(cowplot)
library(broom)
source(here::here("analysis/00-metadata.R"))

communities <- read_csv(paste0(folder_data, "temp/00c-communities.csv"), show_col_types = F)
isolates <- read_csv(paste0(folder_data, "output/isolates.csv"), show_col_types = F)

isolates_rank <- isolates %>%
    group_by(Community) %>%
    #select(Community, Isolate, Rank, RelativeAbundance) %>%
    arrange(Community, desc(RelativeAbundance)) %>%
    drop_na() %>%
    mutate(RankRelativeAbundance = 1:n()) %>%
    left_join(communities) %>%
    group_by(Community) %>%
    mutate(RankNorm = Rank / CommunitySize, RankRelativeAbundanceNorm = RankRelativeAbundance / CommunitySize)


# Check that number of datapoint. It should be 47, the number of ESVs matched to Sanger
nrow(isolates_rank)

# Plot the correlation
p <- isolates_rank %>%
    ggplot(aes(x = RankRelativeAbundance, y = Rank)) +
    geom_point(shape = 21, size = 2, stroke = 1, position = position_jitter(width = 0.15, height = 0.15)) +
    geom_abline(slope = 1, intercept = 0, linetype = 2, color = "red") +
    scale_x_continuous(breaks = 1:12, limits = c(1,12)) +
    scale_y_continuous(breaks = 1:12, limits = c(1,12)) +
    theme_classic() +
    labs()

ggsave(paste0(folder_data, "temp/33-01-abundance_vs_rank.png"), p, width = 4, height = 4)

cor.test(isolates_rank$RankRelativeAbundance, isolates_rank$Rank,
         method = "spearman", alternative = "two.sided", exact = FALSE) %>%
    tidy()

# Check if the abundance is indeed is ranked
p <- isolates_rank %>%
    ggplot() +
    geom_point(aes(x = RankRelativeAbundance, y = RelativeAbundance, color = Family), shape = 21, size = 2, stroke = 1) +
    scale_color_manual(values = RColorBrewer::brewer.pal(7, "Set2")) +
    facet_wrap(~Community, ncol = 3) +
    theme_bw() +
    theme() +
    guides() +
    labs()
ggsave(paste0(folder_data, "temp/33-02-rankAbundance_vs_abundance_comm.png"), p, width = 6, height = 8)


# Check if the more abundant strains are E and P
p <- isolates_rank %>%
    ggplot() +
    geom_point(aes(x = RankRelativeAbundance, y = Rank, color = Family), shape = 21, size = 2, stroke = 1,
               position = position_jitter(width = 0.05, height = 0.05)) +
    scale_color_manual(values = RColorBrewer::brewer.pal(7, "Set2")) +
    theme_bw() +
    theme() +
    guides() +
    labs()
ggsave(paste0(folder_data, "temp/33-03-rankAbundance_vs_abundance_family.png"), p, width = 5, height = 3)

# Normalize the rank by the number of strains in a community
p <- isolates_rank %>%
    ggplot() +
    geom_abline(slope = 1, intercept = 0, linetype = 2, color = "red") +
    geom_point(aes(x = RankRelativeAbundanceNorm, y = RankNorm, color = Family), shape = 21, size = 2, stroke = 1) +
    scale_color_manual(values = RColorBrewer::brewer.pal(7, "Set2")) +
    theme_classic() +
    theme() +
    guides() +
    labs()

ggsave(paste0(folder_data, "temp/33-04-rankAbundance_vs_abundance_family_norm.png"), p, width = 5, height = 3)

cor.test(isolates_rank$RankRelativeAbundanceNorm, isolates_rank$RankNorm,
         method = "pearson", alternative = "two.sided", exact = FALSE) %>%
    tidy()




