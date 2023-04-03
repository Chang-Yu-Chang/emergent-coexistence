library(tidyverse)
library(cowplot)
library(broom)
source(here::here("analysis/00-metadata.R"))

isolates <- read_csv(paste0(folder_data, "output/isolates.csv"), show_col_types = F)

isolates_rank <- isolates %>%
    group_by(Community) %>%
    arrange(Community, desc(RelativeAbundance)) %>%
    drop_na() %>%
    mutate(RankRelativeAbundance = 1:n()) %>%
    left_join(communities) %>%
    group_by(Community) %>%
    mutate(RankNorm = Rank / CommunitySize, RankRelativeAbundanceNorm = RankRelativeAbundance / CommunitySize)

p <- isolates_rank %>%
    ggplot(aes(x = RankRelativeAbundance, y = Rank)) +
    geom_point(shape = 21, size = 2, stroke = 1, position = position_jitter(width = 0.15, height = 0.15)) +
    geom_abline(slope = 1, intercept = 0, linetype = 2, color = "red") +
    scale_x_continuous(breaks = 1:12, limits = c(1,12)) +
    scale_y_continuous(breaks = 1:12, limits = c(1,12)) +
    theme_classic() +
    #labs(x = "Ranked ESV abundance", y = "Competitive rank")
    labs()
ggsave(here::here("plots/FigS13-abundance_vs_rank.png"), p, width = 4, height = 4)

cor.test(isolates_rank$RankRelativeAbundance, isolates_rank$Rank,
         method = "spearman", alternative = "two.sided", exact = FALSE) %>%
    tidy()
