library(tidyverse)
library(cowplot)
library(broom)
source(here::here("analysis/00-metadata.R"))

isolates <- read_csv(paste0(folder_data, "output/isolates_remained.csv"), show_col_types = F)
communities <- read_csv(paste0(folder_data, "output/communities_remained.csv"), show_col_types = F)
isolates_growth_syl <- read_csv(paste0(folder_data, "raw/growth_rate/Estrela_2021_isolates_grmax.csv"), col_types = cols()) %>%
    filter(cs == "glucose") %>%
    rename(ID = SangerID)

# Ranked glucose growth rate
isolates_rank <- isolates %>%
    left_join(isolates_growth_syl) %>%
    group_by(Community) %>%
    # Rank glucose maximum growth rate
    drop_na(gr_max) %>%
    select(ExpID, Community, Isolate, Family, Genus, Game, Win, Lose, Score, Rank, gr_max) %>%
    mutate(Rank_intersect = rank(Rank, ties.method = "min")) %>%
    mutate(rank_gr_max = rank(-gr_max)) # Top growth rate is rank 1


# Growth vs. rank, normalized within community
isolates_norm <- isolates_rank %>%
    group_by(Community) %>%
    mutate(RankNorm = Rank_intersect / n(), rank_gr_max_norm = rank_gr_max / n(), nn = n())

p <- isolates_norm %>%
    ggplot(aes(x = rank_gr_max_norm, y = RankNorm)) +
    geom_abline(slope = 1, intercept = 0, color = "red", linetype = 2) +
    geom_point(shape = 21, size = 2, stroke = 1) +
    scale_x_continuous(limits = c(0,1)) +
    scale_y_continuous(limits = c(0,1)) +
    theme_classic() +
    theme() +
    guides() +
    labs(x = "ranked growth rate (normalized)", y = "competition rank (normalized)")

ggsave(here::here("plots/FigS15-growth_vs_rank.png"), p, width = 4, height = 4)


cor.test(isolates_norm$RankNorm, isolates_norm$rank_gr_max_norm, method = "spearman", alternative = "two.sided", exact = F) %>%
    tidy()



