library(tidyverse)
library(cowplot)
library(broom)
source(here::here("analysis/00-metadata.R"))

isolates <- read_csv(paste0(folder_data, "output/isolates.csv"), show_col_types = F)
communities <- read_csv(paste0(folder_data, "temp/00c-communities.csv"), show_col_types = F)
isolates_growth_syl <- read_csv(paste0(folder_data, "raw/growth_rate/Estrela_2021_isolates_grmax.csv"), col_types = cols()) %>%
    filter(cs == "glucose") %>%
    mutate(ID = as.character(SangerID))

# Ranked glucose growth rate
isolates <- isolates %>%
    left_join(isolates_growth_syl) %>%
    #left_join(isolates_growth_jean) %>%
    group_by(Community) %>%
    # Rank glucose maximum growth rate
    drop_na(gr_max) %>%
    mutate(rank_gr_max = rank(-gr_max)) # Top growth rate is rank 1

# Growth vs. rank, normalized within community
isolates_norm <- isolates %>%
    group_by(Community) %>%
    left_join(communities) %>%
    mutate(RankNorm = Rank / CommunitySize, rank_gr_max_norm = rank_gr_max / CommunitySize)

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

ggsave(here::here("plots/FigS12-growth_vs_rank.png"), p, width = 4, height = 4)


cor.test(isolates_norm$RankNorm, isolates_norm$rank_gr_max_norm,
         method = "pearson", alternative = "two.sided") %>%
    tidy()




