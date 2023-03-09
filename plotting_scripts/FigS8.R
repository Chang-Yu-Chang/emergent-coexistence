library(tidyverse)
library(cowplot)
library(broom)
source(here::here("analysis/00-metadata.R"))

isolates <- read_csv(paste0(folder_data, "output/isolates.csv"), show_col_types = F)

isolates_rank <- isolates %>%
    group_by(Community) %>%
    select(Community, Isolate, Rank, RelativeAbundance) %>%
    arrange(Community, desc(RelativeAbundance)) %>%
    drop_na() %>%
    mutate(RankRelativeAbundance = 1:n())

cor.test(isolates_rank$RankRelativeAbundance, isolates_rank$Rank, method = c("pearson")) %>%
    broom::tidy()

p <- isolates_rank %>%
    ggplot(aes(x = RankRelativeAbundance, y = Rank)) +
    geom_point(shape = 21, size = 2, stroke = 0.5, position = position_jitter(width = 0.15, height = 0.15)) +
    scale_x_continuous(breaks = 1:12) +
    scale_y_continuous(breaks = 1:12) +
    theme_classic() +
    labs(x = "Ranked ESV abundance", y = "Competitive rank")

ggsave(here::here("plots/FigS8-abundance_vs_rank.png"), p, width = 3, height = 3)

