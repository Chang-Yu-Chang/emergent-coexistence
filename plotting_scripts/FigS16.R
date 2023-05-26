library(tidyverse)
library(cowplot)
library(broom)
source(here::here("processing_scripts/00-metadata.R"))

isolates <- read_csv(paste0(folder_data, "output/isolates_remained.csv"), show_col_types = F)
#isolates_rank <- read_csv("~/Downloads/isolates_rank.csv")

k1 <- isolates %>% filter(Fermenter) %>% pull(Rank)
mean(k1) # 2.54
k2 <- isolates %>% filter(!Fermenter) %>% pull(Rank)
mean(k2) # 4.72

# k2 <- isolates_rank %>% filter(!Fermenter) %>% pull(Rank)
# mean(k2) # 4.36

isolates %>%
    mutate(Rank_Abundance = rank(-RelativeAbundance, ties.method = "average")) %>%
    wilcox.test(Rank ~ Fermenter, data = ., exact = F) # p=0.000167

p <- isolates %>%
    mutate(Fermenter = ifelse(Fermenter, "respiro-fermenter", "respirator")) %>%
    mutate(Fermenter = factor(Fermenter, c("respiro-fermenter", "respirator"))) %>%
    ggplot() +
    geom_boxplot(aes(x = Fermenter, y = Rank), outlier.shape = NA) +
    geom_jitter(aes(x = Fermenter, y = Rank), width = 0.3, shape = 21, size = 2) +
    scale_y_continuous(breaks = 1:10) +
    theme_classic() +
    theme() +
    guides() +
    labs(x = "", y = "competitive rank")

ggsave(here::here("plots/FigS16.png"), p, width = 4, height = 4)

table(isolates$Fermenter)

