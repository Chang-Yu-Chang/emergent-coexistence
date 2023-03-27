library(tidyverse)
library(cowplot)
library(broom)
source(here::here("analysis/00-metadata.R"))

pairs <- read_csv(paste0(folder_data, "output/pairs.csv"), show_col_types = F)
communities <- read_csv(paste0(folder_data, "temp/00c-communities.csv"), show_col_types = F)
pairs_freq <- read_csv(paste0(folder_data, "temp/93-pairs_freq.csv"), show_col_types = F)


# OD_B2 <- read_csv(paste0(folder_data, "raw/OD/OD_B2.csv"), show_col_types = F)
# OD_C <- read_csv(paste0(folder_data, "raw/OD/OD_C.csv"), show_col_types = F)
# OD_C2 <- read_csv(paste0(folder_data, "raw/OD/OD_C2.csv"), show_col_types = F)
OD_D <- read_csv(paste0(folder_data, "raw/OD/OD_D.csv"), show_col_types = F)

p <- OD_D %>%
    filter(Wavelength == 620) %>%
    filter(Isolate1 != "blank") %>%
    filter((Isolate1Freq == 95 & MixPlate == 2) | (Isolate1Freq == 50 & MixPlate == 1)) %>%
    mutate(Isolate1 = as.numeric(Isolate1), Isolate2 = as.numeric(Isolate2)) %>%
    left_join(select(pairs, PairID, Community, Isolate1, Isolate2)) %>%
    mutate(UniqeWell = paste0(Layout, Well, Isolate1Freq)) %>%
    ggplot(aes(x = Transfer, y = Abs, group = UniqeWell)) +
    geom_point(shape = 21) +
    geom_line(linewidth = 0.1) +
    scale_x_continuous(breaks = 1:8) +
    theme_classic() +
    guides(color = "none") +
    labs(x = "transfer", y = expression(OD[620]))
ggsave(here::here("plots/FigS92-pairs_ID.png"), p, width = 6, height = 4)
