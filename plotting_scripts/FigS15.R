library(tidyverse)
library(cowplot)
library(broom)
source(here::here("processing_scripts/00-metadata.R"))

isolates <- read_csv(paste0(folder_data, "output/isolates_remained.csv"), show_col_types = F)
communities <- read_csv(paste0(folder_data, "output/communities_remained.csv"), show_col_types = F)
isolates_growth_syl <- read_csv(paste0(folder_data, "raw/growth_rate/Estrela_2021_isolates_grmax.csv"), col_types = cols()) %>%
    filter(cs == "glucose") %>%
    rename(ID = SangerID)


# Ranked glucose growth rate
isolates_rank <- isolates %>%
    left_join(isolates_growth_syl) %>%
    group_by(Community) %>%
    # Remove isolates that do not have growth rate data
    drop_na(gr_max) %>%
    select(ID, Community, Isolate, Rank, gr_max)

nrow(isolates_rank) # 56 strains have growth rate data

# Six strains without growth rate data
isolates %>%
    left_join(isolates_growth_syl) %>%
    group_by(Community) %>%
    filter(is.na(gr_max)) %>%
    select(ExpID, ID, Community, Isolate, Family, Genus)

# growth rate
p1 <- isolates_rank %>%
    ggplot(aes(x = gr_max)) +
    geom_histogram(color = "black", fill = "white", binwidth = 0.1) +
    theme_classic() +
    theme() +
    guides() +
    labs(x = expression(growth~rate(h^-1)))

# competitive rank vs growth rate
p2 <- isolates_rank %>%
    ggplot(aes(x = gr_max, y = Rank)) +
    geom_point(shape = 21, size = 2, stroke = 1) +
    geom_smooth(method = "lm", color = "red", linetype = 2, se = F) +
    scale_x_continuous(breaks = seq(0,1.2,0.4), limits = c(-0.05, 1.25)) +
    scale_y_continuous(breaks = 1:9, limits = c(0.5,9.5)) +
    theme_classic() +
    theme() +
    guides() +
    labs(x = expression(growth~rate(h^-1)), y = "competition rank")

p <- plot_grid(p1, p2, nrow = 1, scale = 0.9, labels = LETTERS[1:2]) + paint_white_background()
ggsave(here::here("plots/FigS15.png"), p, width = 8, height = 4)
cor.test(isolates_rank$Rank, isolates_rank$gr_max, method = "pearson", alternative = "two.sided", exact = F) %>% tidy() #  r = -0.320  P=0.016
lm(Rank ~ gr_max, data = isolates_rank) # y = -3.25x + 5.43

