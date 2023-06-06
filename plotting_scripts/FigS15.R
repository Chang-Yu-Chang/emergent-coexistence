library(tidyverse)
library(cowplot)
library(broom)
source(here::here("processing_scripts/00-metadata.R"))

isolates <- read_csv(paste0(folder_data, "output/isolates_remained.csv"), show_col_types = F)
communities <- read_csv(paste0(folder_data, "output/communities_remained.csv"), show_col_types = F)

# growth rate
p1 <- isolates %>%
    ggplot(aes(x = r)) +
    geom_histogram(color = "black", fill = "white", binwidth = 0.1) +
    scale_x_continuous(breaks = seq(0,1.2,0.4), limits = c(-0.05, 1.25)) +
    scale_y_continuous(breaks = 1:12, limits = c(0, 13)) +
    theme_classic() +
    theme() +
    guides() +
    labs(x = expression(growth~rate(h^-1)))

# competitive rank vs growth rate
p2 <- isolates %>%
    ggplot(aes(x = r, y = Rank)) +
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
cor.test(isolates$Rank, isolates$r, method = "pearson", alternative = "two.sided", exact = F) %>%
    tidy() #  r = -0.314  P=0.0129
lm(Rank ~ r, data = isolates) # y = -3.177x + 5.389

