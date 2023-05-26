library(tidyverse)
library(cowplot)
library(broom)
source(here::here("processing_scripts/00-metadata.R"))

isolates <- read_csv(paste0(folder_data, "output/isolates_remained.csv"), show_col_types = F)
communities <- read_csv(paste0(folder_data, "output/communities_remained.csv"), show_col_types = F)
isolates_growth_syl <- read_csv(paste0(folder_data, "raw/growth_rate/Estrela_2021_isolates_grmax.csv"), col_types = cols()) %>%
    filter(cs == "glucose") %>%
    rename(ID = SangerID)
#isolates_rank <- read_csv("~/Downloads/isolates_rank.csv")

# isolates %>%
#     select(ID, Community, Isolate, Rank) %>%
#     head(20)
# isolates_rank %>%
#     select(ID, Community, Isolate, Rank) %>%
#     head(20)

# Ranked glucose growth rate
isolates_rank <- isolates %>%
    left_join(isolates_growth_syl) %>%
    group_by(Community) %>%
    # Rank glucose maximum growth rate
    drop_na(gr_max) %>%
    select(Community, Isolate, Rank, gr_max)
    #select(ExpID, Community, Isolate, Family, Genus, Game, Win, Lose, Score, Rank, gr_max) %>%
    #mutate(Rank_intersect = rank(Rank, ties.method = "average")) %>%
    #mutate(Rank = rank(Rank, ))
    #mutate(rank_gr_max = rank(-gr_max, ties.method = "average")) # Top growth rate is rank 1

nrow(isolates_rank) # 56 strains have growth rate data
cor.test(isolates_rank$Rank, -isolates_rank$gr_max, method = "pearson", alternative = "two.sided", exact = F)

# growth rate
p1 <- isolates_rank %>%
    ggplot(aes(x = gr_max)) +
    geom_histogram(color = "black", fill = "white", binwidth = 0.1) +
    theme_classic() +
    theme() +
    guides() +
    labs(x = expression(growth~rate(h^-1)))

# comeptitive rank vs grwoth rate rank
p2 <- isolates_rank %>%
    ggplot(aes(x = rank_gr_max, y = Rank)) +
    geom_abline(slope = 1, intercept = 0, color = "red", linetype = 2) +
    #geom_point(shape = 21, size = 2, stroke = 1) +
    geom_point(shape = 21, size = 2, stroke = 1, position = position_jitter(width = 0.1, height = 0.1)) +
    scale_x_continuous(breaks = 1:10, limits = c(0.5,10.5)) +
    scale_y_continuous(breaks = 1:10, limits = c(0.5,10.5)) +
    theme_classic() +
    theme() +
    guides() +
    labs(x = "ranked growth rate", y = "competition rank")

p <- plot_grid(p1, p2, nrow = 1, scale = 0.9, labels = LETTERS[1:2]) + paint_white_background()
ggsave(here::here("plots/FigS15.png"), p, width = 8, height = 4)
cor.test(isolates_rank$Rank, isolates_rank$rank_gr_max, method = "spearman", alternative = "two.sided", exact = F)


if (FALSE) {
# Sylvie's data show all but one
isolates_rank %>%
    ggplot(aes(x = Family, y = gr_max)) +
    geom_boxplot() +
    geom_jitter() +
    theme_classic() +
    theme() +
    guides() +
    labs()

# THe rest 6 strains without growth rate data
isolates %>%
    anti_join(isolates_growth_syl) %>%
    view

isolates_rank %>%
    ggplot() +
    geom_point(aes(x = rank_gr_max, y = gr_max, color = Family)) +
    facet_wrap(~Community, scales = "free_x") +
    scale_x_continuous(breaks = 1:10) +
    theme_classic() +
    theme() +
    guides() +
    labs()


}
