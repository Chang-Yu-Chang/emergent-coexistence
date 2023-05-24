library(tidyverse)
library(cowplot)
library(broom)
source(here::here("processing_scripts/00-metadata.R"))

isolates <- read_csv(paste0(folder_data, "output/isolates_remained.csv"), show_col_types = F)

isolates_rank <- isolates %>%
    group_by(Community) %>%
    select(ExpID, Community, Isolate, CommunityESVID, Rank, RelativeAbundance)

isol <- isolates_rank %>% mutate(Rank_Abundance = rank(-RelativeAbundance, ties.method = "average"))
cor.test(isol$Rank_Abundance, isol$Rank, method = "spearman", alternative = "two.sided", exact = FALSE) %>%
    tidy() # rho = 0.425, p = 0.000583

# 1000 resampling
set.seed(1)
list_rho <- rep(NA, 1000)
for (i in 1:1000) {
    isol <- isolates_rank %>% mutate(Rank_Abundance = rank(-RelativeAbundance, ties.method='random'))
    list_rho[i] <- cor.test(isol$Rank_Abundance, isol$Rank, method = "spearman", alternative = "two.sided", exact = FALSE) %>%
        tidy() %>%
        `[`("estimate") %>% unlist
}

p1 <- tibble(rho = list_rho, BootStrapID = 1:1000) %>%
    ggplot() +
    geom_histogram(aes(x = rho), color = "black", fill = "white") +
    scale_x_continuous(limits = c(0,1)) +
    theme_classic() +
    theme() +
    guides() +
    labs(x = expression(rho), y = "count (# of bootstrap samples)")

# Small communities (n<=7)
isol1 <- isol %>% filter(n() <= 7)
p2 <- isol1 %>%
    ggplot(aes(x = Rank_Abundance, y = Rank)) +
    geom_point(shape = 21, size = 2, stroke = 1, position = position_jitter(width = 0.1, height = 0.1)) +
    geom_abline(slope = 1, intercept = 0, linetype = 2, color = "red") +
    scale_x_continuous(breaks = 1:7, limits = c(0.5,7.5)) +
    scale_y_continuous(breaks = 1:7, limits = c(0.5,7.5)) +
    theme_classic() +
    labs(x = "ranked matched ESV abundance", y = "competition rank")

count(isol1)
cor.test(isol1$Rank_Abundance, isol1$Rank, method = "spearman", alternative = "two.sided", exact = FALSE) %>%
    tidy() # rho = 0.316; p = 0.0392
nrow(isol1) # 43

# Large communities (n= 9, 10)
isol2 <- isol %>% filter(n() > 7)
p3 <- isol2 %>%
    ggplot() +
    geom_jitter(aes(x = Rank_Abundance, y = Rank), shape = 21, size = 2, stroke = 1, width = 0.1, height = 0.1) +
    geom_abline(slope = 1, intercept = 0, linetype = 2, color = "red") +
    scale_x_continuous(breaks = 1:10, limits = c(0.5,10.5)) +
    scale_y_continuous(breaks = 1:10, limits = c(0.5,10.5)) +
    theme_classic() +
    labs(x = "ranked matched ESV abundance", y = "competition rank")
isol2 <- isol %>% filter(n() >= 9)
count(isol2)
cor.test(isol2$Rank_Abundance, isol2$Rank, method = "spearman", alternative = "two.sided", exact = FALSE) %>%
    tidy() # rho = 0.0848; p = 0.73
nrow(isol2) # 19

p <- plot_grid(p1, p2, p3, nrow = 2, scale = 0.9, labels = LETTERS[1:3]) + paint_white_background()
ggsave(here::here("plots/FigS14.png"), p, width = 6, height = 6)


# Range of rho
range(list_rho)

# Number of strains that match to the same ESV
isolates_rank %>%
    group_by(Community, CommunityESVID) %>%
    summarize(n_strains = n()) %>%
    arrange(desc(n_strains)) %>%
    pull(n_strains) %>%
    table()










