library(tidyverse)
library(cowplot)
library(broom)
source(here::here("processing_scripts/00-metadata.R"))

isolates <- read_csv(paste0(folder_data, "output/isolates_remained.csv"), show_col_types = F)

isolates_rank <- isolates %>%
    group_by(Community) %>%
    select(ExpID, Community, Isolate, Genus, Rank, CommunityESVID, RelativeAbundance)

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
range(list_rho) # rho = [0.32, 0.52]

# all 62 isolates
set.seed(1)
isol <- isolates_rank %>% mutate(Rank_Abundance = rank(-RelativeAbundance, ties.method = "random"))
plot_cor <- function (x) {
    x %>%
    group_by(Rank, Rank_Abundance) %>%
        count(name = "number") %>%
        ggplot(aes(x = Rank_Abundance, y = Rank)) +
        geom_point(aes(size = number), fill = "grey", shape = 21, stroke = 0) +
        geom_abline(slope = 1, intercept = 0, linetype = 2, color = "red") +
        scale_size_continuous(limits = c(1,9), breaks = 1:6, range = c(3, 8)) +
        theme_classic() +
        guides(size = guide_legend(title = "")) +
        labs(x = "ranked matched ESV abundance", y = "competition rank")
}
p2 <- isol %>%
    plot_cor() +
    scale_x_continuous(breaks = 1:10, limits = c(0.5,10.5)) +
    scale_y_continuous(breaks = 1:10, limits = c(0.5,10.5))


count(isol)
cor.test(isol$Rank_Abundance, isol$Rank, method = "spearman", alternative = "two.sided", exact = FALSE) %>%
    tidy() # rho = 0.437, p = 0.000380
nrow(isol) # 62

# Small communities (n<=7)
isol1 <- isol %>% filter(n() <= 7)
p3 <- isol1 %>%
    plot_cor() +
    scale_x_continuous(breaks = 1:7, limits = c(0.5,7.5)) +
    scale_y_continuous(breaks = 1:7, limits = c(0.5,7.5))

count(isol1)
cor.test(isol1$Rank_Abundance, isol1$Rank, method = "spearman", alternative = "two.sided", exact = FALSE) %>%
    tidy() # rho = 0.395; p = 0.00867
nrow(isol1) # 43

# Large communities (n= 9, 10)
isol2 <- isol %>% filter(n() > 7)
p4 <- isol2 %>%
    plot_cor() +
    scale_x_continuous(breaks = 1:10, limits = c(0.5,10.5)) +
    scale_y_continuous(breaks = 1:10, limits = c(0.5,10.5))
isol2 <- isol %>% filter(n() >= 9)
count(isol2)
cor.test(isol2$Rank_Abundance, isol2$Rank, method = "spearman", alternative = "two.sided", exact = FALSE) %>%
    tidy() # rho = 0.0248; p = 0.920
nrow(isol2) # 19

p <- plot_grid(p1, p2, p3 + guides(size = "none"), p4 + guides(size = "none"), nrow = 2, scale = 0.9, labels = LETTERS[1:4], axis = "lr", align = "v") + paint_white_background()
ggsave(here::here("plots/FigS14.png"), p, width = 8, height = 6)


# Number of strains that match to the same ESV
isolates_rank %>%
    group_by(Community, CommunityESVID) %>%
    summarize(n_strains = n()) %>%
    arrange(desc(n_strains)) %>%
    pull(n_strains) %>%
    table()










