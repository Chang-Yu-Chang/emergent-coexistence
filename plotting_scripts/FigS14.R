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
    tidy()

#
set.seed(1)
list_rho <- rep(NA, 1000)
for (i in 1:1000) {
    isol <- isolates_rank %>% mutate(Rank_Abundance = rank(-RelativeAbundance, ties.method='random'))
    # isol <- isol %>%
    #     group_by(Community) %>%
    #     mutate(Rank_Norm = Rank / n(), Rank_Abundance_Norm = Rank_Abundance / n())

    list_rho[i] <- cor.test(isol$Rank_Abundance_Norm, isol$Rank_Norm, method = "spearman", alternative = "two.sided", exact = FALSE) %>%
        tidy() %>%
        `[`("estimate") %>% unlist
}


# isol %>%
#     ggplot() +
#     geom_jitter(aes(x = Rank_Abundance_Norm, y = Rank_Norm), shape = 21, width = 0.1, height = 0.1) +
#     theme_classic()


isolates_rank <- isolates %>%
    group_by(Community) %>%
    mutate(Rank_Abundance = rank(-RelativeAbundance, ties.method='random')) %>%
    mutate(Rank_Norm = Rank / n(), Rank_Abundance_Norm = Rank_Abundance / n())



# p1 <- isolates_rank %>%
#     ggplot(aes(x = Rank_Abundance_Norm, y = Rank_Norm)) +
#     geom_point(shape = 21, size = 2, stroke = 1, position = position_jitter(width = 0.15, height = 0.15)) +
#     geom_abline(slope = 1, intercept = 0, linetype = 2, color = "red") +
#     scale_x_continuous(limits = c(0,1)) +
#     scale_y_continuous(limits = c(0,1)) +
#     theme_classic() +
#     labs(x = "ranked matched ESV abundance (normalized)", y = "competition rank (normalized)")

p1 <- tibble(rho = list_rho, BootStrapID = 1:1000) %>%
    ggplot() +
    geom_histogram(aes(x = rho), color = "black", fill = "white") +
    scale_x_continuous(limits = c(0,1)) +
    theme_classic() +
    theme() +
    guides() +
    labs(x = expression(rho), y = "count (# of bootstrap samples)")

p2 <-
    isolates_rank %>%
    ggplot(aes(x = Rank_Abundance_Norm, y = Rank_Norm)) +
    geom_point(shape = 21, size = 2, stroke = 1, position = position_jitter(width = 0.15, height = 0.15)) +
    geom_abline(slope = 1, intercept = 0, linetype = 2, color = "red") +
    scale_x_continuous(limits = c(0,1)) +
    scale_y_continuous(limits = c(0,1)) +
    theme_classic() +
    labs(x = "ranked matched ESV abundance", y = "competition rank")

isolates_rank %>%
    mutate(Rank_Abundance = rank(-RelativeAbundance, ties.method = "average")) %>%
    filter(n() >= 6)
cor.test(isol$Rank_Abundance, isol$Rank, method = "spearman", alternative = "two.sided", exact = FALSE) %>%
    tidy()


#p <- plot_grid(p1, p2, nrow = 1, scale = 0.85, labels = c("A", "B")) + paint_white_background()
p <- p1


ggsave(here::here("plots/FigS14.png"), p, width = 4, height = 4)


# Range of rho
range(list_rho)

# NUmber of strains that match to the same ESV
isolates_rank %>%
    group_by(Community, CommunityESVID) %>%
    summarize(n_strains = n()) %>%
    arrange(desc(n_strains)) %>%
    pull(n_strains) %>%
    table()










