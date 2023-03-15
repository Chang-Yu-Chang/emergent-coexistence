library(tidyverse)
library(cowplot)
source(here::here("analysis/00-metadata.R"))
source(here::here("plotting_scripts/FigS10.R"))


# Calculate Malthusian fitness
communities_abundance_fitness <- communities_abundance_sp %>%
    filter(Transfer %in% c(1, 9:12)) %>%
    mutate(Time = case_when(
        Transfer == 1 ~ "init",
        Transfer %in% 9:12 ~ "end"
    )) %>%
    group_by(Community, ESV_ID, Time) %>%
    summarize(Relative_Abundance = mean(Relative_Abundance, na.rm = T)) %>%
    filter(Community %in% c("C1R4", "C2R6", "C2R8", "C8R4")) %>%
    arrange(Community, ESV_ID, Time) %>%
    pivot_wider(names_from = Time, names_prefix = "T", values_from = Relative_Abundance) %>%
    drop_na() %>%
    mutate(Fitness = log(Tend/Tinit)) %>%
    select(Community, ESV_ID, Fitness)

# x_T1 vs. log(x_T12/x_T1)
temp <- communities_abundance_T1 %>%
    select(Community, ESV_ID, Relative_Abundance) %>%
    filter(Community %in% c("C1R4", "C2R6", "C2R8", "C8R4")) %>%
    group_by(Community, ESV_ID) %>%
    summarize(Relative_Abundance = sum(Relative_Abundance))

p <- communities_abundance_fitness %>%
    left_join(temp) %>%
    ggplot(aes(x = Relative_Abundance, y = Fitness)) +
    geom_point(shape = 21, size = 3, stroke = 1) +
    geom_hline(yintercept = 0, linetype = 2) +
    theme_classic() +
    theme(axis.text = element_text(size = 15),
          axis.title = element_text(size = 15)) +
    labs(x = expression(x[init]), y = expression(log(x[end]/x[init])))

ggsave(here::here("plots/FigS11-species_abundance.png"), p, width = 4, height = 4)

if (FALSE) {
    # log(x[Ti+1]/x[i]) over transfers
    communities_abundance_time <- communities_abundance_sp %>%
        select(ESV_ID, Transfer, Relative_Abundance) %>%
        arrange(ESV_ID, Transfer) %>%
        group_by(ESV_ID) %>%
        # time step change
        mutate(Fitness = log(lead(Relative_Abundance) / Relative_Abundance))

    communities_abundance_time %>%
        ggplot(aes(x = Relative_Abundance, y = Fitness)) +
        geom_point(shape = 21, size = 3, stroke = 1) +
        geom_hline(yintercept = 0, linetype = 2) +
        #facet_wrap(~ESV_ID) +
        #facet_wrap(~Transfer, scales = "free") +
        theme_classic() +
        theme(panel.border = element_rect(color = 1, fill = NA)) +
        labs(x = "x_T", y = "log(x_{Ti+1}/x_{Ti})")

    p3 <- communities_abundance_time %>%
        filter(Transfer != 12) %>%
        ggplot() +
        geom_boxplot(aes(x = Transfer, y = Fitness, group = Transfer), outlier.color = NA) +
        geom_hline(yintercept = 0, linetype = 2) +
        geom_point(aes(x = Transfer, y = Fitness), position = position_jitter(width = 0.1, height = 0),
                   size = 2, shape = 21, stroke = 1) +
        scale_x_continuous(breaks = 1:12) +
        theme_classic() +
        theme() +
        labs(x = "Transfer", y = expression(log(x[Ti+1]/x[i])))

}
