library(tidyverse)
library(cowplot)
library(broom)
source(here::here("analysis/00-metadata.R"))
source(here::here("plotting_scripts/FigS11.R"))


#
communities_abundance_fitness_equilibrium <- communities_abundance_fitness %>%
    filter(Transfer %in% 8:12) %>%
    drop_na() %>%
    mutate(CommunityESV = paste0(Community, ESV_ID)) %>%
    filter(CommunityESV %in% ESV_stable$CommunityESV)

# Linear regression
tb_lm <- communities_abundance_fitness_equilibrium %>%
    nest(data = c(-Community, -ESV_ID)) %>%
    mutate(fit = map(data, ~ lm(Fitness ~ Relative_Abundance, data = .x)),
           tidied = map(fit, tidy),
           r.squared = map(fit, function(x) summary(x)$adj.r.squared)) %>%
    unnest(r.squared) %>%
    unnest(tidied)

tb_lm %>%
    filter(term == "Relative_Abundance") %>%
    select(Community, ESV_ID, estimate, std.error, p.value, r.squared)

ESV_sig <- tb_lm %>%
    filter(term == "Relative_Abundance") %>%
    select(Community, ESV_ID, estimate, std.error, p.value, r.squared) %>%
    filter(p.value < 0.05, estimate < 0, estimate > -100) %>%
    distinct(Community, ESV_ID) %>%
    mutate(CommunityESV = paste0(Community, ESV_ID))



# Abundance vs. fitness for ESVs in final communities, T8-12 ----
p <- communities_abundance_fitness_equilibrium %>%
    ggplot() +
    geom_smooth(data = filter(communities_abundance_fitness_equilibrium, CommunityESV %in% ESV_sig$CommunityESV),
                aes(x = Relative_Abundance, y = Fitness), method = "lm", formula = y~x, se = F) +
    geom_point(aes(x = Relative_Abundance, y = Fitness), shape = 21) +
    geom_hline(yintercept = 0, linetype = 2) +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 3)) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 3)) +
    facet_wrap(Community~ESV_ID, scales = "free", ncol = 9) +
    theme_classic() +
    theme(axis.text = element_text(size = 8, angle = 30, hjust = 1),
          axis.title = element_text(size = 15),
          strip.text = element_text(size = 8),
          panel.border = element_rect(color = 1, fill = NA)) +
    labs(x = expression(x[i]), y = expression(log(x[i+1]/x[i])))

ggsave(here::here("plots/FigS12-species_fitness_equilibrium.png"), p, width = 12, height = 15)

# ESVs that are ephemeral
communities_abundance_fitness %>%
    filter(Transfer %in% 8:12) %>%
    drop_na() %>%
    mutate(CommunityESV = paste0(Community, ESV_ID)) %>%
    filter(!(CommunityESV %in% ESV_stable$CommunityESV)) %>%
    nest(data = c(-Community, -ESV_ID)) %>%
    mutate(fit = map(data, ~ lm(Fitness ~ Relative_Abundance, data = .x)),
           tidied = map(fit, tidy),
           r.squared = map(fit, function(x) summary(x)$adj.r.squared)) %>%
    unnest(r.squared) %>%
    unnest(tidied) %>%
    filter(term == "Relative_Abundance") %>%
    select(Community, ESV_ID, estimate, std.error, p.value, r.squared) %>%
    filter(p.value < 0.05, estimate < 0, estimate > -100)



