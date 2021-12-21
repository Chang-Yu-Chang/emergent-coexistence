# Isolate
library(tidyverse)
library(cowplot)
library(ggsci)

# Data
isolates <- read_csv(here::here("data/output/isolates.csv"))
pairs_meta <- read_csv(here::here("data/output/pairs_meta.csv"))
isolates_growth_syl <- read_csv("~/Dropbox/lab/invasion-network/data/raw/growth_rate/Estrela_2021_isolates_grmax.csv")

# boxplot leakiness
p <- isolates %>%
    filter(!is.na(Fermenter)) %>%
    mutate(Fermenter = ifelse(Fermenter, "fermenter", "respirator")) %>%
    ggplot() +
    geom_boxplot(aes(x = Fermenter, y = leakiness_16hr, color = Fermenter)) +
    geom_jitter(aes(x = Fermenter, y = leakiness_16hr, color = Fermenter), shape = 1, width = 0.5) +
    scale_color_npg() +
    theme_classic() +
    guides(color = "none") +
    labs(x = "", y = "leakiness")
ggsave(here::here("plots/Fig_isolates-leakiness.png"), p, width = 3, height = 4)


isolates %>%
    filter(!is.na(Fermenter)) %>%
    mutate(Fermenter = ifelse(Fermenter, "fermenter", "respirator")) %>%
    group_by(Fermenter) %>%
    summarize(LeakinessMean = mean(leakiness_16hr, na.rm = T), LeakinessSd = sd(leakiness_16hr, na.rm = T))

# r_glu vs. leakiness
p <- isolates %>%
    filter(!is.na(Fermenter)) %>%
    mutate(Fermenter = ifelse(Fermenter, "fermenter", "respirator")) %>%
    ggplot() +
    geom_point(aes(x = r_glucose, y = leakiness_16hr, color = Fermenter), shape = 1, size = 2) +
    scale_y_continuous(limits = c(0, 1)) +
    scale_color_npg() +
    theme_classic() +
    theme(legend.title = element_blank(), legend.position = "top") +
    labs(x = expression(r[glu]), y = "leakiness")
ggsave(here::here("plots/Fig_isolates-r_glu_leakiness.png"), p, width = 4, height = 4)

# r_glu vs. acetate
## with Jean's r_glu
p1 <- isolates %>%
    filter(!is.na(Fermenter)) %>%
    filter(!is.na(X_acetate_16hr)) %>%
    mutate(Fermenter = ifelse(Fermenter, "fermenter", "respirator")) %>%
    ggplot(aes(x = r_glucose, y = X_acetate_16hr, color = Fermenter)) +
    geom_point(shape = 1, size = 2) +
    geom_smooth(method = "lm") +
    scale_color_npg() +
    theme_classic() +
    theme(legend.title = element_blank(), legend.position = "top") +
    labs(x = expression(r[glu]), y = "acetate (mM) at 16hr")
p1
#ggsave(here::here("plots/Fig_isolates-r_glu_acetate.png"), p1, width = 4, height = 4)

## With Sylvie's data
p2 <- isolates %>%
    left_join(isolates_growth_syl %>% filter(cs == "glucose") %>% select(ID = SangerID, gr_max)) %>%
    filter(!is.na(Fermenter)) %>%
    filter(!is.na(X_acetate_16hr)) %>%
    mutate(Fermenter = ifelse(Fermenter, "fermenter", "respirator")) %>%
    ggplot(aes(x = gr_max, y = X_acetate_16hr, color = Fermenter)) +
    geom_point(shape = 1, size = 2) +
    geom_smooth(method = "lm") +
    scale_color_npg() +
    theme_classic() +
    theme(legend.title = element_blank(), legend.position = "top") +
    labs(x = expression(r[glu]), y = "acetate (mM) at 16hr")
p2
#ggsave(here::here("plots/Fig_isolates-r_glu_acetate.png"), p2, width = 4, height = 4)

p <- plot_grid(p1, p2, nrow = 1, axis = "lrtb", align = "hv", labels = c("Jean's r_mid", "Sylvie's r_max"))
ggsave(here::here("plots/Fig_isolates-r_glu_acetate.png"), p, width = 8, height = 4)





