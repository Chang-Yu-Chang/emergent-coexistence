library(tidyverse)
library(cowplot)
source(here::here("processing_scripts/00-metadata.R"))

factorize_communities <- function (x) x %>% mutate(Community = factor(Community, paste0("C", rep(1:12, each = 8), "R", rep(1:8, 12))))
fitness_transient <- read_csv(paste0(folder_data, "temp/15-fitness_transient.csv"), show_col_types = F) %>% factorize_communities
eq_freq_transient <- read_csv(paste0(folder_data, "temp/15-eq_freq_transient.csv"), show_col_types = F) %>% factorize_communities
fitness_transient <- fitness_transient %>% left_join(eq_freq_transient) %>% replace_na(list(Significance = "p>=0.05"))

p <- fitness_transient %>%
    ggplot() +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_point(aes(x = Relative_Abundance, y = Fitness), shape = 21) +
    geom_smooth(aes(x = Relative_Abundance, y = Fitness, color = Significance), method = stats::lm, formula = y ~ x, se = F) +
    scale_color_manual(values = c("p<0.05" = "pink", "p>=0.05" = grey(0.8))) +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 3)) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 3)) +
    facet_wrap(Community~ESV_ID, scales = "free", ncol = 8) +
    theme_classic() +
    theme(
        axis.text = element_text(size = 8, angle = 30, hjust = 1),
        axis.title = element_text(size = 15),
        strip.text = element_text(size = 8),
        panel.border = element_rect(color = 1, fill = NA),
        legend.position = "top"
    ) +
    guides(color = "none") +
    labs(x = expression(x[i]), y = expression("F="~log(x[i]/x[i-1])))

ggsave(here::here("plots/FigS4.png"), p, width = 12, height = 10)

#
nrow(eq_freq_transient) # 46 ESVs
table(eq_freq_transient$Slope < 0) # 41 ESVs with slope <0
table(eq_freq_transient$Significance, useNA = "always") # 10 ESVs have significant correlation
