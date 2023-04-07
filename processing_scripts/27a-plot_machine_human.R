library(tidyverse)
library(cowplot)
source(here::here("analysis/00-metadata.R"))

pairs_freq_machine_human <- read_csv(paste0(folder_data, "temp/27-pairs_freq_machine_human.csv"), show_col_types = F)

# 1. Total count ----
p1 <- pairs_freq_machine_human %>%
    ggplot() +
    geom_abline(slope = 1, intercept = 0, color = "red", linetype = 2) +
    geom_point(aes(x = TotalCount_human, y = TotalCount_machine), shape = 21, size = 2) +
    theme_classic()

p2 <- pairs_freq_machine_human %>%
    ggplot() +
    geom_abline(slope = 1, intercept = 0, color = "red", linetype = 2) +
    geom_point(aes(x = TotalCount_human, y = TotalCount_machine), shape = 21, size = 2) +
    scale_x_log10() +
    scale_y_log10() +
    theme_classic() +
    labs(x = "log(TotalCount_human)", y = "log(TotalCount_machine)")
p <- plot_grid(p1, p2, nrow = 1, axis = "tblr", align = "h", scale = .9) +
    theme(plot.background = element_rect(fill = "white", color = NA))

ggsave(paste0(folder_data, "temp/27a-01-comparison-total_count.png"), p, width = 8, height = 4)


# 2. Frequency ----
p <- pairs_freq_machine_human %>%
    ggplot() +
    geom_abline(slope = 1, intercept = 0, color = "red", linetype = 2) +
    geom_hline(yintercept = c(0,1), color = gray(.8), linetype = 2) +
    geom_vline(xintercept = c(0,1), color = gray(.8), linetype = 2) +
    geom_point(aes(x = Isolate1CFUFreq_human, y = Isolate1CFUFreq_machine),
               shape = 21, size = 2, stroke = .4) +
    scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 21)) +
    theme_classic() +
    ggtitle("")

ggsave(paste0(folder_data, "temp/27a-02-comparison-coculture_frequency.png"), p, width = 4, height = 4)


lm(Isolate1CFUFreq_machine ~ Isolate1CFUFreq_human, data = pairs_freq_machine_human) %>%
    summary()

# 3. freq facets by duplicate pairs ----
p <- pairs_freq_machine_human %>%
    ggplot() +
    geom_abline(slope = 1, intercept = 0, color = "red", linetype = 2) +
    geom_hline(yintercept = c(0,1), color = gray(.8), linetype = 2) +
    geom_vline(xintercept = c(0,1), color = gray(.8), linetype = 2) +
    #geom_smooth(aes(x = Isolate1CFUFreq_human, y = Isolate1CFUFreq_machine), method = "lm") +
    geom_point(aes(x = Isolate1CFUFreq_human, y = Isolate1CFUFreq_machine),
               shape = 21, size = 2, stroke = .4) +
    # geom_segment(aes(x = Isolate1CFUFreq_human, xend = Isolate1CFUFreq_human,
    #                  y = Isolate1CFUFreq_machine + 1* Isolate1CFUFreqSd_machine,
    #                  yend = Isolate1CFUFreq_machine - 1* Isolate1CFUFreqSd_machine),
                 # size = .1) +
    scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 21)) +
    facet_grid(.~PairType) +
    theme_classic() +
    ggtitle("")
ggsave(paste0(folder_data, "temp/27a-03-comparison-coculture_frequency_facet.png"), p, width = 10, height = 4)

