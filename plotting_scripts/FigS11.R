library(tidyverse)
library(cowplot)
library(broom)
source(here::here("analysis/00-metadata.R"))

isolates <- read_csv(paste0(folder_data, "output/isolates.csv"), show_col_types = F)
pairs <- read_csv(paste0(folder_data, "output/pairs.csv"), show_col_types = F)
pairs_freq <- read_csv(paste0(folder_data, "temp/93-pairs_freq.csv"), show_col_types = F)

# Figure S11 ESV abundance vs. CFU frequencies in pairwise coexistence
pairs_freq_ESV <- pairs_freq %>%
    left_join(pairs, by = join_by(Community, Isolate1, Isolate2)) %>%
    filter(InteractionType == "coexistence") %>%
    filter(Time == "T8") %>%
    group_by(PairID, ExpID1, ExpID2) %>%
    # Mean CFU frequency
    summarize(meanIsolate1CFUFreqMean = mean(Isolate1CFUFreqMean),
              sdIsolate1CFUFreqMean = sd(Isolate1CFUFreqMean)) %>%
    left_join(rename_with(select(isolates, ExpID, RelativeAbundance), ~paste0(.x, "1")), by = "ExpID1") %>%
    left_join(rename_with(select(isolates, ExpID, RelativeAbundance), ~paste0(.x, "2")), by = "ExpID2") %>%
    mutate(RelativeRelativeAbundance1 = RelativeAbundance1 / (RelativeAbundance1 + RelativeAbundance2)) %>%
    mutate(RelativeRelativeAbundance2 = RelativeAbundance2 / (RelativeAbundance1 + RelativeAbundance2)) %>%
    filter(!is.na(RelativeAbundance1), !is.na(RelativeAbundance2))

p <- pairs_freq_ESV %>%
    ggplot() +
    geom_point(aes(x = RelativeRelativeAbundance1, y = meanIsolate1CFUFreqMean),
               shape = 21, size = 2, stroke = .5) +
    geom_segment(aes(x = RelativeRelativeAbundance1, xend = RelativeRelativeAbundance1,
                     y = meanIsolate1CFUFreqMean + sdIsolate1CFUFreqMean,
                     yend = meanIsolate1CFUFreqMean - sdIsolate1CFUFreqMean),
                 color = 1) +
#    geom_smooth(aes(x = RelativeRelativeAbundance1, y = meanIsolate1CFUFreqMean), formula = y~x, method = "lm") +
    scale_x_continuous(limits = c(0,1)) +
    scale_y_continuous(limits = c(0,1)) +
    theme_classic() +
    labs(x = "Relative ESV abundance", y = "CFU frequency")

#p <- plot_grid(p1, p2, nrow = 1, labels = LETTERS[1:2], scale = 0.9, align = "hv") + paint_white_background()
ggsave(here::here("plots/FigS11-abundance_vs_pair_frequency.png"), p, width = 3, height = 3)


lm(meanIsolate1CFUFreqMean ~ RelativeRelativeAbundance1, data = pairs_freq_ESV) %>%
    summary()
    tidy()
