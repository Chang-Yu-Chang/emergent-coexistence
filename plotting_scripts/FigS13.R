library(tidyverse)
library(cowplot)
library(broom)
library(infer) # For tidyverse statistics
source(here::here("analysis/00-metadata.R"))

communities <- read_csv(paste0(folder_data, "output/communities_remained.csv"), show_col_types = F)
isolates <- read_csv(paste0(folder_data, "output/isolates_remained.csv"), show_col_types = F)
pairs <- read_csv(paste0(folder_data, "output/pairs_remained.csv"), show_col_types = F)

# 16S mismatch versus probability of coexistence ----
# Data with mismatch
pairs_mismatch <- pairs %>%
    filter(outcome != "5-inconclusive") %>%
    #filter(!is.na(Mismatch)) %>%
    mutate(MismatchGroup = case_when(
        Mismatch < 90 ~ "<90",
        Mismatch >= 90 ~ ">=90"
    )) %>%
    mutate(outcome = case_when(
        outcome %in% c("1-exclusion", "2-exclusion") ~ "exclusion",
        outcome %in% c("3-coexistence", "4-coexistence") ~ "coexistence"
    )) %>%
    mutate(outcomeBinary = ifelse(outcome == "coexistence", 1, 0))

pairs_mismatch %>%
    group_by(outcome) %>%
    count()

pairs_mismatch_group <- pairs_mismatch %>%
    group_by(MismatchGroup, outcome) %>%
    summarize(Count = n()) %>%
    group_by(MismatchGroup) %>%
    mutate(Fraction = Count / sum(Count), TotalCount = sum(Count))


# boxplot by mismatch
p1 <- pairs_mismatch %>%
    ggplot(aes(x = outcome, y = Mismatch, fill = outcome)) +
    geom_boxplot(linewidth = 1, width = 0.4) +
    geom_point(size = 2, stroke = .5, shape = 21, position = position_jitter(width = 0.2)) +
    scale_fill_manual(values = interaction_color) +
    coord_flip() +
    theme_classic() +
    theme(axis.title.y = element_blank()) +
    guides(color = "none", fill = "none") +
    labs(y = "# of nucleotide difference in 16S", y = "")

t.test(
    pairs_mismatch %>% filter(outcome == "coexistence") %>% pull(Mismatch),
    pairs_mismatch %>% filter(outcome == "exclusion") %>% pull(Mismatch)
)

p2 <- pairs_mismatch %>%
    mutate(outcome = factor(outcome, c("coexistence", "exclusion"))) %>%
    group_by(MismatchGroup, outcome) %>%
    summarize(Count = n()) %>%
    mutate(Fraction = Count / sum(Count), TotalCount = sum(Count)) %>%
    ggplot() +
    geom_col(aes(x = outcome, y = Fraction, fill = outcome), color = 1, width = 0.7, position = position_dodge()) +
    geom_text(aes(x = outcome, y = Fraction, label = Count), vjust = 2, position = position_dodge(width = 0.9)) +
    scale_fill_manual(values = interaction_color) +
    facet_grid(.~MismatchGroup) +
    #scale_y_continuous(expand = c(0,0), breaks = seq(0, 1, 0.2), limits = c(0,0.8)) +
    theme_classic() +
    theme(panel.border = element_rect(color = 1, fill = NA),
          axis.text.x = element_text(angle = 20, hjust = 1)) +
    guides(fill = "none") +
    labs(x = "")
matrix(
    pairs_mismatch_group %>% filter(MismatchGroup == "<90") %>% pull(Count),
    pairs_mismatch_group %>% filter(MismatchGroup == ">=90") %>% pull(Count),
    ncol = 2
) %>%
    chisq.test

p <- plot_grid(p1, p2, nrow = 1, labels = LETTERS[1:2], scale = 0.9, axis = "tb", align = "h", rel_widths = c(1,1)) + paint_white_background()
ggsave(here::here("plots/FigS13-mismatch_vs_coexistence.png"), p, width = 8, height = 4)

