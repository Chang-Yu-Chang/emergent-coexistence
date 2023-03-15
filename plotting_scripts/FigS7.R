library(tidyverse)
library(cowplot)
library(broom)
library(infer) # For tidyverse statistics
source(here::here("analysis/00-metadata.R"))

pairs <- read_csv(paste0(folder_data, "output/pairs.csv"), show_col_types = F)
communities <- read_csv(paste0(folder_data, "temp/00c-communities.csv"), show_col_types = F)

# Clean up the pairs data ----
pairs <- pairs %>%
    # Remove no-colony pairs
    unite(col = "Pair", Community, Isolate1, Isolate2, sep = "_", remove = F) %>%
    filter(!(Pair %in% pairs_no_colony)) %>%
    # Remove low-accuracy model pairs
    filter(AccuracyMean > 0.9)

# 16S mismatch versus probability of coexistence ----
# Data with mismatch
pairs_mismatch <- pairs %>%
    filter(!is.na(Mismatch)) %>%
    mutate(MismatchGroup = case_when(
        Mismatch < 90 ~ "<90",
        Mismatch >= 90 ~ ">=90"
    )) %>%
    filter(InteractionType != "unknown") %>%
    mutate(InteractionType = factor(InteractionType, c("exclusion", "coexistence"))) %>%
    mutate(InteractionTypeBinary = ifelse(InteractionType == "coexistence", 1, 0))

pairs_mismatch %>%
    group_by(InteractionType) %>%
    count()

pairs_mismatch_group <- pairs_mismatch %>%
    group_by(MismatchGroup, InteractionType) %>%
    summarize(Count = n()) %>%
    group_by(MismatchGroup) %>%
    mutate(Fraction = Count / sum(Count), TotalCount = sum(Count))


# Statistics
# extract_statistics_equation <- function (model) {
#     p_value <- summary(model)$coefficients[2,4]
#     # For glm model
#     r_square <- with(summary(model), 1 - deviance/null.deviance)
#     paste0(
#         "\n", case_when(
#             p_value < 0.001 ~ "p<0.001",
#             p_value < 0.0001 ~ "p<0.0001",
#             TRUE ~ paste0("p=", as.character(round(p_value, 3)))
#         ),
#         "\nR-squared=", round(r_square,4))
# }
# model1 <- glm(InteractionType ~ Mismatch, family = "binomial", data = pairs_mismatch)


# scatterplot, categorical lm
# p1 <- pairs_mismatch %>%
#     ggplot(aes(x = Mismatch, y = InteractionTypeBinary)) +
#     geom_point(shape = 21, size = 3) +
#     geom_smooth(method = "glm", method.args = list(family = "binomial")) +
#     annotate("text", x = 100, y = 0.8, label = extract_statistics_equation(model1), hjust = 0) +
#     scale_y_continuous(labels = c("0" = "exclusion", "1" = "coexistence"), breaks = c(0,1), expand = c(0,.2)) +
#     theme_classic() +
#     labs(x = "# of nucleotide difference in 16S", y = "")

# boxplot by mismatch
p1 <- pairs_mismatch %>%
    ggplot(aes(x = InteractionType, y = Mismatch, fill = InteractionType)) +
    geom_boxplot(linewidth = 1, width = 0.4) +
    geom_point(size = 2, stroke = .5, shape = 21, position = position_jitter(width = 0.2)) +
    scale_fill_manual(values = interaction_color) +
    coord_flip() +
    theme_classic() +
    theme(axis.title.y = element_blank()) +
    guides(color = "none", fill = "none") +
    labs(y = "# of nucleotide difference in 16S", y = "")

t.test(
    pairs_mismatch %>% filter(InteractionType == "coexistence") %>% pull(Mismatch),
    pairs_mismatch %>% filter(InteractionType == "exclusion") %>% pull(Mismatch)
)

p2 <- pairs_mismatch %>%
    mutate(InteractionType = factor(InteractionType, c("coexistence", "exclusion"))) %>%
    group_by(MismatchGroup, InteractionType) %>%
    summarize(Count = n()) %>%
    mutate(Fraction = Count / sum(Count), TotalCount = sum(Count)) %>%
    ggplot() +
    geom_col(aes(x = MismatchGroup, y = Fraction, fill = InteractionType), color = 1, width = 0.7) +
    geom_text(aes(x = MismatchGroup, label = paste0("n=", TotalCount)), y = 0.9) +
    scale_fill_manual(values = interaction_color) +
    scale_y_continuous(expand = c(0,0), breaks = seq(0, 1, 0.2)) +
    theme_classic() +
    theme(legend.position = "right") +
    guides(fill = guide_legend(title = "")) +
    labs(x = "# of nucleotide difference in 16S")

matrix(
    pairs_mismatch_group %>% filter(MismatchGroup == "<90") %>% pull(Count),
    pairs_mismatch_group %>% filter(MismatchGroup == ">=90") %>% pull(Count),
    ncol = 2
) %>%
    chisq.test

p <- plot_grid(p1, p2, nrow = 1, labels = LETTERS[1:2], scale = 0.9, axis = "tb", align = "h") + paint_white_background()
ggsave(here::here("plots/FigS7-mismatch_vs_coexistence.png"), p, width = 8, height = 4)

