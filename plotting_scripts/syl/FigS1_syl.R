# Figure S1

library(tidyverse)
library(data.table)
library(cowplot)

#
sequences_abundance <- fread("data/temp/sequences_abundance.csv")
communities_name <- c(paste0("C", 1:12, "Rpool"), paste0("C", rep(1:12, each = 8), "R", rep(1:8, 12)))

# Color sets from Sylvie
family_color <- tibble(Color = c("yellow", "deepskyblue3", "blue", "darkorchid2", "firebrick", "orange2", "grey"),
                       Family = c("Aeromonadaceae", "Enterobacteriaceae", "Moraxellaceae", "Pseudomonadaceae","Comamonadaceae","Alcaligenaceae", "Sphingobacteriaceae"))


plot_abundance <- function(df, label_x = "Community", label_y = "RelativeAbundance", fill = "CommunityESVID") {
    ggplot(df) +
        geom_bar(aes_string(x = label_x, y = label_y, fill = fill),
                 position = "stack", stat = "identity", col = 1) +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 90))
}

p <- sequences_abundance %>%
    filter(AlignmentType == "local") %>%
    filter(AllowMismatch == Inf) %>%
    filter(BasePairMismatch <= 4) %>%
    mutate(Community = ordered(Community,  communities_name)) %>%
    plot_abundance(label_x = "Community", label_y = " RelativeAbundance", fill = "Family") +
    labs(x = "", y = "Relative abundance") +
    scale_fill_manual(values = setNames(family_color$Color, family_color$Family)) +
    #scale_fill_manual(values = family_color) +
    scale_x_discrete(expand=c(0,0)) +
    scale_y_continuous(expand=c(0,0), limits = c(0,1), breaks = seq(0,1, .25)) +
    NULL

ggsave("FigS1_syl.png", plot = p, width = 6, height = 4)


# Averaged coverage
sequences_abundance %>%
    filter(AlignmentType == "local") %>%
    filter(AllowMismatch == Inf) %>%
    filter(BasePairMismatch <= 4) %>%
    mutate(Community = ordered(Community,  communities_name)) %>%
    group_by(Community) %>%
    summarize(TotalRelativeAbundance = sum(RelativeAbundance)) %>%
    pull(TotalRelativeAbundance) %>% mean





