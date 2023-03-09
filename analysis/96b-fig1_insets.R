#' Script for figures

library(tidyverse)
library(cowplot)
library(tidygraph)
library(ggraph)
library(gridExtra)
source(here::here("analysis/00-metadata.R"))

isolates <- read_csv(paste0(folder_data, "output/isolates.csv"), show_col_types = F)
pairs <- read_csv(paste0(folder_data, "output/pairs.csv"), show_col_types = F)
communities <- read_csv(paste0(folder_data, "temp/00c-communities.csv"), show_col_types = F)
pairs_freq <- read_csv(paste0(folder_data, "temp/93-pairs_freq.csv"), show_col_types = F)
load(paste0(folder_data, "temp/95-communities_network.Rdata"))
communities_hierarchy <- read_csv(paste0(folder_data, "temp/95-communities_hierarchy.csv"), show_col_types = F)

communities <- communities %>%
    mutate(Community = factor(Community, Community))

#
family_colors <- c(
    Others = grey(0.5),
    Enterobacteriaceae = "#397eb8",
    Pseudomonadaceae = "#e21e26",
    Aeromonadaceae = "#4fb148",
    Sphingobacteriaceae = "#984e9e",
    Moraxellaceae = "firebrick",
    Comamonadaceae = "yellow",
    Alcaligenaceae = "darkorchid2"
)

genus_colors <- c(
    Others = grey(0.5),
    Enterobactor1 = "#225ea8",
    Klebsiella1 = "#3eb6c5",
    Raoultella1 = "#a3d6b2",
    Citrobacter1 = "#fcf8cf",
    Pseudomonas1 = "#7d1517",
    Pseudomonas2 = "#b31e24",
    Pseudomonas3 = "#d63226",
    Pseudomonas4 = "#e44b34",
    Pseudomonas5 = "#ec6448",
    Pseudomonas6 = "#f68d5c",
    Aeromonas1 = "#8fd1c6"
)


# Figure 1 inset ESC abundance over time ----
community_abundance <- read_csv(paste0(folder_data, "raw/community_ESV/CSS_Emergent_Simplicity_Timeseries.csv"), show_col_types = F)
p1 <- community_abundance %>%
    filter(Inoculum == 8, Replicate == 4) %>%
    mutate(Family = ifelse(Relative_Abundance < 0.001 | !(Family %in% names(family_colors)), "Others", Family)) %>%
    mutate(Family = factor(Family, names(family_colors))) %>%
    arrange(Transfer, Family) %>%
    ggplot() +
    geom_col(aes(x = Transfer, y = Relative_Abundance, fill = Family, alpha = ESV), size = .3) +
    scale_fill_manual(values = family_colors) +
    scale_alpha_discrete(range = c(0.5, 1)) +
    scale_x_continuous(breaks = 0:12, expand = c(0,.3)) +
    scale_y_continuous(expand = c(0, 0), breaks = c(0, .5, 1)) +
    theme_classic() +
    theme(panel.grid = element_blank(),
          axis.title = element_text(size = 12),
          axis.text = element_text(color = 1, size = 12),
          panel.border = element_rect(color = 1, size = 1, fill = NA)) +
    guides(alpha = "none") +
    labs(x = "Generation", y = "Relative abundance")


community_abundance %>%
    distinct(Inoculum, Replicate) %>%
    mutate(Community = paste0("C", Inoculum, "R", Replicate)) %>%
    filter(Community %in% )
    distinct(Inoculum, Replicate, Transfer) %>%
    arrange(Inoculum, Replicate, Transfer) %>%
    view
    group_by(Inoculum, Replicate) %>%
    summarize(Count = n()) %>%
    view

# Figure 1 inset isolate abundance in community ----
isolates_abundance <- isolates %>%
    mutate(Community = factor(Community, communities$Community)) %>%
    filter(!is.na(RelativeAbundance)) %>%
    arrange(Community) %>%
    group_by(Community) %>%
    select(Community, Family, Genus, RelativeAbundance) %>%
    ungroup() %>%
    arrange(Community, Family, Genus) %>%
    # Unique genus name
    group_by(Community, Family, Genus) %>%
    mutate(Genus = paste0(Genus, 1:n()))

isolates_abundance <- isolates_abundance %>%
    split.data.frame(.$Community) %>%
    map(function (x) {
        MatchedRelativeAbundance = sum(x$RelativeAbundance)
        bind_rows(x, tibble(Community = unique(x$Community), Family = "Others", Genus = "Others", RelativeAbundance = 1-MatchedRelativeAbundance))
    }) %>%
    bind_rows()


p2 <- isolates_abundance %>%
    mutate(Family = factor(Family, names(family_colors))) %>%
    left_join(communities, by = "Community") %>%
    ggplot() +
    geom_bar(aes(x = CommunityLabel, y = RelativeAbundance, fill = Family, alpha = Genus), size = .3, position = "stack", stat = "identity") +
    theme_bw() +
    scale_fill_manual(values = family_colors, breaks = names(family_colors)[names(family_colors) != "Others"]) +
    scale_alpha_discrete(range = c(1, 0.8)) +
    scale_x_continuous(breaks = 1:13, expand = c(0,.3)) +
    scale_y_continuous(breaks = c(0, .5, 1), expand = c(0,0), limits = c(0, 1)) +
    theme(panel.grid = element_blank(),
          axis.title = element_text(size = 12),
          axis.text = element_text(color = 1, size = 12),
          panel.border = element_rect(color = 1, size = 1)) +
    guides(alpha = "none") +
    labs(x = "Community", y = "Relative abundance")

## Stats
isolates %>%
    group_by(Community) %>%
    summarize(Total = sum(RelativeAbundance, na.rm = T)) %>%
    summarize(Mean = mean(Total))


p <- plot_grid(p1, p2 + guides(fill = "none"), nrow = 1, rel_widths = c(1,1.1), align = "vh", axis = "lrtb")
ggsave(here::here("plots/cartoons/FigS-isolate_abundance.pdf"), p, width = 10, height = 2.5)






