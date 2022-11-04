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
    arrange(CommunitySize) %>%
    mutate(Community = factor(Community, Community))

# Figure SXX isolate abundance in community ----
family_colors <- c(
    Rest = grey(0.5),
    Enterobacteriaceae = "#397eb8",
    Pseudomonadaceae = "#e21e26",
    Aeromonadaceae = "#4fb148",
    Sphingobacteriaceae = "#984e9e",
    Moraxellaceae = "firebrick",
    Comamonadaceae = "yellow",
    Alcaligenaceae = "darkorchid2"
)

genus_colors <- c(
    Rest = grey(0.5),
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
        bind_rows(x, tibble(Community = unique(x$Community), Family = "Rest", Genus = "Rest", RelativeAbundance = 1-MatchedRelativeAbundance))
    }) %>%
    bind_rows()


p <- isolates_abundance %>%
    mutate(Family = factor(Family, names(family_colors))) %>%
    ggplot() +
    geom_bar(aes(x = Community, y = RelativeAbundance, fill = Family), size = .3, position = "stack", stat = "identity") +
    theme_bw() +
    scale_fill_manual(values = family_colors, breaks = names(family_colors)[names(family_colors) != "Rest"]) +
    #scale_x_discrete(labels = 1:13) +
    scale_y_continuous(breaks = c(0, .5, 1), expand = c(0,0), limits = c(0, 1)) +
    theme(panel.grid = element_blank(),
          axis.title = element_text(size = 12),
          axis.text = element_text(color = 1, size = 12),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          panel.border = element_rect(color = 1, size = 1)) +
    labs(y = "Relative abundance")

## Stats
isolates %>%
    group_by(Community) %>%
    summarize(Total = sum(RelativeAbundance, na.rm = T)) %>%
    summarize(Mean = mean(Total))

ggsave(here::here("plots/FigS-isolate_abundance.png"), p, width = 5, height = 3)
ggsave(here::here("plots/FigS-isolate_abundance.pdf"), p, width = 5, height = 3)



# Figure SXX ESC abundance in community ----
#community_abundance <- read_csv(paste0(folder_data, "raw/community_ESV/Emergent_Simplicity_Equilibrium_Data.csv"), show_col_types = F)
community_abundance <- read_csv(paste0(folder_data, "raw/community_ESV/Estrela_2021_16S_data_table_RelAbund_ALL.txt"), show_col_types = F)
community_abundance %>%
    filter(Inoculum == 2)








