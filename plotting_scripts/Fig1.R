library(tidyverse)
library(cowplot)
source(here::here("analysis/00-metadata.R"))

isolates <- read_csv(paste0(folder_data, "output/isolates.csv"), show_col_types = F)
pairs <- read_csv(paste0(folder_data, "output/pairs.csv"), show_col_types = F)
communities <- read_csv(paste0(folder_data, "temp/00c-communities.csv"), show_col_types = F)
community_abundance <- read_csv(paste0(folder_data, "raw/community_ESV/CSS_Emergent_Simplicity_Timeseries.csv"), show_col_types = F)

# Inset1: ESC abundance over time ----
p1 <- community_abundance %>%
    filter(Inoculum == 8, Replicate == 4) %>%
    mutate(Family = ifelse(Relative_Abundance < 0.001 | !(Family %in% names(family_colors)), "Others", Family)) %>%
    mutate(Family = factor(Family, names(family_colors))) %>%
    arrange(Transfer, Family) %>%
    ggplot() +
    geom_col(aes(x = Transfer, y = Relative_Abundance, fill = Family, alpha = ESV), linewidth = .3) +
    annotate("text", x = 0, y = -0.1, label = "inoculum", size = 6, angle = 20, hjust = 0.9, vjust = -0.5) +
    scale_fill_manual(values = family_colors) +
    scale_alpha_discrete(range = c(0.5, 1)) +
    scale_x_continuous(breaks = 0:12, expand = c(0,.3), labels = c("", 1:12)) +
    scale_y_continuous(expand = c(0, 0), breaks = c(0, .5, 1)) +
    coord_cartesian(ylim = c(0, 1), clip = "off") +
    theme_classic() +
    theme(panel.grid = element_blank(),
          axis.title = element_text(size = 17),
          axis.text = element_text(color = 1, size = 17),
          panel.border = element_rect(color = 1, size = 1, fill = NA),
          plot.background = element_blank(),
          panel.background = element_blank()) +
    guides(alpha = "none") +
    labs(x = "Transfer", y = "Relative abundance")

# Inset2: isolate abundance in community ----
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
    coord_cartesian(ylim = c(0, 1), clip = "off") +
    theme(panel.grid = element_blank(),
          axis.title = element_text(size = 17),
          axis.text = element_text(color = 1, size = 17),
          panel.border = element_rect(color = 1, size = 1),
          plot.background = element_blank(),
          panel.background = element_blank()) +
    guides(alpha = "none") +
    labs(x = "Community", y = "Relative abundance")

# Stats
isolates %>%
    group_by(Community) %>%
    summarize(Total = sum(RelativeAbundance, na.rm = T)) %>%
    summarize(Mean = mean(Total))

# Assemble panels
p_legend <- get_legend(p1 + theme(legend.text = element_text(size = 17), legend.title = element_text(size = 17), legend.key.size = unit(1, "cm")))
p <- ggdraw() +
    draw_image(here::here("plots/cartoons/Fig1_cartoon.png")) +
    draw_plot(p1 + guides(fill = "none"), x = 0.54, y = 0.4, width = 0.18, height = 0.26) +
    draw_plot(p2 + guides(fill = "none"), x = 0.54, y = 0.05, width = 0.18, height = 0.26) +
    draw_plot(p_legend, x = 0.76, y = 0.4, width = 0.1, height = 0.1)
ggsave(here::here("plots/Fig1.png"), p, width = 27, height = 15)







