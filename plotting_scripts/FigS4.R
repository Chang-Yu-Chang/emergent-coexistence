library(tidyverse)
library(cowplot)
library(broom)
source(here::here("analysis/00-metadata.R"))
source(here::here("plotting_scripts/Fig1.R"))

isolates <- read_csv(paste0(folder_data, "output/isolates.csv"), show_col_types = F)
pairs <- read_csv(paste0(folder_data, "output/pairs.csv"), show_col_types = F)
communities <- read_csv(paste0(folder_data, "temp/00c-communities.csv"), show_col_types = F)

# ESV-Sanger mismatch
sequences_alignment <- read_csv(paste0(folder_data, "temp/14-sequences_alignment.csv"), show_col_types = F)
isolates_abundance <- read_csv(paste0(folder_data, "temp/14-isolates_abundance.csv"), show_col_types = F) %>%
    left_join(communities) %>%
    mutate(Community = factor(Community, communities$Community))
isolate_ESV_1on1_match <- read_csv(paste0(folder_data, "temp/14-isolate_ESV_1on1_match.csv"), show_col_types = F) %>%
    left_join(communities) %>%
    mutate(Community = factor(Community, communities$Community))

p1 <- sequences_alignment %>%
    left_join(communities) %>%
    mutate(Community = factor(Community, communities$Community)) %>%
    filter(AlignmentType == "local") %>%
    drop_na() %>%
    ggplot() +
    geom_tile(aes(x = ExpID, y = CommunityESVID, fill = BasePairMismatch)) +
    # Each isolate has one ESV match
    geom_tile(data = isolate_ESV_1on1_match, aes(x = ExpID, y = CommunityESVID), fill = NA, color = grey(0.1), linewidth = .5) +
    # Best match
    geom_tile(data = isolates_abundance, aes(x = ExpID, y = CommunityESVID), fill = NA, color = "red", linewidth = .5) +
    facet_wrap(.~CommunityLabel, scales = "free", ncol = 4) +
    scale_fill_gradient(low = grey(0.9), high = "darkblue") +
    scale_x_discrete(expand = c(0,0)) +
    scale_y_discrete(expand = c(0,0)) +
    theme_classic() +
    theme(
        panel.border = element_rect(fill = NA, color = 1, linewidth = 0.5),
        axis.text.x = element_text(angle = 90, vjust = .5, hjust = 1, size = 7),
        axis.text.y = element_text(size = 6),
        strip.background = element_blank(),
        legend.position = c(.6, .05),
        plot.title = element_text(hjust = 0.5)
    ) +
    guides(fill = guide_legend(nrow = 1, title = "bp difference")) +
    labs(x = "", y = "ESV") +
    ggtitle("community")


# Family level figures
family_names <- unique(isolates_abundance$Family)
family_names <- family_names[!is.na(family_names)]
family_colors <- c(
    "Enterobacteriaceae" = rev(RColorBrewer::brewer.pal(11, "RdYlBu"))[2],
    "Pseudomonadaceae" = rev(RColorBrewer::brewer.pal(11, "PRGn"))[2],
    RColorBrewer::brewer.pal(length(family_names)-2, "Set2") %>% setNames(family_names[!(family_names %in% c("Enterobacteriaceae", "Pseudomonadaceae"))]),
    "Other" = "#999999"
)

p3 <- isolates_abundance %>%
    drop_na() %>%
    bind_rows(mutate(isolates_abundance_other, Family = "Other") %>% left_join(communities)) %>%
    mutate(Family = factor(Family, names(family_colors))) %>%
    ggplot() +
    geom_col(aes(x = CommunityLabel, y = RelativeAbundance, fill = Family), position = position_stack(reverse = T)) +
    theme_bw() +
    scale_fill_manual(values = family_colors) +
    scale_x_continuous(breaks = 1:13, expand = c(0,.3)) +
    scale_y_continuous(breaks = c(0, .5, 1), expand = c(0,0), limits = c(0, 1)) +
    coord_cartesian(ylim = c(0, 1), clip = "off") +
    theme(panel.grid = element_blank(),
          axis.title = element_text(size = 17),
          axis.text = element_text(color = 1, size = 17),
          panel.border = element_rect(color = 1, linewidth = 0.5),
          plot.background = element_blank(),
          panel.background = element_blank()) +
    guides(alpha = "none") +
    labs(x = "community", y = "relative abundance")

p5 <- isolates %>%
    group_by(Community) %>%
    summarize(Total = sum(RelativeAbundance, na.rm = T)) %>%
    mutate(Xaxis = 1) %>%
    ggplot(aes(Xaxis, y = Total)) +
    geom_boxplot(width = 0.5, outlier.color = NA) +
    geom_jitter(width = 0.2, shape = 21, stroke = 1, size = 1.5) +
    scale_x_continuous(breaks = 1, limits = c(.5,1.5)) +
    scale_y_continuous(limits = c(0,1), breaks = seq(0,1,.5)) +
    theme_classic() +
    theme(
        panel.grid.major.x = element_line(color = grey(0.9)),
        panel.grid.major.y = element_line(color = grey(0.9)),
        panel.grid.minor.y = element_line(color = grey(0.9)),
        axis.text.x = element_blank(),
        panel.border = element_rect(color = 1, fill = NA),
        axis.text.y = element_text(color = 1, size = 10),
        axis.title = element_text(size = 10)
    ) +
    labs(x = "", y = "total abundance")

p <- plot_grid(
    p1,
    plot_grid(p2 + guides(fill = guide_legend(ncol = 3, title = "ESV")) + theme(
        legend.key.size = unit(5, "mm"),
        legend.title = element_text(size = 10),
        axis.text = element_text(color = 1, size = 10),
        axis.title = element_text(size = 10)),
        NULL,
        nrow = 1, rel_widths = c(40,1)
        ),
    plot_grid(p3 + guides(fill = guide_legend(ncol = 2)) + theme(
        legend.key.size = unit(5, "mm"),
        legend.title = element_text(size = 10),
        axis.text = element_text(color = 1, size = 10),
        axis.title = element_text(size = 10)),
        p5,
        nrow = 1, rel_widths = c(3.9,1), axis = "tb", align = "h", labels = c("", "D"), label_x = -.2, label_y = 1.1
        ),
    ncol = 1, axis = "lr", scale = c(.95, .9, .9), labels = LETTERS[1:3], rel_heights = c(4,1,1)
) + paint_white_background()

ggsave(here::here("plots/FigS4-isolate_composition.png"), p, width = 8, height = 10)

# ESV richness per community
communities_abundance <- read_csv(paste0(folder_data, "raw/community_ESV/Emergent_Comunity_Data.csv"), show_col_types = F) %>%
    filter(Carbon_Source == "Glucose" | Carbon_Source == "Original") %>%
    mutate(Community = factor(paste0("C", Inoculum, "R", Replicate), paste0("C", rep(1:12, each = 8), "R", rep(1:8, 12))))%>%
    arrange(Community, Family, Transfer, ESV)

communities_abundance %>%
    filter(Community %in% communities$Community) %>%
    filter(Transfer == 12) %>%
    #filter(Relative_Abundance > 0.01) %>%
    group_by(Community) %>%
    count() %>%
    left_join(communities) %>%
    arrange(CommunityLabel)

# Matched ESV richness per community
isolates_abundance %>%
    drop_na() %>%
    group_by(Community) %>%
    count()

# Total abundance
isolates %>%
    group_by(Community) %>%
    summarize(Total = sum(RelativeAbundance, na.rm = T)) %>%
    summarize(Mean = mean(Total))



















