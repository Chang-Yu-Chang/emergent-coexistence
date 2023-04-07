library(tidyverse)
library(cowplot)
library(broom)
source(here::here("processing_scripts/00-metadata.R"))

communities_abundance_T0 <- read_csv(paste0(folder_data, "temp/13-communities_abundance_T0.csv"), show_col_types = F)
rarefaction <- read_csv(paste0(folder_data, "temp/13-rarefaction.csv"), show_col_types = F)

# Panel A: T0 soil sample composition ----
communities_abundance_T0_order <- communities_abundance_T0 %>%
    group_by(Inoculum, Order) %>%
    summarize(Relative_Abundance = sum(Relative_Abundance)) %>%
    filter(Relative_Abundance > 0.01) %>%
    ungroup() %>%
    distinct(Order) %>%
    arrange(Order)

order_colors <- c(RColorBrewer::brewer.pal(9, "Set1"), RColorBrewer::brewer.pal(8, "Set3"), RColorBrewer::brewer.pal(8, "Dark2"))[1:(length(communities_abundance_T0_order$Order))]
order_colors <- c(setNames(order_colors, communities_abundance_T0_order$Order), "Others" = "snow")

p1 <- communities_abundance_T0 %>%
    mutate(Order = ifelse(Order %in% communities_abundance_T0_order$Order, Order, "Others")) %>%
    mutate(Order = factor(Order, c(communities_abundance_T0_order$Order, "Others"))) %>%
    mutate(Inoculum = factor(Inoculum, c(1:10, 12))) %>%
    ggplot() +
    geom_col(aes(x = Inoculum, y = Relative_Abundance, color = Order, fill = Order), width = .8, linewidth = .05, position = position_stack(reverse = T)) +
    #scale_x_continuous(breaks = c(1:10, 12), expand = c(0,0.1)) +
    scale_y_continuous(breaks = seq(0, 1, 0.2), limits = c(0,1), expand = c(0,0)) +
    scale_fill_manual(values = order_colors) +
    scale_color_manual(values = order_colors) +
    theme_classic() +
    theme(
        strip.background = element_rect(fill = NA, color = NA),
        panel.border = element_rect(fill = NA, color = 1, linewidth = 1),
        axis.title = element_text(size = 10),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        legend.key.width = unit(5, "mm"),
        legend.key.height = unit(2, "mm"),
    ) +
    guides(fill = guide_legend(ncol = 1, override.aes = list(color = 1, linewidth = .5))) +
    labs(x = "inoculum", y = "relative abundance")

# Panel B: Rarefaction ----
p2 <- rarefaction %>%
    mutate(Inoculum = factor(Inoculum, 1:12)) %>%
    ggplot() +
    scale_y_log10(
        breaks = scales::trans_breaks("log10", function(x) 10^x)(c(1,1e3)),
        labels = scales::trans_format("log10", scales::math_format(10^.x))
    ) +
    annotation_logticks(sides = "lr") +
    geom_line(aes(x = SamplingSize, y = MeanRichness, color = Inoculum), linewidth = .5) +
    scale_x_continuous(expand = c(0,0), limits = c(0, 10000), breaks = seq(0, 10000, 2000)) +
    scale_color_manual(values = c(RColorBrewer::brewer.pal(10, "Paired"), "salmon4")) +
    theme_classic() +
    theme(
        panel.border = element_rect(color = 1, fill = NA),
        axis.title = element_text(size = 10),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.key.size = unit(3, "mm"),
        legend.background = element_blank(),
        legend.position = c(0.7, 0.3)
    ) +
    guides(color = guide_legend(ncol = 2)) +
    labs(x = "sampling size [number of reads]", y = "richness")

p <- plot_grid(p1, p2, nrow = 1, labels = LETTERS[1:2], scale = .9, rel_widths = c(1.3, 1)) + paint_white_background()
ggsave(here::here("plots/FigS1-inoculum_diversity.png"), p, width = 8, height = 3)

# T0 samples
aux = read_csv(paste0(folder_data, 'raw/community_ESV/metadata.csv'), show_col_types = F)
Taxonomy_Data = read_csv(paste0(folder_data, 'raw/community_ESV/taxonomy.csv'), show_col_types = F)
Abundance_Data = read_csv(paste0(folder_data, 'raw/community_ESV/otu_table.csv'), show_col_types = F) #Data is actually at an ESV level

aux = aux %>% filter(Transfer == 0)
Abundance_Data = Abundance_Data[,aux$ID]
Abundance_Data = as.matrix(Abundance_Data[,1:ncol(Abundance_Data)])
range(colSums(Abundance_Data)) # Number of reads
range(apply(Abundance_Data, 2, function(x) x>0) %>% colSums()) # number of unique ESVs



