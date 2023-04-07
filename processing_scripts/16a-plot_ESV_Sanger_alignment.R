#' This script plots the heatmap for ESV Sange alignnebt

library(tidyverse)
library(cowplot)
source(here::here("analysis/00-metadata.R"))

communities <- read_csv(paste0(folder_data, "temp/00c-communities.csv"), show_col_types = F)
sequences_alignment <- read_csv(paste0(folder_data, "temp/16-sequences_alignment.csv"), show_col_types = F)
isolates_abundance <- read_csv(paste0(folder_data, "temp/16-isolates_abundance.csv"), show_col_types = F)


# 1. Heatmap, ESV and sanger alignment ----
isolates_algn_bp_mismatch <- isolates_abundance %>%
    select(Community, ExpID, Genus) %>%
    left_join(algn_Sanger_ESV1) %>%
    mutate(Community = factor(Community, communities$Community)) %>%
    mutate(ExpIDGenus = paste0(ExpID, "-", Genus))

p <- sequences_alignment %>%
    mutate(Community = factor(Community, communities$Community)) %>%
    mutate(ExpIDGenus = paste0(ExpID, "-", Genus)) %>%
    ggplot() +
    geom_tile(aes(x = ExpIDGenus, y = CommunityESVID, fill = BasePairMismatch)) +
    # Each Sanger has one ESV match
    geom_tile(data = isolates_algn_bp_mismatch, aes(x = ExpIDGenus, y = CommunityESVID), fill = NA, color = "gold", linewidth = 1) +
    facet_wrap(.~Community, scales = "free", ncol = 4) +
    scale_fill_gradient(low = RColorBrewer::brewer.pal(n = 9, "Blues")[9], high = grey(0.9)) +
    scale_x_discrete(expand = c(0,0)) +
    scale_y_discrete(expand = c(0,0)) +
    theme_classic() +
    theme(
        panel.border = element_rect(fill = NA, color = 1, linewidth = 0.5),
        axis.text.x = element_text(angle = 30, vjust = 1.1, hjust = 1, size = 7),
        axis.text.y = element_text(size = 7),
        legend.position = c(.6, .05)
    ) +
    guides(fill = guide_legend(nrow = 1)) +
    labs()

ggsave(paste0(folder_data, "temp/16a-01-ESV_Sanger_bp_mismatch.png"), p, width = 12, height = 12)


# 2. Heatmap, ESV and sanger alignment, binned by mismatch ----
p <- sequences_alignment %>%
    mutate(Community = factor(Community, communities$Community)) %>%
    mutate(ExpIDGenus = paste0(ExpID, "-", Genus)) %>%
    mutate(BasePairMismatchCategory = case_when(
        BasePairMismatch == 0 ~ "0",
        BasePairMismatch == 1 ~ "1",
        BasePairMismatch == 2 ~ "2",
        BasePairMismatch == 3 ~ "3",
        BasePairMismatch == 4 ~ "4",
        # BasePairMismatch >= 4 & BasePairMismatch <= 10 ~ "4-10",
        # BasePairMismatch >= 11 & BasePairMismatch <= 29 ~ "11-29",
        # BasePairMismatch >= 30 ~ ">=30"
        BasePairMismatch >= 5 ~ ">=5"
        #) %>% ordered(c(0:3, "4-10", "11-29", ">=30"))
    ) %>% ordered(c(0:4, ">=5"))
    ) %>%
    ggplot() +
    geom_tile(aes(x = ExpIDGenus, y = CommunityESVID, fill = BasePairMismatchCategory)) +
    # Each Sanger has one ESV match
    geom_tile(data = isolates_algn_bp_mismatch, aes(x = ExpIDGenus, y = CommunityESVID), fill = NA, color = "gold", linewidth = 1) +
    facet_wrap(.~Community, scales = "free", ncol = 4) +
    scale_fill_manual(values = rev(RColorBrewer::brewer.pal(n = 6, "Blues"))) +
    scale_x_discrete(expand = c(0,0)) +
    scale_y_discrete(expand = c(0,0)) +
    theme_classic() +
    theme(
        panel.border = element_rect(fill = NA, color = 1, linewidth = 0.5),
        axis.text.x = element_text(angle = 30, vjust = 1.1, hjust = 1, size = 7),
        axis.text.y = element_text(size = 7),
        legend.position = c(.6, .05)
    ) +
    guides(fill = guide_legend(nrow = 1)) +
    labs()

ggsave(paste0(folder_data, "temp/16a-02-ESV_Sanger_bp_mismatch_binned.png"), p, width = 12, height = 12)

# 3. Concensus length ----
p <- sequences_alignment %>%
    mutate(Community = factor(Community, communities$Community)) %>%
    mutate(ExpIDGenus = paste0(ExpID, "-", Genus)) %>%
    ggplot() +
    geom_tile(aes(x = ExpIDGenus, y = CommunityESVID, fill = ConsensusLength)) +
    # Each Sanger has one ESV match
    geom_tile(data = isolates_algn_bp_mismatch, aes(x = ExpIDGenus, y = CommunityESVID), fill = NA, color = "gold", linewidth = 1) +
    facet_wrap(.~Community, scales = "free", ncol = 4) +
    scale_fill_gradient(low = RColorBrewer::brewer.pal(n = 9, "Blues")[9], high = grey(0.9)) +
    scale_x_discrete(expand = c(0,0)) +
    scale_y_discrete(expand = c(0,0)) +
    theme_classic() +
    theme(
        panel.border = element_rect(fill = NA, color = 1, linewidth = 0.5),
        axis.text.x = element_text(angle = 30, vjust = 1.1, hjust = 1, size = 7),
        axis.text.y = element_text(size = 7),
        legend.position = c(.6, .05)
    ) +
    guides(fill = guide_legend(nrow = 1)) +
    labs()

ggsave(paste0(folder_data, "temp/16a-03-ESV_Sanger_bp_concensus_length.png"), p, width = 12, height = 12)


# 4. Alignment score ----
algn_Sanger_ESV1 <- sequences_alignment %>%
    # For each Sanger, find the ESV which it has least base pair mismatch with
    group_by(ExpID) %>%
    filter(AlignmentScore == max(AlignmentScore)) %>%
    # When a Sanger has two ESV that both has zero, pick the first one
    slice(1) %>%
    ungroup()

isolates_algn_score <- isolates_abundance %>%
    select(Community, ExpID, Genus) %>%
    left_join(algn_Sanger_ESV1) %>%
    mutate(Community = factor(Community, communities$Community)) %>%
    mutate(ExpIDGenus = paste0(ExpID, "-", Genus))

p <- sequences_alignment %>%
    mutate(Community = factor(Community, communities$Community)) %>%
    mutate(ExpIDGenus = paste0(ExpID, "-", Genus)) %>%
    ggplot() +
    geom_tile(aes(x = ExpIDGenus, y = CommunityESVID, fill = AlignmentScore)) +
    # Each Sanger has one ESV match
    geom_tile(data = isolates_algn_score, aes(x = ExpIDGenus, y = CommunityESVID), fill = NA, color = "gold", linewidth = 1) +
    facet_wrap(.~Community, scales = "free", ncol = 4) +
    scale_fill_gradient(low = RColorBrewer::brewer.pal(n = 9, "Blues")[9], high = grey(0.9)) +
    scale_x_discrete(expand = c(0,0)) +
    scale_y_discrete(expand = c(0,0)) +
    theme_classic() +
    theme(
        panel.border = element_rect(fill = NA, color = 1, linewidth = 0.5),
        axis.text.x = element_text(angle = 30, vjust = 1.1, hjust = 1, size = 7),
        axis.text.y = element_text(size = 7),
        legend.position = c(.6, .05)
    ) +
    guides(fill = guide_legend(nrow = 1)) +
    labs()

ggsave(paste0(folder_data, "temp/16a-04-ESV_Sanger_score.png"), p, width = 12, height = 12)

# 6. plot the isolate abundance ----
# Check if the number is correct. Facet by community
isolates_abundance %>%
    mutate(Community = factor(Community, communities$Community)) %>%
    group_by(Community) %>% count()

p <- isolates_abundance %>%
    #filter(Community == "C8R4") %>%
    distinct(Community, CommunityESVID, .keep_all = T) %>%
    mutate(Community = factor(Community, communities$Community)) %>%
    ggplot() +
    geom_col(aes(x = CommunityESVID, y = RelativeAbundance, fill = ESVGenus)) +
    scale_fill_manual(values = c(RColorBrewer::brewer.pal(8, "Set1"), RColorBrewer::brewer.pal(8, "Set2"))) +
    facet_wrap(~Community, ncol = 3, scales = "free_x") +
    theme_classic() +
    theme(
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)
    ) +
    guides() +
    labs()

ggsave(paste0(folder_data, "temp/16a-06-matched_ESV_abundance_comm.png"), p, width = 12, height = 12)

# 7. plot total abundance by community ----
p1 <- isolates_abundance %>%
    distinct(Community, CommunityESVID, .keep_all = T) %>%
    mutate(Community = factor(Community, paste0("C", rep(1:12, each = 8), "R", rep(1:8, 12)))) %>%
    #mutate(Community = factor(Community, communities_remained$Community)) %>%
    ggplot() +
    geom_col(aes(x = Community, y = RelativeAbundance, fill = ESVGenus, group = CommunityESVID), color = "black", position = position_stack(reverse = T)) +
    scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
    scale_fill_manual(values = c(RColorBrewer::brewer.pal(8, "Set1"), RColorBrewer::brewer.pal(8, "Set2"))) +
    #facet_wrap(~Community, ncol = 3, scales = "free_x") +
    theme_classic() +
    theme(
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        panel.border = element_rect(color = 1, fill = NA, linewidth = 0.5)
    ) +
    guides() +
    labs() +
    ggtitle("ESV genus")

p2 <- isolates_abundance %>%
    distinct(Community, CommunityESVID, .keep_all = T) %>%
    mutate(Community = factor(Community, paste0("C", rep(1:12, each = 8), "R", rep(1:8, 12)))) %>%
    #mutate(Community = factor(Community, communities_remained$Community)) %>%
    ggplot() +
    geom_col(aes(x = Community, y = RelativeAbundance, fill = Genus, group = CommunityESVID), color = "black", position = position_stack(reverse = T)) +
    scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
    scale_fill_manual(values = c(RColorBrewer::brewer.pal(8, "Set1"), RColorBrewer::brewer.pal(8, "Set2"))) +
    #facet_wrap(~Community, ncol = 3, scales = "free_x") +
    theme_classic() +
    theme(
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        panel.border = element_rect(color = 1, fill = NA, linewidth = 0.5)
    ) +
    guides() +
    labs() +
    ggtitle("Sanger genus")
p <- plot_grid(p1, p2, ncol = 1, align = "v", axis = "lr")
ggsave(paste0(folder_data, "temp/16a-07-matched_ESV_abundance_color.png"), p, width = 6, height = 8)


# 8. plot total abundance by family ----
total_abundance <- isolates_abundance %>%
    distinct(Community, CommunityESVID, .keep_all = T) %>%
    group_by(Community) %>%
    summarize(TotalAbundance = round(sum(RelativeAbundance),2)) %>%
    mutate(Community = factor(Community, paste0("C", rep(1:12, each = 8), "R", rep(1:8, 12))))
total_abundance
mean(total_abundance$TotalAbundance)

p <- isolates_abundance %>%
    distinct(Community, CommunityESVID, .keep_all = T) %>%
    mutate(CommunityESVID = factor(CommunityESVID)) %>%
    group_by(Community) %>%
    mutate(ESVLabel = 1:n()) %>%
    mutate(Community = factor(Community, paste0("C", rep(1:12, each = 8), "R", rep(1:8, 12)))) %>%
    ggplot() +
    geom_col(aes(x = Community, y = RelativeAbundance, fill = Family, group = ESVLabel), color = "black") +
    geom_text(data = total_abundance, aes(x =  Community, label = TotalAbundance), y = 0.95, size = 2) +
    scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
    scale_fill_manual(values = c(RColorBrewer::brewer.pal(8, "Set1"), RColorBrewer::brewer.pal(8, "Set2"))) +
    #facet_wrap(~Community, ncol = 3, scales = "free_x") +
    theme_classic() +
    theme(
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        panel.border = element_rect(color = 1, fill = NA, linewidth = 0.5)
    ) +
    guides() +
    labs()

ggsave(paste0(folder_data, "temp/16a-08-matched_ESV_abundance_family.png"), p, width = 6, height = 4)

