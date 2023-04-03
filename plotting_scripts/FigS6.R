library(tidyverse)
library(cowplot)
library(broom)
source(here::here("analysis/00-metadata.R"))

communities <- read_csv(paste0(folder_data, "output/communities_remained.csv"), show_col_types = F) %>%
    mutate(CommunityLabel2 = paste0("#", CommunityLabel, ": ", Community))
isolates <- read_csv(paste0(folder_data, "output/isolates.csv"), show_col_types = F) # 68 isolates
sequences_alignment <- read_csv(paste0(folder_data, "temp/32-sequences_alignment.csv"), show_col_types = F) %>%
    arrange(Community, Isolate)

# Number of isolates with good match
sequences_alignment %>%
    filter(ConsensusLength >= 200, BasePairMismatch <= 4) %>%
    distinct(Community, Isolate) %>%
    nrow()

#
algn_Sanger_ESV <- sequences_alignment %>%
    filter(ConsensusLength >= 200, BasePairMismatch <= 4) %>%
    # For each Sanger, find the ESV which it has least base pair mismatch with
    group_by(ExpID) %>%
    filter(BasePairMismatch == min(BasePairMismatch)) %>%
    # When a Sanger has two ESV that both has zero, pick the first one
    slice(1) %>%
    ungroup()

nrow(algn_Sanger_ESV)

algn_Sanger_ESV %>% # Number of mismatches for each isolate Sanger
    group_by(BasePairMismatch) %>%
    count()

algn_Sanger_ESV1 <- algn_Sanger_ESV %>%
    group_by(Community, CommunityESVID) %>%
    arrange(BasePairMismatch) %>%
    slice(1)

table(algn_Sanger_ESV1$BasePairMismatch) # table of base pair mismatch
range(algn_Sanger_ESV1$ConsensusLength) # range of consensus length
nrow(algn_Sanger_ESV1) # Numer of unique matches

clean_isolate_names <- function (x) {
    y <- x %>%
        group_by(Community) %>%
        mutate(ExpIDGenus = paste0(ExpID, "-", Genus)) %>%
        mutate(IsolateGenus = paste0(Isolate, "-", Genus))
    y %>%
        mutate(IsolateGenus = factor(IsolateGenus, unique(y$IsolateGenus)))
}

isolates_algn_bp_mismatch <- isolates %>%
    select(Community, Isolate, ExpID, Genus) %>%
    left_join(algn_Sanger_ESV) %>%
    left_join(communities) %>%
    mutate(CommunityLabel2 = factor(CommunityLabel2, communities$CommunityLabel2)) %>%
    clean_isolate_names

isolates_algn_bp_mismatch1 <- isolates %>%
    select(Community, Isolate, ExpID, Genus) %>%
    left_join(algn_Sanger_ESV1) %>%
    left_join(communities) %>%
    mutate(Community = factor(Community, communities$Community)) %>%
    clean_isolate_names

# Heatmap, ESV and sanger alignment, binned by mismatch
p <- sequences_alignment %>%
    left_join(communities) %>%
    mutate(CommunityLabel2 = factor(CommunityLabel2, communities$CommunityLabel2)) %>%
    #mutate(Community = factor(Community, communities$Community)) %>%
    clean_isolate_names %>%
    mutate(BasePairMismatchCategory = case_when(
        BasePairMismatch == 0 ~ "0",
        BasePairMismatch == 1 ~ "1",
        BasePairMismatch == 2 ~ "2",
        BasePairMismatch == 3 ~ "3",
        BasePairMismatch == 4 ~ "4",
        BasePairMismatch >= 5 ~ ">=5"
    ) %>% ordered(c(0:4, ">=5"))
    ) %>%
    ggplot() +
    geom_tile(aes(x = IsolateGenus, y = CommunityESVID, fill = BasePairMismatchCategory)) +
    # Each Sanger has one ESV match
    geom_tile(data = isolates_algn_bp_mismatch, aes(x = IsolateGenus, y = CommunityESVID), fill = NA, color = "gold", linewidth = 1) +
    scale_fill_manual(values = rev(RColorBrewer::brewer.pal(n = 6, "Blues"))) +
    scale_x_discrete(expand = c(0,0)) +
    scale_y_discrete(expand = c(0,0)) +
    facet_wrap(.~CommunityLabel2, scales = "free", ncol = 4) +
    theme_classic() +
    theme(
        panel.border = element_rect(fill = NA, color = 1, linewidth = 0.5),
        axis.text.x = element_text(angle = 45, vjust = 1.05, hjust = 1, size = 6),
        axis.text.y = element_text(size = 6),
        legend.position = c(.6, .05),
        legend.spacing.x = unit(4, "mm"),
        strip.background = element_blank(),
        plot.title = element_text(hjust = 0.5)
    ) +
    guides(fill = guide_legend(title = "base pair mismatch", nrow = 1, override.aes = list(color = "black"))) +
    labs(x = "isolate", y = "ESV") +
    ggtitle("community")


ggsave(here::here("plots/FigS6-Sanger_ESV_alignment.png"), p, width = 10, height = 10)


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



















