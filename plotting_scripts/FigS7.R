library(tidyverse)
library(cowplot)
library(broom)
source(here::here("processing_scripts/00-metadata.R"))

factorize_communities <- function (x) x %>% mutate(Community = factor(Community, paste0("C", rep(1:12, each = 8), "R", rep(1:8, 12))))
communities <- read_csv(paste0(folder_data, "output/communities_remained.csv"), show_col_types = F) %>%
    mutate(CommunityLabel2 = paste0("#", CommunityLabel, ": ", Community))
sequences_alignment <- read_csv(paste0(folder_data, "temp/16-sequences_alignment.csv"), show_col_types = F)
isolates <- read_csv(paste0(folder_data, "output/isolates.csv"), show_col_types = F)


clean_isolate_names <- function (x) {
    y <- x %>%
        group_by(Community) %>%
        arrange(Isolate) %>%
        mutate(ExpIDGenus = paste0(ExpID, "-", Genus)) %>%
        mutate(IsolateGenus = paste0(Isolate, "-", Genus))
    y %>%
        mutate(IsolateGenus = factor(IsolateGenus, unique(y$IsolateGenus)))
}

isolates_algn_bp_mismatch <- isolates %>%
    filter(!is.na(BasePairMismatch)) %>%
    select(Community, Isolate, CommunityESVID, ExpID, Genus) %>%
    left_join(communities) %>%
    mutate(CommunityLabel2 = factor(CommunityLabel2, communities$CommunityLabel2)) %>%
    clean_isolate_names

nrow(isolates_algn_bp_mismatch) # 62 isolates that match ESV

p <- sequences_alignment %>%
    left_join(communities) %>%
    mutate(CommunityLabel2 = factor(CommunityLabel2, communities$CommunityLabel2)) %>%
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
    facet_wrap(.~CommunityLabel2, scales = "free", ncol = 3) +
    theme_classic() +
    theme(
        panel.border = element_rect(fill = NA, color = 1, linewidth = 0.5),
        axis.text.x = element_text(angle = 45, vjust = 1.05, hjust = 1, size = 6),
        axis.text.y = element_text(size = 6),
        legend.position = "bottom",
        legend.spacing.x = unit(4, "mm"),
        strip.background = element_blank(),
        plot.title = element_text(hjust = 0.5)
    ) +
    guides(fill = guide_legend(title = "base pair mismatch", nrow = 1, override.aes = list(color = "black"))) +
    labs(x = "isolate", y = "ESV") +
    ggtitle("community")


ggsave(here::here("plots/FigS7.png"), p, width = 10, height = 12, dpi = 500)

#
sequences_alignment %>% distinct(Community, CommunityESVID) %>% nrow() # 102 ESVs
sequences_alignment %>% distinct(Community, ExpID) %>% nrow() # 65 isolates

# 62 isolates match to ESVs
isolates %>%
    drop_na(BasePairMismatch) %>%
    pull(BasePairMismatch) %>%
    table()

# 0  1  2  3  4
# 40 11  6  2  3
# 40 isolates has full match, 11 has one mismatch, 6 has two mismatches, 2 has three mismatches, and 3 has four mismatches

# 44 ESVs match to isolates
isolates %>%
    drop_na(CommunityESVID) %>%
    group_by(Community, CommunityESVID) %>%
    count() %>%
    pull(n) %>%
    table

# 1  2  3  4
# 32  8  2  2
# 32 ESVs match 1 isolate, 8 ESVs match 2 isolates, 2 ESVs match 3 isolates, 2 ESVs match 4 isolates



