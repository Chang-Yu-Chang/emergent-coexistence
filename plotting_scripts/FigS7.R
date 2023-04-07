library(tidyverse)
library(cowplot)
library(broom)
source(here::here("analysis/00-metadata.R"))

factorize_communities <- function (x) x %>% mutate(Community = factor(Community, paste0("C", rep(1:12, each = 8), "R", rep(1:8, 12))))
communities <- read_csv(paste0(folder_data, "output/communities_remained.csv"), show_col_types = F) %>%
    bind_rows(tibble(Community = "C10R2", CommunityLabel = NA)) %>%
    mutate(CommunityLabel2 = paste0("#", CommunityLabel, ": ", Community))
sequences_alignment <- read_csv(paste0(folder_data, "temp/16-sequences_alignment.csv"), show_col_types = F)
isolates <- read_csv(paste0(folder_data, "output/isolates.csv"), show_col_types = F)

clean_isolate_names <- function (x) {
    y <- x %>%
        group_by(Community) %>%
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

nrow(isolates_algn_bp_mismatch) # 64 isolates

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


ggsave(here::here("plots/FigS7-Sanger_ESV_alignment.png"), p, width = 10, height = 12)


# ESV richness per community
communities_abundance <- read_csv(paste0(folder_data, "temp/14-communities_abundance.csv"), show_col_types = F) %>% factorize_communities

communities_abundance %>%
    filter(Community %in% communities$Community) %>%
    filter(Transfer == 12) %>%
    #filter(Relative_Abundance > 0.01) %>%
    group_by(Community) %>%
    count() %>%
    left_join(communities) %>%
    arrange(CommunityLabel)

# Matched ESV richness per community
isolates %>%
    drop_na(BasePairMismatch) %>%
    group_by(Community) %>%
    count()



