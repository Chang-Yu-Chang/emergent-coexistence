library(tidyverse)
source(here::here("analysis/00-metadata.R"))

communities <- read_csv(paste0(folder_data, "output/communities_remained.csv"), show_col_types = F) %>%
    bind_rows(tibble(Community = "C10R2", CommunityLabel = NA)) %>%
    mutate(Community = factor(Community)) %>%
    mutate(CommunityLabel2 = paste0("#", CommunityLabel, ": ", Community))

isolates_RDP <- read_csv(paste0(folder_data, "temp/12-isolates_RDP.csv"), show_col_types = F)
pairs_mismatch <- read_csv(paste0(folder_data, "temp/22-pairs_mismatch.csv"), show_col_types = F) %>%
    mutate(ID1 = as.character(ID1), ID2 = as.character(ID2)) %>%
    left_join(select(isolates_RDP, ID1 = ID, Genus1 = Genus)) %>%
    left_join(select(isolates_RDP, ID2 = ID, Genus2 = Genus))


clean_isolate_names <- function (x) {
    y <- x %>%
        group_by(Community) %>%
        arrange(Isolate1, Isolate2) %>%
        mutate(IsolateGenus1 = paste0(Isolate1, "-", Genus1)) %>%
        mutate(IsolateGenus2 = paste0(Isolate2, "-", Genus2))
    y %>%
        mutate(IsolateGenus1 = factor(IsolateGenus1, unique(y$IsolateGenus1))) %>%
        mutate(IsolateGenus2 = factor(IsolateGenus2, rev(unique(y$IsolateGenus2))))
}


p <- pairs_mismatch %>%
    left_join(communities) %>%
    clean_isolate_names %>%
    mutate(CommunityLabel2 = factor(CommunityLabel2, communities$CommunityLabel2)) %>%
    ggplot() +
    geom_tile(aes(x = IsolateGenus1, y = IsolateGenus2, fill = Mismatch), color = "black", linewidth = 0.1) +
    #scale_fill_manual(values = rev(RColorBrewer::brewer.pal(n = 6, "Blues"))) +
    scale_fill_gradient(low = "dodgerblue3", high = "white") +
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
    guides(fill = guide_legend(direction = "horizontal")) +
    labs()

ggsave(paste0(folder_data, "temp/22a-01-pairs_mismatch.png"), p, width = 12, height = 12)
