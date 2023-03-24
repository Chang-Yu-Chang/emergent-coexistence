library(tidyverse)
library(cowplot)
library(broom)
source(here::here("analysis/00-metadata.R"))

communities <- read_csv(paste0(folder_data, "temp/00c-communities.csv"), show_col_types = F) %>%
    mutate(Community = factor(Community, Community))
communities_abundance <- read_csv(paste0(folder_data, "raw/community_ESV/Emergent_Comunity_Data.csv"), show_col_types = F) %>%
    filter(Carbon_Source == "Glucose" | Carbon_Source == "Original") %>%
    mutate(Community = factor(paste0("C", Inoculum, "R", Replicate), paste0("C", rep(1:12, each = 8), "R", rep(1:8, 12))))%>%
    arrange(Community, Family, Transfer, ESV)


curate_abundant_genus <- function () {
    abundant_genus <- c(
        "Enterobacteriaceae","Klebsiella", "Enterobacter", "Raoultella", "Citrobacter", "Salmonella", "Pantoea", "Yersinia",
        "Pseudomonadaceae", "Pseudomonas", "Azomonas", "Aeromonas", "Acinetobacter", "Enterococcus", "Bordetella",
        "Stenotrophomonas", "Delftia", "Serratia"
    )
    abundant_genus <- paste0(rep(abundant_genus, each = 101),
                             rep(c("", paste0(".", 1:100)), length(abundant_genus)))
    return(abundant_genus)
}
abundant_genus <- curate_abundant_genus()
bin_ESV_names <- function (comm_abundance) {
    comm_abundance %>%
        #mutate(ESV_ID = str_replace(ESV_ID, "Enterobacteriaceae", "Enterobacter")) %>%
        mutate(ESV_ID = case_when(
            str_detect(ESV_ID, "Enterobacteriaceae") ~ str_replace(ESV_ID, "\\.\\d+", ""),
            str_detect(ESV_ID, "Pantoea") ~ str_replace(ESV_ID, "\\.\\d+", ""),
            str_detect(ESV_ID, "Raoultella") ~ str_replace(ESV_ID, "\\.\\d+", ""),
            str_detect(ESV_ID, "Salmonella") ~ str_replace(ESV_ID, "\\.\\d+", ""),
            str_detect(ESV_ID, "Citrobacter") ~ str_replace(ESV_ID, "\\.\\d+", ""),
            str_detect(ESV_ID, "Enterobacter") ~ str_replace(ESV_ID, "\\.\\d+", ""),
            str_detect(ESV_ID, "Klebsiella") ~ str_replace(ESV_ID, "\\.\\d+", ""),
            str_detect(ESV_ID, "Acinetobacter") ~ str_replace(ESV_ID, ".\\d+", ""),
            #str_detect(ESV_ID, "Pseudomonadaceae") ~ str_replace(ESV_ID, "Pseudomonadaceae", "Pseudomonadas"),
            str_detect(ESV_ID, "Pseudomonadaceae") ~ str_replace(ESV_ID, "\\.\\d+", ""),
            str_detect(ESV_ID, "Pseudomonas.1\\d+") ~ str_replace(ESV_ID, "\\.1\\d+", "\\.1"),
            str_detect(ESV_ID, "Pseudomonas.2\\d+") ~ str_replace(ESV_ID, "\\.2\\d+", "\\.2"),
            str_detect(ESV_ID, "Pseudomonas.3\\d+") ~ str_replace(ESV_ID, "\\.3\\d+", "\\.3"),
            str_detect(ESV_ID, "Pseudomonas.4\\d+") ~ str_replace(ESV_ID, "\\.4\\d+", "\\.4"),
            str_detect(ESV_ID, "Pseudomonas.5\\d+") ~ str_replace(ESV_ID, "\\.5\\d+", "\\.5"),
            str_detect(ESV_ID, "Pseudomonas.5\\d+") ~ str_replace(ESV_ID, "\\.5\\d+", "\\.5"),
            str_detect(ESV_ID, "Azomonas") ~ str_replace(ESV_ID, ".\\d+", ""),
            str_detect(ESV_ID, "Stenotrophomonas") ~ str_replace(ESV_ID, ".\\d+", ""),
            str_detect(ESV_ID, "Aeromonas") ~ str_replace(ESV_ID, ".\\d+", ""),
            str_detect(ESV_ID, "Delftia") ~ str_replace(ESV_ID, ".\\d+", ""),
            str_detect(ESV_ID, "Enterococcus") ~ str_replace(ESV_ID, ".\\d+", ""),
            T ~ ESV_ID
        ))
}
clean_ESV_names <- function (comm_abundance) {
    comm_abundance %>%
        mutate(ESV_ID = ifelse(!(ESV_ID %in% abundant_genus), "Other", ESV_ID)) %>%
        return()
}
get_ESV_colors <- function (comm_abundance) {
    communities_abundance_ESV_ID <- comm_abundance %>%
        distinct(Family, ESV_ID) %>%
        mutate(ESV_ID = factor(ESV_ID, abundant_genus)) %>%
        arrange(Family, ESV_ID)

    # Entero
    temp1 <- communities_abundance_ESV_ID %>%
        filter(Family %in% c("Enterobacteriaceae")) %>%
        drop_na() %>%
        #mutate(ESV_color = viridis::viridis_pal(option = "viridis")(n())) %>%
        #mutate(ESV_color = scales::div_gradient_pal(low = "#313797", mid = "#ffffbf", high = "#a50026", space = "Lab")(seq(0, 1, length.out = n()))) %>%
        filter(ESV_ID != "Other")

    if (nrow(temp1) >= 12) {
        temp1 <- temp1 %>% mutate(ESV_color = c(rev(RColorBrewer::brewer.pal(11, "RdYlBu")), rev(RColorBrewer::brewer.pal(9, "OrRd")))[1:n()])
    } else if (nrow(temp1) >= 3) {
        temp1 <- temp1 %>% mutate(ESV_color = rev(RColorBrewer::brewer.pal(n(), "RdYlBu")))
    } else if (nrow(temp1) < 3) {
        temp1 <- temp1 %>% mutate(ESV_color = rev(RColorBrewer::brewer.pal(3, "RdYlBu"))[1:n()])
    }

    # Pseudo
    temp2 <- communities_abundance_ESV_ID %>%
        filter(Family == c("Pseudomonadaceae")) %>%
        drop_na() %>%
        #mutate(ESV_color = scales::seq_gradient_pal(low = "#d73027", high = "#fef1e7", space = "Lab")(seq(0, 1, length.out = n()))) %>%
        #mutate(ESV_color = scales::div_gradient_pal(low = "#00451a", mid = "#f6f6f7", high = "#41004a", space = "Lab")(seq(0, 1, length.out = n()))) %>%
        filter(ESV_ID != "Other")

    if (nrow(temp2) >= 12) {
        temp2 <- temp2 %>% mutate(ESV_color = c(rev(RColorBrewer::brewer.pal(11, "PRGn")), rev(RColorBrewer::brewer.pal(9, "YlGn")))[1:n()])
    } else if (nrow(temp2) >= 3) {
        temp2 <- temp2 %>% mutate(ESV_color = rev(RColorBrewer::brewer.pal(n(), "PRGn")))
    } else if (nrow(temp2) < 3) {
        temp2 <- temp2 %>% mutate(ESV_color =rev(RColorBrewer::brewer.pal(3, "PRGn"))[1:n()])
    }

    # Other
    temp3 <- communities_abundance_ESV_ID %>%
        filter(!(Family %in% c("Enterobacteriaceae", "Pseudomonadaceae")) & ESV_ID != "Other") %>%
        mutate(ESV_color = c(RColorBrewer::brewer.pal(8, "Set2"), RColorBrewer::brewer.pal(8, "Set3"))[1:n()]) %>%
        filter(ESV_ID != "Other")

    ESV_colors <- c(temp1$ESV_color %>% setNames(temp1$ESV_ID),
                    temp2$ESV_color %>% setNames(temp2$ESV_ID),
                    temp3$ESV_color %>% setNames(temp3$ESV_ID),
                    Other = "#999998")


    return(ESV_colors)

}

# Shared color code
ESV_colors <- communities_abundance %>%
    bin_ESV_names() %>%
    clean_ESV_names() %>%
    get_ESV_colors


# T12 composition. Replicate reflect the actual order. Fig. S6 in Goldford2018
communities_abundance_T12 <- communities_abundance %>%
    filter(Transfer == 12) %>%
    bin_ESV_names() %>%
    clean_ESV_names()

p1 <- communities_abundance_T12 %>%
    mutate(ESV_ID = ifelse(ESV_ID %in% names(ESV_colors[-length(ESV_colors)]), ESV_ID, "Other")) %>%
    mutate(ESV_ID = factor(ESV_ID, rev(names(ESV_colors)))) %>%
    #mutate(Community = factor(Community, rev(paste0("C", rep(1:12, each = 8), "R", rep(1:8, 12))))) %>%
    ggplot() +
    geom_col(aes(x = Community, y = Relative_Abundance, fill = ESV_ID), width = 1, color = NA) +
    #scale_x_discrete(position = "top") +
    scale_y_continuous(breaks = seq(0, 1, 0.5), expand = c(0,0)) +
    scale_fill_manual(values = ESV_colors, breaks = names(ESV_colors)) +
    facet_wrap(Inoculum~., scales = "free_x", nrow = 1) +
    theme_classic() +
    theme(
        strip.background = element_rect(fill = NA, color = NA),
        panel.border = element_rect(fill = NA, color = 1, linewidth = 0.5),
        panel.spacing.x = unit(0, "mm"),
        axis.text.x.bottom = element_text(angle = 90, hjust = 0.5),
        axis.text.x.top = element_blank(),
        axis.line.x.top = element_blank(),
        axis.ticks.x.top = element_blank(),
        plot.title = element_text(size = 10, hjust = 0.5)
    ) +
    guides(color = "none", fill = "none", x.sec = guide_axis(title = "inoculum", check.overlap = T)) +
    labs(x = "community", y = "relative abundance")

if (FALSE) {
    # Vertical facet wrap
p1 <- communities_abundance_T12 %>%
    mutate(ESV_ID = ifelse(ESV_ID %in% names(ESV_colors[-length(ESV_colors)]), ESV_ID, "Other")) %>%
    mutate(ESV_ID = factor(ESV_ID, rev(names(ESV_colors)))) %>%
    mutate(Community = factor(Community, rev(paste0("C", rep(1:12, each = 8), "R", rep(1:8, 12))))) %>%
    ggplot() +
    geom_col(aes(x = Community, y = Relative_Abundance, fill = ESV_ID), width = 1, color = NA) +
    scale_x_discrete(position = "top") +
    scale_y_continuous(breaks = seq(0, 1, 0.5), expand = c(0,0)) +
    scale_fill_manual(values = ESV_colors, breaks = names(ESV_colors)) +
    coord_flip() +
    facet_wrap(Inoculum~., scales = "free_y", ncol = 1, strip.position = "left") +
    theme_classic() +
    theme(
        strip.background = element_rect(fill = NA, color = NA),
        panel.border = element_rect(fill = NA, color = 1, linewidth = 0.5),
        panel.spacing.y = unit(0, "mm"),
        axis.text.y.left = element_blank(),
        axis.line.y.left = element_blank(),
        axis.ticks.y.left = element_blank(),
        plot.title = element_text(size = 10, hjust = 0.5)
    ) +
    guides(color = "none", fill = "none", y.sec = guide_axis(title = "inoculum", check.overlap = T)) +
    labs(x = "community", y = "relative abundance")

}


# Temporal dynamics of all 26 communities -----
communities_abundance_temporal <- communities_abundance %>%
    distinct(Community, Transfer) %>%
    arrange(Community, Transfer) %>%
    filter(Transfer == 4)
communities_abundance_temporal_26 <- communities_abundance %>%
    filter(Community %in% communities_abundance_temporal$Community | Transfer == 0) %>%
    filter(!is.na(Community)) %>%
    bin_ESV_names() %>%
    clean_ESV_names()

communities_abundance_T0_for_plot <- bind_rows(
    communities_abundance %>% filter(Transfer == 0) %>% mutate(Community = paste0("C", Inoculum, "R4")),
    communities_abundance %>% filter(Transfer == 0, Inoculum == 2) %>% slice(rep(1:n(), 7)) %>% mutate(Community = rep(paste0("C", rep(2,7), "R", c(1:3, 5:8)), each = n()/7)),
    communities_abundance %>% filter(Transfer == 0, Inoculum == 6) %>% slice(rep(1:n(), 7)) %>% mutate(Community = rep(paste0("C", rep(6,7), "R", c(1:3, 5:8)), each = n()/7))
)


p2 <- communities_abundance_temporal_26 %>%
    bind_rows(communities_abundance_T0_for_plot) %>%
    mutate(Community = factor(Community, c(paste0("C", rep(c(2,6), each = 8), "R", rep(1:8,2)), paste0("C", c(1,3:5, 7:12),"R", 4)))) %>%
    mutate(ESV_ID = ifelse(ESV_ID %in% names(ESV_colors[-length(ESV_colors)]), ESV_ID, "Other")) %>%
    mutate(ESV_ID = factor(ESV_ID, rev(names(ESV_colors)))) %>%
    ggplot() +
    geom_col(aes(x = Transfer, y = Relative_Abundance, fill = ESV_ID), width = 1, linewidth = .1) +
    scale_x_continuous(breaks = seq(0, 12, 2), limits = c(-0.5, 12.5), expand = c(0,0)) +
    scale_y_continuous(breaks = seq(0, 1, 0.5), expand = c(0,0)) +
    scale_fill_manual(values = ESV_colors, breaks = names(ESV_colors)) +
    facet_wrap(~Community, ncol = 4, labeller = labeller(Inoculum = label_both)) +
    theme_classic() +
    theme(
        strip.background = element_rect(fill = NA, color = NA),
        panel.border = element_rect(fill = NA, color = 1, linewidth = 1),
        legend.title = element_blank(),
        legend.key.size = unit(6, "mm"),
        legend.position = "right"
    ) +
    guides(color = "none", fill = guide_legend(ncol = 1)) +
    labs(x = "transfer", y = "relative abundance")

# p <- plot_grid(
#     p1,
#     p2,
#     nrow = 2, labels = LETTERS[1:2], scale = .95, rel_heights = c(1,3)
# ) + paint_white_background()
p <- p2
ggsave(here::here("plots/FigS1-temporal_dynamics.png"), p, width = 8, height = 8)


