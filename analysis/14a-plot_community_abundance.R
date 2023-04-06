#' This script plots the community abundance
library(tidyverse)
source(here::here("analysis/00-metadata.R"))

communities <- read_csv(paste0(folder_data, "temp/00c-communities.csv"), show_col_types = F) %>%
    mutate(Community = factor(Community, Community))

communities_abundance <- read_csv(paste0(folder_data, 'temp/14-communities_abundance.csv'), show_col_types = F)

communities_temporal <- communities_abundance %>%
    distinct(Community, Transfer) %>%
    arrange(Community, Transfer) %>%
    filter(Transfer == 4)

curate_abundant_genus <- function () {
    abundant_genus <- c(
        "Enterobacteriaceae","Klebsiella", "Enterobacter", "Raoultella", "Citrobacter", "Salmonella", "Pantoea", "Yersinia",
        "Pseudomonadaceae", "Pseudomonas", "Azomonas", "Aeromonas", "Acinetobacter", "Enterococcus",
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
    #comm_abundance <- temp
    communities_abundance_ESV_ID <- comm_abundance %>%
        #filter(Transfer == 12) %>%
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
        mutate(ESV_color = RColorBrewer::brewer.pal(9, "Set1")[1:n()]) %>%
        filter(ESV_ID != "Other")

    ESV_colors <- c(temp1$ESV_color %>% setNames(temp1$ESV_ID),
                    temp2$ESV_color %>% setNames(temp2$ESV_ID),
                    temp3$ESV_color %>% setNames(temp3$ESV_ID),
                    Other = "#999999")


    return(ESV_colors)

}


# T12 composition. Replicate reflect the actual order. Fig. S6 in Goldford2018 ----
communities_abundance_T12 <- communities_abundance %>%
    filter(Transfer == 12) %>%
    bin_ESV_names() %>%
    clean_ESV_names()

ESV_colors <- communities_abundance_T12 %>% get_ESV_colors

p <- communities_abundance_T12 %>%
    mutate(ESV_ID = ifelse(ESV_ID %in% names(ESV_colors[-length(ESV_colors)]), ESV_ID, "Other")) %>%
    mutate(ESV_ID = factor(ESV_ID, rev(names(ESV_colors)))) %>%
    ggplot() +
    geom_col(aes(x = Replicate, y = Relative_Abundance, fill = ESV_ID), width = .8, color = NA) +
    scale_x_continuous(breaks = 1:8, limits = c(0.5,8.5), expand = c(0,0)) +
    scale_y_continuous(breaks = seq(0, 1, 0.5), expand = c(0,0)) +
    scale_fill_manual(values = ESV_colors, breaks = names(ESV_colors)) +
    coord_flip() +
    facet_wrap(~Inoculum, labeller = labeller(Inoculum = label_both), ncol = 3, scales = "free_y") +
    theme_classic() +
    theme(strip.background = element_rect(fill = NA, color = NA),
          panel.border = element_rect(fill = NA, color = 1, linewidth = 1)) +
    guides(color = "none") +
    labs(x = "replicate", y = "relative abundance")

ggsave(paste0(folder_data, "temp/14a-01-Goldford2018_FigS6_actual_replicate.png"), p, width = 10, height = 8)

# T12 composition. Replicate reordered according to the ESV abundance Fig.S6 in Goldford2018 ----
communities_abundance_reordered <- communities_abundance %>%
    filter(Transfer == 12, !is.na(Community)) %>%
    select(ESV_ID, Community, Inoculum, Replicate, Relative_Abundance) %>%
    bin_ESV_names() %>%
    clean_ESV_names() %>%
    group_by(ESV_ID, Community, Inoculum, Replicate) %>%
    summarize(Relative_Abundance = sum(Relative_Abundance)) %>%
    ungroup() %>%
    pivot_wider(names_from = ESV_ID, values_from = Relative_Abundance, values_fill = 0) %>%
    # Replicate ordered according to ESV abundance
    arrange(Inoculum, Klebsiella, Enterobacteriaceae, Raoultella, Citrobacter, Pseudomonas) %>%
    group_by(Inoculum) %>%
    mutate(ReplicateReordered = 1:n()) %>%
    select(Community, Inoculum, Replicate, ReplicateReordered)


communities_abundance_T12_ordered <- communities_abundance %>%
    left_join(communities_abundance_reordered) %>%
    filter(Transfer == 12) %>%
    bin_ESV_names() %>%
    clean_ESV_names()

ESV_colors <- communities_abundance_T12_ordered %>% get_ESV_colors

p <- communities_abundance_T12_ordered %>%
    mutate(ESV_ID = ifelse(ESV_ID %in% names(ESV_colors[-length(ESV_colors)]), ESV_ID, "Other")) %>%
    mutate(ESV_ID = factor(ESV_ID, rev(names(ESV_colors)))) %>%
    ggplot() +
    geom_col(aes(x = ReplicateReordered, y = Relative_Abundance, fill = ESV_ID), width = .8, color = NA) +
    scale_x_continuous(breaks = 1:8, limits = c(0.5,8.5), expand = c(0,0)) +
    scale_y_continuous(breaks = seq(0, 1, 0.5), expand = c(0,0)) +
    scale_fill_manual(values = ESV_colors, breaks = names(ESV_colors)) +
    coord_flip() +
    facet_wrap(~Inoculum, labeller = labeller(Inoculum = label_both), ncol = 3, scales = "free_y") +
    theme_classic() +
    theme(strip.background = element_rect(fill = NA, color = NA),
          panel.border = element_rect(fill = NA, color = 1, linewidth = 1)) +
    guides(color = "none") +
    labs(x = "replicate (reordered)", y = "relative abundance")

ggsave(paste0(folder_data, "temp/14a-02-Goldford2018_FigS6.png"), p, width = 10, height = 8)



# Temporal dynamics of all CXXR4 communities. One replicate per inoculum. Fig.S2 in Goldford2018 ----
communities_abundance_R4 <- communities_abundance %>%
    filter(Community %in% communities_temporal$Community | Transfer == 0) %>%
    filter(Replicate == 4 | Transfer == 0) %>%
    bin_ESV_names() %>%
    clean_ESV_names()

ESV_colors <- communities_abundance_R4 %>% get_ESV_colors

p <- communities_abundance_R4 %>%
    mutate(Community = factor(Community, paste0("C", rep(1:12, each = 8), "R", rep(1:8, 12)))) %>%
    mutate(ESV_ID = ifelse(ESV_ID %in% names(ESV_colors[-length(ESV_colors)]), ESV_ID, "Other")) %>%
    mutate(ESV_ID = factor(ESV_ID, rev(names(ESV_colors)))) %>%
    ggplot() +
    geom_col(aes(x = Transfer, y = Relative_Abundance, fill = ESV_ID), width = 1, linewidth = .1) +
    scale_x_continuous(breaks = seq(0, 12, 2), limits = c(-0.5, 12.5), expand = c(0,0)) +
    scale_y_continuous(breaks = seq(0, 1, 0.2), expand = c(0,0)) +
    scale_fill_manual(values = ESV_colors, breaks = names(ESV_colors)) +
    #scale_color_manual(values = ESV_colors, breaks = names(ESV_colors)) +
    facet_wrap(~Inoculum, ncol = 3, scales = "free", labeller = labeller(Inoculum = label_both)) +
    theme_classic() +
    theme(strip.background = element_rect(fill = NA, color = NA),
          panel.border = element_rect(fill = NA, color = 1, linewidth = 1)) +
    guides(color = "none") +
    labs(x = "transfer", y = "relative abundance")
ggsave(paste0(folder_data, "temp/14a-03-Goldford2018_FigS2.png"), p, width = 10, height = 8)


# Temporal dynamics of all C2RXX communities.  ----
communities_abundance_T0 <- communities_abundance %>%
    filter(Transfer == 0, Inoculum == 2) %>%
    slice(rep(1:n(), each = 8)) %>%
    mutate(Community = paste0("C2R", rep(1:8, n()/8)))

communities_abundance_C2 <- communities_abundance %>%
    filter(Community %in% communities_temporal$Community) %>%
    filter(Inoculum == 2) %>%
    bin_ESV_names() %>%
    clean_ESV_names()

ESV_colors <- communities_abundance_C2 %>% get_ESV_colors

p <- communities_abundance_C2 %>%
    bind_rows(communities_abundance_T0) %>%
    mutate(Community = factor(Community, paste0("C", rep(1:12, each = 8), "R", rep(1:8, 12)))) %>%
    mutate(ESV_ID = ifelse(ESV_ID %in% names(ESV_colors[-length(ESV_colors)]), ESV_ID, "Other")) %>%
    mutate(ESV_ID = factor(ESV_ID, rev(names(ESV_colors)))) %>%
    ggplot() +
    geom_col(aes(x = Transfer, y = Relative_Abundance, fill = ESV_ID), width = 1) +
    scale_x_continuous(breaks = seq(0, 12, 2), limits = c(-0.5, 12.5), expand = c(0,0)) +
    scale_y_continuous(breaks = seq(0, 1, 0.2), expand = c(0,0)) +
    scale_fill_manual(values = ESV_colors, breaks = names(ESV_colors)) +
    facet_wrap(~Community, ncol = 2, scales = "free") +
    theme_classic() +
    theme(strip.background = element_rect(fill = NA, color = NA),
          panel.border = element_rect(fill = NA, color = 1, linewidth = 1)) +
    guides(color = "none") +
    labs(x = "transfer", y = "relative abundance")
ggsave(paste0(folder_data, "temp/14a-04-Goldford2018_Inoculum2.png"), p, width = 6, height = 8)


# Temporal dynamics of all C6RXX communities. Fig S7 in Goldford2018 ----
communities_abundance_T0 <- communities_abundance %>%
    filter(Transfer == 0, Inoculum == 6) %>%
    slice(rep(1:n(), each = 8)) %>%
    mutate(Community = paste0("C6R", rep(1:8, n()/8)))

communities_abundance_C6 <- communities_abundance %>%
    filter(Community %in% communities_temporal$Community) %>%
    filter(Inoculum == 6) %>%
    bin_ESV_names() %>%
    clean_ESV_names()

ESV_colors <- communities_abundance_C6 %>% get_ESV_colors

p <- communities_abundance_C6 %>%
    bind_rows(communities_abundance_T0) %>%
    mutate(Community = factor(Community, paste0("C", rep(1:12, each = 8), "R", rep(1:8, 12)))) %>%
    mutate(ESV_ID = ifelse(ESV_ID %in% names(ESV_colors[-length(ESV_colors)]), ESV_ID, "Other")) %>%
    mutate(ESV_ID = factor(ESV_ID, rev(names(ESV_colors)))) %>%
    ggplot() +
    geom_col(aes(x = Transfer, y = Relative_Abundance, fill = ESV_ID), width = 1) +
    scale_x_continuous(breaks = seq(0, 12, 2), limits = c(-0.5, 12.5), expand = c(0,0)) +
    scale_y_continuous(breaks = seq(0, 1, 0.2), expand = c(0,0)) +
    scale_fill_manual(values = ESV_colors, breaks = names(ESV_colors)) +
    facet_wrap(~Community, ncol = 2, scales = "free") +
    theme_classic() +
    theme(strip.background = element_rect(fill = NA, color = NA),
          panel.border = element_rect(fill = NA, color = 1, linewidth = 1)) +
    guides(color = "none") +
    labs(x = "transfer", y = "relative abundance")
ggsave(paste0(folder_data, "temp/14a-05-Goldford2018_FigS7.png"), p, width = 6, height = 8)




# Temporal dynamics of only the four communities used for pairwise competition experiments ----
communities_abundance_T0 <- bind_rows(
    communities_abundance %>% filter(Transfer == 0, Inoculum == 1) %>% mutate(Community = "C1R4"),
    communities_abundance %>% filter(Transfer == 0, Inoculum == 2) %>% mutate(Community = "C2R6"),
    communities_abundance %>% filter(Transfer == 0, Inoculum == 2) %>% mutate(Community = "C2R8"),
    communities_abundance %>% filter(Transfer == 0, Inoculum == 8) %>% mutate(Community = "C8R4"),
)

temp <- communities_abundance %>%
    filter(Community %in% communities_temporal$Community) %>%
    filter(Community %in% communities$Community) %>%
    bin_ESV_names() %>%
    clean_ESV_names()

ESV_colors <- temp %>% get_ESV_colors

p <- temp %>%
    bind_rows(communities_abundance_T0) %>%
    mutate(Community = factor(Community, paste0("C", rep(1:12, each = 8), "R", rep(1:8, 12)))) %>%
    mutate(ESV_ID = ifelse(ESV_ID %in% names(ESV_colors[-length(ESV_colors)]), ESV_ID, "Other")) %>%
    mutate(ESV_ID = factor(ESV_ID, rev(names(ESV_colors)))) %>%
    ggplot() +
    geom_col(aes(x = Transfer, y = Relative_Abundance, fill = ESV_ID)) +
    scale_x_continuous(breaks = 1:12, limits = c(0.5, 12.5), expand = c(0,0)) +
    scale_y_continuous(breaks = seq(0, 1, 0.2), expand = c(0,0)) +
    scale_fill_manual(values = ESV_colors) +
    facet_wrap(~Community, ncol = 1, scales = "free") +
    theme_classic() +
    theme(strip.background = element_rect(fill = NA, color = NA),
          panel.border = element_rect(fill = NA, color = 1, linewidth = 1)) +
    guides(color = "none", fill = guide_legend(ncol = 1)) +
    labs(x = "", y = "replicate")
ggsave(paste0(folder_data, "temp/14a-07-four_communities_temporal.png"), p, width = 5, height = 8)

