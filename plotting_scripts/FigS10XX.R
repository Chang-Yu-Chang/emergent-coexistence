library(tidyverse)
library(cowplot)
source(here::here("analysis/00-metadata.R"))

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
    # Set rare taxa into Other
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

communities <- read_csv(paste0(folder_data, "temp/00c-communities.csv"), show_col_types = F) %>%
    mutate(Community = factor(Community, Community))
communities_abundance <- read_csv(paste0(folder_data, "raw/community_ESV/Emergent_Comunity_Data.csv"), show_col_types = F) %>%
    filter(Carbon_Source == "Glucose" | Carbon_Source == "Original") %>%
    mutate(Community = factor(paste0("C", Inoculum, "R", Replicate), paste0("C", rep(1:12, each = 8), "R", rep(1:8, 12))))%>%
    bin_ESV_names() %>%
    clean_ESV_names() %>%
    arrange(Community, Family, Transfer, ESV)

# Subset T1 and T12
communities_abundance_T1 <- communities_abundance %>%
    filter(Transfer == 1)
communities_abundance_T12 <- communities_abundance %>%
    filter(Transfer == 12) %>%
    filter(Relative_Abundance > 0.01)


# Subset species that are present both at T1 and at T12
communities_abundance_sp <- communities_abundance %>%
    right_join(distinct(communities_abundance_T12, Community, ESV_ID)) %>%
    filter(Transfer != 0)


# Plot community abundance over time ----
# Community relative abundance over time
communities_abundance_temporal <- communities_abundance %>%
    distinct(Community, Transfer) %>%
    arrange(Community, Transfer) %>%
    filter(Transfer == 4)

# T0 inoculum
communities_abundance_T0 <- bind_rows(
    communities_abundance %>% filter(Transfer == 0, Inoculum == 1) %>% mutate(Community = "C1R4"),
    communities_abundance %>% filter(Transfer == 0, Inoculum == 2) %>% mutate(Community = "C2R6"),
    communities_abundance %>% filter(Transfer == 0, Inoculum == 2) %>% mutate(Community = "C2R8"),
    communities_abundance %>% filter(Transfer == 0, Inoculum == 8) %>% mutate(Community = "C8R4"),
)

temp <- communities_abundance %>%
    filter(Community %in% communities_abundance_temporal$Community) %>%
    filter(Community %in% communities$Community)
ESV_colors <- temp %>% get_ESV_colors

p1 <- temp %>%
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
    guides(color = "none") +
    labs(x = "transfer", y = "relative abundance")

# Species relative abundance over time ----
plot_species_abundance <- function (species_abundance) {
    species_abundance %>%
        mutate(ESV_ID = ifelse(ESV_ID %in% names(ESV_colors[-length(ESV_colors)]), ESV_ID, "Other")) %>%
        mutate(ESV_ID = factor(ESV_ID, rev(names(ESV_colors)))) %>%
        ggplot() +
        geom_col(aes(x = Transfer, y = Relative_Abundance, fill = ESV_ID), color = "black", linewidth = .3) +
        geom_hline(yintercept = 0.01, linetype = 2) +
        scale_fill_manual(values = ESV_colors) +
        scale_x_continuous(breaks = seq(0,12,2)) +
        facet_wrap(~ESV_ID, nrow = 1, scales = "free_y") +
        theme_classic() +
        theme(panel.border = element_rect(color = 1, fill = NA),
              strip.background = element_rect(color = NA, fill = NA),
              plot.margin = unit(rep(0.5, 4), "cm")
        ) +
        guides(alpha = "none", fill = "none") +
        labs(x = "transfer", y = "relative abundance")
}

temp_plot_list <- rep(list(NA), 4)
comms <- c("C1R4", "C2R6", "C2R8", "C8R4")
for (i in 1:4) {
    comm <- comms[i]
    temp <- communities_abundance %>%
        filter(Community == comm) %>%
        filter(ESV_ID %in% communities_abundance_T12$ESV_ID[communities_abundance_T12$Community == comm],
               ESV_ID %in% communities_abundance_T1$ESV_ID[communities_abundance_T1$Community == comm])

    # Use the ESV_colors from p1
    temp_plot_list[[i]] <- temp %>%
        group_by(Transfer, ESV_ID) %>%
        summarize(Relative_Abundance = sum(Relative_Abundance)) %>%
        plot_species_abundance()

}


p2 <- plot_grid(
    plot_grid(temp_plot_list[[1]], NULL, rel_widths = c(3.2,1)),
    plot_grid(temp_plot_list[[2]], NULL, rel_widths = c(4,0)),
    plot_grid(temp_plot_list[[3]], NULL, rel_widths = c(4,0)),
    plot_grid(temp_plot_list[[4]], NULL, rel_widths = c(3.2,1)),
    ncol = 1
)

p <- plot_grid(
    plot_grid(p1 + guides(fill = "none"), p2, nrow = 1, rel_widths = c(1, 4), labels = c("A", "B")),
    get_legend(p1 + theme(legend.position = "bottom") +guides(fill = guide_legend(nrow = 3))),
    ncol = 1, rel_heights = c(5, 1)
    ) + paint_white_background()
ggsave(here::here("plots/FigS10-communities_abundance.png"), p, width = 13, height = 10)



