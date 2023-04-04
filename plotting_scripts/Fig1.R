library(tidyverse)
library(cowplot)
library(broom)
library(grid) # For drawing polygon
source(here::here("analysis/00-metadata.R"))


communities <- read_csv(paste0(folder_data, "output/communities_remained.csv"), show_col_types = F)
isolates <- read_csv(paste0(folder_data, "output/isolates_remained.csv"), show_col_types = F)
pairs <- read_csv(paste0(folder_data, "output/pairs_remained.csv"), show_col_types = F)
communities_abundance <- read_csv(paste0(folder_data, "raw/community_ESV/Emergent_Comunity_Data.csv"), show_col_types = F) %>%
    filter(Carbon_Source == "Glucose" | Carbon_Source == "Original") %>%
    mutate(Community = factor(paste0("C", Inoculum, "R", Replicate), paste0("C", rep(1:12, each = 8), "R", rep(1:8, 12))))%>%
    arrange(Community, Family, Transfer, ESV)
isolates_abundance <- read_csv(paste0(folder_data, "temp/32-isolates_abundance.csv"), col_types = cols())

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
            #str_detect(ESV_ID, "Salmonella") ~ str_replace(ESV_ID, "\\.\\d+", ""),
            str_detect(ESV_ID, "Citrobacter") ~ str_replace(ESV_ID, "\\.\\d+", ""),
            str_detect(ESV_ID, "Enterobacter") ~ str_replace(ESV_ID, "\\.\\d+", ""),
            str_detect(ESV_ID, "Klebsiella") ~ str_replace(ESV_ID, "\\.\\d+", ""),
            #str_detect(ESV_ID, "Acinetobacter") ~ str_replace(ESV_ID, ".\\d+", ""),
            #str_detect(ESV_ID, "Pseudomonadaceae") ~ str_replace(ESV_ID, "Pseudomonadaceae", "Pseudomonadas"),
            str_detect(ESV_ID, "Pseudomonadaceae") ~ str_replace(ESV_ID, "\\.\\d+", ""),
            # str_detect(ESV_ID, "Pseudomonas.1\\d+") ~ str_replace(ESV_ID, "\\.1\\d+", "\\.1"),
            # str_detect(ESV_ID, "Pseudomonas.2\\d+") ~ str_replace(ESV_ID, "\\.2\\d+", "\\.2"),
            # str_detect(ESV_ID, "Pseudomonas.3\\d+") ~ str_replace(ESV_ID, "\\.3\\d+", "\\.3"),
            # str_detect(ESV_ID, "Pseudomonas.4\\d+") ~ str_replace(ESV_ID, "\\.4\\d+", "\\.4"),
            # str_detect(ESV_ID, "Pseudomonas.5\\d+") ~ str_replace(ESV_ID, "\\.5\\d+", "\\.5"),
            str_detect(ESV_ID, "Azomonas") ~ str_replace(ESV_ID, ".\\d+", ""),
            #str_detect(ESV_ID, "Stenotrophomonas") ~ str_replace(ESV_ID, ".\\d+", ""),
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


comm = "C8R4"

# Shared color palette for inset 1 and 2 ----
communities_abundance_comm <- communities_abundance %>%
    filter(Community == comm) %>%
    bind_rows(filter(communities_abundance, Inoculum == str_sub(comm, 2, 2), Transfer == 0)) %>%
    select(Transfer, Family, ESV_ID, Relative_Abundance) %>%
    arrange(Transfer, Family, ESV_ID)
isolates_abundance_factor <- isolates_abundance %>%
    distinct(Community, CommunityESVID, .keep_all = T) %>%
    drop_na() %>%
    rename(ESV_ID = CommunityESVID)
ESV_comm <- bind_rows(
    communities_abundance_comm %>% filter(Relative_Abundance > 0.02) %>% distinct(Family, ESV_ID),
    distinct(isolates_abundance_factor, Family, ESV_ID)
) %>%
    mutate(Family = factor(Family, c("Enterobacteriaceae", "Pseudomonadaceae", "Comamonadaceae", "Aeromonadaceae", "Moraxellaceae", "Alcaligenaceae"))) %>%
    arrange(Family, ESV_ID) %>%
    drop_na() %>%
    distinct(Family, ESV_ID) %>%
    mutate(ESV_colors = c(rev(RColorBrewer::brewer.pal(8, "YlGnBu")),
                          rev(RColorBrewer::brewer.pal(8, "OrRd")),
                          rev(RColorBrewer::brewer.pal(4, "Purples")),
                          "gold2","forestgreen", "hotpink2", "salmon4")
    )

ESV_colors1 <- ESV_comm$ESV_colors %>% setNames(ESV_comm$ESV_ID)
ESV_colors1 <- c(ESV_colors1, Other = "snow", `Not isolated` = "#999998")


# Panel B Inset 1: temporal dynamics of C8R4 ----
pB1 <- communities_abundance_comm %>%
    mutate(ESV_ID = ifelse(ESV_ID %in% names(ESV_colors1), ESV_ID, "Other")) %>%
    mutate(ESV_ID = factor(ESV_ID, names(ESV_colors1))) %>%
    #drop_na() %>%
    ggplot() +
    geom_col(aes(x = Transfer, y = Relative_Abundance, fill = ESV_ID, color = ESV_ID), linewidth = .3, position = position_stack(reverse = T)) +
    annotate("text", x = 0, y = -0.1, label = "inoculum", size = 4, angle = 20, hjust = 0.9, vjust = -0.5) +
    scale_fill_manual(values = ESV_colors1, breaks = names(ESV_colors1), drop = F) +
    scale_color_manual(values = ESV_colors1) +
    scale_x_continuous(breaks = 0:12, expand = c(0,.3), labels = c("", 1:12)) +
    scale_y_continuous(breaks = seq(0,1,0.25), expand = c(0, 0)) +
    coord_cartesian(ylim = c(0, 1), clip = "off") +
    theme_classic() +
    theme(
        panel.grid = element_blank(),
        panel.border = element_rect(color = 1, linewidth = .5, fill = NA),
        panel.background = element_rect(color = 1, fill = NA),
        axis.title = element_text(size = 10),
        axis.text = element_text(color = 1, size = 10),
        plot.background = element_rect(color = NA, fill = NA),
        legend.background = element_rect(color = NA, fill = NA),
        legend.key.width = unit(5, "mm"),
        legend.key.height = unit(2, "mm"),
        legend.text = element_text(size = 10),
        legend.spacing.y = unit(5, "mm"),
        legend.title = element_blank(),
    ) +
    guides(alpha = "none", color = "none", fill = guide_legend(ncol = 1, byrow = F)) +
    labs(x = "transfer", y = "relative abundance")

# Panel B Inset 2: isolate abundance in community ----
isolates_abundance_other <- isolates_abundance_factor %>%
    group_by(Community) %>%
    summarize(TotalAbundance = sum(RelativeAbundance)) %>%
    mutate(RelativeAbundance = 1-TotalAbundance) %>%
    mutate(ESV_ID = "Not isolated")


pB2 <- isolates_abundance_factor %>%
    bind_rows(isolates_abundance_other) %>%
    mutate(ESV_ID = ifelse(ESV_ID %in% names(ESV_colors1), ESV_ID, "Other")) %>%
    mutate(ESV_ID = factor(ESV_ID, names(ESV_colors1))) %>%
    left_join(communities, by = "Community") %>%
    ggplot() +
    geom_col(aes(x = CommunityLabel, y = RelativeAbundance, fill = ESV_ID, group = ESV_ID), position = position_stack(reverse = T)) +
    theme_bw() +
    scale_fill_manual(values = ESV_colors1, breaks = names(ESV_colors1)) +
    #scale_fill_manual(values = c(RColorBrewer::brewer.pal(8, "Set1"), c(RColorBrewer::brewer.pal(8, "Set2")))) +
    scale_x_continuous(breaks = 1:13, expand = c(0,.3)) +
    scale_y_continuous(breaks = seq(0,1,0.25), expand = c(0,0), limits = c(0, 1)) +
    coord_cartesian(ylim = c(0, 1), clip = "off") +
    theme(
        panel.grid = element_blank(),
        axis.title = element_text(size = 10),
        axis.text = element_text(color = 1, size = 10),
        panel.border = element_rect(color = 1, linewidth = .5),
        plot.background = element_blank(),
        panel.background = element_blank()
    ) +
    guides(alpha = "none") +
    labs(x = "community", y = "relative abundance")


# Panel C ----
communities_abundance_temporal <- communities_abundance %>%
    filter(Transfer != 0) %>%
    # Filter for those that has temporal data
    filter(Inoculum %in% c(2,6) | Replicate == 4) %>%
    select(Community, Transfer, ESV_ID, Relative_Abundance) %>%
    arrange(Community, Transfer, ESV_ID)

# A complete list of ESV-Transfer
communities_abundance_temporal_complete <- communities_abundance_temporal %>%
    distinct(Community, ESV_ID) %>%
    slice(rep(1:n(), each = 12)) %>%
    mutate(Transfer = rep(1:12, n()/12))
nrow(distinct(communities_abundance_temporal_complete, Community, ESV_ID)) # 755 unique ESVs
nrow(communities_abundance_temporal_complete) # 755 ESVs * 12 transfers = 9060 rows

# Calculate fitness
communities_abundance_fitness <- communities_abundance_temporal %>%
    right_join(communities_abundance_temporal_complete) %>%
    group_by(Community, ESV_ID) %>%
    arrange(Community, ESV_ID, Transfer) %>%
    mutate(Fitness = log(lead(Relative_Abundance) / Relative_Abundance))
ESV_stable <- communities_abundance_temporal %>%
    filter(Transfer %in% c(9:12)) %>%
    pivot_wider(id_cols = c(Community, ESV_ID), names_from = Transfer, names_prefix = "T", values_from = Relative_Abundance) %>%
    filter(
        # Species with complete presence at T8-12
        (!is.na(T9) & !is.na(T10) & !is.na(T11) & !is.na(T12) ) |
            # Community C4R4, C9R4 have missing time points at T11. For these communities, add back the species
            (Community %in% c("C4R4", "C9R4") & (!is.na(T9) & !is.na(T10) & !is.na(T12))) |
            # Community C5R4 have missing time points at T9. For these communities, add back the species
            (Community %in% c("C5R4") & (!is.na(T10) & !is.na(T11) & !is.na(T12)))
    ) %>%
    distinct(Community, ESV_ID) %>%
    mutate(CommunityESV = paste0(Community, ESV_ID))

communities_abundance_fitness_stable <- communities_abundance_fitness %>%
    drop_na() %>%
    mutate(CommunityESV = paste0(Community, ESV_ID)) %>%
    filter(CommunityESV %in% ESV_stable$CommunityESV)

# Correlation
cor_test_long <- function (long_data) cor.test(long_data$Relative_Abundance, long_data$Fitness, method = "spearman", exact = F, alternative = c("less"))
tb_cor_stable <- communities_abundance_fitness_stable %>%
    nest(data = c(-Community, -ESV_ID)) %>%
    mutate(fit = map(data, ~ cor_test_long(.x)),
           tidied = map(fit, tidy)) %>%
    unnest(tidied)
ESV_sig_stable <- tb_cor_stable %>%
    select(Community, ESV_ID, estimate, p.value) %>%
    filter(p.value < 0.05, estimate < 0) %>%
    mutate(CommunityESV = paste0(Community, ESV_ID)) %>%
    mutate(Significance = "p<0.05")

nrow(tb_cor_stable) # 99 stable ESVs
nrow(ESV_sig_stable) # 52 shows significant negative correlation
range(ESV_sig_stable$estimate) # range of correlation [-0.95, -0.53]
tb_cor_stable %>% filter(Community %in% c("C4R4", "C5R4", "C9R4")) %>% nrow() # 9 ESVs from C4R4, C5R4, C9R4 have only three time points

# Empirical equilibrium frequency
eq_freq_stable <- communities_abundance_temporal %>%
    mutate(CommunityESV = paste0(Community, ESV_ID)) %>%
    filter(Transfer %in% 9:12) %>%
    filter(CommunityESV %in% ESV_stable$CommunityESV) %>%
    arrange(Community, ESV_ID, CommunityESV, Transfer) %>%  # 387 rows. 99 * 4 - 9 = 387. 9 ESVs from C4R4, C5R4, C9R4 have only three time points
    group_by(Community, ESV_ID, CommunityESV) %>%
    summarize(EquilibriumAbundance = mean(Relative_Abundance)) # 99 ESVs

# Linear model predicted equilibrium frequency
xintercept_stable <- communities_abundance_fitness_stable %>%
    nest(data = c(-Community, -ESV_ID, -CommunityESV)) %>%
    mutate(fit = map(data, ~ lm(Fitness ~ Relative_Abundance, data = .x)),
           tidied = map(fit, tidy)) %>%
    unnest(tidied) %>%
    select(Community, ESV_ID, CommunityESV, term, estimate) %>%
    pivot_wider(names_from = term, values_from = estimate) %>%
    # X intercept
    mutate(Xintercept = -`(Intercept)` / Relative_Abundance) %>%
    select(Community, ESV_ID, CommunityESV, Xintercept, Slope = Relative_Abundance, Yintercept = `(Intercept)`) %>%
    mutate(ESVType = "stable") %>%
    # Mark significant correlation
    left_join(select(tb_cor_stable, Community, ESV_ID, estimate, p.value))


eq_freq_stable_both <- eq_freq_stable %>%
    left_join(xintercept_stable) %>%
    drop_na(Xintercept)

eq_freq_stable_comm <- eq_freq_stable_both %>%
    mutate(InThisStudy = case_when(
        Community %in% communities$Community ~ "four communities\nin current study",
        T ~ "other 22 communities"
    )) %>%
    arrange(desc(InThisStudy)) # 99 ESVs

eq_freq_stable_comm_filtered <- eq_freq_stable_comm %>%
    # Negative slope
    filter(Slope < 0)  # 95 ESVs

eq_freq_stable_comm_filtered %>%
    group_by(InThisStudy, Community) %>%
    count() %>%
    pull(InThisStudy) %>%
    table() # 4 communities in current study, 22 other communities


pC <- eq_freq_stable_comm_filtered %>%
    ggplot() +
    geom_abline(intercept = 0, slope = 1, linetype = 2, color = "black") +
    geom_point(aes(x = EquilibriumAbundance, y = Xintercept, color = InThisStudy), shape = 21, size = 2, stroke = 1) +
    scale_color_manual(values = c("four communities\nin current study" = "red", "other 22 communities" = grey(0.8))) +
    scale_x_continuous(limits = c(0,1), breaks = seq(0,1,0.25)) +
    scale_y_continuous(limits = c(0,1), breaks = seq(0,1,0.25)) +
    theme_classic() +
    theme(
        panel.border = element_rect(color = 1, fill = NA, linewidth = 1),
        axis.text = element_text(color = 1),
        axis.title = element_text(size = 10),
        #legend.position = c(0.7, 0.15),
        legend.position = "top",
        legend.title = element_blank(),
        legend.margin = margin(0,3,1,1, unit = "mm"),
        legend.box.margin = margin(0,-1,-3,0, unit = "mm"),
        legend.text = element_text(size = 8),
        legend.key.size = unit(3, "mm"),
        legend.spacing.x = unit(1, "mm"),
        legend.spacing.y = unit(0, "mm"),
        legend.background = element_rect(color = NA, fill = NA),
        plot.background = element_blank()
    ) +
    guides(color = guide_legend(nrow = 1)) +
    labs(x = "empirical equilibrium abundance", y = "equilibrium abundance\npredicted from\nassembly dynamics")

cor.test(eq_freq_stable_comm$EquilibriumAbundance, eq_freq_stable_comm$Xintercept, method = "pearson") %>%
    tidy()


# Panel D: Negative frequency dependent selection of 2 ESVs from C8R4 ----
pD <- communities_abundance_fitness_stable %>%
    filter(Community == "C8R4") %>%
    mutate(ESV_ID = factor(ESV_ID, c("Pseudomonas.10", "Klebsiella"))) %>%
    ggplot() +
    geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, fill = ESV_ID), alpha = 0.1) +
    geom_smooth(aes(x = Relative_Abundance, y = Fitness), method = stats::lm, formula = y ~ x, se = T, color = "gold2") +
    geom_point(aes(x = Relative_Abundance, y = Fitness), shape = 21, stroke = 1, size = 2) +
    geom_hline(yintercept = 0, linetype = 2) +
    scale_fill_manual(values = ESV_colors1) +
    scale_x_continuous(limits = c(0,1), breaks = seq(0,1,0.25)) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 3)) +
    facet_wrap(~ESV_ID, ncol = 1) +
    theme_classic() +
    theme(
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 10, color = 1),
        strip.text = element_text(size = 10),
        strip.background = element_blank(),
        panel.border = element_rect(color = 1, linewidth = 1, fill = NA),
        plot.background = element_blank()
    ) +
    guides(color = "none", fill = "none") +
    labs(x = "relative abundance", y = "invasion fitness")


# Assemble panels ----
p_legend_ESV <- {pB1 +
        theme(
            legend.text = element_text(size = 10),
            legend.title = element_blank(),
            legend.spacing.y = unit(2, "mm"),
            #legend.title = element_text(size = 12),
            legend.key.size = unit(5, "mm")) +
        guides(fill = guide_legend(ncol = 2, override.aes = list(color = 1)))
} %>% get_legend()

make_polygon <- function (x1, x2, x3, x4, y1, y2, y3, y4) {
    polygonGrob(x = c(x1, x2, x3, x4),
                y = c(y1, y2, y3, y4),
                gp = gpar(fill = "grey", alpha = 0.3, col = NA))

}
zoom_polygon1 <- make_polygon(0.74, 0.74, 0.825, 0.825,
                              0.606, 0.653, 0.648, 0.558)
zoom_polygon2 <- make_polygon(0.74, 0.74, 0.825, 0.825,
                              0.437, 0.598, 0.53, 0.442)

p <- ggdraw() +
    draw_image(here::here("plots/cartoons/Fig1_cartoon.png")) +
    draw_grob(zoom_polygon1) +
    draw_grob(zoom_polygon2) +
    draw_plot(pB1 + guides(fill = "none"), x = 0.50, y = 0.40, width = 0.25, height = 0.26) +
    draw_plot(pB2 + guides(fill = "none"), x = 0.50, y = 0.08, width = 0.25, height = 0.26) +
    draw_plot(plot_grid(pC, labels = "C", label_size = 20, scale = .95), x = 0.76, y = 0.69, width = 0.20, height = 0.25) +
    draw_plot(plot_grid(pD, labels = "D", label_size = 20, scale = .9), x = 0.79, y = 0.39, width = 0.18, height = 0.3) +
    draw_plot(p_legend_ESV, x = 0.85, y = 0.2, width = 0.08, height = 0.1)

ggsave(here::here("plots/Fig1.png"), p, width = 15, height = 12)


# Stats
isolates %>%
    distinct(Community, CommunityESVID, .keep_all = T) %>%
    group_by(Community) %>%
    summarize(Total = sum(RelativeAbundance, na.rm = T)) %>%
    summarize(Mean = mean(Total))








