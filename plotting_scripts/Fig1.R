library(tidyverse)
library(cowplot)
library(broom)
library(grid) # For drawing polygon
source(here::here("analysis/00-metadata.R"))


isolates <- read_csv(paste0(folder_data, "output/isolates.csv"), show_col_types = F)
pairs <- read_csv(paste0(folder_data, "output/pairs.csv"), show_col_types = F)
communities <- read_csv(paste0(folder_data, "temp/00c-communities.csv"), show_col_types = F)
communities_abundance <- read_csv(paste0(folder_data, "raw/community_ESV/Emergent_Comunity_Data.csv"), show_col_types = F) %>%
    filter(Carbon_Source == "Glucose" | Carbon_Source == "Original") %>%
    mutate(Community = factor(paste0("C", Inoculum, "R", Replicate), paste0("C", rep(1:12, each = 8), "R", rep(1:8, 12))))%>%
    arrange(Community, Family, Transfer, ESV)
isolates_abundance <- read_csv(paste0(folder_data, "temp/14-isolates_abundance.csv"), col_types = cols())

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
            #str_detect(ESV_ID, "Pseudomonas.1\\d+") ~ str_replace(ESV_ID, "\\.1\\d+", "\\.1"),
            str_detect(ESV_ID, "Pseudomonas.2\\d+") ~ str_replace(ESV_ID, "\\.2\\d+", "\\.2"),
            str_detect(ESV_ID, "Pseudomonas.3\\d+") ~ str_replace(ESV_ID, "\\.3\\d+", "\\.3"),
            str_detect(ESV_ID, "Pseudomonas.4\\d+") ~ str_replace(ESV_ID, "\\.4\\d+", "\\.4"),
            str_detect(ESV_ID, "Pseudomonas.5\\d+") ~ str_replace(ESV_ID, "\\.5\\d+", "\\.5"),
            str_detect(ESV_ID, "Pseudomonas.5\\d+") ~ str_replace(ESV_ID, "\\.5\\d+", "\\.5"),
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


comm = "C2R8"
# Color code shared by both insets I and II
temp <- communities_abundance %>%
    #filter((Inoculum == 8 & Replicate == 4) | (Transfer == 12)) %>%
    filter(Community == comm | (Transfer == 12)) %>%
    bind_rows(filter(communities_abundance, Inoculum == str_sub(comm, 2, 2), Transfer == 0)) %>%
    bind_rows(select(drop_na(isolates_abundance), Family, ESV_ID = CommunityESVID)) %>%
    #distinct(Family, ESV_ID) %>%
    bin_ESV_names() %>%
    clean_ESV_names() %>%
    mutate(Family =  ifelse(Family %in% unique(isolates_abundance$Family), Family, "Others")) %>%
    select(Transfer, Family, ESV_ID, Relative_Abundance) %>%
    arrange(Transfer, Family, ESV_ID) %>%
    distinct(Family, ESV_ID)

ESV_colors <- temp %>% get_ESV_colors()


# Inset 1: ESC abundance over time ----
communities_abundance_comm <- communities_abundance %>%
    filter(Community == comm) %>%
    bind_rows(filter(communities_abundance, Inoculum == str_sub(comm, 2, 2), Transfer == 0)) %>%
    bin_ESV_names() %>%
    clean_ESV_names() %>%
    mutate(Family =  ifelse(Family %in% c("Enterobacteriaceae", "Pseudomonadaceae"), Family, "Others")) %>%
    select(Transfer, Family, ESV_ID, Relative_Abundance) %>%
    arrange(Transfer, Family, ESV_ID)

ESV_colors <- communities_abundance_comm %>% get_ESV_colors()

p1 <- communities_abundance_comm %>%
    mutate(ESV_ID = ifelse(ESV_ID %in% names(ESV_colors[-length(ESV_colors)]), ESV_ID, "Other")) %>%
    mutate(ESV_ID = factor(ESV_ID, names(ESV_colors))) %>%
    ggplot() +
    geom_col(aes(x = Transfer, y = Relative_Abundance, fill = ESV_ID, color = ESV_ID), linewidth = .3, position = position_stack(reverse = T)) +
    annotate("text", x = 0, y = -0.1, label = "inoculum", size = 6, angle = 20, hjust = 0.9, vjust = -0.5) +
    scale_fill_manual(values = ESV_colors) +
    scale_color_manual(values = ESV_colors) +
    scale_x_continuous(breaks = 0:12, expand = c(0,.3), labels = c("", 1:12)) +
    scale_y_continuous(expand = c(0, 0), breaks = c(0, .5, 1)) +
    coord_cartesian(ylim = c(0, 1), clip = "off") +
    theme_classic() +
    theme(panel.grid = element_blank(),
          panel.border = element_rect(color = 1, linewidth = .5, fill = NA),
          panel.background = element_rect(color = 1, fill = "#FBF2E4"),
          axis.title = element_text(size = 17),
          axis.text = element_text(color = 1, size = 17),
          plot.background = element_rect(color = NA, fill = NA),
          legend.background = element_rect(color = NA, fill = NA),
          legend.key.size = unit(10, "mm"),
          legend.text = element_text(size = 10),
          legend.title = element_blank(),
    ) +
    guides(alpha = "none", color = "none", fill = guide_legend(ncol = 2)) +
    labs(x = "transfer", y = "relative abundance")


# Inset 2: isolate abundance in community ----
isolates_abundance_factor <- isolates_abundance %>%
    drop_na() %>%
    rename(ESV_ID = CommunityESVID) %>%
    mutate(ESV_ID = ifelse(ESV_ID %in% names(ESV_colors[-length(ESV_colors)]), ESV_ID, "Other")) %>%
    mutate(ESV_ID = factor(ESV_ID, names(ESV_colors)))

isolates_abundance_other <- isolates_abundance_factor %>%
    group_by(Community) %>%
    summarize(TotalAbundance = sum(RelativeAbundance)) %>%
    mutate(RelativeAbundance = 1-TotalAbundance) %>%
    mutate(ESV_ID = "haha")

ESV_colors <- isolates_abundance_factor %>% get_ESV_colors()

#isolates_RDP <- read_csv(paste0(folder_data, "temp/12-isolates_RDP.csv"), show_col_types = F)

p2 <- isolates_abundance_factor %>%
    #bind_rows(isolates_abundance_other) %>%
    mutate(ESV_ID = factor(ESV_ID, names(ESV_colors))) %>%
    left_join(communities, by = "Community") %>%
    ggplot() +
    geom_col(aes(x = CommunityLabel, y = RelativeAbundance, fill = Genus), position = position_stack(reverse = T)) +
    theme_bw() +
    #scale_fill_manual(values = ESV_colors) +
    scale_fill_manual(values = c(RColorBrewer::brewer.pal(8, "Set1"), c(RColorBrewer::brewer.pal(5, "Set2")))) +
    scale_x_continuous(breaks = 1:13, expand = c(0,.3)) +
    scale_y_continuous(breaks = c(0, .5, 1), expand = c(0,0), limits = c(0, 1)) +
    coord_cartesian(ylim = c(0, 1), clip = "off") +
    theme(panel.grid = element_blank(),
          axis.title = element_text(size = 17),
          axis.text = element_text(color = 1, size = 17),
          panel.border = element_rect(color = 1, linewidth = .5),
          plot.background = element_blank(),
          panel.background = element_blank()) +
    guides(alpha = "none") +
    labs(x = "community", y = "relative abundance")
p2

# Inset 3: Negative frequency dependent selection of three ESVs from C1R4 ----
communities_abundance_temporal <- communities_abundance %>%
    filter(Transfer != 0) %>%
    # Filter for those that has temporal data
    filter(Inoculum %in% c(2,6) | Replicate == 4) %>%
    select(Community, Transfer, ESV_ID, Relative_Abundance) %>%
    arrange(Community, Transfer, ESV_ID)
communities_abundance_temporal_complete <- communities_abundance_temporal %>%
    distinct(Community, ESV_ID) %>%
    slice(rep(1:n(), each = 12)) %>%
    mutate(Transfer = rep(1:12, n()/12))

# Calculate fitness
communities_abundance_fitness <- communities_abundance_temporal %>%
    right_join(communities_abundance_temporal_complete) %>%
    group_by(Community, ESV_ID) %>%
    arrange(Community, ESV_ID, Transfer) %>%
    mutate(Fitness = log(lead(Relative_Abundance) / Relative_Abundance))

communities_abundance_comm <- communities_abundance %>%
    filter(Community == comm, Transfer == 12, Relative_Abundance > 0.02)

communities_abundance_fitness_comm <- communities_abundance_fitness %>%
    filter(Community == comm) %>%
    filter(ESV_ID %in% communities_abundance_comm$ESV_ID) %>%
    mutate(CommunityESV = paste0(Community, ESV_ID))

cor_test_long <- function (long_data) cor.test(long_data$Relative_Abundance, long_data$Fitness, method = "spearman", exact = F, alternative = c("less"))
tb_cor <- communities_abundance_fitness_comm %>%
    nest(data = c(-Community, -ESV_ID)) %>%
    mutate(fit = map(data, ~ cor_test_long(.x)),
           tidied = map(fit, tidy)) %>%
    unnest(tidied)

ESV_sig <- tb_cor %>%
    #filter(term == "Relative_Abundance") %>%
    select(Community, ESV_ID, estimate, p.value) %>%
    filter(p.value < 0.05, estimate < 0) %>%
    #distinct(Community, ESV_ID) %>%
    mutate(CommunityESV = paste0(Community, ESV_ID)) %>%
    mutate(Significance = "p<0.05")

p3 <- communities_abundance_fitness_comm %>%
    mutate(ESV_ID = factor(ESV_ID, c("Aeromonas", "Pseudomonas.2", "Citrobacter", "Raoultella", "Enterobacteriaceae"))) %>%
    ggplot() +
    #geom_smooth(aes(x = Relative_Abundance, y = Fitness), method = "lm", formula = y~x, se = F) +
    geom_smooth(data = communities_abundance_fitness_comm %>% left_join(ESV_sig) %>% replace_na(list(Significance = "p>=0.05")),
                aes(x = Relative_Abundance, y = Fitness, color = Significance),
                method = stats::lm, formula = y ~ x, se = F) +
    geom_point(aes(x = Relative_Abundance, y = Fitness), shape = 21, stroke = 1, size = 3) +
    geom_hline(yintercept = 0, linetype = 2) +
    scale_color_manual(values = c("p<0.05" = "pink", "p>=0.05" = grey(0.8))) +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 3), limits = c(0,1)) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 3)) +
    facet_wrap(~ESV_ID, ncol = 1) +
    theme_classic() +
    theme(
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 15),
        strip.text = element_text(size = 13),
        strip.background = element_blank(),
        panel.border = element_rect(color = 1, linewidth = 0.5, fill = NA),
        plot.background = element_blank()
    ) +
    guides(color = "none") +
    labs(x = "relative abundance", y = "fitness")

communities_abundance_fitness_comm %>%
    nest(data = c(-Community, -ESV_ID)) %>%
    mutate(fit = map(data, ~ lm(Fitness ~ Relative_Abundance, data = .x)),
           tidied = map(fit, tidy),
           r.squared = map(fit, function(x) summary(x)$adj.r.squared)) %>%
    unnest(r.squared) %>%
    unnest(tidied) %>%
    filter(term == "Relative_Abundance")

# Assemble panels
p_legend_ESV <- {p1 +
        theme(
            legend.text = element_text(size = 13),
            legend.title = element_blank(),
            #legend.title = element_text(size = 12),
            legend.key.size = unit(5, "mm")) +
        guides(fill = guide_legend(ncol = 2))
    } %>% get_legend()


make_polygon <- function (x1, x2, x3, x4, y1, y2, y3, y4) {
    polygonGrob(x = c(x1, x2, x3, x4),
                y = c(y1, y2, y3, y4),
                gp = gpar(fill = "grey", alpha = 0.3, col = NA))

}

zoom_polygon1 <- make_polygon(0.76, 0.76, 0.865, 0.865,
                              0.635, 0.65, 0.755, 0.68)
zoom_polygon2 <- make_polygon(0.76, 0.76, 0.87, 0.87,
                              0.58, 0.64, 0.67, 0.59)
zoom_polygon3 <- make_polygon(0.76, 0.76, 0.87, 0.87,
                              0.5, 0.57, 0.58, 0.52)
zoom_polygon4 <- make_polygon(0.76, 0.76, 0.87, 0.87,
                              0.5, 0.57, 0.58, 0.52)
zoom_polygon5 <- make_polygon(0.76, 0.76, 0.87, 0.87,
                              0.5, 0.57, 0.58, 0.52)


p <- ggdraw() +
    draw_image(here::here("plots/cartoons/Fig1_cartoon.png")) +
    draw_grob(zoom_polygon1) +
    draw_grob(zoom_polygon2) +
    draw_grob(zoom_polygon3) +
    draw_grob(zoom_polygon4) +
    draw_grob(zoom_polygon5) +
    draw_plot(p1 + guides(fill = "none"), x = 0.52, y = 0.40, width = 0.25, height = 0.26) +
    draw_plot(p2 + guides(fill = "none"), x = 0.52, y = 0.05, width = 0.25, height = 0.26) +
    draw_plot(plot_grid(p3 + guides(fill = "none"), labels = "C", label_size = 25), x = 0.83, y = 0.28, width = 0.12, height = 0.5) +
    draw_plot(p_legend_ESV, x = 0.85, y = 0.1, width = 0.08, height = 0.1)

ggsave(here::here("plots/Fig1.png"), p, width = 20, height = 15)


# Stats
isolates %>%
    group_by(Community) %>%
    summarize(Total = sum(RelativeAbundance, na.rm = T)) %>%
    summarize(Mean = mean(Total))


isolates_abundance %>%
    filter(Genus != str_replace(CommunityESVID, "\\.\\d+", "")) %>%
    view




