library(tidyverse)
library(cowplot)
library(broom)
library(grid) # For drawing polygon
source(here::here("processing_scripts/00-metadata.R"))


communities <- read_csv(paste0(folder_data, "output/communities_remained.csv"), show_col_types = F)
isolates <- read_csv(paste0(folder_data, "output/isolates_remained.csv"), show_col_types = F)
communities_abundance <- read_csv(paste0(folder_data, "temp/14-communities_abundance.csv"), show_col_types = F)
eq_freq_stable <- read_csv(paste0(folder_data, "temp/15-eq_freq_stable.csv"), show_col_types = F)
fitness_stable <- read_csv(paste0(folder_data, "temp/15-fitness_stable.csv"), show_col_types = F)

comm = "C8R4"

communities_abundance_comm <- communities_abundance %>%
    filter(Community == comm) %>%
    bind_rows(filter(communities_abundance, Inoculum == str_sub(comm, 2, 2), Transfer == 0)) %>%
    select(Transfer, Family, ESV_ID, Relative_Abundance) %>%
    arrange(Transfer, Family, ESV_ID)
isolates_abundance_factor <- isolates %>%
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
                          rev(RColorBrewer::brewer.pal(5, "Purples")),
                          "gold2","forestgreen", "hotpink2", "salmon4")
    )

ESV_colors1 <- ESV_comm$ESV_colors %>% setNames(ESV_comm$ESV_ID)
ESV_colors1 <- c(ESV_colors1, Other = "snow", `Not isolated` = "#999998")


# Panel B Inset 1: temporal dynamics of C8R4 ----
pB1 <- communities_abundance_comm %>%
    mutate(ESV_ID = ifelse(ESV_ID %in% names(ESV_colors1), ESV_ID, "Other")) %>%
    mutate(ESV_ID = factor(ESV_ID, names(ESV_colors1))) %>%
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
    labs(x = "transfer", y = "frequency")

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
    labs(x = "community", y = "frequency")

# Mean ESV abundance isolated
isolates_abundance <- isolates %>% select(Community, Isolate, CommunityESVID, RelativeAbundance)
isolates_abundance %>%
    distinct(Community, CommunityESVID, .keep_all = T) %>%
    group_by(Community) %>%
    summarize(Total = sum(RelativeAbundance, na.rm = T)) %>%
    summarize(Mean = mean(Total)) # mean abundance is 0.894

# ESV richness
xx <- communities_abundance %>%
    filter(Community %in% communities$Community, Transfer == 12) %>%
    group_by(Community) %>%
    count()
range(xx$n) # ESV abundance is 5-13


# Panel C ----
eq_freq_stable_comm_filtered <- eq_freq_stable %>%
    mutate(InThisStudy = case_when(
        Community %in% communities$Community ~ "four communities\nin current study",
        T ~ "other 22 communities"
    )) %>%
    arrange(desc(InThisStudy)) %>% # 99 ESVs
    # Negative slope
    filter(Slope < 0) # 95 ESVs

pC <- eq_freq_stable_comm_filtered %>%
    ggplot() +
    geom_abline(intercept = 0, slope = 1, linetype = 2, color = "black") +
    geom_point(aes(x = EmpiricalEqAbundance, y = PredictedEqAbundance, color = InThisStudy), shape = 21, size = 2, stroke = 1) +
    scale_color_manual(values = c("four communities\nin current study" = "red", "other 22 communities" = grey(0.8))) +
    scale_x_continuous(limits = c(0,1), breaks = seq(0,1,0.25)) +
    scale_y_continuous(limits = c(0,1), breaks = seq(0,1,0.25)) +
    theme_classic() +
    theme(
        panel.border = element_rect(color = 1, fill = NA, linewidth = 1),
        axis.text = element_text(color = 1),
        axis.title = element_text(size = 10),
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
    labs(x = "average frequency over the last four transfers", y = "predicted from frequency dependent selection (x*)") +
    ggtitle("ESV equilibrium frequency")

#
model <- lm(PredictedEqAbundance ~ EmpiricalEqAbundance, data = eq_freq_stable_comm_filtered)
# R2
summary(model)$r.squared # R2 = 0.91
# RMSD; root mean square deviation
sqrt(mean(model$residuals^2)) # RMSD=0.089


# Panel D: Negative frequency dependent selection of 2 ESVs from C8R4 ----
eq_comm <- eq_freq_stable %>%
    filter(Community == "C8R4") %>%
    mutate(ESV_ID = factor(ESV_ID, c("Pseudomonas.10", "Klebsiella"))) %>%
    mutate(hjust = c(1.1,-0.1))
eq_comm$PredictedEqAbundance[eq_comm$ESV_ID == "Pseudomonas.10"] # predicted eq abundance Pseudo x*=0.173
eq_comm$PredictedEqAbundance[eq_comm$ESV_ID == "Klebsiella"] # predicted eq abundance Klebs x*=0.762

pD <- fitness_stable %>%
    filter(Community == "C8R4") %>%
    mutate(ESV_ID = factor(ESV_ID, c("Pseudomonas.10", "Klebsiella"))) %>%
    ggplot() +
    geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, fill = ESV_ID), alpha = 0.1) +
    geom_smooth(aes(x = Relative_Abundance, y = Fitness), method = stats::lm, formula = y ~ x, se = T, color = "gold2") +
    geom_point(aes(x = Relative_Abundance, y = Fitness), shape = 21, stroke = 1, size = 2) +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_vline(data = eq_comm, aes(xintercept = PredictedEqAbundance), linetype = 2) +
    geom_text(data = eq_comm, aes(x = PredictedEqAbundance, y = Inf, label = paste0("x*=", round(PredictedEqAbundance,3)), hjust = hjust), vjust = 2) +
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
    labs(x = "frequency", y = "invasion fitness")

eq_freq_stable_comm_filtered %>% filter(Community == "C8R4")
# Pseudo: y = -11.8x + 2.05, p<0.001, R^2= 0.9212
# Klebsi: y = -1.81x + 1.38, p<0.001, R^2 = 0.7089
fit <- fitness_stable %>%
    nest(data = c(-Community, -ESV_ID, -CommunityESV)) %>%
    mutate(fit = map(data, ~ lm(Fitness ~ Relative_Abundance, data = .x)),
           tidied = map(fit, tidy)) %>%
    unnest(tidied) %>%
    filter(Community == "C8R4")
summary(fit$fit[3][[1]]) # R^2= 0.9291
summary(fit$fit[1][[1]]) # R^2 = 0.7089


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


# Save vector based
p <- ggdraw() +
    draw_image(here::here("plots/cartoons/Fig1_cartoon.png")) +
    draw_grob(zoom_polygon1) +
    draw_grob(zoom_polygon2) +
    draw_plot(pB1 + guides(fill = "none"), x = 0.50, y = 0.40, width = 0.25, height = 0.26) +
    draw_plot(pB2 + guides(fill = "none"), x = 0.50, y = 0.08, width = 0.25, height = 0.26) +
    draw_plot(plot_grid(pC, labels = "C", label_size = 20, scale = .95), x = 0.76, y = 0.69, width = 0.20, height = 0.25) +
    draw_plot(plot_grid(pD, labels = "D", label_size = 20, scale = .9), x = 0.79, y = 0.39, width = 0.18, height = 0.3) +
    draw_plot(p_legend_ESV, x = 0.85, y = 0.2, width = 0.08, height = 0.1)

ggsave(here::here("plots/Fig1.pdf"), p, width = 15, height = 12)







