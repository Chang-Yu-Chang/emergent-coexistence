#' This script reads the sampled modeling matrices and draws heatmaps
#'
#' N is the total number of species
#' R is the total number of resources
#'
#' 1. c matrix: NxR matrix describing the consumption rate of species i on resource j
#' 2. D matrix: NxN stoicheometic matrix depicting species metabolism
#' 3. l matrix: NXR matrix describing the leakiness rate of species i on producing resource j

library(tidyverse)
library(cowplot)
source(here::here("analysis/00-metadata.R"))

input_parameters <- read_csv(here::here("simulation/01-input_parameters.csv"), col_types = cols())
input_monocultures <- read_csv(here::here("simulation/02a-input_monocultures.csv"), col_types = cols())
output_dir <- paste0(folder_simulation, "02a-monocultures/")
Dm_F0 <- read_csv(paste0(output_dir, "00-D_S0.csv"), skip = 1, col_types = cols()) # D matrix
Dm_F1 <- read_csv(paste0(output_dir, "00-D_S500.csv"), skip = 1, col_types = cols()) # D matrix
cm <- read_csv(paste0(output_dir, "00-c.csv"), skip = 1, col_types = cols()) # c matrix
lm <- read_csv(paste0(output_dir, "00-l.csv"), skip = 1, col_types = cols()) # l matrix

# Generate family-species and class-resource matching tibble
sa <- input_monocultures$sa[1]
ma <- input_monocultures$ma[1]
sal <- tibble(Family = paste0("F", c(rep(0, sa), rep(1, sa))), Species = paste0("S", 0:(sa * 2 - 1))) %>% mutate(SpeciesID = 1:n())
mal <- tibble(Class = paste0("T", c(rep(0, ma), rep(1, ma))), Resource = paste0("R", 0:(ma * 2 - 1))) %>% mutate(ResourceID = 1:n())

# 1. c matrix ----
cml <- cm %>% # c matrix longer
    pivot_longer(cols = starts_with("R"), names_to = "Resource", values_to = "ConsumptionRate") %>%
    rename(Family = ...1, Species = ...2) %>%
    left_join(mal) %>%
    left_join(sal) %>%
    select(Family, Species, SpeciesID, Class, Resource, ResourceID, ConsumptionRate)


p2 <- cml %>%
    ggplot() +
    geom_tile(aes(x = ResourceID, y = SpeciesID, fill = ConsumptionRate)) +
    # Color bar
    geom_segment(aes(color = "sugar"), x = 0.5, xend = ma+0.5, y = -10, yend = -10, linewidth = 2) +
    geom_segment(aes(color = "acid"), x = ma+0.5, xend = ma*2+0.5, y = -10, yend = -10, linewidth = 2) +
    geom_segment(aes(color = "fermenter"), x = 0, xend = 0, y = 0.5, yend = sa+0.5, linewidth = 2) +
    geom_segment(aes(color = "respirator"), x = 0, xend = 0, y = sa+0.5, yend = sa*2+0.5, linewidth = 2) +
    # Axis label
    annotate("text", x = ma*0.5, y = -.5, label = "sugar", hjust = 0.5, vjust = -1, color = category_color["sugar"], fontface = "bold") +
    annotate("text", x = ma*1.5, y = -.5, label = "acid", hjust = 0.5, vjust = -1, color = category_color["acid"], fontface = "bold") +
    annotate("text", x = -.5, y = sa*0.5, label = "fermenter", angle = 90, vjust = -1, color = category_color["fermenter"], fontface = "bold") +
    annotate("text", x = -.5, y = sa*1.5, label = "respirator", angle = 90, vjust = -1, color = category_color["respirator"], fontface = "bold") +
    scale_fill_gradient(high = "black", low = "white") +
    scale_color_manual(values = category_color) +
    scale_x_continuous(position = "top", breaks = 1:(ma*2)) +
    coord_cartesian(xlim = c(1, ma*2), ylim = c(sa*2, 1), clip = "off") +
    theme_classic() +
    theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank()) +
    guides(fill = guide_colorbar(title = expression(c[i][alpha])), color = "none") +
    labs(x = "", y = "") +
    ggtitle("Uptake rate (c)")

ggsave(here::here("simulation/plots/21-matrix2_c.png"), p2, width = 4, height = 6)


# 2. D matrix ----
Dm_F0l <- Dm_F0 %>% # D matrix longer
    pivot_longer(cols = starts_with("R"), names_to = "Resource1", values_to = "SecretionFlux") %>%
    rename(Class2 = ...1, Resource2 = ...2) %>%
    select(-Class2) %>%
    left_join(rename_with(mal, ~ paste0(., "1"), everything()), by = join_by(Resource1)) %>%
    left_join(rename_with(mal, ~ paste0(., "2"), everything()), by = join_by(Resource2)) %>%
    select(Class1, Resource1, ResourceID1, Class2, Resource2, ResourceID2, SecretionFlux)

Dm_F1l <- Dm_F1 %>% # D matrix longer
    pivot_longer(cols = starts_with("R"), names_to = "Resource1", values_to = "SecretionFlux") %>%
    rename(Class2 = ...1, Resource2 = ...2) %>%
    select(-Class2) %>%
    left_join(rename_with(mal, ~ paste0(., "1"), everything()), by = join_by(Resource1)) %>%
    left_join(rename_with(mal, ~ paste0(., "2"), everything()), by = join_by(Resource2)) %>%
    select(Class1, Resource1, ResourceID1, Class2, Resource2, ResourceID2, SecretionFlux)

p1_1 <- Dm_F0l %>%
    ggplot() +
    geom_tile(aes(x = ResourceID1, y = ResourceID2, fill = SecretionFlux)) +
    # Color bar
    geom_segment(aes(color = "sugar"), x = 0.5, xend = ma+0.5, y = 0, yend = 0, linewidth = 2) +
    geom_segment(aes(color = "acid"), x = ma+0.5, xend = ma*2+0.5, y = 0, yend = 0, linewidth = 2) +
    geom_segment(aes(color = "sugar"), x = 0, xend = 0, y = 0.5, yend = ma+0.5, linewidth = 2) +
    geom_segment(aes(color = "acid"), x = 0, xend = 0, y = ma+0.5, yend = ma*2+0.5, linewidth = 2) +
    # Axis label
    annotate("text", x = ma*0.5, y = -.5, label = "sugar", hjust = 0.5, color = category_color["sugar"], fontface = "bold") +
    annotate("text", x = ma*1.5, y = -.5, label = "acid", hjust = 0.5, color = category_color["acid"], fontface = "bold") +
    annotate("text", x = -.5, y = ma*0.5, label = "sugar", angle = 90, vjust = 0, color = category_color["sugar"], fontface = "bold") +
    annotate("text", x = -.5, y = ma*1.5, label = "acid", angle = 90, vjust = 0, color = category_color["acid"], fontface = "bold") +
    scale_fill_gradient(high = "black", low = "white") +
    scale_color_manual(values = category_color, breaks = c("sugar", "acid")) +
    scale_x_continuous(position = "top", breaks = 1:(ma*2)) +
    coord_cartesian(xlim = c(1, ma*2), ylim = c(ma*2, 1), clip = "off") +
    theme_classic() +
    theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank()) +
    guides(fill = guide_colorbar(title = expression(D[i][beta][alpha])),
           color = "none")+
    labs(x = expression(R[alpha]), y = expression(R[beta])) +
    ggtitle("Metabolic matrix (D) of a fermenter")

p1_2 <- Dm_F1l %>%
    ggplot() +
    geom_tile(aes(x = ResourceID1, y = ResourceID2, fill = SecretionFlux)) +
    # Color bar
    geom_segment(aes(color = "sugar"), x = 0.5, xend = ma+0.5, y = 0, yend = 0, linewidth = 2) +
    geom_segment(aes(color = "acid"), x = ma+0.5, xend = ma*2+0.5, y = 0, yend = 0, linewidth = 2) +
    geom_segment(aes(color = "sugar"), x = 0, xend = 0, y = 0.5, yend = ma+0.5, linewidth = 2) +
    geom_segment(aes(color = "acid"), x = 0, xend = 0, y = ma+0.5, yend = ma*2+0.5, linewidth = 2) +
    # Axis label
    annotate("text", x = ma*0.5, y = -.5, label = "sugar", hjust = 0.5, color = category_color["sugar"], fontface = "bold") +
    annotate("text", x = ma*1.5, y = -.5, label = "acid", hjust = 0.5, color = category_color["acid"], fontface = "bold") +
    annotate("text", x = -.5, y = ma*0.5, label = "sugar", angle = 90, vjust = 0, color = category_color["sugar"], fontface = "bold") +
    annotate("text", x = -.5, y = ma*1.5, label = "acid", angle = 90, vjust = 0, color = category_color["acid"], fontface = "bold") +
    scale_fill_gradient(high = "black", low = "white") +
    scale_color_manual(values = category_color, breaks = c("sugar", "acid")) +
    scale_x_continuous(position = "top", breaks = 1:(ma*2)) +
    coord_cartesian(xlim = c(1, ma*2), ylim = c(ma*2, 1), clip = "off") +
    theme_classic() +
    theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank()) +
    guides(fill = guide_colorbar(title = expression(D[i][beta][alpha])), color = "none")+
    labs(x = expression(R[alpha]), y = expression(R[beta])) +
    ggtitle("Metabolic matrix (D) of a respirator")

p1 <- plot_grid(p1_1, p1_2, nrow = 1, align = "h")

ggsave(here::here("simulation/plots/21-matrix1_D.png"), p1, width = 8, height = 6)


# 3. l matrix ----
lml <- lm %>% # l matrix longer
    pivot_longer(cols = starts_with("R"), names_to = "Resource", values_to = "Leakiness") %>%
    rename(Family = ...1, Species = ...2) %>%
    left_join(mal) %>%
    left_join(sal) %>%
    select(Family, Species, SpeciesID, Class, Resource, ResourceID, Leakiness)
p3 <- lml %>%
    ggplot() +
    geom_tile(aes(x = ResourceID, y = SpeciesID, fill = Leakiness)) +
    # Color bar
    geom_segment(aes(color = "sugar"), x = 0.5, xend = ma+0.5, y = -10, yend = -10, linewidth = 2) +
    geom_segment(aes(color = "acid"), x = ma+0.5, xend = ma*2+0.5, y = -10, yend = -10, linewidth = 2) +
    geom_segment(aes(color = "fermenter"), x = 0, xend = 0, y = 0.5, yend = sa+0.5, linewidth = 2) +
    geom_segment(aes(color = "respirator"), x = 0, xend = 0, y = sa+0.5, yend = sa*2+0.5, linewidth = 2) +
    # Axis label
    annotate("text", x = ma*0.5, y = -.5, label = "sugar", hjust = 0.5, vjust = -1, color = category_color["sugar"], fontface = "bold") +
    annotate("text", x = ma*1.5, y = -.5, label = "acid", hjust = 0.5, vjust = -1, color = category_color["acid"], fontface = "bold") +
    annotate("text", x = -.5, y = sa*0.5, label = "fermenter", angle = 90, vjust = -1, color = category_color["fermenter"], fontface = "bold") +
    annotate("text", x = -.5, y = sa*1.5, label = "respirator", angle = 90, vjust = -1, color = category_color["respirator"], fontface = "bold") +
    scale_fill_gradient(high = "black", low = "white") +
    scale_color_manual(values = category_color) +
    scale_x_continuous(position = "top", breaks = 1:(ma*2)) +
    coord_cartesian(xlim = c(1, ma*2), ylim = c(sa*2, 1), clip = "off") +
    theme_classic() +
    theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank()) +
    guides(fill = guide_colorbar(title = expression(l[i][alpha])), color = "none") +
    labs(x = "", y = "") +
    ggtitle("Leakiness (l)")

ggsave(here::here("simulation/plots/21-matrix3_l.png"), p3, width = 4, height = 6)
























