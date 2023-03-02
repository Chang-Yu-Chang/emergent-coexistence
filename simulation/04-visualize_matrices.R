#' This script reads the sampled modeling matrices and draws heatmaps
#'
#' N is the total number of species
#' R is the total number of resources
#'
#' 1. D matrix: NxN stoicheometic matrix depicting species metabolism
#' 2. c matrix: NxR matrix describing the consumption rate of species i on resource j
#' 3. l matrix: NXR matrix describing the leakiness rate of species i on producing resource j

library(tidyverse)
library(cowplot)
source(here::here("analysis/00-metadata.R"))

#
#input_parameters <- read_csv(here::here("simulation/01-input_parameters.csv"), col_types = cols())
input_monocultures <- read_csv(here::here("simulation/01a-input_monocultures.csv"), col_types = cols())
output_dir <- input_monocultures$output_dir[1]
Dm_S0 <- read_csv(paste0(output_dir, "D_S0_seed1.csv"), skip = 1, col_types = cols()) # D matrix
Dm_S500 <- read_csv(paste0(output_dir, "D_S500_seed1.csv"), skip = 1, col_types = cols()) # D matrix
cm <- read_csv(paste0(output_dir, "c_seed1.csv"), skip = 1, col_types = cols()) # c matrix
lm <- read_csv(paste0(output_dir, "l_seed1.csv"), skip = 1, col_types = cols()) # l matrix
category_colors <- c(sugar = "#ED6A5A", acid = "#03CEA4", fermenter = "#8A89C0", respirator = "#FFCB77")

# Generate family-species and class-resource matching tibble
sa <- input_monocultures$sa[1]
ma <- input_monocultures$ma[1]
sal <- tibble(Family = paste0("F", c(rep(0, sa), rep(1, sa))), Species = paste0("S", 0:(sa * 2 - 1)))
mal <- tibble(Class = paste0("T", c(rep(0, ma), rep(1, ma))), Resource = paste0("R", 0:(ma * 2 - 1)))

# D matrix
Dm_S0l <- Dm_S0 %>% # D matrix longer
    pivot_longer(cols = starts_with("R"), names_to = "Resource1", values_to = "SecretionFlux") %>%
    rename(Class2 = ...1, Resource2 = ...2) %>%
    left_join(rename_with(mal, ~ paste0(., "1"), everything()), by = join_by(Resource1)) %>%
    mutate(Resource1 = ordered(Resource1, mal$Resource), Resource2 = ordered(Resource2, rev(mal$Resource))) %>%
    select(Class1, Resource1, Class2, Resource2, SecretionFlux)

Dm_S500l <- Dm_S500 %>% # D matrix longer
    pivot_longer(cols = starts_with("R"), names_to = "Resource1", values_to = "SecretionFlux") %>%
    rename(Class2 = ...1, Resource2 = ...2) %>%
    left_join(rename_with(mal, ~ paste0(., "1"), everything()), by = join_by(Resource1)) %>%
    mutate(Resource1 = ordered(Resource1, mal$Resource), Resource2 = ordered(Resource2, rev(mal$Resource))) %>%
    select(Class1, Resource1, Class2, Resource2, SecretionFlux)

p1_1 <- Dm_S0l %>%
    ggplot() +
    geom_tile(aes(x = Resource1, y = Resource2, fill = SecretionFlux)) +
    # Color bar
    geom_segment(aes(color = "sugar"), x = "R0", xend = "R19", y = "R0", yend = "R0", lwd = 2) +
    geom_segment(aes(color = "acid"), x = "R20", xend = "R39", y = "R0", yend = "R0", lwd = 2) +
    geom_segment(aes(color = "sugar"), x = "R0", xend = "R0", y = "R0", yend = "R19", lwd = 2) +
    geom_segment(aes(color = "acid"), x = "R0", xend = "R0", y = "R20", yend = "R39", lwd = 2) +
    # Axis label
    annotate("text", x = "R10", y = 41, label = "sugar", hjust = 0.5, color = category_colors["sugar"], fontface = "bold") +
    annotate("text", x = "R30", y = 41, label = "acid", hjust = 0.5, color = category_colors["acid"], fontface = "bold") +
    annotate("text", x = -0.1, y = "R10", label = "sugar", angle = 90, vjust = 0, color = category_colors["sugar"], fontface = "bold") +
    annotate("text", x = -0.1, y = "R30", label = "acid", angle = 90, vjust = 0, color = category_colors["acid"], fontface = "bold") +
    scale_fill_gradient(high = "black", low = "white") +
    scale_color_manual(values = category_colors, breaks = c("sugar", "acid")) +
    scale_x_discrete(limits = mal$Resource, position = "top") +
    #scale_y_discrete(limits = mal$Resource) +
    coord_cartesian(xlim = c(0, 40), ylim = c(0, 40), clip = "off") +
    theme_classic() +
    theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank()) +
    guides(fill = guide_colorbar(title = expression(D[i][beta][alpha])),
           color = "none")+
    labs(x = expression(R[alpha]), y = expression(R[beta])) +
    ggtitle("fermenter")

p1_2 <- Dm_S500l %>%
    ggplot() +
    geom_tile(aes(x = Resource1, y = Resource2, fill = SecretionFlux)) +
    # Color bar
    geom_segment(aes(color = "sugar"), x = "R0", xend = "R19", y = "R0", yend = "R0", lwd = 2) +
    geom_segment(aes(color = "acid"), x = "R20", xend = "R39", y = "R0", yend = "R0", lwd = 2) +
    geom_segment(aes(color = "sugar"), x = "R0", xend = "R0", y = "R0", yend = "R19", lwd = 2) +
    geom_segment(aes(color = "acid"), x = "R0", xend = "R0", y = "R20", yend = "R39", lwd = 2) +
    # Axis label
    annotate("text", x = "R10", y = 41, label = "sugar", hjust = 0.5, color = category_colors["sugar"], fontface = "bold") +
    annotate("text", x = "R30", y = 41, label = "acid", hjust = 0.5, color = category_colors["acid"], fontface = "bold") +
    annotate("text", x = -0.1, y = "R10", label = "sugar", angle = 90, vjust = 0, color = category_colors["sugar"], fontface = "bold") +
    annotate("text", x = -0.1, y = "R30", label = "acid", angle = 90, vjust = 0, color = category_colors["acid"], fontface = "bold") +
    scale_fill_gradient(high = "black", low = "white") +
    scale_color_manual(values = category_colors, breaks = c("sugar", "acid")) +
    scale_x_discrete(limits = mal$Resource, position = "top") +
    coord_cartesian(xlim = c(0, 40), ylim = c(0, 40), clip = "off") +
    theme_classic() +
    theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank()) +
    guides(fill = guide_colorbar(title = expression(D[i][beta][alpha])), color = "none")+
    labs(x = expression(R[alpha]), y = expression(R[beta])) +
    ggtitle("respirator")

p1 <- plot_grid(p1_1, p1_2, nrow = 1, align = "h")

ggsave(here::here("simulation/plots/matrix1_D.png"), p1, width = 8, height = 6)


# c matrix
cml <- cm %>% # c matrix longer
    pivot_longer(cols = starts_with("R"), names_to = "Resource", values_to = "ConsumptionRate") %>%
    rename(Family = ...1, Species = ...2) %>%
    left_join(mal) %>%
    mutate(Species = ordered(Species, rev(sal$Species)), Resource = ordered(Resource, mal$Resource)) %>%
    select(Family, Species, Class, Resource, ConsumptionRate)

p2 <- cml %>%
    ggplot() +
    geom_tile(aes(x = Resource, y = Species, fill = ConsumptionRate)) +
    # Color bar
    geom_segment(aes(color = "sugar"), x = "R0", xend = "R19", y = "S0", yend = "S0", lwd = 2) +
    geom_segment(aes(color = "acid"), x = "R20", xend = "R39", y = "S0", yend = "S0", lwd = 2) +
    geom_segment(aes(color = "fermenter"), x = "R0", xend = "R0", y = "S0", yend = "S499", lwd = 2) +
    geom_segment(aes(color = "respirator"), x = "R0", xend = "R0", y = "S500", yend = "S999", lwd = 2) +
    # Axis label
    annotate("text", x = "R10", y = 1020, label = "sugar", hjust = 0.5, color = category_colors["sugar"], fontface = "bold") +
    annotate("text", x = "R30", y = 1020, label = "acid", hjust = 0.5, color = category_colors["acid"], fontface = "bold") +
    annotate("text", x = -0.1, y = "S250", label = "fermenter", angle = 90, vjust = 0, color = category_colors["fermenter"], fontface = "bold") +
    annotate("text", x = -0.1, y = "S750", label = "respirator", angle = 90, vjust = 0, color = category_colors["respirator"], fontface = "bold") +
    scale_fill_gradient(high = "black", low = "white") +
    scale_color_manual(values = category_colors) +
    scale_x_discrete(limits = mal$Resource, position = "top") +
    coord_cartesian(xlim = c(0, 40), ylim = c(0, 1000), clip = "off") +
    theme_classic() +
    theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank()) +
    guides(fill = guide_colorbar(title = expression(c[i][alpha])),
           color = "none") +
    labs()

ggsave(here::here("simulation/plots/matrix2_c.png"), p2, width = 4, height = 6)


# l matrix
lml <- lm %>% # l matrix longer
    pivot_longer(cols = starts_with("R"), names_to = "Resource", values_to = "Leakiness") %>%
    rename(Family = ...1, Species = ...2) %>%
    left_join(mal) %>%
    mutate(Species = ordered(Species, rev(sal$Species)), Resource = ordered(Resource, mal$Resource)) %>%
    select(Family, Species, Class, Resource, Leakiness)
p3 <- lml %>%
    ggplot() +
    geom_tile(aes(x = Resource, y = Species, fill = Leakiness)) +
    # Color bar
    geom_segment(aes(color = "sugar"), x = "R0", xend = "R19", y = "S0", yend = "S0", lwd = 2) +
    geom_segment(aes(color = "acid"), x = "R20", xend = "R39", y = "S0", yend = "S0", lwd = 2) +
    geom_segment(aes(color = "fermenter"), x = "R0", xend = "R0", y = "S0", yend = "S499", lwd = 2) +
    geom_segment(aes(color = "respirator"), x = "R0", xend = "R0", y = "S500", yend = "S999", lwd = 2) +
    # Axis label
    annotate("text", x = "R10", y = 1020, label = "sugar", hjust = 0.5, color = category_colors["sugar"], fontface = "bold") +
    annotate("text", x = "R30", y = 1020, label = "acid", hjust = 0.5, color = category_colors["acid"], fontface = "bold") +
    annotate("text", x = -0.1, y = "S250", label = "fermenter", angle = 90, vjust = 0, color = category_colors["fermenter"], fontface = "bold") +
    annotate("text", x = -0.1, y = "S750", label = "respirator", angle = 90, vjust = 0, color = category_colors["respirator"], fontface = "bold") +
    scale_fill_gradient(high = "black", low = "white") +
    scale_color_manual(values = category_colors) +
    scale_x_discrete(limits = mal$Resource, position = "top") +
    coord_cartesian(xlim = c(0, 40), ylim = c(0, 1000), clip = "off") +
    theme_classic() +
    theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank()) +
    guides(fill = guide_colorbar(title = expression(l[i][alpha])),
           color = "none") +
    labs()
ggsave(here::here("simulation/plots/matrix3_l.png"), p3, width = 4, height = 6)
























