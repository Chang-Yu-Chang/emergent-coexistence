library(tidyverse)
library(cowplot)
source(here::here("analysis/00-metadata.R"))
source(here::here("simulation/01-generate_input.R"))

# 0. parameters ----
input_parameters <- read_csv(here::here("simulation/01-input_parameters.csv"), col_types = cols())
input_monocultures <- read_csv(here::here("simulation/02a-input_monocultures.csv"), col_types = cols())
output_dir <- paste0(folder_simulation, "02a-monocultures/")
Dm <- read_csv(paste0(output_dir, "00-D.csv"), skip = 1, col_types = cols()) # D matrix
cm <- read_csv(paste0(output_dir, "00-c.csv"), skip = 1, col_types = cols()) # c matrix
lm <- read_csv(paste0(output_dir, "00-l.csv"), skip = 1, col_types = cols()) # l matrix

mcrm_family_colors <- RColorBrewer::brewer.pal(10, "Set3") %>% setNames(paste0("F", 0:9))
mcrm_resource_colors <- c(RColorBrewer::brewer.pal(9, "Set1"), RColorBrewer::brewer.pal(3, "Set2")[2]) %>% setNames(paste0("R", 0:9))

family_names <- paste0("family ", 1:10) %>% setNames(paste0("F", 0:9))
resource_names <- LETTERS[1:10] %>% setNames(paste0("R", 0:9))


# 1. c matrix ----
cml <- cm %>% # c matrix longer
    pivot_longer(cols = starts_with("R"), names_to = "Resource", values_to = "ConsumptionRate") %>%
    rename(Family = ...1, Species = ...2) %>%
    left_join(mal) %>%
    left_join(sal) %>%
    select(Family, Species, SpeciesID, Class, Resource, ResourceID, ConsumptionRate)


# Stacked bar plot

p1 <- cml %>%
    #filter(Species %in% paste0("S", c(0:10))) %>%
    filter(Species %in% (sal %>% group_by(Family) %>% slice(1:10) %>% pull(Species))) %>%
    mutate(Species = factor(Species, sal$Species)) %>%
    #mutate(Resource = factor(Resource, rev(names(mcrm_resource_colors)))) %>%
    # summarize(sumCR = sum(ConsumptionRate)) %>%
    ggplot()+
    geom_col(aes(x = Species, y = ConsumptionRate, fill = Resource), width = 1) +
    scale_fill_manual(values = mcrm_resource_colors, label = resource_names) +
    scale_y_continuous(expand = c(0,0), limits = c(0, 1.05), breaks = c(0,0.5, 1)) +
    facet_grid(.~Family, scales = "free_x", labeller = as_labeller(family_names)) +
    theme_classic() +
    theme(panel.border = element_rect(color = 1, fill = NA),
          strip.background = element_rect(color = NA, fill = NA),
          panel.spacing.x = unit(0, "mm"),
          axis.text.x = element_blank()) +
    guides(fill = guide_legend(ncol = 2)) +
    labs(x = "species", y = "uptake rate")

# 2. D matrix
Dml <- Dm %>% #
    pivot_longer(cols = starts_with("R"), names_to = "Resource1", values_to = "SecretionFlux") %>%
    rename(Class2 = ...1, Resource2 = ...2) %>%
    select(-Class2) %>%
    left_join(rename_with(mal, ~ paste0(., "1"), everything()), by = join_by(Resource1)) %>%
    left_join(rename_with(mal, ~ paste0(., "2"), everything()), by = join_by(Resource2)) %>%
    select(Class1, Resource1, ResourceID1, Class2, Resource2, ResourceID2, SecretionFlux)

p2 <- Dml %>%
    ggplot() +
    geom_tile(aes(x = ResourceID1, y = ResourceID2, fill = SecretionFlux)) +
    # Color bar
    # geom_segment(aes(color = "sugar"), x = 0.5, xend = ma+0.5, y = 0, yend = 0, linewidth = 2) +
    # geom_segment(aes(color = "acid"), x = ma+0.5, xend = ma*2+0.5, y = 0, yend = 0, linewidth = 2) +
    # geom_segment(aes(color = "sugar"), x = 0, xend = 0, y = 0.5, yend = ma+0.5, linewidth = 2) +
    # geom_segment(aes(color = "acid"), x = 0, xend = 0, y = ma+0.5, yend = ma*2+0.5, linewidth = 2) +
    # # Axis label
    # annotate("text", x = ma*0.5, y = -.5, label = "sugar", hjust = 0.5, color = category_color["sugar"], fontface = "bold") +
    # annotate("text", x = ma*1.5, y = -.5, label = "acid", hjust = 0.5, color = category_color["acid"], fontface = "bold") +
    # annotate("text", x = -.5, y = ma*0.5, label = "sugar", angle = 90, vjust = 0, color = category_color["sugar"], fontface = "bold") +
    # annotate("text", x = -.5, y = ma*1.5, label = "acid", angle = 90, vjust = 0, color = category_color["acid"], fontface = "bold") +
    scale_fill_gradient(high = "black", low = "white") +
    scale_color_manual(values = category_color, breaks = c("sugar", "acid")) +
    scale_x_continuous(position = "top", breaks = 1:(ma*fa)) +
    coord_cartesian(xlim = c(1, ma*fa), ylim = c(ma*fa, 1), clip = "off") +
    theme_classic() +
    theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank()) +
    guides(fill = guide_colorbar(title = expression(D[i][beta][alpha])),
           color = "none") +
    labs(x = expression(R[alpha]), y = expression(R[beta]))

# Assemble panels ----
p <- plot_grid(
    p1,
    #plot_grid(p3 + guides(fill = "none"), p2 + guides(fill = "none"), nrow = 1, rel_widths = c(1, 3), scale = c(.9, .95), labels = c("B", "C")),
    plot_grid(p2, NULL, rel_widths = c(1,1)),
    ncol = 1, rel_heights = c(1,1), labels = c("A", "B"), scale = c(.95, .95)) +
    paint_white_background()
ggsave(here::here("plots/FigS15-mcrm_parameters.png"), p, width = 8, height = 6)










