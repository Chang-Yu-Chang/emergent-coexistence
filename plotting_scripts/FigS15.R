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

family_names <- paste0("", 1:10) %>% setNames(paste0("F", 0:9))
resource_names <- LETTERS[1:10] %>% setNames(paste0("R", 0:9))

p1 <- ggdraw()

# 1. c matrix ----
cml <- cm %>% # c matrix longer
    pivot_longer(cols = starts_with("R"), names_to = "Resource", values_to = "ConsumptionRate") %>%
    rename(Family = ...1, Species = ...2) %>%
    left_join(mal) %>%
    left_join(sal) %>%
    select(Family, Species, SpeciesID, Class, Resource, ResourceID, ConsumptionRate)


# Stacked bar plot
p2 <- cml %>%
    filter(Species %in% (sal %>% group_by(Family) %>% slice(1:10) %>% pull(Species))) %>%
    mutate(Species = factor(Species, sal$Species)) %>%
    ggplot()+
    geom_col(aes(x = Species, y = ConsumptionRate, fill = Resource), width = 1) +
    scale_fill_manual(values = mcrm_resource_colors, label = resource_names) +
    scale_y_continuous(expand = c(0,0), limits = c(0, 1.05), breaks = c(0,0.5, 1)) +
    facet_grid(.~Family, scales = "free_x", labeller = as_labeller(family_names)) +
    theme_classic() +
    theme(panel.border = element_rect(color = 1, fill = NA),
          strip.background = element_rect(color = NA, fill = NA),
          panel.spacing.x = unit(0, "mm"),
          plot.title = element_text(hjust = 0.5, vjust = -1, size = 10),
          axis.text.x = element_blank()) +
    guides(fill = guide_legend(ncol = 2)) +
    labs(x = "species", y = "uptake rate") +
    ggtitle("family")

# 2. D matrix
Dml <- Dm %>% #
    pivot_longer(cols = starts_with("R"), names_to = "Resource1", values_to = "SecretionFlux") %>%
    rename(Class2 = ...1, Resource2 = ...2) %>%
    select(-Class2) %>%
    left_join(rename_with(mal, ~ paste0(., "1"), everything()), by = join_by(Resource1)) %>%
    left_join(rename_with(mal, ~ paste0(., "2"), everything()), by = join_by(Resource2)) %>%
    select(Class1, Resource1, ResourceID1, Class2, Resource2, ResourceID2, SecretionFlux)

p3 <- Dml %>%
    ggplot() +
    geom_tile(aes(x = ResourceID1, y = ResourceID2, fill = SecretionFlux)) +
    scale_fill_gradient(high = "black", low = "white") +
    scale_color_manual(values = category_color, breaks = c("sugar", "acid")) +
    scale_x_continuous(position = "top", breaks = 1:(ma*fa)) +
    coord_cartesian(xlim = c(1, ma*fa), ylim = c(ma*fa, 1), clip = "off") +
    theme_classic() +
    theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank()) +
    guides(fill = guide_colorbar(title = expression(D[beta][alpha])),
           color = "none") +
    labs(x = expression(R[alpha]), y = expression(R[beta]))

# Assemble panels ----
p <- plot_grid(
    p1,
    plot_grid(
        p2,
        plot_grid(p3, NULL, rel_widths = c(1,1)),
        ncol = 1, rel_heights = c(1,1), labels = c("B", "C"), scale = c(.95, .95)),
    nrow = 1, labels = c("A", ""), rel_widths = c(1,2)
) +
    paint_white_background()
ggsave(here::here("plots/FigS15-mcrm_parameters.png"), p, width = 8, height = 4)










