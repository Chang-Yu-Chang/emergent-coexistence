library(tidyverse)
library(cowplot)
library(latex2exp)
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
resource_names <- 1:10 %>% setNames(paste0("R", 0:9))

cml <- cm %>% # c matrix longer
    pivot_longer(cols = starts_with("R"), names_to = "Resource", values_to = "ConsumptionRate") %>%
    rename(Family = ...1, Species = ...2) %>%
    left_join(mal) %>%
    left_join(sal) %>%
    select(Family, Species, SpeciesID, Class, Resource, ResourceID, ConsumptionRate)

#
annotate_title <- function (x, y, label) annotate("text", x = y, y = x, hjust = 0, size = 3, label = label)
annotate_equation <- function (x, y, equation) annotate("text", x = y, y = x, hjust = 0, size = 2.5, parse = TRUE, label = TeX(equation, output = 'character'))

p1 <- ggplot() +
    annotate_title(x = 0.1, y = 0.9, label = "1. Generate Dirichlet concentration parameter") +
    geom_col(data = filter(cml, Species == "S0") %>% mutate(ConsumptionRate = ConsumptionRate / sum(ConsumptionRate)) %>% mutate(x = 0.75),
             aes(x = x, y = ConsumptionRate, fill = Resource), width = 0.05, position = position_stack(reverse = T)) +
    annotate("segment", y = 0.23, x = 0.8, yend = 0.23, xend = 0.75, arrow = arrow(length = unit(2, "mm"))) +
    annotate_equation(x = 0.11, y = 0.83, equation = r'($ \theta_{\alpha=f} \sim Gaussian (\mu_f, \sigma^2_f) $)') +
    annotate_equation(x = 0.5, y = 0.65, equation = r'($ \theta'_{\alpha \neq f} = Uniform (0,1) $)') +
    annotate_equation(x = 0.5, y = 0.62, equation = r'($ \theta_{\alpha \neq f} = (1-\theta_{\alpha=f}) (\theta'_{\alpha \neq f})/(\Sigma_{z \neq f} \theta'_z) $)') +
    annotate("segment", y = 0.4, x = 0.72, yend = 1, xend = 0.72, arrow = arrow(ends = "both", angle = 90, length = unit(1, "mm"))) +
    annotate("segment", y = 0.7, x = 0.67, yend = 0.7, xend = 0.7, arrow = arrow(length = unit(2, "mm"))) +
    annotate_title(x = 0.1, y = 0.55, label = "2. Sample proportion of uptake rates") +
    annotate_equation(x = 0.2, y = 0.5, equation = r'($c'_{i,1}, c'_{i,2}, ..., c'_{i,M} \sim Dirichlet (\theta_1 \Omega_f, \theta_2 \Omega_f, ..., \theta_M \Omega_f)$)') +
    annotate_title(x = 0.1, y = 0.45, label = "3. Sample total uptake capacity") +
    annotate_equation(x = 0.2, y = 0.41, equation = r'($T_i \sim Gaussian (\mu_T, \sigma^2_T)$)') +
    annotate_title(x = 0.1, y = 0.36, label = "4. Calculate uptake rates") +
    annotate_equation(x = 0.2, y = 0.32, equation = r'($c_{i\alpha} = T_i\;c'_{i\alpha}$)') +
    annotate_title(x = 0.1, y = 0.25, label = "5. Sample stoichiometric matrix") +
    annotate_equation(x = 0.2, y = 0.21, equation = r'($D_{\beta\alpha} \sim Uniform (1,\;1/M) $)') +
    scale_x_continuous(limits = c(0, 1)) +
    scale_y_continuous(limits = c(0, 1)) +
    scale_fill_manual(values = mcrm_resource_colors) +
    coord_flip() +
    #theme_classic() +
    theme_nothing() +
    theme(plot.margin = margin(0,0,0,0,"mm")) +
    labs()


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
          panel.spacing.x = unit(0, "mm"),
          strip.background = element_rect(color = NA, fill = NA),
          axis.text.x = element_blank(),
          plot.title = element_text(hjust = 0.5, vjust = -1, size = 10),
          legend.title = element_text(size = 8),
          legend.text = element_text(size = 8),
          legend.key.size = unit(3, "mm")) +
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
        plot_grid(p3, NULL, rel_widths = c(2.2,1)),
        ncol = 1, rel_heights = c(1,1), labels = c("B", "C"), scale = c(.95, .95)),
    nrow = 1, labels = c("A", ""), rel_widths = c(1,1.5)
) +
    paint_white_background()
ggsave(here::here("plots/FigS15-mcrm_parameters.png"), p, width = 8, height = 5)










