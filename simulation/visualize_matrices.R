library(tidyverse)
library(cowplot)

#
input_parameters <- read_csv(here::here("simulation/input_parameters.csv"))
output_dir <- input_parameters$output_dir[1]
Dm <- read_csv(paste0(output_dir, "D_seed1.csv"), skip = 1) # D matrix
cm <- read_csv(paste0(output_dir, "c_seed1.csv"), skip = 1) # c matrix
lm <- read_csv(paste0(output_dir, "l_seed1.csv"), skip = 1) # l matrix
category_colors <- c(sugar = "#ED6A5A", acid = "#03CEA4", waste = "#51513D", fermenter = "#8A89C0", respirator = "#FFCB77")

# Generate family-species and class-resource matching tibble
sa <- input_parameters$sa[1]
ma <- input_parameters$ma[1]
sal <- tibble(Family = paste0("F", c(rep(0, sa), rep(1, sa))), Species = paste0("S", 0:(sa * 2 - 1)))
mal <- tibble(Class = paste0("T", c(rep(0, ma), rep(1, ma), rep(2, ma))), Resource = paste0("R", 0:(ma * 3 - 1)))

# D matrix
Dml <- Dm %>% # D matrix longer
    pivot_longer(cols = starts_with("R"), names_to = "Resource1", values_to = "SecretionFlux") %>%
    rename(Class2 = ...1, Resource2 = ...2) %>%
    left_join(rename_with(mal, ~ paste0(., "1"), everything())) %>%
    mutate(Resource1 = ordered(Resource1, mal$Resource), Resource2 = ordered(Resource2, mal$Resource)) %>%
    select(Class1, Resource1, Class2, Resource2, SecretionFlux)
p1 <- Dml %>%
    ggplot() +
    geom_tile(aes(x = Resource1, y = Resource2, fill = SecretionFlux)) +
    geom_segment(aes(color = "sugar"), x = "R0", xend = paste0("R", ma-1), y = "spacer", yend = "spacer", lwd = 2) +
    geom_segment(aes(color = "acid"), x = paste0("R", ma), xend = paste0("R", 2*ma-1), y = "spacer", yend = "spacer", lwd = 2) +
    geom_segment(aes(color = "waste"), x = paste0("R", 2*ma), xend = paste0("R", 3*ma-1), y = "spacer", yend = "spacer", lwd = 2) +
    geom_segment(aes(color = "sugar"), x = "spacer", xend = "spacer", y = "R0", yend = paste0("R", ma-1), lwd = 2) +
    geom_segment(aes(color = "acid"), x = "spacer", xend = "spacer", y = paste0("R", ma), yend = paste0("R", 2*ma-1), lwd = 2) +
    geom_segment(aes(color = "waste"), x = "spacer", xend = "spacer", y = paste0("R", 2*ma), yend = paste0("R", 3*ma-1), lwd = 2) +
    scale_fill_gradient(high = "black", low = "white") +
    scale_color_manual(values = category_colors, breaks = c("sugar", "acid", "waste")) +
    scale_x_discrete(limits = c("spacer", mal$Resource)) +
    scale_y_discrete(limits = c("spacer", mal$Resource)) +
    theme_classic() +
    theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank()) +
    guides(fill = guide_colorbar(title = expression(D[i][beta][alpha])),
           color = guide_legend(title = "")) +
    labs(x = expression(R[alpha]), y = expression(R[beta]))
ggsave(here::here("simulation/plots/matrix1_D.png"), p1, width = 6, height = 6)


# c matrix
cml <- cm %>% # c matrix longer
    pivot_longer(cols = starts_with("R"), names_to = "Resource", values_to = "ConsumptionRate") %>%
    rename(Family = ...1, Species = ...2) %>%
    left_join(mal) %>%
    mutate(Species = ordered(Species, sal$Species), Resource = ordered(Resource, mal$Resource)) %>%
    select(Family, Species, Class, Resource, ConsumptionRate)
p2 <- cml %>%
    ggplot() +
    geom_tile(aes(x = Resource, y = Species, fill = ConsumptionRate)) +
    geom_segment(aes(color = "sugar"), x = "R0", xend = paste0("R", ma-1), y = "spacer", yend = "spacer", lwd = 2) +
    geom_segment(aes(color = "acid"), x = paste0("R", ma), xend = paste0("R", 2*ma-1), y = "spacer", yend = "spacer", lwd = 2) +
    geom_segment(aes(color = "waste"), x = paste0("R", 2*ma), xend = paste0("R", 3*ma-1), y = "spacer", yend = "spacer", lwd = 2) +
    geom_segment(aes(color = "fermenter"), x = "spacer", xend = "spacer", y = "S0", yend = paste0("S", sa-1), lwd = 2) +
    geom_segment(aes(color = "respirator"), x = "spacer", xend = "spacer", y = paste0("S", sa), yend = paste0("S", 2*sa-1), lwd = 2) +
    scale_fill_gradient(high = "black", low = "white") +
    scale_color_manual(values = category_colors) +
    scale_x_discrete(limits = c("spacer", mal$Resource)) +
    scale_y_discrete(limits = c("spacer", sal$Species)) +
    theme_classic() +
    theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank()) +
    guides(fill = guide_colorbar(title = expression(c[i][alpha])),
           color = guide_legend(title = "")) +
    labs()
ggsave(here::here("simulation/plots/matrix2_c.png"), p2, width = 6, height = 6)


# l matrix
lml <- lm %>% # l matrix longer
    pivot_longer(cols = starts_with("R"), names_to = "Resource", values_to = "Leakiness") %>%
    rename(Family = ...1, Species = ...2) %>%
    left_join(mal) %>%
    mutate(Species = ordered(Species, sal$Species), Resource = ordered(Resource, mal$Resource)) %>%
    select(Family, Species, Class, Resource, Leakiness)
p3 <- lml %>%
    ggplot() +
    geom_tile(aes(x = Resource, y = Species, fill = Leakiness)) +
    geom_segment(aes(color = "sugar"), x = "R0", xend = paste0("R", ma-1), y = "spacer", yend = "spacer", lwd = 2) +
    geom_segment(aes(color = "acid"), x = paste0("R", ma), xend = paste0("R", 2*ma-1), y = "spacer", yend = "spacer", lwd = 2) +
    geom_segment(aes(color = "waste"), x = paste0("R", 2*ma), xend = paste0("R", 3*ma-1), y = "spacer", yend = "spacer", lwd = 2) +
    geom_segment(aes(color = "fermenter"), x = "spacer", xend = "spacer", y = "S0", yend = paste0("S", sa-1), lwd = 2) +
    geom_segment(aes(color = "respirator"), x = "spacer", xend = "spacer", y = paste0("S", sa), yend = paste0("S", 2*sa-1), lwd = 2) +
    scale_fill_gradient(high = "black", low = "white") +
    scale_color_manual(values = category_colors) +
    scale_x_discrete(limits = c("spacer", mal$Resource)) +
    scale_y_discrete(limits = c("spacer", sal$Species)) +
    theme_classic() +
    theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank()) +
    guides(fill = guide_colorbar(title = expression(l[i][alpha])),
           color = guide_legend(title = "")) +
    labs()
ggsave(here::here("simulation/plots/matrix3_l.png"), p3, width = 6, height = 6)
























