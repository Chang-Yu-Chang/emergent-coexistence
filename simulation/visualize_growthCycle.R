# Visualize a test growth cycle
library(tidyverse)
library(cowplot)

#
input_parameters <- read_csv(here::here("simulation/input_parameters.csv"))
output_dir <- input_parameters$output_dir[1]
category_colors <- c(sugar = "#ED6A5A", acid = "#03CEA4", waste = "#51513D", fermenter = "#8A89C0", respirator = "#FFCB77")
input_row <- input_parameters[2,]

# Generate family-species and class-resource matching tibble
sa <- input_parameters$sa[1]
ma <- input_parameters$ma[1]
sal <- tibble(Family = paste0("F", c(rep(0, sa), rep(1, sa))), Species = paste0("S", 0:(sa * 2 - 1)))
mal <- tibble(Class = paste0("T", c(rep(0, ma), rep(1, ma), rep(2, ma))), Resource = paste0("R", 0:(ma * 3 - 1)))

#
read_wide_file <- function(x) {
    read_csv(x, col_types = cols()) %>%
        pivot_longer(cols = starts_with("W"), names_to = "Well", values_to = "Abundance") %>%
        rename(Family = ...1, Species = ...2)
}
N_long <- list.files(output_dir, pattern = "id1_") %>%
    paste0(output_dir, .) %>%
    #`[`(1:2) %>%
    lapply(function(x) read_wide_file(x) %>% mutate(Time = str_replace(x, paste0(output_dir, "id1_t"), "") %>% str_replace(".csv", "") %>% as.numeric)) %>%
    bind_rows %>%
    mutate(Well = ordered(Well, paste0("W", 0:(input_row$n_wells-1)))) %>%
    arrange(Well)


# Line plot
p1 <- N_long %>%
    filter(Well == "W1") %>%
    filter(Abundance != 0 ) %>%
    ggplot(aes(x = Time, y = Abundance, color = Species)) +
    geom_line(aes(linetype = Family), lwd = 1) +
    geom_point(aes(shape = Family), size = 2, stroke = 2, fill = "white") +
    scale_shape_manual(values = c("F0" = 16, "F1" = 21)) +
    scale_linetype_manual(values = c("F0" = 1, "F1" = 2)) +
    scale_y_log10() +
    theme_classic() +
    guides(color = "none") +
    labs()

# Barplot over time
p2 <- N_long %>%
    filter(Abundance != 0) %>%
    ggplot(aes(x = Time, y = Abundance, fill = Family, color = Species)) +
    geom_col() +
    scale_fill_manual(values = c("F0" = "#8A89C0", "F1" = "#FFCB77")) +
    scale_color_manual(values = rep("black", length(sal$Species))) +
    facet_wrap(.~Well, ncol = 5) +
    theme_classic() +
    guides(color = "none", fill = "none") +
    labs()
ggsave(here::here("simulation/plots/02-barplot_time.png"), p2, width = 10, height = 10)

# Barplot over time, standard
p3 <- N_long %>%
    filter(Abundance != 0) %>%
    ggplot(aes(x = Time, y = Abundance, fill = Family, color = Species)) +
    geom_col(position = "fill") +
    scale_fill_manual(values = c("F0" = "#8A89C0", "F1" = "#FFCB77")) +
    scale_color_manual(values = rep("black", length(sal$Species))) +
    facet_wrap(.~Well, ncol = 5) +
    theme_classic() +
    guides(color = "none", fill = "none") +
    labs()
ggsave(here::here("simulation/plots/03-barplot_time_standard.png"), p3, width = 10, height = 10)



# Barplot final time point
temp <- N_long %>%
    filter(Time == max(Time)) %>%
    filter(Abundance != 0) %>%
    group_by(Well) %>%
    summarize(TotalAbundance = round(sum(Abundance), 0))
p4 <- N_long %>%
    filter(Time == 24) %>%
    filter(Abundance != 0) %>%
    ggplot() +
    geom_col(aes(x = Well, y = Abundance, fill = Family), position = "fill", color = 1) +
    # Total abundance (biomass)
    geom_text(data = temp, aes(x = Well, label = TotalAbundance), y = Inf, vjust = 1) +
    scale_fill_manual(values = c("F0" = "#8A89C0", "F1" = "#FFCB77")) +
    theme_classic() +
    guides(color = "none") +
    labs()

ggsave(here::here("simulation/plots/04-barplot_final.png"), p4, width = 10, height = 5)












