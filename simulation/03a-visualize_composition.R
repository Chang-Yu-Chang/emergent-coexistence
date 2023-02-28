#' This scripts generate the barplot for N and R compositions
#' N: consumer abundance
#' R: resource abundance

library(tidyverse)
library(cowplot)
source(here::here("analysis/00-metadata.R"))

# 0. parameters ----
# Read parameters
input_parameters <- read_csv(here::here("simulation/01-input_parameters.csv"), col_types = cols())
input_monocultures <- read_csv(here::here("simulation/01a-input_monocultures.csv"), col_types = cols())
input_communities <- read_csv(here::here("simulation/01b-input_communities.csv"), col_types = cols())
category_colors <- c(sugar = "#ED6A5A", acid = "#03CEA4", waste = "#51513D", fermenter = "#8A89C0", respirator = "#FFCB77")
input_row <- input_parameters[1,]

# Generate family-species and class-resource tibble for matching
sa <- input_parameters$sa[1]
ma <- input_parameters$ma[1]
sal <- tibble(Family = paste0("F", c(rep(0, sa), rep(1, sa))), Species = paste0("S", 0:(sa * 2 - 1)))
mal <- tibble(Class = paste0("T", c(rep(0, ma), rep(1, ma))), Resource = paste0("R", 0:(ma * 2 - 1)))

# 1. Self-assembled community composition ----
# Read simulation output files
read_wide_file <- function(x, type = "N") {
    temp <- read_csv(x, col_types = cols()) %>%
        pivot_longer(cols = starts_with("W"), names_to = "Well", values_to = "Abundance")
    if ("...1" %in% colnames(temp)) {
        if (type == "N") temp <- temp %>% rename(Family = ...1, Species = ...2)
        if (type == "R") temp <- temp %>% rename(Class = ...1, Resource = ...2)
    }

    return(temp)
}
N_long <- list.files(input_communities$output_dir[1], pattern = "selfAssembly-1-N") %>%
    map(c("N_T\\d+\\.csv", "init"), str_subset, string = .) %>%
    reduce(union) %>%
    paste0(input_communities$output_dir[1], .) %>%
    lapply(function(x) {
        time_point <- str_split_i(x, "-", -1) %>%
            str_replace("N_", "") %>%
            str_replace(".csv", "")
        read_wide_file(x) %>%
            mutate(Time = time_point) %>%
            return()
    }) %>%
    bind_rows %>%
    mutate(Well = ordered(Well, paste0("W", 0:(input_row$n_wells-1)))) %>%
    mutate(Time = ordered(Time, c("init", paste0("T", 1:20)))) %>%
    arrange(Well, Time)


# Line plot
p1 <- N_long %>%
    filter(Abundance != 0 ) %>%
    ggplot(aes(x = Time, y = Abundance, color = Family, group = Species)) +
    geom_line(linewidth = .2) +
    geom_point(size = 1, shape = 21) +
    scale_color_manual(values = c("F0" = "#8A89C0", "F1" = "#FFCB77")) +
    scale_linetype_manual(values = c("F0" = 1, "F1" = 2)) +
    facet_wrap(.~Well, ncol = 5) +
    theme_classic() +
    theme(panel.border = element_rect(color = 1, fill = NA)) +
    guides(alpha = "none") +
    labs()
ggsave(here::here("simulation/plots/03a-community_line.png"), p1, width = 12, height = 10)

# Barplot over time
p2 <- N_long %>%
    filter(Abundance != 0) %>%
    ggplot(aes(x = Time, y = Abundance, fill = Family, color = Species)) +
    geom_col() +
    scale_fill_manual(values = c("F0" = "#8A89C0", "F1" = "#FFCB77")) +
    scale_color_manual(values = rep("black", length(sal$Species))) +
    facet_wrap(.~Well, ncol = 5) +
    theme_classic() +
    guides(color = "none") +
    labs()

ggsave(here::here("simulation/plots/03a-community_bar.png"), p2, width = 12, height = 10)

# Barplot over time, standard
p3 <- N_long %>%
    filter(Abundance != 0) %>%
    ggplot(aes(x = Time, y = Abundance, fill = Family, color = Species)) +
    geom_col(position = "fill") +
    scale_fill_manual(values = c("F0" = "#8A89C0", "F1" = "#FFCB77")) +
    scale_color_manual(values = rep("black", length(sal$Species))) +
    facet_wrap(.~Well, ncol = 5) +
    theme_classic() +
    guides(color = "none") +
    labs()
ggsave(here::here("simulation/plots/03a-community_bar_standard.png"), p3, width = 12, height = 10)



# Barplot final time point
N_long_abundant <- N_long %>%
    filter(Time == "T5") %>%
    filter(Abundance != 0) %>%
    group_by(Well) %>%
    mutate(TotalAbundance = sum(Abundance)) %>%
    #filter(Well == "W0") %>%
    filter(Abundance > 0.01 * sum(Abundance))

N_summmary <- N_long_abundant %>%
    summarize(Richness = n())
p4 <- N_long_abundant %>%
    filter(Time == "T5") %>%
    filter(Abundance != 0) %>%
    group_by(Well) %>%
    filter(Abundance > 0.01 * sum(Abundance)) %>%
    mutate(RelativeAbundance = Abundance / sum(Abundance)) %>%
    ggplot() +
    geom_col(aes(x = Well, y = RelativeAbundance, fill = Family), color = 1) +
    geom_text(data = N_summmary, aes(x = Well, label = Richness), y = 1.1) +
    scale_fill_manual(values = c("F0" = "#8A89C0", "F1" = "#FFCB77")) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 1.1)) +
    coord_cartesian(ylim = c(0, 1), clip = "off") +
    theme_classic() +
    theme(plot.margin = unit(c(1,.5,.5,.5), "cm")) +
    guides(color = "none") +
    labs()
ggsave(here::here("simulation/plots/03a-community_bar_final.png"), p4, width = 9, height = 3)




# 2. Monoculture ----
# Read monoculture data
N_long <- list.files(input_monocultures$output_dir[1], pattern = "monoculture-1-N") %>%
    map(c("N_T\\d+\\.csv", "init"), str_subset, string = .) %>%
    reduce(union) %>%
    paste0(input_monocultures$output_dir[1], .) %>%
    lapply(function(x) {
        time_point <- str_split_i(x, "-", -1) %>%
            str_replace("N_", "") %>%
            str_replace(".csv", "")
        read_wide_file(x) %>%
            mutate(Time = time_point) %>%
            return()
    }) %>%
    bind_rows %>%
    mutate(Well = ordered(Well, paste0("W", 0:(input_row$n_wells-1)))) %>%
    mutate(Time = ordered(Time, c("init", paste0("T", 1:20)))) %>%
    arrange(Well, Time)


# Line plot
p1 <- N_long %>%
    #filter(Well == "W1") %>%
    filter(Abundance != 0) %>%
    ggplot(aes(x = Time, y = Abundance, color = Family, group = Species)) +
    geom_line(linewidth = .2) +
    geom_point(size = 1, shape = 21) +
    scale_color_manual(values = c("F0" = "#8A89C0", "F1" = "#FFCB77")) +
    scale_linetype_manual(values = c("F0" = 1, "F1" = 2)) +
    #scale_y_log10() +
    theme_classic() +
    theme(panel.border = element_rect(color = 1, fill = NA)) +
    guides(alpha = "none") +
    labs()
ggsave(here::here("simulation/plots/03a-monoculture_line.png"), p1, width = 5, height = 4)





























