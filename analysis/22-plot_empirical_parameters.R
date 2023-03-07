#' This script read the cleaned the empirical data for empirical parameters for MCRM models
#'
#' 1. c matrix from the growth curves
#' 2. D matrix from the metabolomics
#' 3. l matrix from the metabolites

library(tidyverse)
library(cowplot)
source(here::here("analysis/00-metadata.R"))

# 0. carbon source list ----
csl <- read_csv(paste0(folder_data, "temp/21-csl.csv"), col_types = cols())
stl <- read_csv(paste0(folder_data, "temp/21-stl.csv"), col_types = cols())


# 1. c matrix from the growth curves ----
isolates_u <- read_csv(paste0(folder_data, "temp/21-isolates_u.csv"), col_types = cols())
p1 <- isolates_u %>%
    left_join(csl) %>%
    filter(Type == "sugar") %>%
    ggplot() +
    geom_histogram(aes(x = u, fill = Fermenter), color = 1, position = "identity", alpha = 0.5) +
    facet_wrap(Source~., ncol = 1) +
    scale_fill_manual(values = category_color, breaks = c("fermenter", "respirator")) +
    theme_classic() +
    theme(legend.title = element_blank()) +
    guides(fill = "none") +
    ggtitle("Sugar")
p2 <- isolates_u %>%
    left_join(csl) %>%
    filter(Type == "acid") %>%
    ggplot() +
    geom_histogram(aes(x = u, fill = Fermenter), color = 1, position = "identity", alpha = 0.5) +
    facet_wrap(Source~., nrow = 5) +
    scale_fill_manual(values = category_color, breaks = c("fermenter", "respirator")) +
    theme_classic() +
    theme(legend.title = element_blank(), legend.position = "right") +
    ggtitle("Organic acids")
p <- plot_grid(p1, p2, nrow = 1, rel_widths = c(1, 3), axis = "tb", align = "h")

ggsave(paste0(folder_data, "temp/22-update_rate.png"), p, width = 10, height = 10)


##
isolates_u %>%
    left_join(select(csl, Source, SourceType = Type)) %>%
    drop_na() %>%
    group_by(Fermenter, ExpID, SourceType) %>%
    # Name
    mutate(SourceType = factor(SourceType, c("sugar", "acid"))) %>%
    mutate(Fermenter = factor(Fermenter, c("fermenter", "respirator"))) %>%
    summarize(uSum = sum(u)) %>%
    group_by(Fermenter, SourceType) %>%
    summarize(uSumMean = mean(uSum, na.rm = T), SumSd = sd(uSum, na.rm = T))


# 2. D matrix from the metabolomics ----
metabolomics <- read_csv(paste0(folder_data, "temp/21-metabolomics.csv"), col_types = cols())

## Plot the raw data
p <- metabolomics %>%
    mutate(Source = ordered(Source, csl$Source)) %>%
    mutate(Metabolite = ordered(Metabolite, csl$Source)) %>%
    mutate(MetaboliteConc = log(MetaboliteConc)) %>%
    ggplot() +
    geom_tile(aes(x = Source, y = Metabolite, fill = MetaboliteConc)) +
    facet_grid(.~Strain, scales = c("free_x")) +
    scale_fill_gradient(low = "white", high = "blue") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5))

ggsave(paste0(folder_data, "temp/22-metabolomics.png"), p, width = 15, height = 5)


## Calculate the fraction of conversion
metabolomics %>%
    mutate(SourceType = factor(SourceType, c("sugar", "acid"))) %>%
    mutate(MetaboliteType = factor(MetaboliteType, c("sugar", "acid"))) %>%
    # Use only two strains. Check if treating the E and P isolates as replicates
    #filter(Strain %in% c("Ecoli", "Pputida")) %>%
    # Average for group-CS-metabolite
    group_by(Fermenter, SourceType, Source, MetaboliteType) %>%
    summarize(MetaboliteConc = sum(MetaboliteConc)) %>%
    mutate(RelativeMetaboliteConc = MetaboliteConc / sum(MetaboliteConc)) %>%
    # Average for group-CS
    group_by(Fermenter, SourceType, MetaboliteType) %>%
    summarize(RelativeMetaboliteConc = mean(RelativeMetaboliteConc)) %>%
    mutate(Abbreviation = paste0("f", str_sub(Fermenter, 1, 1), str_sub(SourceType, 1, 1), str_sub(MetaboliteType, 1, 1)))




# 3. l matrix from the metabolites ----

isolates_leakiness <- read_csv(paste0(folder_data, "temp/21-isolates_leakiness.csv"), col_types = cols())
isolates_leakiness %>%
    mutate(Fermenter = ifelse(Fermenter, "fermenter", "respirator")) %>%
    group_by(Fermenter) %>%
    summarize(leakiness_16hr_mean = mean(leakiness_16hr),
              leakiness_16hr_sd = sd(leakiness_16hr))














