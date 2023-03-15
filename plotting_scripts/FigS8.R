library(tidyverse)
library(cowplot)
library(broom)
source(here::here("analysis/00-metadata.R"))

isolates <- read_csv(paste0(folder_data, "output/isolates.csv"), show_col_types = F)
pairs <- read_csv(paste0(folder_data, "output/pairs.csv"), show_col_types = F)
communities <- read_csv(paste0(folder_data, "temp/00c-communities.csv"), show_col_types = F)

# 0.1 isolate taxa ----
isolates_RDP <- read_csv(paste0(folder_data, "temp/12-isolates_RDP.csv"), col_types = cols())
isolates_ID <- read_csv(paste0(folder_data, "temp/00c-isolates_ID.csv"), col_types = cols()) %>%
    left_join(select(isolates_RDP, ExpID, Family, Genus, Fermenter))


# 0. carbons source list ----
# Metabolite types
csl <- tibble(Source = c("glucose", "fructose", "galactose", "ribose", "arabinose",

                         "malate", "acetate", "glycerol", "pyruvate", "succinate", "fumarate",

                         "glycine", "acetyl-ornithine", "alanine", "succinate", "lactate",
                         "serine", "valine", "pyruvate", "alpha-ketoglutarate", "putrescine",
                         "valerate", "asparagine", "butyrate", "methionine", "hippurate",
                         "propionate", "citrate", "acetate", "beta-hydroxybutyrate", "fumarate",
                         "gluconate", "glutamine", "ketogluconate", "leucine", "oxoglutarate", "phenylalanine"),
              Type = c(rep("sugar", 5), rep("acid", 6+20+6)),
              Carbon = c(6, 6, 6, 5, 5,
                         4, 2, 3, 3, 4, 4,
                         2, 7, 3, 4, 3,
                         3, 5, 3, 5, 4,
                         5, 4, 4, 5, 9,
                         3, 6, 2, 4, 4,
                         6, 5, 6, 6, 5, 9)) %>%
    distinct(Source, Carbon, Type)

# Strain type
stl <- tibble(Strain = c("none", "Ecoli", "Enterobacter", "Pputida", "Pseudomonas"),
              Fermenter = c("none", "fermenter", "fermenter", "respirator", "respirator"))

write_csv(csl, paste0(folder_data, "temp/21-csl.csv"))
write_csv(stl, paste0(folder_data, "temp/21-stl.csv"))


# Data from Jean
isolates_curves <- read_csv(paste0(folder_data, "raw/growth_rate/20GC_Data.csv"), col_types = cols()) %>%
    select(ID = SangerID, CS, Time, OD620) %>%
    filter(ID %in% isolates_ID$ID) %>%
    # Clean the Carbon Source names
    mutate(CS = tolower(CS)) %>%
    mutate(CS = str_replace(CS, "d-", "") %>% str_replace("l-", "") %>% str_replace("2-", "")) %>%
    # Clean the Time points
    mutate(Time = cut_width(Time, 1, boundary = 0, labels = F))  %>%
    distinct(ID, CS, Time, .keep_all = T) %>%
    filter(CS != "ddh20")

# Growth rates using the time points 12, 16, 28 hr -----
calculate_r <- function(N0, N1, T0, T1) (log10(N1)-log10(N0)) / (T1 - T0)
isolates_curves_T0 <- isolates_curves %>%
    group_by(ID, CS) %>%
    filter(Time == min(Time)) %>%
    mutate(T0 = ifelse(Time == min(Time), 0, Time), N0 = OD620) %>%
    # Assign a minimum OD value to prevent error in log(0)
    mutate(N0 = ifelse(N0 == 0, 0.001, N0)) %>%
    select(-Time, -OD620)
isolates_growth <- isolates_curves %>%
    filter(Time %in% c(12, 16, 28)) %>%
    group_by(ID, CS) %>%
    arrange(ID, CS, Time) %>%
    mutate(T1 = Time, N1 = OD620) %>%
    select(-Time, -OD620) %>%
    left_join(isolates_curves_T0) %>%
    # Remove negative OD reads
    filter(N1 > 0) %>%
    # Calculate r
    mutate(r = calculate_r(N0, N1, T0, T1)) %>%
    select(ID, CS, Time = T1, r) %>%
    # Remove contamination
    filter(r>0) %>%
    # Average
    group_by(ID, CS, Time) %>%
    summarize(r = mean(r)) %>%
    pivot_wider(names_from = c(Time, CS), values_from = r, names_glue = "r_{CS}_{Time}hr") %>%
    ungroup()

isolates_OD <- isolates_curves %>%
    filter(Time %in% c(12, 16, 28)) %>%
    group_by(ID, CS) %>%
    pivot_wider(names_from = c(Time, CS), values_from = OD620, names_glue = "OD_{CS}_{Time}hr")


# ----
isolates <- isolates %>%
    left_join(isolates_growth) %>%
    group_by(Community) %>%
    mutate(rank_r_glucose_16hr = rank(r_glucose_16hr))

# Figure -----
p <- isolates %>%
    ggplot(aes(x = rank_r_glucose_16hr, y = Rank, group = ExpID)) +
    geom_jitter(shape = 21, size = 2, stroke = 1, width = 0.1, height = 0.1) +
    scale_fill_manual(values = interaction_color) +
    scale_x_continuous(breaks = 1:12) +
    scale_y_continuous(breaks = 1:12) +
    theme_classic() +
    guides(fill = guide_legend(title = "")) +
    labs(x = expression(ranked~growth~rate), y = "competitive rank")

ggsave(here::here("plots/FigS8-growth_vs_rank.png"), p, width = 4, height = 4)

#
cor.test(isolates$Rank, isolates$rank_r_glucose_16hr,
         method = "spearman", alternative = "two.sided", exact = FALSE)

