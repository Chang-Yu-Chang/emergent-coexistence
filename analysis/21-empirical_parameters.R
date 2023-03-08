#' This script cleans up the empirical data used to find the empirical parameters for MCRM models
#'
#' 1. c matrix from the growth curves
#' 2. D matrix from the metabolomics
#' 3. l matrix from the metabolites

library(tidyverse)
library(janitor)
source(here::here("analysis/00-metadata.R"))

# 0.1 isolate taxa ----
isolates_RDP <- read_csv(paste0(folder_data, "temp/12-isolates_RDP.csv"), col_types = cols())
isolates_ID <- read_csv(paste0(folder_data, "temp/00c-isolates_ID.csv"), col_types = cols()) %>%
    left_join(select(isolates_RDP, ExpID, Family, Genus, Fermenter))

# 0.2 carbons source list ----

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



# 1. c matrix from the growth curves ----

# Growth curve from Jean
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


# Estimate the uptake rates
isolates_u <- isolates_ID %>%
    select(ID, ExpID, Community, Fermenter) %>%
    left_join(isolates_growth) %>%
    left_join(isolates_OD) %>%
    #drop_na() %>%
    mutate(Fermenter = ifelse(Fermenter, "fermenter", "respirator")) %>%
    distinct(ExpID, .keep_all = T) %>%
    select(Fermenter, ExpID, starts_with("r") & ends_with("16hr"), starts_with("OD") & ends_with("16hr")) %>%
    pivot_longer(cols = ends_with("16hr")) %>%
    separate(name, into = c("name", "CS", "Time")) %>%
    mutate(Time = str_replace(Time, "hr", "") %>% as.numeric) %>%
    pivot_wider(names_from = name, values_from = value) %>%
    # Uptake rates
    mutate(u = r * OD) %>%
    arrange(Fermenter, CS) %>%
    rename(Source = CS)

write_csv(isolates_u, paste0(folder_data, "temp/21-isolates_u.csv"))


# 2. D matrix from the metabolomics ----
mb <- read_csv(paste0(folder_data, "raw/metabolomics/TMIC_LCMS.csv"), col_types = cols())

# Clean up names
mb <- mb %>%
    rename_with(~str_replace(., "_", "")) %>%
    #remove_empty() %>%
    mutate(Source = tolower(CarbonSource)) %>%
    mutate(Metabolite = tolower(Metabolite)) %>%
    mutate(Strain = ifelse(Strain == "None", "none", Strain)) %>%
    mutate(Metabolite = str_replace(Metabolite, "d-\\|l-", "")) %>%
    mutate(Source = str_replace(Source, "d-|l-", "")) %>%
    filter(Source != "all") %>%
    mutate(Metabolite = case_when(
        Metabolite == "acetic acid" ~ "acetate",
        Metabolite == "succinic acid" ~ "succinate",
        Metabolite == "lactic acid" ~ "lactate",
        Metabolite == "citric acid" ~ "citrate",
        Metabolite == "pyruvic acid" ~ "pyruvate",
        Metabolite == "butyric acid" ~ "butyrate",
        Metabolite == "hippuric acid" ~ "hippurate",
        Metabolite == "propionic acid" ~ "propionate",
        Metabolite == "fumaric acid" ~ "fumarate",
        Metabolite == "valeric acid" ~ "valerate",
        Metabolite == "beta-hydroxybutyric acid" ~ "beta-hydroxybutyrate",
        Metabolite == "alpha-ketoglutaric acid" ~ "alpha-ketoglutarate",
        TRUE ~ Metabolite
    ))

mb <- mb %>%
    left_join(rename(csl, SourceType = Type, SourceCarbon = Carbon)) %>%
    left_join(rename(csl, Metabolite = Source, MetaboliteType = Type, MetaboliteCarbon = Carbon)) %>%
    left_join(stl)


# Subtract by blank. Use one single blank (no strain, no CS) for all
mb_blank <- mb %>%
    filter(Source == "none") %>%
    filter(Strain == "none") %>%
    select(Metabolite, BlankMetaboliteConc = MetaboliteConc)

metabolomics <- mb %>%
    left_join(mb_blank) %>%
    filter(Source != "none") %>%
    mutate(MetaboliteConc = MetaboliteConc - BlankMetaboliteConc) %>%
    # Set negative values to 0
    mutate(MetaboliteConc = ifelse(MetaboliteConc <0, 0, MetaboliteConc)) %>%
    # Also only use time point 48hr
    filter(Timepoint == 48) %>%
    select(Fermenter, Strain, SourceType, Source, SourceCarbon,
           MetaboliteType, Metabolite, MetaboliteCarbon, Replicate, MetaboliteConc)

write_csv(metabolomics, paste0(folder_data, "temp/21-metabolomics.csv"))

# 3. l matrix from the metabolites ----
# Byproduct measurement on glucose. Data from Sylvie
isolates_byproduct_time <- read_csv(paste0(folder_data, "raw/growth_rate/Estrela_2021_isolates_ph_OAs.csv"), col_types = cols()) %>%
    select(ID = SangerID, Time = time_hours, OD620, pH, Glucose_perc, acetate_mM, succinate_mM, lactate_mM)

# Leakiness
isolates_leakiness <- isolates_byproduct_time %>%
    select(ID, Time, gluConc = Glucose_perc, X_acetate = acetate_mM, X_succinate = succinate_mM, X_lactate = lactate_mM) %>%
    mutate(X_sum = 2 * X_acetate + 3 * X_lactate + 4 * X_succinate) %>%
    # Convert the glucose from mass fraction to mM per carbon
    mutate(gluConsumption = (0.2 - gluConc) / 180 / 0.1 * 1000 * 6) %>%
    mutate(leakiness = X_sum / gluConsumption) %>%
    filter(Time != 0) %>%
    select(ID, Time, leakiness) %>%
    pivot_wider(names_from = Time, values_from = leakiness, names_glue = "leakiness_{Time}hr") %>%
    # Replace the leakiness of 0 glucose consumption time point by NA
    mutate_all(function(x) ifelse(is.infinite(x), NA, x))

isolates_leakiness <- isolates_ID %>%
    left_join(isolates_leakiness) %>%
    drop_na()

write_csv(isolates_leakiness, paste0(folder_data, "temp/21-isolates_leakiness.csv"))













