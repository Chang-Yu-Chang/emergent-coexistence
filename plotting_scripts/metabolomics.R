library(tidyverse)
library(janitor)
source(here::here("plotting_scripts/network_functions.R"))
mb <- read_csv("~/Dropbox/lab/invasion-network/data/raw/metabolomics/TMIC_LCMS.csv")

# Clean up names
mb <- mb %>%
    rename_with(~str_replace(., "_", "")) %>%
    remove_empty() %>%
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

# Metabolite types
csl <- tibble(Source = c("glucose", "fructose", "galactose", "ribose", "arabinose",
                         "gluconate", "ketogluconate",
                         "malate", "acetate", "glycerol", "pyruvate", "succinate", "fumarate",
                         "glycine", "acetyl-ornithine", "alanine", "succinate", "lactate",
                         "serine", "valine", "pyruvate", "alpha-ketoglutarate", "putrescine",
                         "valerate", "asparagine", "butyrate", "methionine", "hippurate",
                         "propionate", "citrate", "acetate", "beta-hydroxybutyrate", "fumarate",
                         "glutamine", "leucine", "oxoglutarate", "phenylalanine"),
              Type = c(rep("sugar", 5+2), rep("acid", 6+20+4)),
              Carbon = c(6, 6, 6, 5, 5,
                         6, 6,
                         4, 2, 3, 3, 4, 4,
                         2, 7, 3, 4, 3,
                         3, 5, 3, 5, 4,
                         5, 4, 4, 5, 9,
                         3, 6, 2, 4, 4,
                         5, 6, 5, 9)) %>%
    distinct(Source, Carbon, Type)
# Strain type
stl <- tibble(Strain = c("none", "Ecoli", "Enterobacter", "Pputida", "Pseudomonas"),
              Fermenter = c("none", "fermenter", "fermenter", "respirator", "respirator"))

mb <- mb %>%
    left_join(rename(csl, SourceType = Type, SourceCarbon = Carbon)) %>%
    left_join(rename(csl, Metabolite = Source, MetaboliteType = Type, MetaboliteCarbon = Carbon)) %>%
    left_join(stl)


# Subtract by blank. Use one single blank (no strain, no CS) for all
mb_blank <- mb %>%
    filter(Source == "none") %>%
    filter(Strain == "none") %>%
    select(Metabolite, BlankMetaboliteConc = MetaboliteConc)

mb <- mb %>%
    left_join(mb_blank) %>%
    filter(Source != "none") %>%
    mutate(MetaboliteConc = MetaboliteConc - BlankMetaboliteConc) %>%
    # Set negative values to 0
    mutate(MetaboliteConc = ifelse(MetaboliteConc <0, 0, MetaboliteConc)) %>%
    # Also only use time point 48hr
    filter(Timepoint == 48) %>%
    select(Fermenter, Strain, SourceType, Source, SourceCarbon,
           MetaboliteType, Metabolite, MetaboliteCarbon, Replicate, MetaboliteConc)

write_csv(mb, here::here("data/output/metabolomics.csv"))



#
p <- mb %>%
    #filter(Strain == "Ecoli") %>%
    mutate(Source = ordered(Source, csl$Source)) %>%
    mutate(Metabolite = ordered(Metabolite, csl$Source)) %>%
    mutate(MetaboliteConc = log(MetaboliteConc)) %>%
    ggplot() +
    geom_tile(aes(x = Source, y = Metabolite, fill = MetaboliteConc)) +
    facet_grid(.~Strain, scales = c("free_x")) +
    #facet_grid(SourceType~Strain, scales = c("free_x")) +
    #facet_grid(SourceType~MetaboliteType, scales = c("free_x")) +
    scale_fill_gradient(low = "white", high = "blue") +
    #scale_fill_gradient(low = "lightblue", high = "red") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5))

ggsave(here::here("plots/FigS6-metabolomics.png"), p, width = 15, height = 5)


## Calculate the fraction of conversion
mb %>%
    # Use only two strains. Check if treating the E and P isolates as replicates
    filter(Strain %in% c("Ecoli", "Pputida")) %>%
    # For one fermenter using one CS, what is the fraction of metabolic fluxes into acids and the fluxes into sugars ?
    group_by(Fermenter, SourceType, Source, MetaboliteType) %>%
    summarize(MetaboliteConc = sum(MetaboliteConc)) %>%
    mutate(RelativeMetaboliteConc = MetaboliteConc / sum(MetaboliteConc)) %>%
    # Take the average within acid-type CS, and within sugar-type CS
    group_by(Fermenter, SourceType, MetaboliteType) %>%
    summarize(RelativeMetaboliteConc = mean(RelativeMetaboliteConc))




# Estimate the uptake rates
isolates <- read_csv(here::here("data/output/isolates.csv"))

isolates_u <- isolates %>%
    drop_na() %>%
    mutate(Fermenter = ifelse(Fermenter, "fermenter", "respirator")) %>%
    distinct(ExpID, .keep_all = T) %>%
    select(Fermenter, ExpID, starts_with("r") & ends_with("16hr"), starts_with("OD") & ends_with("16hr")) %>%
    pivot_longer(cols = ends_with("16hr")) %>%
    separate(name, into = c("name", "CS", "Time")) %>%
    mutate(Time = str_replace(Time, "hr", "") %>% as.numeric) %>%
    pivot_wider(names_from = name, values_from = value) %>%
    # Calculate uptake rate
    #mutate(w = OD / Time, u = r / w) %>%
    mutate(u = r * OD) %>%
    arrange(Fermenter, CS) %>%
    rename(Source = CS) %>%
    left_join(csl) %>%
    filter(Source != "ddh20")


p <- isolates_u %>%
    ggplot() +
    #geom_histogram(aes(x = u, fill = Fermenter), color = 1) +
    geom_histogram(aes(x = u, fill = Fermenter), color = 1, position="identity", alpha = 0.5) +
    facet_grid(Fermenter~Type) +
    scale_fill_manual(values = category_color, breaks = c("fermenter", "respirator")) +
    theme_classic() +
    theme(legend.title = element_blank(), panel.border = element_rect(color = 1, fill = NA))

ggsave(here::here("plots/FigS7-uptake_rate.png"), p, width = 6, height = 6)


##
isolates_u %>%
    group_by(Fermenter, Type) %>%
    summarize(uMean = mean(u), uSd = sd(u))



# Leakiness
isolates %>%
    drop_na() %>%
    mutate(Fermenter = ifelse(Fermenter, "fermenter", "respirator")) %>%
    distinct(ExpID, .keep_all = T) %>%
    select(Fermenter, ExpID, starts_with("leakiness") & ends_with("16hr")) %>%
    group_by(Fermenter) %>%
    summarize(lMean = mean(leakiness_16hr), lSd = sd(leakiness_16hr))
























