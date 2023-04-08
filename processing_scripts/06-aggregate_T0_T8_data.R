#' This script convert the T0 OD to CFU and read the image processing T8 data
#' 1. Calculate epsilon for T0, and convert T0 OD to CFU. Output temp/06-isolate_epsilon.csv and temp/06-pairs_T0.csv
#' 2. Aggregate T8 CFU frequency from random forest outputs and clean column names. Output temp/06-pairs_T8.csv

library(tidyverse)
library(cowplot)
source(here::here("processing_scripts/00-metadata.R"))

pairs_freq_ID <- read_csv(paste0(folder_data, "temp/00c-pairs_freq_ID.csv"), show_col_types = F)

# 1. Compute epsilon to convert T0 OD to CFU frequencies ----
# 1.1 Read CFU data ----
# The CFU data is from the machine counted images
temp <- rep(list(NA), length(batch_names))
for (j in 1:length(batch_names)) {
    image_names <- list.files(paste0(folder_pipeline, "images/", batch_names[j], "-07-feature/green")) %>%
        # Remove folders
        str_subset(".csv") %>%
        #str_subset("T8") %>%
        str_replace(".csv", "") %>%
        # Remove all that contain _-, which was a naming convention for the different dilution factors on the same plate
        str_subset("^((?!_-).)*$")
    name_length <- image_names %>% str_split("_") %>% sapply(length)
    index_monoculture <- name_length == 4
    image_names <- image_names[index_monoculture]

    temp[[j]] <- tibble(image_name = image_names) %>%
        rowwise() %>%
        mutate(ColonyCount = paste0(folder_pipeline, "images/", batch_names[j], "-07-feature/green/", image_name, ".csv") %>%
                   read_csv(show_col_types = F) %>%
                   nrow()) %>%
        ungroup()
}
isolates_CFU <- bind_rows(temp) %>%
    separate(image_name, into = c("Batch", "Time", "Community", "Isolate"), sep = "_", remove = F)

isolates_CFU <- isolates_CFU %>%
    # This image is used to train the model because it has many colonies, but it does not have OD data
    filter(image_name != "C_stock_C11R1_1") %>%
    # To use the OD data, add the image C_T0_C11R1_1 that has OD data but does not have high resolution. Manually count colony
    bind_rows(tibble(image_name = "C_T0_C11R1_1", Batch = "C", Time = "T0", Community = "C11R1", Isolate = "1", ColonyCount = 8))


# 1.2 Read OD data ----
# Remove the contaminated data: B2 C11R1 isolate1 and C2 C11R2 isolate 13 (streplococcus contamination)
OD_B2 <- read_csv(paste0(folder_data, "raw/OD/OD_B2.csv"), col_types = cols()) %>%
    mutate(Batch = "B2", Layout = as.character(Layout), Isolate1Freq = as.character(Isolate1Freq), Isolate2Freq = as.character(Isolate2Freq)) %>%
    filter(!((Community == "C11R1" & Isolate1 == 1) | (Community == "C11R1" & Isolate2 == 1)))
OD_C <- read_csv(paste0(folder_data, "raw/OD/OD_C.csv"), col_types = cols()) %>%
    mutate(Batch = "C", Isolate1 = as.character(Isolate1), Isolate2 = as.character(Isolate2)) %>%
    mutate(Isolate1Freq = as.character(Isolate1Freq), Isolate2Freq = as.character(Isolate2Freq))
OD_C2 <- read_csv(paste0(folder_data, "raw/OD/OD_C2.csv"), col_types = cols()) %>%
    mutate(Batch = "C2", Isolate1 = as.character(Isolate1), Isolate2 = as.character(Isolate2)) %>%
    filter(!((Community == "C11R2" & Isolate1 == 13) | (Community == "C11R2" & Isolate2 == 13)))
OD_D <- read_csv(paste0(folder_data, "raw/OD/OD_D.csv"), col_types = cols()) %>%
    mutate(Batch = "D", Isolate1 = as.character(Isolate1), Isolate2 = as.character(Isolate2))


# Join OD data from different batches
OD <- bind_rows(OD_B2, OD_C, OD_C2, OD_D) %>%
    # Select essential variables for analysis. Remove rows of experimental blanks
    select(Batch, Well, Community, Time = Transfer, Isolate1, Isolate2, Isolate1Freq, Isolate2Freq, MixIsolate, MixPlate, Wavelength, Abs) %>%
    filter(Isolate1 != "blank") %>%
    mutate(Community = factor(Community, paste0("C", rep(1:12, each = 8), "R", rep(1:8, 12)))) %>%
    arrange(Community, Time, Isolate1, Isolate2, Isolate1Freq) %>%
    filter(Wavelength == 620) %>%
    rename(OD620 = Abs)

# Subset single culture that are not mixed
isolates_OD <- OD %>%
    filter(Isolate1 == Isolate2, MixIsolate == F) %>%
    select(Batch, Time, Community, Isolate = Isolate1, OD620) %>%
    group_by(Batch, Time, Community, Isolate) %>%
    summarize(OD620 = mean(OD620)) %>%
    ungroup() %>%
    arrange(Community) %>%
    mutate(Time = paste0("T", Time)) %>%
    unite(col = "image_name", Batch, Time, Community, Isolate, sep = "_") %>%
    # Set negative OD values to 0 or it causes issues in division
    mutate(OD620 = ifelse(OD620 < 0, 0, OD620))

# 1.3 Match the monoculture OD and CFU data ----
# Manually key in dilution factors
isolates_OD_CFU <- isolates_CFU %>%
    left_join(isolates_OD, by = "image_name") %>%
    # For C11R1, use the better plate in either batch B2 and C
    filter(!(Batch == "C" & Community == "C11R1" & Time == "T8")) %>%
    filter(!(image_name %in% c("B2_T0_C11R1_2", "B2_T1_C11R1_8", "B2_T0_C11R1_9"))) %>%
    mutate(Community = factor(Community, paste0("C", rep(1:12, each = 8), "R", rep(1:8, 12)))) %>%
    arrange(Community, Isolate) %>%
    # For T0 monoculture (inoculum), OD was standardized to OD=0.1
    mutate(OD620 = ifelse(Time == "T0", 0.1, OD620)) %>%
    # Dilution factor for all plating practice is 10^-5
    mutate(DilutionFactor = 10^-5)

# 1.4 Calculate the OD-CFU conversion coefficient epsilon for each isolate ----
#' \epsilon_A = \frac{CFU_A}{OD_A DF_A V_A}, where
#' CFU is the colony count at a given dilution factor. Obtained from image processing
#' DF is the dilution factor for plating (10^-5)
#' OD is the optical density either measured for T1 and T8, or set to 0.1 for T0
#' V is the plating volume (20 uL)
isolates_epsilon <- isolates_OD_CFU %>%
    mutate(
        CFU = ColonyCount,
        OD = OD620,
        V = 20, # Plating volume is 20 uL for all plating practice
        Epsilon = CFU/(OD*DilutionFactor*V)
    ) %>%
    mutate(Isolate = as.numeric(Isolate))
write_csv(isolates_epsilon, paste0(folder_data, "temp/06-isolates_epsilon.csv"))

# 1.5 Convert T0 OD to CFU ----
pairs_T0 <- pairs_freq_ID %>%
    left_join(select(isolates_epsilon, Community, Isolate1 = Isolate, Epsilon1 = Epsilon), by = c("Community", "Isolate1")) %>%
    left_join(select(isolates_epsilon, Community, Isolate2 = Isolate, Epsilon2 = Epsilon), by = c("Community", "Isolate2"))

#' The idea is to convert OD (0.1*Isolate1InitialODFreq) of type A to CFU (cfu_A),
#' and use cfu_A as the mean the parameterize poisson. Draw n_A ~ Pois(cfu_A).
#' Repeat it for type B and obtain n_B. The bootstrapped freq_A = n_A / (n_A + n_B)
#' here I am using dilution factor 10^-4 such that the CFU count is no smaller than 10
pairs_T0 <- pairs_T0 %>%
    mutate(cfu_A = 0.1 * (Isolate1InitialODFreq / 100) * 10^(-4) * 20 * Epsilon1,
           cfu_B = 0.1 * (Isolate2InitialODFreq / 100) * 10^(-4) * 20 * Epsilon2) %>%
    select(Batch, Community, Isolate1, Isolate2, Isolate1InitialODFreq, Isolate2InitialODFreq, cfu_A, cfu_B)

write_csv(pairs_T0, paste0(folder_data, "temp/06-pairs_T0.csv"))



# 2. Aggregate T8 predicted results ----
temp <- rep(list(NA), nrow(pairs_freq_ID))
for (i in 1:nrow(pairs_freq_ID)) {
    # Skip images with no colony
    if (pairs_freq_ID$image_name_pair[i] %in% plates_no_colony) {cat("\nno colony\t", pairs_freq_ID$image_name_pair[i]); next}
    # Skip cocultures with no plate images
    if (is.na(pairs_freq_ID$image_name_pair[i])) {cat("\nno colony\t", pairs_freq_ID$image_name_pair[i]); next}
    temp[[i]] <- read_csv(paste0(folder_pipeline, "images/", pairs_freq_ID$Batch[i],"-08-random_forest/", pairs_freq_ID$image_name_pair[i], ".csv"), show_col_types = F)
    cat("\n", pairs_freq_ID$image_name_pair[i], i, "/", nrow(pairs_freq_ID))
}


pairs_T8 <- bind_rows(temp[which(!is.na(temp))]) %>%
    mutate(Group = factor(Group, c("predicted isolate1", "predicted isolate2"))) %>%
    group_by(image_name, Group, .drop = F) %>%
    count(name = "Count") %>%
    group_by(image_name) %>%
    pivot_wider(names_from = Group, values_from = Count) %>%
    rename(Isolate1Count = `predicted isolate1`, Isolate2Count = `predicted isolate2`) %>%
    mutate(TotalCount = Isolate1Count + Isolate2Count,
           Isolate1CFUFreq = Isolate1Count/TotalCount) %>%
    rename(image_name_pair = image_name) %>%
    ungroup() %>%
    # Correct the isolate order
    separate(col = image_name_pair, into = c("Batch", "Time", "Community", "Isolate2InitialODFreq", "Isolate1InitialODFreq", "Isolate1", "Isolate2"), remove = F, convert = T) %>%
    rowwise() %>%
    mutate(Isolate1InitialODFreq = ifelse(Isolate1 > Isolate2, 5, Isolate1InitialODFreq),
           Isolate2InitialODFreq = ifelse(Isolate1 > Isolate2, 95, Isolate2InitialODFreq),
           Isolate1Count = ifelse(Isolate1 > Isolate2, TotalCount-Isolate1Count, Isolate1Count),
           Isolate1CFUFreq = ifelse(Isolate1 > Isolate2, 1-Isolate1CFUFreq, Isolate1CFUFreq),
           FlipOrder = ifelse(Isolate1 > Isolate2, T, F)
    ) %>%
    mutate(temp = min(Isolate1,Isolate2), Isolate2 = max(Isolate1, Isolate2), Isolate1 = temp) %>%
    select(-temp) %>%
    mutate(Community = factor(Community, paste0("C", rep(1:12, each = 8), "R", rep(1:8, 12)))) %>%
    select(image_name_pair, Batch, Community, Isolate1, Isolate2,
           Isolate1InitialODFreq, Isolate2InitialODFreq, Time,
           Isolate1Count, TotalCount, Isolate1CFUFreq) %>%
    arrange(Community, Isolate1, Isolate2, Isolate1InitialODFreq) %>%
    ungroup()

write_csv(pairs_T8, paste0(folder_data, "temp/06-pairs_T8.csv"))




















