#' This script convert the T0 OD to CFU and bootstrapping both T0 andT8 frequencies
#' 1. Calculate epsilon for T0, and convert T0 OD to CFU
#' 2. Bootstrap T0 CFU frequency. Output meta/06-pairs_T0_boots.csv
#' 3. Aggregate T8 CFU frequency from random forest outputs and clean column names. Output meta/06-pairs_T8.csv
#' 4. Bootstrap T8 CFU frequencies. Output meta/06-pairs_T8_boots.csv

library(tidyverse)
library(cowplot)
source(here::here("analysis/00-metadata.R"))

pairs_freq_ID <- read_csv(paste0(folder_data, "temp/00c-pairs_freq_ID.csv"), show_col_types = F)

# 1. Compute epsilon to convert T0 OD to CFU frequencies ----
# 1.1 Read CFU data ----
#' The CFU data is from the machine counted images. Then match them to the OD data
#' For T0 monoculture plates, there is no OD measurement data. Instead the OD should be standardized to OD=0.1
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

# Manually C2 C11R1 isolate 1
isolates_CFU <- isolates_CFU %>%
    filter(image_name != "C_stock_C11R1_1") %>%
    bind_rows(tibble(image_name = "C_T0_C11R1_1", Batch = "C", Time = "T0", Community = "C11R1", Isolate = "1", ColonyCount = 8))


# 1.2 Read OD data ----
# Remove the contaminated data: B2 community C11R1 isolate1 and C2 C11R2 isolate 13 (streplococcus contamination)
OD_B2 <- read_csv(paste0(folder_data, "raw/OD/OD_B2.csv"), col_types = cols()) %>%
    mutate(Batch = "B2", Layout = as.character(Layout), Isolate1Freq = as.character(Isolate1Freq), Isolate2Freq = as.character(Isolate2Freq)) %>%
    filter(!((Community == "C11R1" & Isolate1 == 1) | (Community == "C11R1" & Isolate2 == 1)))
OD_C <- read_csv(paste0(folder_data, "raw/OD/OD_C.csv"), col_types = cols()) %>% # Data from batch C C11R1 plate is also included
    mutate(Batch = "C", Isolate1 = as.character(Isolate1), Isolate2 = as.character(Isolate2)) %>%
    mutate(Isolate1Freq = as.character(Isolate1Freq), Isolate2Freq = as.character(Isolate2Freq))
OD_C2 <- read_csv(paste0(folder_data, "raw/OD/OD_C2.csv"), col_types = cols()) %>% # Data from batch C C11R1 plate is also included
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
    #' I do not have the T8 OD data for batch C. Remove them
    filter(!(Batch == "C" & Community == "C11R1" & Isolate != 1)) %>%
    filter(!(Batch == "B2" & Community == "C11R1" & Isolate == 1)) %>%
    # Remove staph contamination
    filter(!(Community == "C11R2" & Isolate == 13)) %>%
    mutate(Community = factor(Community, paste0("C", rep(1:12, each = 8), "R", rep(1:8, 12)))) %>%
    mutate(Isolate = factor(Isolate, 1:13)) %>%
    arrange(Community, Isolate) %>%
    # For T0 monoculture, I do not have OD data, but instead they are standardized to OD=0.1, so assign 0.1 to them
    mutate(OD620 = ifelse(Time == "T0", 0.1, OD620)) %>%
    # Dilution factor for these plates are 5
    mutate(DilutionFactor = 5)

# 1.4 Calculate the OD-CFU conversion coefficient epsilon for each isolate ----
#' $\epsilon_A = \frac{OD_A DF_A v_A}{CFU_A}$, where
#' D is the dilution factor for plating (for instance, $10^{5}$),
#' CFU is the colony count under given dilution factor,
#' OD is the optical density measured at transfer 8 after 48 hours,
#' $v$ is the plating volume (20 uL)
isolates_epsilon <- isolates_OD_CFU %>%
    mutate(
        CFU = ColonyCount, OD = OD620, V=20, # Plating volume
        V1 = 10, V2 = 90, # Volume used in serial dilution
        ErrorCFU = sqrt(CFU), ErrorOD = 0.001, ErrorV = 0.4,
        ErrorV1 = 0.4, ErrorV2 = 2
    ) %>%
    mutate( # Variable uncertainties
        n = DilutionFactor, # 4 or 5. shortened to n for convenience
        PartialV1 = n * V1^(n-1) * V2 / (V1+V2)^(n+1), #-n*(V1+V2)^(n-1)*V2/(V1^(n+1)),
        PartialV2 = -n * V1^(n-1) / (V1+V2)^(n+1), #n*(V1+V2)^(n-1)/(V1^n),
        DF = (V1/(V1+V2))^n,
        ErrorDF = sqrt((PartialV1*ErrorV1)^2 + (PartialV2*ErrorV2)^2)
    ) %>%
    select(-n, -V1, -V2, -ErrorV1, -ErrorV2, -PartialV1, -PartialV2) %>% # Remove transient variables
    mutate( # Partial derivatives
        PartialCFU = 1/(OD*DF*V),
        PartialOD = -CFU/(OD^2*DF*V),
        PartialDF = -CFU/(OD*DF^2*V),
        PartialV = -CFU/(OD*CFU*V^2)
    ) %>%
    mutate(
        Epsilon = CFU/(OD*DF*V),
        ErrorEpsilon = sqrt((PartialDF*ErrorDF)^2 + (PartialCFU*ErrorCFU)^2 + (PartialOD*ErrorOD)^2 + (PartialV*ErrorV)^2)
    ) %>%
    # Some filters
    # Batch B2 C11R1 Is contaminated, so set its epsilon to NA
    mutate(Epsilon = ifelse((Batch == "B2" & Community == "C11R1" & Isolate == 1), NA, Epsilon),
           ErrorEpsilon = ifelse((Batch == "B2" & Community == "C11R1" & Isolate == 1), NA, ErrorEpsilon)) %>%
    # Set epsilon to NA if OD <= 0 or CFU <= 10
    rowwise() %>%
    mutate(Epsilon = ifelse(OD <= 0 | CFU <= 0, NA, Epsilon),
           ErrorEpsilon = ifelse(OD <= 0 | CFU <= 0, NA, ErrorEpsilon)) %>%
    # Set epsilon to NA if colony counts <= 10
    mutate(Epsilon = ifelse(OD <= 0 | CFU <= 0, NA, Epsilon),
           ErrorEpsilon = ifelse(OD <= 0 | CFU <= 0, NA, ErrorEpsilon)) %>%
    #select(image_name, Batch, Time, Community, Isolate, ColonyCount, Epsilon, ErrorEpsilon) %>%
    mutate(Isolate = as.numeric(Isolate))
write_csv(isolates_epsilon, paste0(folder_data, "temp/06-isolates_epsilon.csv"))


# 1.5 Convert T0 OD to CFU ----
pairs_epsilon <- pairs_freq_ID %>%
    left_join(select(isolates_epsilon, Community, Isolate1 = Isolate, Epsilon1 = Epsilon, ErrorEpsilon1 = ErrorEpsilon), by = c("Community", "Isolate1")) %>%
    left_join(select(isolates_epsilon, Community, Isolate2 = Isolate, Epsilon2 = Epsilon, ErrorEpsilon2 = ErrorEpsilon), by = c("Community", "Isolate2"))

#' The idea is to convert OD (0.1*Isolate1InitialODFreq) of type A to CFU (cfu_A),
#' and use cfu_A as the mean the parameterize poisson. Draw n_A ~ Pois(cfu_A).
#' Repeat it for type B and obtain n_B. The bootstrapped freq_A = n_A / (n_A + n_B)
#' here I am using dilution factor 10^-4 such that the CFU count is no smaller than 10
pairs_epsilon <- pairs_epsilon %>%
    mutate(cfu_A = 0.1 * (Isolate1InitialODFreq / 100) * 10^(-4) * 20 * Epsilon1,
           cfu_B = 0.1 * (Isolate2InitialODFreq / 100) * 10^(-4) * 20 * Epsilon2)

# 2. Bootstrap T0 freq_A from Poisson ----
n_bootstraps = 1000
set.seed(1)
pairs_T0_boots <- pairs_epsilon %>%
    mutate(Time = "T0", RawDataType = "ODtoCFU") %>%
    rowwise() %>%
    mutate(bootstrap = list(
        tibble(BootstrapID = 1:n_bootstraps,
               n_A = rpois(n_bootstraps, cfu_A),
               n_B = rpois(n_bootstraps, cfu_B),
               Isolate1CFUFreq = n_A / (n_A + n_B))
    )) %>%
    unnest(cols = c(bootstrap))

write_csv(pairs_T0_boots, paste0(folder_data, "temp/06-pairs_T0_boots.csv"))


# 3. Aggregate T8 predicted results ----
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
    #separate(col = Freqs, sep = "-", into = c("Isolate2InitialODFreq", "Isolate1InitialODFreq"), convert = T) %>%
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

# 4. Bootstrap T8 freq_A from Poisson -----
temp <- rep(list(NA), nrow(pairs_T8))
for (i in 1:nrow(pairs_T8)) {
    bootstrap_freq <- function (total_count, cfu_freq) sum(runif(total_count, 0, 1) < cfu_freq)/total_count
    cfu_freqs <- NULL
    for (k in 1:n_bootstraps)  cfu_freqs[k] <- bootstrap_freq(pairs_T8$TotalCount[i], pairs_T8$Isolate1CFUFreq[i])

    temp[[i]] <- tibble(
        Batch = pairs_T8$Batch[i],
        Community = pairs_T8$Community[i],
        Isolate1 = pairs_T8$Isolate1[i],
        Isolate2 = pairs_T8$Isolate2[i],
        Isolate1InitialODFreq = pairs_T8$Isolate1InitialODFreq[i],
        Isolate2InitialODFreq = pairs_T8$Isolate2InitialODFreq[i],
        Time = "T8",
        RawDataType = "CFU",
        BootstrapID = 1:n_bootstraps,
        Isolate1CFUFreq = cfu_freqs
    )
    cat("\t", i)
}

pairs_T8_boots <- bind_rows(temp)
write_csv(pairs_T8_boots, paste0(folder_data, "temp/06-pairs_T8_boots.csv"))




















