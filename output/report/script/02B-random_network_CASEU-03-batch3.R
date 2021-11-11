#' Read the CASEU result of random network
#' In this CASEU RN batch 3 (CASEU_RN3), there are three plates sent for sequencing
#' The raw Sanger data are in three different folders
#' - P1: T0 C P2
#' - P2: T0 AD P2
#' - P3: T3 AD P2
#'
#'
#' 1. Read the plate layout
#' 2. Read the isolates trace matrices from the raw Sanger sequences
#' 3. Read the mixture trace matrices. Map the plate layout
#' 4. CASEU predicts the relative abundance
library(tidyverse)
library(data.table)
library(CASEU)
`%notin%` <- Negate(`%in%`)
source(here::here("output/report/script/misc.R"))

# Read plate layout ----
plates_random <- fread(here::here("data/output/plates_random.csv")) %>% as_tibble()
plates_C_P2 <- plates_random %>% filter(PlateLayout == "C", MixPlate == "P2") %>% mutate(Sample = 1:96)
plates_AD_P2 <- plates_random %>% filter(PlateLayout == "AD", MixPlate == "P2") %>% mutate(Sample = 1:96)

# Read isolates trace matrices ----
#' Object traces_isolate is a R list with all 32 isolates used in
#' random networks. The CASEU data is obtained from CASEU_RN2
load(here::here("data/temp/CASEU_RN_isolates_trace.Rdata"))

# Read mixtures trace matrices ----
## Folder of raw Sanger sequence data
folder_directory_RN3_P1 <- "data/raw/Sanger/CASEU_RN3/30-333649143/Caseu_RN_3_30-333649143_ab1/"
folder_directory_RN3_P2 <- "data/raw/Sanger/CASEU_RN3/30-333649365/Caseu_RN_3_30-333649365_ab1/"
folder_directory_RN3_P3 <- "data/raw/Sanger/CASEU_RN3/30-333649517/Caseu_RN_3_30-333649517_ab1/"

## Function for reading trace matrices in different folders
read_trace_temp <- function(folder, name) {
  paste0(folder, name, "-27F.ab1") %>%
    here::here() %>%
    sangerseqR::read.abif() %>%
    CASEU::extractElectropherogram()
}

## Create a R list for reading muxture traces
traces_mixture <- rep(list(NA), 96 * 3)
names(traces_mixture) <- 1:(96*3)

for (i in 1:(96*3)){
  if (i <= 96) traces_mixture[[i]] <- read_trace_temp(folder_directory_RN3_P1, i)
  if (i > 96 & i <= 192) traces_mixture[[i]] <- read_trace_temp(folder_directory_RN3_P2, i)
  if (i > 192 & i <= 288) traces_mixture[[i]] <- read_trace_temp(folder_directory_RN3_P3, i)
}


## Read the mixture layout on the plate
temp_df1 <- plates_C_P2 %>%
  mutate(Isolate1Freq = Isolate1Freq * 100, Isolate2Freq = Isolate2Freq * 100) %>%
  mutate(Time = "T0") %>%
  unite("Mixture", c("Community", "Isolate1Freq", "Isolate2Freq", "Isolate1", "Isolate2")) %>%
  unite("Plate", c("Time", "PlateLayout", "MixPlate")) %>%
  select(Sample, Plate, Mixture)

temp_df2 <- plates_AD_P2 %>%
  mutate(Isolate1Freq = Isolate1Freq * 100, Isolate2Freq = Isolate2Freq * 100) %>%
  mutate(Time = "T0") %>%
  unite("Mixture", c("Community", "Isolate1Freq", "Isolate2Freq", "Isolate1", "Isolate2")) %>%
  unite("Plate", c("Time", "PlateLayout", "MixPlate")) %>%
  mutate(Sample = Sample + 96) %>%
  select(Sample, Plate, Mixture)

temp_df3 <- plates_AD_P2 %>%
  mutate(Isolate1Freq = Isolate1Freq * 100, Isolate2Freq = Isolate2Freq * 100) %>%
  mutate(Time = "T3") %>%
  unite("Mixture", c("Community", "Isolate1Freq", "Isolate2Freq", "Isolate1", "Isolate2")) %>%
  unite("Plate", c("Time", "PlateLayout", "MixPlate")) %>%
  mutate(Sample = Sample + 192) %>%
  select(Sample, Plate, Mixture)

mixture_list <- bind_rows(temp_df1, temp_df2, temp_df3)
names(traces_mixture) <- mixture_list$Mixture

## Subset for only mixtures. sample 65-96 on plate T0 C P2 are isolates
mixture_list <- filter(mixture_list, Sample != 65:96)
traces_mixture <- traces_mixture[-(65:96)]


# CASEU prediction ----
## Separate the mixture
CASEU_RN3 <- mixture_list %>%
  # Separate the mixture column
  tidyr::separate(col = Mixture, sep = "_", remove = F, into = c("Community", "Isolate1Freq", "Isolate2Freq", "Isolate1", "Isolate2") ) %>%
  # Separate the Plate column
  tidyr::separate(col = Plate, sep = "_", remove = T, into = c("Time", "PlateLayout", "MixPlate")) %>%
  mutate(Isolate1 = as.numeric(as.character(Isolate1)),
    Isolate2 = as.numeric(as.character(Isolate2)),
    Isolate1Freq = as.numeric(as.character(Isolate1Freq))/100,
    Isolate2Freq = as.numeric(as.character(Isolate2Freq))/100) %>%
  select(Sample, Mixture, Time, Community, Isolate1, Isolate2, Isolate1Freq, Isolate2Freq)


# CASEU prediction. Loop through all experimental pairs/mixtures
## This may take ~45 minutes
caseu_prediction <- rep(list(NA), nrow(mixture_list))
names(caseu_prediction) <- mixture_list$Mixture
tt <- proc.time()

for (i in 1:length(caseu_prediction)) {
  # Community and isolates in a pair
  community <- CASEU_RN3$Community[i]
  isolate1 <- CASEU_RN3$Isolate1[i]
  isolate2 <- CASEU_RN3$Isolate2[i]

  # Trace matrices
  trace_mixture <- traces_mixture[[i]]
  trace_isolate1 <- traces_isolate[[match(paste0(community, "_", isolate1), names(traces_isolate))]]
  trace_isolate2 <- traces_isolate[[match(paste0(community, "_", isolate2), names(traces_isolate))]]

  if (isolate1 == isolate2) next # Skip single culture
  if (is.na(trace_mixture)) next # Skip if there is no trace matrix data available for the mixture
  if (is.na(trace_isolate1) || is.na(trace_isolate2)) next # Skip if either of the isolate trace is missed
  caseu_prediction[[i]] <- CASEU::fitSangerMixture( # CASEU package function
    mixture = trace_mixture, # Mixture trace
    knots = seq(1500, min(c(nrow(trace_mixture), nrow(trace_isolate1), nrow(trace_isolate2))), by=1500),
    components = list(trace_isolate1, trace_isolate2)
  )

  print(i)
  cat((proc.time() - tt)[3], "seconds\n")
}

### Read CASEU prediction outcomes
temp <- caseu_prediction %>%
  lapply(function(x) {
    data.frame(
      Isolate1FreqPredicted = unlist(x["frac"])[1],
      Isolate2FreqPredicted = unlist(x["frac"])[2],
      RSquare = unlist(x["r2"])[1])
  }) %>%
  rbindlist(idcol = "Mixture") %>%
  mutate(Sample = c(1:64, 97:288))

### Format the CASEU df
CASEU_RN3 <- CASEU_RN3 %>%
  left_join(temp, by = c("Sample", "Mixture")) %>%
  filter(Isolate1 != Isolate2) %>%
#  select(-Mixture) %>%
#  select(Time, Community, Isolate1, Isolate2, Isolate1Freq, Isolate2Freq, Isolate1FreqPredicted, Isolate2FreqPredicted) %>%
#  arrange(Isolate1, Isolate2, Isolate1Freq, Isolate2Freq)
  # Switch isolates so that isolate2 is always greater than isolate1
  switch_pairwise_column() %>%
  arrange(Sample) %>%
{.}



# Save raw CASEU output
CASEU_RN3_raw_output <- caseu_prediction # Raw CASEU output
CASEU_RN3_mixture_list <- mixture_list # list of the mixture
#save(CASEU_RN3_raw_output, CASEU_RN3_mixture_list, file = here::here("data/temp/CASEU_RN3_raw_output.Rdata"))

fwrite(CASEU_RN3, here::here("data/temp/CASEU_RN3.csv"))
