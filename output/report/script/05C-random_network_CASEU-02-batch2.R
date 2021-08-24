#' Read the CASEU result of random network
#' Batch 2: plate layout T3 C P2
#'
#' This script involve three parts
#' 1. Read the plate layout
#' 2. Read the trace matrices from the raw Sanger sequences
#' 3. CASEU predicts the relative abundance
library(tidyverse)
library(data.table)
library(CASEU)
`%notin%` <- Negate(`%in%`)
source(here::here("output/report/script/misc.R"))

# Read plate layout ----
plates_random <- fread(here::here("data/output/plates_random.csv")) %>% as_tibble()
plates_C_P2 <- plates_random %>% filter(PlateLayout == "C", MixPlate == "P2") %>% mutate(Sample = 1:96)

# Folder of raw Sanger sequence data
folder_directory_RN2 <- "data/raw/Sanger/CASEU_RN2/CASEU_RN_2_30-323299648_ab1/"
available_sanger <- list.files(here::here(paste0(folder_directory_RN2))) %>% gsub("-27F.ab1", "", .) %>% factor(1:96) %>%  sort() %>%  as.character()

# Read trace matrices for isolates and mixtures from Sanger sequences ----
## Isolates
isolates_list <- plates_C_P2 %>% filter(Isolate1 == Isolate2, is.na(Isolate1Freq)) %>%
  unite("Isolate", c("Community", "Isolate1")) %>%
  select(Sample, Isolate) %>%
  mutate(DataAvailability = Sample %in% available_sanger)

traces_isolate <- rep(list(NA), nrow(isolates_list))
names(traces_isolate) <- isolates_list$Isolate

for (i in 1:nrow(isolates_list)) {
  # Skip the non-available data
  if (isolates_list$Sample[i] %notin% available_sanger) next
  # Read trace
  traces_isolate[[i]] <-
    paste0(folder_directory_RN2, isolates_list$Sample[i], "-27F.ab1") %>%
    here::here() %>%
    sangerseqR::read.abif() %>%
    CASEU::extractElectropherogram()
}

## Mixtures
mixture_list <- plates_C_P2 %>% filter(MixIsolate == T) %>%
  mutate(Isolate1Freq = Isolate1Freq * 100, Isolate2Freq = Isolate2Freq * 100) %>%
  unite("Mixture", c("Community", "Isolate1Freq", "Isolate2Freq", "Isolate1", "Isolate2")) %>%
  select(Sample, Mixture) %>%
  mutate(DataAvailability = Sample %in% available_sanger)

traces_mixture <- rep(list(NA), nrow(mixture_list))
names(traces_mixture) <- mixture_list$Mixture

for (i in 1:nrow(mixture_list)) {
  # Skip the non-available data
  if (mixture_list$Sample[i] %notin% available_sanger) next
  # Read trace
  traces_mixture[[i]] <-
    paste0(folder_directory_RN2, mixture_list$Sample[i], "-27F.ab1") %>%
    here::here() %>%
    sangerseqR::read.abif() %>%
    CASEU::extractElectropherogram()
}



# CASEU prediction ----
### Parse mixture column
CASEU_RN2 <- mixture_list %>%
  tidyr::separate(col = Mixture, sep = "_", remove = F, into = c("Community", "Isolate1Freq", "Isolate2Freq", "Isolate1", "Isolate2") ) %>%
  mutate(
    Isolate1 = as.numeric(as.character(Isolate1)),
    Isolate2 = as.numeric(as.character(Isolate2)),
    Isolate1Freq = as.numeric(as.character(Isolate1Freq))/100,
    Isolate2Freq = as.numeric(as.character(Isolate2Freq))/100) %>%
  # Transfer
  mutate(Transfer = "T3")

#### Available sanger seuqencing result
#CASEU_RN2$DataAvailability <- CASEU_RN2$Community == "RanAss1" & !(CASEU_RN2$Isolate1 %in% 5 | CASEU_RN2$Isolate2 %in% 5 | CASEU_RN2$Isolate1 == CASEU_RN2$Isolate2)
CASEU_RN2$DataAvailability <- T

## Fit Sanger sequnece electorpheogram of mixture by using CASEU packages. This may take a few minutes.
### Keep raw CASEU outputs
caseu_prediction <- rep(list(NA), nrow(mixture_list))
names(caseu_prediction) <- mixture_list$Mixture
tt <- proc.time()

for (i in 1:length(caseu_prediction)) {
#load(file = here::here("data/temp/CASEU_RN2_raw_output.Rdata"))
#caseu_prediction <- CASEU_RN2_raw_output
#for (i in which(is.na(CASEU_RN2_raw_output) & CASEU_RN2$Isolate1 != CASEU_RN2$Isolate2)) {
  community <- CASEU_RN2$Community[i]
  isolate1 <- CASEU_RN2$Isolate1[i]
  isolate2 <- CASEU_RN2$Isolate2[i]

  trace_mixture <- traces_mixture[[i]]
  trace_isolate1 <- traces_isolate[[match(paste0(community, "_", isolate1), names(traces_isolate))]]
  trace_isolate2 <- traces_isolate[[match(paste0(community, "_", isolate2), names(traces_isolate))]]

  if (isolate1 == isolate2) next
  if (is.na(traces_mixture[[i]])) next
  if (CASEU_RN2$DataAvailability[i]) {
    caseu_prediction[[i]] <- CASEU::fitSangerMixture( # CASEU package function
      mixture = trace_mixture, # Mixture trace
      knots = seq(1500, min(c(nrow(trace_mixture), nrow(trace_isolate1), nrow(trace_isolate2))), by=1500),
      components = list(trace_isolate1, trace_isolate2)
    )
  }

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
  rbindlist(idcol = "Mixture")

### Format the CASEU df
CASEU_RN2 <- CASEU_RN2 %>%
  left_join(temp) %>%
  select(-Mixture, -DataAvailability, -Sample) %>%
  filter(!is.na(Isolate1FreqPredicted)) %>%
  select(Community, Isolate1, Isolate2, Isolate1Freq, Isolate2Freq, Isolate1FreqPredicted, Isolate2FreqPredicted) %>%
  arrange(Isolate1, Isolate2, Isolate1Freq, Isolate2Freq) %>%
  # Switch isolates so that isolate2 is always greater than isolate1
  switch_pairwise_column()


# Save raw CASEU output
save(traces_isolate, file = here::here("data/temp/CASEU_RN_isolates_trace.Rdata"))
CASEU_RN2_raw_output <- caseu_prediction # Raw CASEU output
CASEU_RN2_mixture_list <- mixture_list # list of the mixture
save(CASEU_RN2_raw_output, CASEU_RN2_mixture_list, file = here::here("data/temp/CASEU_RN2_raw_output.Rdata"))

fwrite(CASEU_RN2, file = here::here("data/temp/CASEU_RN2.csv"))

