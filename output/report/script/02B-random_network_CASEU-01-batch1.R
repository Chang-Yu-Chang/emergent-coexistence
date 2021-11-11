#' Read the CASEU result of random network
#' Batch 1: plate layout T3 C P2
#'
#' This script involve three parts
#' 1. Read the plate layout
#' 2. Read the trace matrices from the raw Sanger sequences
#' 3. CASEU predicts the relative abundance
library(tidyverse)
library(data.table)
library(CASEU)
`%notin%` <- Negate(`%in%`)
"
Important information before starting: The sample 3, 4, 5, 12, 13, 27, 29, 33, 73, 85, 86, 91, and 94
arrived Genewiz in dry tube, so the samples are not included in the batch.
"

# Read plate layout ----
plates_random <- fread(here::here("data/output/plates_random.csv")) %>% as_tibble()
plates_C_P2 <- plates_random %>% filter(PlateLayout == "C", MixPlate == "P2") %>% mutate(Sample = 1:96)

# Folder of raw Sanger sequence data
folder_directory_RN1 <- "data/raw/Sanger/CASEU_RN1/Caseu_RN_1_30-322411127_ab1/"
available_sanger <- list.files(here::here(paste0(folder_directory_RN1))) %>% gsub("-27F.ab1", "", .)

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
    paste0(folder_directory_RN1, isolates_list$Sample[i], "-27F.ab1") %>%
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
    paste0(folder_directory_RN1, mixture_list$Sample[i], "-27F.ab1") %>%
    here::here() %>%
    sangerseqR::read.abif() %>%
    CASEU::extractElectropherogram()
}


# CASEU prediction ----
### Parse mixture column
CASEU_RN1 <- mixture_list %>%
  tidyr::separate(col = Mixture, sep = "_", remove = F, into = c("Community", "Isolate1Freq", "Isolate2Freq", "Isolate1", "Isolate2") ) %>%
  mutate(
    Isolate1 = as.numeric(as.character(Isolate1)),
    Isolate2 = as.numeric(as.character(Isolate2)),
    Isolate1Freq = as.numeric(as.character(Isolate1Freq))/100,
    Isolate2Freq = as.numeric(as.character(Isolate2Freq))/100) %>%
  # Transfer
  mutate(Transfer = "T3")

#### Available sanger seuqencing result
CASEU_RN1$DataAvailability <- CASEU_RN1$Community == "RanAss1" & !(CASEU_RN1$Isolate1 %in% 5:6 | CASEU_RN1$Isolate2 %in% 5:6 | CASEU_RN1$Isolate1 == CASEU_RN1$Isolate2)

## Fit Sanger sequnece electorpheogram of mixture by using CASEU packages. This may take a few minutes.
### Keep raw CASEU outputs
caseu_prediction <- rep(list(NA), nrow(mixture_list))
names(caseu_prediction) <- mixture_list$Mixture
tt <- proc.time()

for (i in 1:length(caseu_prediction)) {
#for (i in 1:2) {
  community <- CASEU_RN1$Community[i]
  isolate1 <- CASEU_RN1$Isolate1[i]
  isolate2 <- CASEU_RN1$Isolate2[i]

  trace_mixture <- traces_mixture[[i]]
  trace_isolate1 <- traces_isolate[[match(paste0(community, "_", isolate1), names(traces_isolate))]]
  trace_isolate2 <- traces_isolate[[match(paste0(community, "_", isolate2), names(traces_isolate))]]

  if (isolate1 == isolate2) next
  if (is.na(traces_mixture[[i]])) next
  if (CASEU_RN1$DataAvailability[i]) {
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
switch_pairwise_column <- function (df, bypair = T) {
  if (any(is.factor(df$Isolate1))) df$Isolate1 <- as.numeric(df$Isolate1); df$Isolate2 <- as.numeric(df$Isolate2)
  if ("Isolate1FreqPredicted" %in% colnames(df)) {
    if (bypair == T) {
      temp_index <- df$Isolate1 > df$Isolate2
      df[temp_index, c("Isolate1", "Isolate2", "Isolate1Freq", "Isolate2Freq", "Isolate1FreqPredicted", "Isolate2FreqPredicted")] <-
        df[temp_index, c("Isolate2", "Isolate1", "Isolate2Freq", "Isolate1Freq", "Isolate2FreqPredicted", "Isolate1FreqPredicted")]

      df %>% arrange(Isolate1, Isolate2, Isolate1Freq) %>% return()
    } else if (bypair == F) {
      temp_index <- df$Isolate1Freq == 5
      df[temp_index, c("Isolate1", "Isolate2", "Isolate1Freq", "Isolate2Freq", "Isolate1FreqPredicted", "Isolate2FreqPredicted")] <-
        df[temp_index, c("Isolate2", "Isolate1", "Isolate2Freq", "Isolate1Freq", "Isolate2FreqPredicted", "Isolate1FreqPredicted")]

      df %>% arrange(Isolate1Freq, Isolate1, Isolate2) %>% return()
    }
  } else {

    if (bypair == T) {
      temp_index <- df$Isolate1 > df$Isolate2
      df[temp_index, c("Isolate1", "Isolate2", "Isolate1Freq", "Isolate2Freq")] <-
        df[temp_index, c("Isolate2", "Isolate1", "Isolate2Freq", "Isolate1Freq")]

      df %>% arrange(Isolate1, Isolate2, Isolate1Freq) %>% return()
    } else if (bypair == F) {
      temp_index <- df$Isolate1Freq == 5
      df[temp_index, c("Isolate1", "Isolate2", "Isolate1Freq", "Isolate2Freq")] <-
        df[temp_index, c("Isolate2", "Isolate1", "Isolate2Freq", "Isolate1Freq")]

      df %>% arrange(Isolate1Freq, Isolate1, Isolate2) %>% return()
    }
  }
}
CASEU_RN1 <- CASEU_RN1 %>%
  left_join(temp) %>%
  select(-Mixture, -DataAvailability, -Sample) %>%
  filter(!is.na(Isolate1FreqPredicted)) %>%
  select(Community, Isolate1, Isolate2, Isolate1Freq, Isolate2Freq, Isolate1FreqPredicted, Isolate2FreqPredicted) %>%
  arrange(Isolate1, Isolate2, Isolate1Freq, Isolate2Freq) %>%
  # Switch isolates so that isolate2 is always greater than isolate1
  switch_pairwise_column()


# Save raw CASEU output
save(traces_isolate, file = here::here("data/temp/CASEU_RN_isolates_trace.Rdata"))
CASEU_RN1_raw_output <- caseu_prediction # Raw CASEU output
CASEU_RN1_mixture_list <- mixture_list # list of the mixture
#save(CASEU_RN1_raw_output, CASEU_RN1_mixture_list, file = here::here("data/temp/CASEU_RN1_raw_output.Rdata"))


fwrite(CASEU_RN1, file = here::here("data/temp/CASEU_RN1.csv"))
