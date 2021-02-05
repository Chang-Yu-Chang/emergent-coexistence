#' This script reads Sylvie's pilot2 Sanger sequences (n=20) and applies CASEU
#' method to calculate the relative abundances of strains in mixtures.
library(tidyverse)
library(data.table)
library(CASEU)
read_trace_matrix <- function(abif.file) {
  x <- sangerseqR::read.abif(abif.file) %>% sangerseqR::sangerseq()
  return(x@traceMatrix)
}

# Read name list
folder_directory <- "data/raw/Sanger/CASEU_pilot2/30-280311837_ab1/" # Replace it with the correct folder directory
genewiz_pilot2_sylvie <- fread(here::here("output/protocol/tab_fig/caseu_pilot_2_samples_ID_SE.csv"))


# Read trace matrices for isolates and mixtures from Sanger sequences ----
list_seq <- paste0(genewiz_pilot2_sylvie$`DNA Name`, "-27F")

## Isolates
isolates_names <- c("K+", "K-", "P", "A")
list_isolates <- paste0("SE_", 1:4)
trace_isolates <- rep(list(NA), length(list_isolates)) %>% setNames(isolates_names)
for (i in 1:length(list_isolates)) trace_isolates[[i]] <- read_trace_matrix(here::here(paste0(folder_directory, list_isolates[i], "-27F.ab1")))


## Mixtures
list_mix <- paste0("SE_", 5:16)
trace_mix <- rep(list(NA), length(list_mix))
for (i in 1:length(list_mix)) trace_mix[[i]] <- read_trace_matrix(here::here(paste0(folder_directory, list_mix[i], "-27F.ab1")))


# Make dataframe for CASEU ----
CASEU_result_pilot2_sylvie <- data.frame(
  MixtureName = paste0("SE_", 5:16),
  Isolate1 = c("K+", "K+", "K+", "K-", "K-", "A", "K+", "K+", "K-", "K+", "K+", "K-"),
  Isolate2 = c("A", "P", "K-", "A", "P", "P", "A", "P", "K+", "A", "P", "K+"),
  Isolate1Freq = c(rep(50, 6), rep(95, 6)),
  Isolate2Freq = c(rep(50, 6), rep(5, 6))
)


# Fit mixture Sanger electorpheogram by using CASEU packages. This may take a few minutes. ----
sanger_mixture_sylvie <- rep(list(NA), length(list_mix))
for (i in 1:length(list_mix)) { # 1:length(list_mix)
  sanger_mixture_sylvie[[i]] <- CASEU::fitSangerMixture( # CASEU package function
    mixture = trace_mix[[i]], # Mixture data
    components = list(
      trace_isolates[[match(CASEU_result_pilot2_sylvie$Isolate1[i], isolates_names)]], # Isolate 1
      trace_isolates[[match(CASEU_result_pilot2_sylvie$Isolate2[i], isolates_names)]]) # Isolate 2
  )
  print(i)
}


# Read CASEU predicted fraction ----
for (i in 1:length(list_mix)) CASEU_result_pilot2_sylvie[i,c("Isolate1_freq_predicted", "Isolate2_freq_predicted")] <- sanger_mixture_sylvie[[i]]$frac
CASEU_result_pilot2_sylvie <- CASEU_result_pilot2_sylvie %>% arrange(Isolate1, Isolate2, Treatment)


CASEU_result_pilot2_sylvie %>%
  filter(!is.na(Isolate1_freq_predicted)) %>%
  fwrite(here::here("data/temp/CASEU_pilot2_sylvie.csv"))
