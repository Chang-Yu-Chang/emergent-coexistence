#' This script runs the analysus of CASEU pilot1 Sanger sequences raw data
library(tidyverse)
library(data.table)
library(CASEU)
read_trace_matrix <- function(abif.file) {
  x <- sangerseqR::read.abif(abif.file) %>% sangerseqR::sangerseq()
  return(x@traceMatrix)
}

folder_directory <- "data/raw/Sanger/CASEU_pilot2/30-280311837_ab1/" # Replace it with the correct folder directory
genewiz_pilot2 <- fread(here::here("output/protocol/tab_fig/protocol_20190910_Sanger_seq_prep-genewiz_table_CYC.csv"))

# Read trace matrices for isolates and mixtures from Sanger sequences ----
list_seq <- paste0(genewiz_pilot2$`DNA Name`)

## Isolates
isolates_names <- grep("single", list_seq, value =T) %>% gsub("-27F", "", .) %>% strsplit("_") %>% unlist() %>% matrix(ncol = 4, byrow = T) %>% `[`(,4)
list_isolates <- grep("single", list_seq, value = T)
trace_isolates <- rep(list(NA), length(list_isolates)) %>% setNames(isolates_names)
for (i in 1:length(list_isolates)) trace_isolates[[i]] <- read_trace_matrix(here::here(paste0(folder_directory, list_isolates[i], "-27F.ab1")))


## Mixtures
list_mix <- list_seq[!(list_seq %in% list_isolates)]
trace_mix <- rep(list(NA), length(list_mix))
for (i in 1:length(list_mix)) trace_mix[[i]] <- read_trace_matrix(here::here(paste0(folder_directory, list_mix[i], "-27F.ab1")))


# Make dataframe for CASEU ----
CASEU_pilot2 <- data.frame(MixtureName = list_mix)
temp <- CASEU_pilot2$MixtureName %>%
  gsub("-27F|-27F_R", "", .) %>%
  strsplit("_") %>%
  unlist() %>%
  matrix(ncol=5, byrow=T) %>%
  as.data.frame() %>%
  setNames(c("User", "Treatment", "Isolate1Freq", "Isolate2Freq", "Isolates")) %>%
  mutate(Isolate1 = substr(Isolates, 1, 1), Isolate2 = substr(Isolates, 2, 2),
    Isolate1Freq = as.numeric(as.character(Isolate1Freq))/100,
    Isolate2Freq = as.numeric(as.character(Isolate2Freq))/100) %>%
  select(-Isolates)
CASEU_pilot2 <- cbind(CASEU_pilot2, temp)


# Fit mixture Sanger electorpheogram by using CASEU packages. This may take a few minutes. ----
sanger_mixture <- rep(list(NA), length(list_mix))
for (i in 1:length(list_mix)) {
  sanger_mixture[[i]] <- CASEU::fitSangerMixture( # CASEU package function
    mixture = trace_mix[[i]], # Mixture data
    components = list(
      trace_isolates[[match(CASEU_pilot2$Isolate1[i], isolates_names)]], # Isolate 1
      trace_isolates[[match(CASEU_pilot2$Isolate2[i], isolates_names)]]) # Isolate 2
  )
  print(i)
}


# Read CASEU predicted fraction ----
for (i in 1:length(list_mix)) CASEU_pilot2[i,c("Isolate1FreqPredicted", "Isolate2FreqPredicted")] <- sanger_mixture[[i]]$frac
fwrite(CASEU_pilot2, here::here("data/temp/CASEU_pilot2.csv"))

# Save raw CASEU output
# CASEU_pilot2_raw_output <- sanger_mixture # Raw CASEU output
# CASEU_pilot2_list_mix <- list_mix # list of the mixture
# save(CASEU_pilot2_raw_output, CASEU_pilot2_list_mix, file = here::here("data/temp/CASEU_pilot2_raw_output.Rdata"))

