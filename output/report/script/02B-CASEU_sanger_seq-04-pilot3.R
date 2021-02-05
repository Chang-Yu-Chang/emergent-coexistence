#' This script runs the analysus of CASEU pilot3 Sanger sequences raw data
library(tidyverse)
library(data.table)
library(CASEU)
read_trace_matrix <- function(abif.file) {
  x <- sangerseqR::read.abif(abif.file) %>% sangerseqR::sangerseq()
  return(x@traceMatrix)
}
folder_directory <- "data/raw/Sanger/CASEU_pilot3/30-286536939_ab1/" # Replace it with the correct folder directory
genewiz_pilot3 <- fread(here::here("output/protocol/tab_fig/protocol_20190924_Sanger_seq_prep-genewiz_table_CYC.csv"))

# Read trace matrices for isolates and mixtures from Sanger sequences ----
list_seq <- paste0(genewiz_pilot3$`DNA Name`)

## Isolates
isolates_names <- list_seq[43:63]
trace_isolates <- rep(list(NA), length(isolates_names)) %>% setNames(isolates_names)
#for (i in 1:length(isolates_names)) trace_isolates[[i]] <- invnet::read_trace_matrix(here::here(paste0(folder_directory, isolates_names[i], "-27F.ab1")))
for (i in 1:length(isolates_names)) trace_isolates[[i]] <- CASEU::extractElectropherogram(sangerseqR::read.abif(here::here(paste0(folder_directory, isolates_names[i], "-27F.ab1"))))

## Mixtures
list_mix <- list_seq[1:42]
trace_mix <- rep(list(NA), length(list_mix))
#for (i in 1:length(list_mix)) trace_mix[[i]] <- invnet::read_trace_matrix(here::here(paste0(folder_directory, list_mix[i], "-27F.ab1")))
for (i in 1:length(list_mix)) trace_mix[[i]] <- CASEU::extractElectropherogram(sangerseqR::read.abif(here::here(paste0(folder_directory, list_mix[i], "-27F.ab1"))))



# Make dataframe for CASEU ----
CASEU_pilot3 <- data.frame(MixtureName = list_mix)
temp <- CASEU_pilot3$MixtureName %>%
  gsub("-27F|-27F_R", "", .) %>%
  strsplit("_") %>%
  unlist() %>%
  matrix(ncol=5, byrow=T) %>%
  as.data.frame() %>%
  setNames(c("Community", "Isolate1Freq", "Isolate2Freq", "Isolate1", "Isolate2")) %>%
  mutate(
    Isolate1 = as.numeric(as.character(Isolate1)),
    Isolate2 = as.numeric(as.character(Isolate2)),
    Isolate1Freq = as.numeric(as.character(Isolate1Freq))/100,
    Isolate2Freq = as.numeric(as.character(Isolate2Freq))/100)
CASEU_pilot3 <- cbind(CASEU_pilot3, temp)


# Fit mixture Sanger electorpheogram by using CASEU packages. This may take a few minutes. ----
sanger_mixture <- rep(list(NA), length(list_mix))
#for (i in c(20, 22, 24, 28, 37, 38, 40)) {
for (i in 1:length(sanger_mixture)) {
  sanger_mixture[[i]] <- CASEU::fitSangerMixture( # CASEU package function
    mixture = trace_mix[[i]], # Mixture data
    components = list(
      trace_isolates[[match(paste0(CASEU_pilot3$Community[i], "_", CASEU_pilot3$Isolate1[i]), isolates_names)]], # Isolate 1
      trace_isolates[[match(paste0(CASEU_pilot3$Community[i], "_", CASEU_pilot3$Isolate2[i]), isolates_names)]]), # Isolate 2
    tol = 0.01
  )
  print(i)
}


# Read CASEU predicted fraction ----
for (i in 1:length(list_mix)) CASEU_pilot3[i,c("Isolate1FreqPredicted", "Isolate2FreqPredicted", "RSquare")] <- c(sanger_mixture[[i]]$frac, sanger_mixture[[i]]$r2)

# Save raw CASEU output
CASEU_pilot3_raw_output <- sanger_mixture # Raw CASEU output
CASEU_pilot3_list_mix <- list_mix # list of the mixture
save(CASEU_pilot3_raw_output, CASEU_pilot3_list_mix, file = here::here("data/temp/CASEU_pilot3_raw_output.Rdata"))

#
CASEU_pilot3 <- CASEU_pilot3 %>% arrange(Community, Isolate1, Isolate2, Isolate1Freq)
fwrite(CASEU_pilot3, here::here("data/temp/CASEU_pilot3.csv"))



