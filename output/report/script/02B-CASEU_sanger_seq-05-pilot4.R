#' This script analyszes CASEU pilot4 Sanger sequences raw data
#' For the isolate sequences, I need to use isolate sequences from CASEU pilot3
#'
library(tidyverse)
library(data.table)
library(CASEU)
read_trace_matrix <- function(abif.file) {
  x <- sangerseqR::read.abif(abif.file) %>% sangerseqR::sangerseq()
  return(x@traceMatrix)
}

folder_directory_pilot3 <- "data/raw/Sanger/CASEU_pilot3/30-286536939_ab1/"
folder_directory_pilot4 <- "data/raw/Sanger/CASEU_pilot4/Caseu_4_30-292390814_ab1/"
genewiz_pilot3 <- fread(here::here("output/protocol/tab_fig/protocol_20190924_Sanger_seq_prep-genewiz_table_CYC.csv"))
genewiz_pilot4 <- fread(here::here("output/protocol/tab_fig/protocol_20191007_Sanger_seq_prep-genewiz_table_CYC.csv"))

# Read trace matrices for isolates and mixtures from Sanger sequences ----
list_seq <- paste0(genewiz_pilot4$`DNA Name`)

## Isolates
isolates_names <- c(list_seq[31:40], genewiz_pilot3$`DNA Name`[43:63]) # Some isolates are sequenced in CASEU pilot3
trace_isolates <- rep(list(NA), length(isolates_names)) %>% setNames(isolates_names)
#for (i in 1:length(isolates_names)) trace_isolates[[i]] <- invnet::read_trace_matrix(here::here(paste0(folder_directory, isolates_names[i], "-27F.ab1")))
for (i in 1:length(isolates_names)) {
  if (i <= 10) {
    trace_isolates[[i]] <- # Isolates sequenced in CASEU pilot4
      paste0(folder_directory_pilot4, isolates_names[i], "-27F.ab1") %>%
      here::here() %>%
      sangerseqR::read.abif() %>%
      CASEU::extractElectropherogram()
  } else if (i > 10) { # Isolates sequenced in CASEU pilot3
    trace_isolates[[i]] <-
      paste0(folder_directory_pilot3, isolates_names[i], "-27F.ab1") %>%
      here::here() %>%
      sangerseqR::read.abif() %>%
      CASEU::extractElectropherogram()
  }
}

## Mixtures
list_mix <- list_seq[1:30]
trace_mix <- rep(list(NA), length(list_mix))
for (i in 1:length(list_mix)) {
  trace_mix[[i]] <-
    paste0(folder_directory_pilot4, list_mix[i], "-27F.ab1") %>%
    here::here() %>%
    sangerseqR::read.abif() %>%
    CASEU::extractElectropherogram()
}

# Make dataframe for CASEU ----
CASEU_pilot4 <- data.frame(MixtureName = list_mix)
temp <- CASEU_pilot4$MixtureName %>%
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
CASEU_pilot4 <- cbind(CASEU_pilot4, temp)


# Fit mixture Sanger electorpheogram by using CASEU packages. This may take a few minutes. ----
sanger_mixture <- rep(list(NA), length(list_mix))
tt = proc.time()
for (i in 1:length(sanger_mixture)) {
#for (i in 7) {
  sanger_mixture[[i]] <- CASEU::fitSangerMixture( # CASEU package function
    mixture = trace_mix[[i]], # Mixture data
    components = list(
      trace_isolates[[match(paste0(CASEU_pilot4$Community[i], "_", CASEU_pilot4$Isolate1[i]), isolates_names)]], # Isolate 1
      trace_isolates[[match(paste0(CASEU_pilot4$Community[i], "_", CASEU_pilot4$Isolate2[i]), isolates_names)]]), # Isolate 2
    knots = seq(1500, 8000, by=1500),
    tol = 0.1 # Tolerace for R2
  )
  print(i)
  cat((proc.time() - tt)[3], "seconds\n")
}

# Read CASEU predicted fraction ----
for (i in 1:length(list_mix)) {
  if (is.list(sanger_mixture[[i]])) CASEU_pilot4[i,c("Isolate1FreqPredicted", "Isolate2FreqPredicted", "RSquare")] <- c(sanger_mixture[[i]]$frac, sanger_mixture[[i]]$r2)
}

# Save raw CASEU output
CASEU_pilot4_raw_output <- sanger_mixture # Raw CASEU output
CASEU_pilot4_list_mix <- list_mix # list of the mixture
save(CASEU_pilot4_raw_output, CASEU_pilot4_list_mix, file = here::here("data/temp/CASEU_pilot4_raw_output.Rdata"))

CASEU_pilot4 <- CASEU_pilot4 %>% arrange(Community, Isolate1, Isolate2, Isolate1Freq)
fwrite(CASEU_pilot4, here::here("data/temp/CASEU_pilot4.csv"))


