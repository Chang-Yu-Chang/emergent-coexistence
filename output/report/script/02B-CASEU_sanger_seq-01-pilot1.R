#' This script runs the analysis of CASEU pilot1 Sanger sequences raw data
library(tidyverse)
library(data.table)
library(CASEU)
read_trace_matrix <- function(abif.file) {
  x <- sangerseqR::read.abif(abif.file) %>% sangerseqR::sangerseq()
  return(x@traceMatrix)
}

# Read trace matrices for isolates and mixtures from Sanger sequences ----
## Isolates
trace_isolate <- rep(list(NA), 4) %>% setNames(LETTERS[1:4])
for (i in 1:4) {
  trace_isolate[[i]] <- read_trace_matrix(here::here(paste0("data/raw/Sanger/CASEU_pilot/Caseu_pilot_30-265554768_ab1/", LETTERS[i],"-27F.ab1")))
}
trace_isolate[[1]] <- read_trace_matrix(here::here(paste0("data/raw/Sanger/CASEU_pilot/[Repeats]Caseu_pilot_30-269695307_ab1/A-27F_R.ab1"))) ## Repeat file by Genewiz

## Mixtures
list_mix <- c("AM_1_99_CD-27F", "AM_5_95_CD-27F", "OD_1_99_CD-27F", "OD_5_95_CD-27F",
  "AM_5_95_AD-27F", "OD_1_99_AD-27F", "OD_5_95_AC-27F", "OD_50_50_AD-27F", "AM_1_99_AD-27F", "AM_5_95_AD-27F", "OD_50_50_AC-27F_R",
  "OD_1_99_AB-27F_R", "OD_1_99_AC-27F_R", "OD_5_95_AB-27F_R", "OD_5_95_AD-27F_R", "AM_1_99_AC-27F_R")

trace_mix <- rep(list(NA), length(list_mix))

for (i in 1:length(list_mix)) {
  if (list_mix[i] == "OD_50_50_AC-27F_R") {
    trace_mix[[i]] <- read_trace_matrix(here::here(paste0("data/raw/Sanger/CASEU_pilot/[Repeats]Caseu_pilot_30-269695307_ab1/", list_mix[i], ".ab1")))
  } else if (grepl("_R$", list_mix[i])) { # Repeat
    trace_mix[[i]] <- read_trace_matrix(here::here(paste0("data/raw/Sanger/CASEU_pilot/[Repeats]Caseu_pilot_30-270218908_ab1/", list_mix[i], ".ab1")))
  } else {
    trace_mix[[i]] <- read_trace_matrix(here::here(paste0("data/raw/Sanger/CASEU_pilot/Caseu_pilot_30-265554768_ab1/", list_mix[i],".ab1")))
  }
}


# Make dataframe for CASEU ----
CASEU_pilot1 <- data.frame(MixtureName = list_mix)
temp <- CASEU_pilot1$MixtureName %>%
  gsub("-27F|-27F_R", "", .) %>%
  strsplit("_") %>%
  unlist() %>%
  matrix(ncol=4, byrow=T) %>%
  as.data.frame() %>%
  setNames(c("Treatment", "Isolate1Freq", "Isolate2Freq", "Isolates")) %>%
  mutate(Isolate1 = substr(Isolates, 1, 1), Isolate2 = substr(Isolates, 2, 2),
    Isolate1Freq = as.numeric(as.character(Isolate1Freq))/100,
    Isolate2Freq = as.numeric(as.character(Isolate2Freq))/100) %>%
  select(-Isolates)
CASEU_pilot1 <- cbind(CASEU_pilot1, temp)


# Fit mixture Sanger electorpheogram by using CASEU packages. This may take a few minutes. ----
sanger_mixture <- rep(list(NA), length(list_mix))
for (i in 1:length(list_mix)) {
  sanger_mixture[[i]] <- CASEU::fitSangerMixture(
    mixture = trace_mix[[i]],
    components = list(trace_isolate[[match(CASEU_pilot1$Isolate1[i], LETTERS[1:4])]],
      trace_isolate[[match(CASEU_pilot1$Isolate2[i], LETTERS[1:4])]])
  )
  print(i)
}


# Read CASEU predicted fraction ----
for (i in 1:length(list_mix)) CASEU_pilot1[i,c("Isolate1FreqPredicted", "Isolate2FreqPredicted")] <- sanger_mixture[[i]]$frac
CASEU_pilot1 <- CASEU_pilot1 %>% arrange(Isolate1, Isolate2, Treatment)

fwrite(CASEU_pilot1, here::here("data/temp/CASEU_pilot1.csv"))







