#' Read the CASEU result of random network
#' In this CASEU RN batch 4 (CASEU_RN4), one plate was sent for Sanger sequencing
#' T3 BD P3

library(tidyverse)
library(CASEU)

read_trace_matrix <- function(abif.file) {
    x <- sangerseqR::read.abif(abif.file) %>% sangerseqR::sangerseq()
    return(x@traceMatrix)
}

folder_directory_RN4 <- "~/Dropbox/lab/emergent-coexistence/data/raw/Sanger/CASEU_RN4/CASEU_RN4_30-481538036_ab1/"
if (FALSE) {
plates_random <- read_csv("~/Dropbox/lab/emergent-coexistence/data/output/plates_random.csv", col_types = cols())
plates_RN4 <- plates_random %>%
    filter(PlateLayout == "BD", MixPlate == "P3") %>%
    mutate(Well = factor(Well, paste0(rep(LETTERS[1:8], 12), sprintf("%02d", rep(1:12, each = 8))))) %>%
    arrange(Well) %>%
    mutate(Sample = 1:96) %>%
    mutate(Time = "T3") %>%
    select(Sample, everything())

write_csv(plates_RN4, "~/Dropbox/lab/emergent-coexistence/data/raw/Sanger/CASEU_RN4/plates_RN4.csv")

}
plates_RN4 <- read_csv("~/Dropbox/lab/emergent-coexistence/data/raw/Sanger/CASEU_RN4/plates_RN4.csv", col_types = cols())


# Isolate trace: load trace_isolate from RN2
load("~/Dropbox/lab/emergent-coexistence/data/temp/CASEU_RN2.Rdata")
cat("\nStart reading isolate traces from batch RN2")
trace_isolate <- CASEU_RN2_trace_isolate
cat("\nFinish reading", nrow(trace_isolate), "isolate traces")

# Mixture trace
cat("\nStart reading mixture traces")
trace_mixture <- plates_RN4 %>%
    mutate(FileName = paste0(folder_directory_RN4, Sample, "-27F.ab1")) %>%
    filter(MixIsolate == T, Isolate1 != Isolate2) %>%
    select(FileName, Sample, Community, Isolate1, Isolate2, Isolate1Freq, Isolate2Freq, Time) %>%
    rowwise() %>%
    mutate(Trace = list(read_trace_matrix(FileName))) %>%
    mutate(TraceLength = nrow(Trace)) %>%
    arrange(Community, Isolate1, Isolate2) %>%
    mutate(CASEU = NA, Isolate1FreqPredicted = NA)
cat("\nFinish reading", nrow(trace_mixture), "mixture traces")


# Fit mixture Sanger electropherogram  using CASEU packages. This may take a few minutes.
cat("\nStart caseu fitting", nrow(trace_mixture), "pairs")
for (i in 1:nrow(trace_mixture)) {
    minimal_trace_length <- min(c(trace_mixture$TraceLength, trace_isolate$TraceLength))
    isolate1_index <- which(trace_isolate$Community == trace_mixture$Community[i] & trace_isolate$Isolate == trace_mixture$Isolate1[i])
    isolate2_index <- which(trace_isolate$Community == trace_mixture$Community[i] & trace_isolate$Isolate == trace_mixture$Isolate2[i])

    # CASEU output
    trace_mixture$CASEU[i] <- CASEU::fitSangerMixture(
        mixture = trace_mixture$Trace[[i]],
        components = list(trace_isolate$Trace[[isolate1_index]], trace_isolate$Trace[[isolate2_index]]),
        knots = seq(1500, minimal_trace_length, by = 1500),
        tol = 0.01
    ) %>% list

    # Extract predicted fraction
    trace_mixture$Isolate1FreqPredicted[i] <- trace_mixture$CASEU[[i]]$frac[1]

    # Output intermediate
    CASEU_RN4 <- trace_mixture %>% select(FileName,Community, Isolate1, Isolate2, Isolate1Freq, Isolate2Freq, Time, Isolate1FreqPredicted)
    write_csv(CASEU_RN4, "~/Dropbox/lab/emergent-coexistence/data/temp/CASEU_RN4.csv")
    cat("\nrow", i, "out of", nrow(trace_mixture), "is finished")
}

CASEU_RN4_trace_isolate <- trace_isolate
CASEU_RN4_trace_mixture <- trace_mixture
save(CASEU_RN4_trace_isolate, CASEU_RN4_trace_mixture, file = "~/Dropbox/lab/emergent-coexistence/data/temp/CASEU_RN4.Rdata")











if (FALSE) {

# Read plate layout ----
plates_random <- fread(here::here("data/output/plates_random.csv")) %>% as_tibble()
plates_BD_P3 <- plates_random %>% filter(PlateLayout == "BD", MixPlate == "P3") %>% mutate(Sample = 1:96)

# Read isolates trace matrices ----
#' Object traces_isolate is a R list with all 32 isolates used in
#' random networks. The CASEU data is obtained from CASEU_RN2
load(here::here("data/temp/CASEU_RN_isolates_trace.Rdata"))

# Read mixtures trace matrices ----
## Folder of raw Sanger sequence data
folder_directory_RN4_P1 <- "data/raw/Sanger/CASEU_RN4/CASEU_RN4_30-481538036_ab1/"

## Function for reading trace matrices in different folders
read_trace_temp <- function(folder, name) {
  paste0(folder, name, "-27F.ab1") %>%
    here::here() %>%
    sangerseqR::read.abif() %>%
    CASEU::extractElectropherogram()
}

## Create a R list for reading muxture traces
traces_mixture <- rep(list(NA), 96)
names(traces_mixture) <- 1:96

for (i in 1:96){
  if (i <= 96) traces_mixture[[i]] <- read_trace_temp(folder_directory_RN4_P1, i)
}


## Read the mixture layout on the plate
temp_df1 <- plates_BD_P3 %>%
  mutate(Isolate1Freq = Isolate1Freq * 100, Isolate2Freq = Isolate2Freq * 100) %>%
  mutate(Time = "T3") %>%
  unite("Mixture", c("Community", "Isolate1Freq", "Isolate2Freq", "Isolate1", "Isolate2")) %>%
  unite("Plate", c("Time", "PlateLayout", "MixPlate")) %>%
  select(Sample, Plate, Mixture)

mixture_list <- bind_rows(temp_df1)
names(traces_mixture) <- mixture_list$Mixture


# CASEU prediction ----
## Separate the mixture
CASEU_RN4 <- mixture_list %>%
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
## This may take ~20 minutes
caseu_prediction <- rep(list(NA), nrow(mixture_list))
names(caseu_prediction) <- mixture_list$Mixture
tt <- proc.time()

for (i in 1:length(caseu_prediction)) {
  # Community and isolates in a pair
  community <- CASEU_RN4$Community[i]
  isolate1 <- CASEU_RN4$Isolate1[i]
  isolate2 <- CASEU_RN4$Isolate2[i]

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
  mutate(Sample = c(1:96))

### Format the CASEU df
CASEU_RN4 <- CASEU_RN4 %>%
  left_join(temp, by = c("Sample", "Mixture")) %>%
  filter(Isolate1 != Isolate2) %>%
  switch_pairwise_column() %>%
  arrange(Sample) %>%
{.}



# Save raw CASEU output
CASEU_RN4_raw_output <- caseu_prediction # Raw CASEU output
CASEU_RN4_mixture_list <- mixture_list # list of the mixture
#save(CASEU_RN4_raw_output, CASEU_RN4_mixture_list, file = here::here("data/temp/CASEU_RN4_raw_output.Rdata"))

fwrite(CASEU_RN4, here::here("data/temp/CASEU_RN4.csv"))
}





