#' This script runs the analysis of CASEU pilot4 Sanger sequences raw data
#' For the isolate sequences, I need to use isolate sequences from CASEU pilot3
library(tidyverse)
library(CASEU)

read_trace_matrix <- function(abif.file) {
    x <- sangerseqR::read.abif(abif.file) %>% sangerseqR::sangerseq()
    return(x@traceMatrix)
}


folder_directory_pilot3 <- "~/Dropbox/lab/emergent-coexistence/data/raw/Sanger/CASEU_pilot3/30-286536939_ab1/"
folder_directory_pilot4 <- "~/Dropbox/lab/emergent-coexistence/data/raw/Sanger/CASEU_pilot4/Caseu_4_30-292390814_ab1/"
genewiz_pilot3 <- read_csv("~/Dropbox/lab/emergent-coexistence/data/raw/Sanger/CASEU_pilot3/protocol_20190924_Sanger_seq_prep-genewiz_table_CYC.csv", col_types = cols())
genewiz_pilot4 <- read_csv("~/Dropbox/lab/emergent-coexistence/data/raw/Sanger/CASEU_pilot4/protocol_20191007_Sanger_seq_prep-genewiz_table_CYC.csv", col_types = cols())

# Isolate trace
cat("\nStart reading isolate traces")
trace_isolate <- tibble(FileName = c(paste0(folder_directory_pilot3, genewiz_pilot3$`DNA Name`, "-27F.ab1"),
                                     paste0(folder_directory_pilot4, genewiz_pilot4$`DNA Name`, "-27F.ab1"))) %>%
    # Choose isolates
    filter(!str_detect(FileName, "50_50"), !str_detect(FileName, "5_95"), !str_detect(FileName, "95_5"), !str_detect(FileName, "mock")) %>%
    rowwise() %>%
    mutate(temp = str_split(FileName, "/", simplify = T) %>% last() %>% str_replace("(-27F.ab1)|(-27F_R.ab1)", "")) %>%
    separate(temp, into = c("Community", "Isolate"), sep = "_", remove = F) %>%
    mutate(Isolate = str_replace(Isolate, "-27F.ab1", "")) %>%
    select(FileName, Community, Isolate) %>%
    rowwise() %>%
    mutate(Trace = list(read_trace_matrix(FileName))) %>%
    mutate(TraceLength = nrow(Trace))
cat("\nFinish reading", nrow(trace_isolate), "isolate traces")

# Mixture trace
cat("\nStart reading mixture traces")
trace_mixture <- tibble(FileName = paste0(folder_directory_pilot4, genewiz_pilot4$`DNA Name`, "-27F.ab1")) %>%
    # Remove isolates
    filter(str_detect(FileName, "50_50") |str_detect(FileName, "5_95") | str_detect(FileName, "95_5"), !str_detect(FileName, "mock")) %>%
    rowwise() %>%
    mutate(temp = str_split(FileName, "/", simplify = T) %>% last() %>% str_replace("(-27F.ab1)|(-27F_R.ab1)", "")) %>%
    separate(temp, into = c("Community", "Isolate1Freq", "Isolate2Freq", "Isolate1", "Isolate2"), sep = "_", remove = F) %>%
    mutate(Isolate2 = str_replace(Isolate2, "-27F.ab1", "")) %>%
    select(FileName, Community, Isolate1, Isolate2, Isolate1Freq, Isolate2Freq) %>%
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
    CASEU_pilot4 <- trace_mixture %>% select(FileName,Community, Isolate1, Isolate2, Isolate1Freq, Isolate2Freq, Isolate1FreqPredicted)
    write_csv(CASEU_pilot4, "~/Dropbox/lab/emergent-coexistence/data/temp/CASEU_pilot4.csv")
    cat("\nrow", i, "out of", nrow(trace_mixture), "is finished")
}


CASEU_pilot4_trace_isolate <- trace_isolate
CASEU_pilot4_trace_mixture <- trace_mixture
save(CASEU_pilot4_trace_isolate, CASEU_pilot4_trace_mixture, file = "~/Dropbox/lab/emergent-coexistence/data/temp/CASEU_pilot4.Rdata")
















if (FALSE) {


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
}



