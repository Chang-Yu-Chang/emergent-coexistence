#' This script runs the analysis of CASEU pilot3 Sanger sequences raw data
library(tidyverse)
library(CASEU)

read_trace_matrix <- function(abif.file) {
    x <- sangerseqR::read.abif(abif.file) %>% sangerseqR::sangerseq()
    return(x@traceMatrix)
}

folder_directory <- "~/Dropbox/lab/emergent-coexistence/data/raw/Sanger/CASEU_pilot3/30-286536939_ab1/" # Replace it with the correct folder directory
genewiz_pilot3 <- read_csv("~/Dropbox/lab/emergent-coexistence/data/raw/Sanger/CASEU_pilot3/protocol_20190924_Sanger_seq_prep-genewiz_table_CYC.csv", col_types = cols())

# Isolate trace
trace_isolate <- tibble(FileName = genewiz_pilot3$`DNA Name`) %>%
    mutate(FileName = paste0(folder_directory, FileName, "-27F.ab1")) %>%
    # Choose isolates
    filter(!str_detect(FileName, "50_50"), !str_detect(FileName, "5_95"), !str_detect(FileName, "95_5"), !str_detect(FileName, "mock")) %>%
    rowwise() %>%
    mutate(temp = str_split(FileName, "/", simplify = T) %>% last() %>% str_replace("(-27F.ab1)|(-27F_R.ab1)", "")) %>%
    separate(temp, into = c("Community", "Isolate"), sep = "_", remove = F) %>%
    mutate(Isolate = str_replace(Isolate, "-27F.ab1", "")) %>%
    select(FileName, Community, Isolate) %>%
    mutate(Trace = NA)
cat("\nStart reading", nrow(trace_isolate), "isolate traces")
for (i in 1:nrow(trace_isolate)) trace_isolate$Trace[i] <- read_trace_matrix(trace_isolate$FileName[i]) %>% list()


# Mixture trace
trace_mixture <- tibble(FileName = genewiz_pilot3$`DNA Name`) %>%
    mutate(FileName = paste0(folder_directory, FileName, "-27F.ab1")) %>%
    # Remove isolates
    filter(str_detect(FileName, "50_50") |str_detect(FileName, "5_95") | str_detect(FileName, "95_5"), !str_detect(FileName, "mock")) %>%
    rowwise() %>%
    mutate(temp = str_split(FileName, "/", simplify = T) %>% last() %>% str_replace("(-27F.ab1)|(-27F_R.ab1)", "")) %>%
    separate(temp, into = c("Community", "Isolate1Freq", "Isolate2Freq", "Isolate1", "Isolate2"), sep = "_", remove = F) %>%
    mutate(Isolate2 = str_replace(Isolate2, "-27F.ab1", "")) %>%
    select(FileName, Community, Isolate1, Isolate2, Isolate1Freq, Isolate2Freq) %>%
    mutate(Trace = NA) %>%
    arrange(Community, Isolate1, Isolate2)
cat("\nStart reading", nrow(trace_mixture), "mixture traces")
for (i in 1:nrow(trace_mixture)) trace_mixture$Trace[i] <- read_trace_matrix(trace_mixture$FileName[i]) %>% list


## Mixture trace file length
trace_mixture <- trace_mixture %>%
    rowwise() %>%
    mutate(TraceLength = nrow(Trace)) %>%
    mutate(CASEU = NA, Isolate1FreqPredicted = NA)

# Fit mixture Sanger electropherogram  using CASEU packages. This may take a few minutes.
cat("\nStart caseu fitting", nrow(trace_mixture), "pairs")
for (i in 1:nrow(trace_mixture)) {
    isolate1_index <- which(trace_isolate$Community == trace_mixture$Community[i] & trace_isolate$Isolate == trace_mixture$Isolate1[i])
    isolate2_index <- which(trace_isolate$Community == trace_mixture$Community[i] & trace_isolate$Isolate == trace_mixture$Isolate2[i])
    # CASEU output
    trace_mixture$CASEU[i] <- CASEU::fitSangerMixture(
        mixture = trace_mixture$Trace[[i]],
        components = list(trace_isolate$Trace[[isolate1_index]], trace_isolate$Trace[[isolate2_index]]),
        tol = 0.01
    ) %>% list

    # Extract predicted fraction
    trace_mixture$Isolate1FreqPredicted[i] <- trace_mixture$CASEU[[i]]$frac[1]

    # Output intermediate
    CASEU_pilot3 <- trace_mixture %>% select(FileName,Community, Isolate1, Isolate2, Isolate1Freq, Isolate2Freq, Isolate1FreqPredicted)
    write_csv(CASEU_pilot3, "~/Dropbox/lab/emergent-coexistence/data/temp/CASEU_pilot3.csv")
    cat("\nrow", i, "out of", nrow(trace_mixture), "is finished")
}


CASEU_pilot3_trace_isolate <- trace_isolate
CASEU_pilot3_trace_mixture <- trace_mixture
save(CASEU_pilot3_trace_isolate, CASEU_pilot3_trace_mixture, file = "~/Dropbox/lab/emergent-coexistence/data/temp/CASEU_pilot3.Rdata")


