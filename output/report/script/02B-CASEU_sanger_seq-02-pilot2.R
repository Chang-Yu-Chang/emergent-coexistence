#' This script runs the analysis of CASEU pilot2 Sanger sequences raw data
library(tidyverse)
library(CASEU)

read_trace_matrix <- function(abif.file) {
    x <- sangerseqR::read.abif(abif.file) %>% sangerseqR::sangerseq()
    return(x@traceMatrix)
}

folder_directory <- "~/Dropbox/lab/emergent-coexistence/data/raw/Sanger/CASEU_pilot2/30-280311837_ab1/" # Replace it with the correct folder directory
genewiz_pilot2 <- read_csv("~/Dropbox/lab/emergent-coexistence/data/raw/Sanger/CASEU_pilot2/protocol_20190910_Sanger_seq_prep-genewiz_table_CYC.csv", col_types = cols())

# Isolate trace
trace_isolate <- tibble(FileName = genewiz_pilot2$`DNA Name`) %>%
    mutate(FileName = paste0(folder_directory, FileName, "-27F.ab1")) %>%
    # Choose isolates
    filter(str_detect(FileName, "single")) %>%
    rowwise() %>%
    mutate(temp = str_split(FileName, "/", simplify = T) %>% last() %>% str_replace("(-27F.ab1)|(-27F_R.ab1)", "")) %>%
    separate(temp, into = c("temp1", "Community", "temp2", "Isolate"), sep = "_") %>%
    mutate(Isolate = str_replace(Isolate, "-27F.ab1", "")) %>%
    select(FileName, Community, Isolate) %>%
    mutate(Trace = NA)
cat("\nStart reading", nrow(trace_isolate), "isolate traces")
for (i in 1:nrow(trace_isolate)) trace_isolate$Trace[i] <- read_trace_matrix(trace_isolate$FileName[i]) %>% list()


# Mixture trace
trace_mixture <- tibble(FileName = genewiz_pilot2$`DNA Name`) %>%
    mutate(FileName = paste0(folder_directory, FileName, "-27F.ab1")) %>%
    # Remove isolates
    filter(!str_detect(FileName, "single")) %>%
    rowwise() %>%
    mutate(temp = str_split(FileName, "/", simplify = T) %>% last() %>% str_replace("(-27F.ab1)|(-27F_R.ab1)", "")) %>%
    separate(temp, into = c("temp1", "Community", "Isolate1Freq", "Isolate2Freq", "Isolate"), sep = "_") %>%
    mutate(Isolate1 = str_sub(Isolate, 1, 1), Isolate2 = str_sub(Isolate, 2,2)) %>%
    select(FileName, Community, Isolate1, Isolate2, Isolate1Freq, Isolate2Freq) %>%
    mutate(Trace = NA)
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
        components = list(trace_isolate$Trace[[isolate1_index]], trace_isolate$Trace[[isolate2_index]])
    ) %>% list

    # Extract predicted fraction
    trace_mixture$Isolate1FreqPredicted[i] <- trace_mixture$CASEU[[i]]$frac[1]

    # Output intermediate
    CASEU_pilot2 <- trace_mixture %>% select(FileName,Community, Isolate1, Isolate2, Isolate1Freq, Isolate2Freq, Isolate1FreqPredicted)
    write_csv(CASEU_pilot2, "~/Dropbox/lab/emergent-coexistence/data/temp/CASEU_pilot2.csv")
    cat("\nrow", i, "out of", nrow(trace_mixture), "is finished")
}


CASEU_pilot2_trace_isolate <- trace_isolate
CASEU_pilot2_trace_mixture <- trace_mixture
save(CASEU_pilot2_trace_isolate, CASEU_pilot2_trace_mixture, file = "~/Dropbox/lab/emergent-coexistence/data/temp/CASEU_pilot2.Rdata")


