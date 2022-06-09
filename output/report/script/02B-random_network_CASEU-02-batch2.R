#' Read the CASEU result of random network
#' Batch 2: plate layout T3 C P2

library(tidyverse)
library(CASEU)

read_trace_matrix <- function(abif.file) {
    x <- sangerseqR::read.abif(abif.file) %>% sangerseqR::sangerseq()
    return(x@traceMatrix)
}

folder_directory_RN2 <- "~/Dropbox/lab/emergent-coexistence/data/raw/Sanger/CASEU_RN2/CASEU_RN_2_30-323299648_ab1/"
plates_RN2 <- read_csv("~/Dropbox/lab/emergent-coexistence/data/raw/Sanger/CASEU_RN2/plates_RN2.csv", col_types = cols())

# Isolate trace
cat("\nStart reading isolate traces")
trace_isolate <- plates_RN2 %>%
    mutate(FileName = paste0(folder_directory_RN2, Sample, "-27F.ab1")) %>%
    filter(Isolate1 == Isolate2, MixIsolate == F) %>%
    select(FileName, Sample, Community, Isolate = Isolate1, MixIsolate) %>%
    rowwise() %>%
    mutate(Trace = list(read_trace_matrix(FileName))) %>%
    mutate(TraceLength = nrow(Trace))
cat("\nFinish reading", nrow(trace_isolate), "isolate traces")

# Mixture trace
cat("\nStart reading mixture traces")
trace_mixture <- plates_RN2 %>%
    filter(MixIsolate == T, Isolate1 != Isolate2) %>%
    mutate(FileName = paste0(folder_directory_RN2, Sample, "-27F.ab1")) %>%
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
    CASEU_RN2 <- trace_mixture %>% select(FileName,Community, Isolate1, Isolate2, Isolate1Freq, Isolate2Freq, Time, Isolate1FreqPredicted)
    write_csv(CASEU_RN2, "~/Dropbox/lab/emergent-coexistence/data/temp/CASEU_RN2.csv")
    cat("\nrow", i, "out of", nrow(trace_mixture), "is finished")
}

CASEU_RN2_trace_isolate <- trace_isolate
CASEU_RN2_trace_mixture <- trace_mixture
save(CASEU_RN2_trace_isolate, CASEU_RN2_trace_mixture, file = "~/Dropbox/lab/emergent-coexistence/data/temp/CASEU_RN2.Rdata")

