#' This script runs the analysis of CASEU pilot1 Sanger sequences raw data
#' Genewiz repeat the sanger twice so the one isolate and a couple of pairs have repeated sequence data
library(tidyverse)
library(CASEU)
library(sangerseqR)

# Read trace matrices for isolates and mixtures from Sanger sequences
read_trace_matrix <- function(abif.file) {
    x <- sangerseqR::read.abif(abif.file) %>% sangerseqR::sangerseq()
    return(x@traceMatrix)
}

# Isolates
trace_isolate <- tibble(FileName = paste0(LETTERS[1:4], "-27F.ab1"), Isolate = LETTERS[1:4], Trace = NA) %>%
    mutate(FileName = paste0("~/Dropbox/lab/emergent-coexistence/data/raw/Sanger/CASEU_pilot1/Caseu_pilot_30-265554768_ab1/", FileName))
cat("\nStart reading", nrow(trace_isolate), "isolate traces")
for (i in 1:4) trace_isolate$Trace[i] <- read_trace_matrix(trace_isolate$FileName[i]) %>% list()
trace_isolate$Trace[1] <- read_trace_matrix(paste0("~/Dropbox/lab/emergent-coexistence/data/raw/Sanger/CASEU_pilot1/[Repeats]Caseu_pilot_30-269695307_ab1/A-27F_R.ab1")) %>% list ## Repeat file by Genewiz

# Mixtures full list
list_mixture <- read_csv("~/Dropbox/lab/emergent-coexistence/data/raw/Sanger/CASEU_pilot1/protocol_20190813_Sanger_seq_prep-genewiz_table.csv", col_types = cols()) %>%
    pull(`DNA Name`) %>%
    str_replace("  ", "") %>%
    str_subset("_") %>%
    # Amplicon treatment actually uses 3:99 because the small volume. Replace it by 3 such that the file names match
    str_replace("3", "1") %>%
    paste0("~/Dropbox/lab/emergent-coexistence/data/raw/Sanger/CASEU_pilot1/Caseu_pilot_30-265554768_ab1/", ., "-27F.ab1")

## List of mixture repeated by Genewiz
temp <- list.files("~/Dropbox/lab/emergent-coexistence/data/raw/Sanger/CASEU_pilot1/[Repeats]Caseu_pilot_30-270218908_ab1", full.names = T) %>%
    c("~/Dropbox/lab/emergent-coexistence/data/raw/Sanger/CASEU_pilot1/[Repeats]Caseu_pilot_30-269695307_ab1/OD_50_50_AC-27F_R.ab1")
list_mixture_repeat <- temp[!str_detect(temp, "/B")]

trace_mixture <- tibble(FileName = c(list_mixture, list_mixture_repeat)) %>%
    mutate(Repeat = ifelse(str_detect(FileName, "_R.ab1"), T, F)) %>%
    rowwise() %>%
    mutate(temp = str_split(FileName, "/", simplify = T) %>% last() %>% str_replace("(-27F.ab1)|(-27F_R.ab1)", "")) %>%
    separate(col = temp, into = c("Treatment", "Isolate1Freq", "Isolate2Freq", "Isolate")) %>%
    mutate(Isolate1 = str_sub(Isolate, 1, 1), Isolate2 = str_sub(Isolate, 2, 2)) %>%
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
    # CASEU output
    trace_mixture$CASEU[i] <- CASEU::fitSangerMixture(
        mixture = trace_mixture$Trace[[i]],
        components = list(trace_isolate$Trace[[which(trace_isolate$Isolate == trace_mixture$Isolate1[i])]],
                          trace_isolate$Trace[[which(trace_isolate$Isolate == trace_mixture$Isolate2[i])]])
    ) %>% list

    # Extract predicted fraction
    trace_mixture$Isolate1FreqPredicted[i] <- trace_mixture$CASEU[[i]]$frac[1]

    # Output intermediate
    CASEU_pilot1 <- trace_mixture %>%
    select(FileName, Repeat, Treatment, Isolate1, Isolate2, Isolate1Freq, Isolate2Freq, Isolate1FreqPredicted)
    write_csv(CASEU_pilot1, "~/Dropbox/lab/emergent-coexistence/data/temp/CASEU_pilot1.csv")
    cat("\nrow", i, "out of", nrow(trace_mixture), "is finished")
}


CASEU_pilot1_trace_isolate <- trace_isolate
CASEU_pilot1_trace_mixture <- trace_mixture
save(CASEU_pilot1_trace_isolate, CASEU_pilot1_trace_mixture, file = "~/Dropbox/lab/emergent-coexistence/data/temp/CASEU_pilot1.Rdata")
















