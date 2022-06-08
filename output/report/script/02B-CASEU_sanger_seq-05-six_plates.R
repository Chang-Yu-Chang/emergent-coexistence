#' Read the CASEU result of six plates
#' In this CASEU batch, six plate were sent for Sanger sequencing
#'   B2 T7 933 P2
#'   B2 T7 444 P2
#'   C2 T7 13A P2
#'   C2 T7 13B P2
#'   D T7 75 P2
#'   D T7 5543 P2
#'   These are all P2, meaning that they contains 95-5 pairs
#'   The samples names were labelled as 1-96 of a plate, so we have to match it to the plate layout
#' 1. Read the plate layout
#' 2. Read the isolates trace matrices from the raw Sanger sequences
#' 3. Read the mixture trace matrices. Map the plate layout
#' 4. CASEU predicts the relative abundance
library(tidyverse)
library(CASEU)

read_trace_matrix <- function(abif.file) {
    x <- sangerseqR::read.abif(abif.file) %>% sangerseqR::sangerseq()
    return(x@traceMatrix)
    # paste0(folder_directory, Seq_ID, "-27F.ab1") %>%
    #     sangerseqR::read.abif() %>%
    #     CASEU::extractElectropherogram()
}

folder_directory <- "~/Dropbox/lab/emergent-coexistence/data/raw/Sanger/CASEU_six_plates/all_ab1/"

# Plate layout: join experiment plates and sequencing plates
plates_experiment <- read_csv("~/Dropbox/lab/emergent-coexistence/data/output/plates.csv", col_types = cols()) %>%
    filter(Plate == "P2", PlateLayout %in% c("933", "444", "13A", "13B", "75", "5543"))
plates_sequencing <- read_csv("~/Dropbox/lab/emergent-coexistence/data/raw/Sanger/CASEU_six_plates/20211202_Sanger_seq_prep-genewiz_table_CYC.csv", col_types = cols()) %>%
    select(Sample, FileName = `DNA Name`) %>%
    separate(col = FileName, remove = F, into = c("Batch", "Time", "PlateLayout", "Plate", "Well_ID"), sep = "_") %>%
    mutate(FileName = paste0(folder_directory, FileName, "-27F.ab1")) %>%
    select(-Well_ID) %>%
    # This is important: caseu label the wells by column
    mutate(Well = paste0(rep(LETTERS[1:8], 12), rep(sprintf("%02d", 1:12), each = 8)) %>% rep(6))

plates <- plates_experiment %>%
    left_join(plates_sequencing) %>%
    select(Sample, FileName, Community, Isolate1, Isolate2, Isolate1Freq, Isolate2Freq, MixIsolate)
    #select(Sample, FileName, Batch, PlateLayout, Community, Isolate1, Isolate2, Isolate1Freq, Isolate2Freq)


# Isolate trace
cat("\nStart reading isolate traces")
trace_isolate <- plates %>%
    filter(Isolate1 == Isolate2, Community != "blank") %>%
    arrange(Community, Isolate1, Isolate2) %>%
    select(Sample, FileName, Community, Isolate = Isolate1) %>%
    rowwise() %>%
    mutate(Trace = list(read_trace_matrix(FileName))) %>%
    mutate(TraceLength = nrow(Trace)) %>%
    arrange(Community, Isolate) %>%
    # Remove the short-trace sample due to bad quality sequencing
    filter(TraceLength > 1500) %>%
    group_by(Community, Isolate) %>%
    slice(1) %>%
    # Remove the contaminant
    filter(!(Community == "C11R2" & Isolate == 13))
cat("\nFinish reading", nrow(trace_isolate), "isolate traces")


# Read mixture trace data. 372 pairs that are not diagonal and not blank
trace_mixture <- plates %>%
    filter(Isolate1 != Isolate2, MixIsolate == T, Community != "blank") %>%
    arrange(Community, Isolate1, Isolate2) %>%
    select(-MixIsolate) %>%
    rowwise() %>%
    mutate(Trace = list(read_trace_matrix(FileName))) %>%
    mutate(TraceLength = nrow(Trace)) %>%
    arrange(Community, Isolate1, Isolate2) %>%
    # Remove the short-trace sample due to bad quality sequencing
    # filter(TraceLength > 1500) %>%
    # Remove the contaminant
    filter(!(Community == "C11R2" & Isolate1 == 13), !(Community == "C11R2" & Isolate2 == 13)) %>%
    mutate(CASEU = NA, Isolate1FreqPredicted = NA)
cat("\nFinish reading", nrow(trace_mixture), "mixture traces")


## 19 pairs have short traces
trace_mixture %>%
    arrange(Community, Isolate1, Isolate2) %>%
    filter(TraceLength < 1500) %>%
    nrow()


# Save a master csv
CASEU_sixplates <- trace_mixture %>%
    select(FileName,Community, Isolate1, Isolate2, Isolate1Freq, Isolate2Freq, Isolate1FreqPredicted)
#write_csv(CASEU_sixplates, "~/Dropbox/lab/emergent-coexistence/data/temp/CASEU_sixplates.csv")


# Fit mixture Sanger electropherogram  using CASEU packages. This may take a few minutes.
cat("\nStart caseu fitting", nrow(trace_mixture), "pairs")
# current_data <- read_csv("~/Dropbox/lab/emergent-coexistence/data/temp/CASEU_sixplates.csv", col_types = cols()) %>% mutate(Isolate1 = as.character(Isolate1), Isolate2 = as.character(Isolate2))
# to_tested <- which(is.na(trace_mixture$Isolate1FreqPredicted))
for (i in 1:nrow(trace_mixture)) {
#for (i in to_tested) {
    minimal_trace_length <- min(c(trace_mixture$TraceLength[trace_mixture$TraceLength > 1000], trace_isolate$TraceLength))
    isolate1_index <- which(trace_isolate$Community == trace_mixture$Community[i] & trace_isolate$Isolate == trace_mixture$Isolate1[i])
    isolate2_index <- which(trace_isolate$Community == trace_mixture$Community[i] & trace_isolate$Isolate == trace_mixture$Isolate2[i])

    # Skip short-trace mixtures
    if (trace_mixture$TraceLength[i] < 1000) {
        cat(paste0("\nPair trace ", trace_mixture$FileName[i], "\n",
                   trace_mixture$Community[i], " Isolate ", trace_mixture$Isolate1[i], " and ", trace_mixture$Isolate2[i], " is too short"))
        next
    }


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
    CASEU_sixplates <- trace_mixture %>% select(FileName,Community, Isolate1, Isolate2, Isolate1Freq, Isolate2Freq, Isolate1FreqPredicted)
    write_csv(CASEU_sixplates, "~/Dropbox/lab/emergent-coexistence/data/temp/CASEU_sixplates.csv")
    cat("\nrow", i, "out of", nrow(trace_mixture), "is finished")
}

CASEU_sixplates_trace_isolate <- trace_isolate
CASEU_sixplates_trace_mixture <- trace_mixture
save(CASEU_sixplates_trace_isolate, CASEU_sixplates_trace_mixture, file = "~/Dropbox/lab/emergent-coexistence/data/temp/CASEU_sixplates.Rdata")



if (FALSE) {

# CASEU prediction
# trace_mixture$CASEU_output <- NA
# trace_mixture$Isolate1FreqPredicted <- NA
# trace_mixture$Isolate2FreqPredicted <- NA
# trace_mixture$RSquare <- NA

current_data <- read_csv("~/Dropbox/lab/invasion-network/data/temp/CASEU_sixplatest.csv") %>% mutate(Isolate1 = as.character(Isolate1), Isolate2 = as.character(Isolate2))
is.na(current_data$Isolate1FreqPredicted) %>% sum
trace_mixture <- trace_mixture %>% left_join(current_data)

to_tested <- which(is.na(trace_mixture$Isolate1FreqPredicted))

tt <- proc.time()
#for (i in 1:nrow(trace_mixture)) {
for (i in to_tested) {
    # Community and isolates in a pair
    community <- trace_mixture$Community[i]
    isolate1 <- trace_mixture$Isolate1[i]
    isolate2 <- trace_mixture$Isolate2[i]

    # Trace matrices
    trace_mixture <- trace_mixture$Trace[[i]]
    trace_isolate1 <- trace_isolate %>% filter(Community == community, Isolate == isolate1) %>% pull(Trace) %>% `[[`(1)
    trace_isolate2 <- trace_isolate %>% filter(Community == community, Isolate == isolate2) %>% pull(Trace) %>% `[[`(1)

    # Skip single culture
    if (isolate1 == isolate2) next
    # Skip if there is no trace matrix data available for the mixture
    if (is.na(trace_mixture)) next
    # Skip if either of the isolate trace is missed
    if (is.na(trace_isolate1) || is.na(trace_isolate2)) {
        cat("\nOne of the isolate trace is missed")
        next
    }
    # Skip short-trace pairs
    if (nrow(trace_mixture) < 1000) {
        cat(paste0("\nPair trace ", trace_mixture$Seq_ID[i], " ",
                   trace_mixture$Community[i], " Isolate ", trace_mixture$Isolate1[i], " and ", trace_mixture$Isolate2[i], " is too short"))
        next
    }

    if (FALSE) {
    # Skip short trace from isolates
    if (nrow(trace_isolate1) < 1000) {
        temp <- trace_isolate %>% filter(Community == community, Isolate == isolate1)
        cat(paste0("\nIsolate trace ", temp$Seq_ID, " ", temp$Community, " Isolate ", temp$Isolate, " is too short"))
        next
    }

    # Skip short trace from isolates
    if (nrow(trace_isolate2) < 1000) {
        temp <- trace_isolate %>% filter(Community == community, Isolate == isolate2)
        cat(paste0("\nIsolate trace ", temp$Seq_ID, " ", temp$Community, " Isolate ", temp$Isolate, " is too short"))
        next
    }

    }


    trace_mixture$CASEU_output[[i]] <- CASEU::fitSangerMixture(
        mixture = trace_mixture,
        knots = seq(1500, min(c(nrow(trace_mixture), nrow(trace_isolate1), nrow(trace_isolate2))), by=1500),
        components = list(trace_isolate1, trace_isolate2)
    ) %>%
        list()


    #
    trace_mixture[i, c("Isolate1FreqPredicted", "Isolate2FreqPredicted", "RSquare")] <-
        trace_mixture$CASEU_output[[i]] %>%
        lapply(function(x) {
            data.frame(
                Isolate1FreqPredicted = unlist(x["frac"])[1],
                Isolate2FreqPredicted = unlist(x["frac"])[2],
                RSquare = unlist(x["r2"])[1])
        }) %>%
        `[[`(1)

    trace_mixture %>%
        select(-Trace) %>%
        write_csv(file = "~/Dropbox/lab/invasion-network/data/temp/CASEU_sixplates.csv")
    #
    cat("\n", i)
    cat("\n", (proc.time() - tt)[3], "seconds")

}

#save(trace_mixture, file = "~/Dropbox/lab/invasion-network/data/temp/CASEU_sixplates.Rdata")
}


