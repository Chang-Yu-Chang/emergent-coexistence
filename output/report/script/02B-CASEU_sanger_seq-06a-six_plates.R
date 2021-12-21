#' Read the CASEU result of six plates
#' In this CASEU batch, six plate were sent for Sanger sequencing
#'   B2 T7 933 P2
#'   B2 T7 444 P2
#'   C2 T7 13A P2
#'   C2 T7 13B P2
#'   D T7 75 P2
#'   D T7 5543 P2
#' 1. Read the plate layout
#' 2. Read the isolates trace matrices from the raw Sanger sequences
#' 3. Read the mixture trace matrices. Map the plate layout
#' 4. CASEU predicts the relative abundance
library(tidyverse)
library(CASEU)
`%notin%` <- Negate(`%in%`)

# Read plate layout
plates <- read_csv(here::here("data/output/plates.csv")) %>%
    filter(Plate == "P2", PlateLayout %in% c("933", "444", "13A", "13B", "75", "5543")) %>%
    mutate(Experiment = str_replace(Experiment, "Transitivity_", ""))

# Read sequencing ID
plates_submitted <- read_csv("~/Dropbox/lab/invasion-network/data/raw/Sanger/CASEU_six_plates/20211202_Sanger_seq_prep-genewiz_table_CYC.csv") %>%
    select(Sample, Seq_ID = `DNA Name`) %>%
    separate(col = Seq_ID, remove = F, into = c("Experiment", "Time", "PlateLayout", "Plate", "Well_ID"), sep = "_", convert = T) %>%
    select(-Well_ID) %>%
    mutate(Well = paste0(rep(LETTERS[1:8], 12), rep(sprintf("%02d", 1:12), each = 8)) %>% rep(6))
plates <- plates_submitted %>% left_join(plates)

# Read isolate trace data. 69 isolates
folder_directory <- "~/Dropbox/lab/invasion-network/data/raw/Sanger/CASEU_six_plates/all_ab1/"
read_trace <- function(Seq_ID) {
    paste0(folder_directory, Seq_ID, "-27F.ab1") %>%
        sangerseqR::read.abif() %>%
        CASEU::extractElectropherogram()
}

isolates_trace <- plates %>%
    filter(Isolate1 == Isolate2, MixIsolate == F, Community != "blank") %>%
    arrange(Sample) %>%
    select(Sample, Seq_ID, Experiment, PlateLayout, Plate, Well, Community, Isolate = Isolate1) %>%
    rowwise() %>%
    mutate(Trace = list(read_trace(Seq_ID)))

# Read mixture trace data. 396 pairs
pairs_trace <- plates %>%
    filter(Isolate1 != Isolate2, MixIsolate == T, Community != "blank") %>%
    arrange(Sample) %>%
    select(-MixIsolate) %>%
    rowwise() %>%
    mutate(Trace = list(read_trace(Seq_ID)))


# CASEU prediction
# pairs_trace$CASEU_output <- NA
# pairs_trace$Isolate1FreqPredicted <- NA
# pairs_trace$Isolate2FreqPredicted <- NA
# pairs_trace$RSquare <- NA

current_data <- read_csv("~/Dropbox/lab/invasion-network/data/temp/CASEU_six_plates_result.csv") %>% mutate(Isolate1 = as.character(Isolate1), Isolate2 = as.character(Isolate2))
pairs_trace <- pairs_trace %>% left_join(current_data)


tt <- proc.time()
#for (i in 1:nrow(pairs_trace)) {
for (i in 289:nrow(pairs_trace)) {
    # Community and isolates in a pair
    community <- pairs_trace$Community[i]
    isolate1 <- pairs_trace$Isolate1[i]
    isolate2 <- pairs_trace$Isolate2[i]

    # Trace matrices
    trace_mixture <- pairs_trace$Trace[[i]]
    trace_isolate1 <- isolates_trace %>% filter(Community == community, Isolate == isolate1) %>% pull(Trace) %>% `[[`(1)
    trace_isolate2 <- isolates_trace %>% filter(Community == community, Isolate == isolate2) %>% pull(Trace) %>% `[[`(1)

    # Skip single culture
    if (isolate1 == isolate2) next
    # Skip if there is no trace matrix data available for the mixture
    if (is.na(trace_mixture)) next
    # Skip if either of the isolate trace is missed
    if (is.na(trace_isolate1) || is.na(trace_isolate2)) {
        cat("\nOne of the isolate trace is missed")
        next
    }
    # Skip short trace
    if (nrow(trace_mixture) < 1000) {
        cat(paste0("\nPair trace ", pairs_trace$Seq_ID[i], " ",
                   pairs_trace$Community[i], " Isolate ", pairs_trace$Isolate1[i], " and ", pairs_trace$Isolate2[i], " is too short"))
        next
    }

    # Skip short trace from isolates
    if (nrow(trace_isolate1) < 1000) {
        temp <- isolates_trace %>% filter(Community == community, Isolate == isolate1)
        cat(paste0("\nIsolate trace ", temp$Seq_ID, " ", temp$Community, " Isolate ", temp$Isolate, " is too short"))
        next
    }

    # Skip short trace from isolates
    if (nrow(trace_isolate2) < 1000) {
        temp <- isolates_trace %>% filter(Community == community, Isolate == isolate2)
        cat(paste0("\nIsolate trace ", temp$Seq_ID, " ", temp$Community, " Isolate ", temp$Isolate, " is too short"))
        next
    }

    pairs_trace$CASEU_output[[i]] <- CASEU::fitSangerMixture(
        mixture = trace_mixture,
        knots = seq(1500, min(c(nrow(trace_mixture), nrow(trace_isolate1), nrow(trace_isolate2))), by=1500),
        components = list(trace_isolate1, trace_isolate2)
    ) %>%
        list()


    #
    pairs_trace[i, c("Isolate1FreqPredicted", "Isolate2FreqPredicted", "RSquare")] <-
        pairs_trace$CASEU_output[[i]] %>%
        lapply(function(x) {
            data.frame(
                Isolate1FreqPredicted = unlist(x["frac"])[1],
                Isolate2FreqPredicted = unlist(x["frac"])[2],
                RSquare = unlist(x["r2"])[1])
        }) %>%
        `[[`(1)

    pairs_trace %>%
        select(-Trace) %>%
        write_csv(file = "~/Dropbox/lab/invasion-network/data/temp/CASEU_six_plates_result.csv")
    #
    cat("\n", i)
    cat("\n", (proc.time() - tt)[3], "seconds")

}

save(pairs_trace, file = "~/Dropbox/lab/invasion-network/data/temp/CASEU_six_plates_result.Rdata")

