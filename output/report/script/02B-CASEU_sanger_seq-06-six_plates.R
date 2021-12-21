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
#library(data.table)
library(CASEU)
`%notin%` <- Negate(`%in%`)
#source(here::here("output/report/script/misc.R"))

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
pairs_trace$CASEU_output <- NA
pairs_trace$Isolate1FreqPredicted <- NA
pairs_trace$Isolate2FreqPredicted <- NA
pairs_trace$RSquare <- NA

tt <- proc.time()
for (i in 1:nrow(pairs_trace)) {
    # Community and isolates in a pair
    community <- pairs_trace$Community[i]
    isolate1 <- pairs_trace$Isolate1[i]
    isolate2 <- pairs_trace$Isolate2[i]

    # Trace matrices
    trace_mixture <- pairs_trace$Trace[[i]]
    trace_isolate1 <- isolates_trace %>% filter(Community == community, Isolate == isolate1) %>% pull(Trace) %>% `[[`(1)
    trace_isolate2 <- isolates_trace %>% filter(Community == community, Isolate == isolate2) %>% pull(Trace) %>% `[[`(1)

    #if (isolate1 == isolate2) next # Skip single culture
    #if (is.na(trace_mixture)) next # Skip if there is no trace matrix data available for the mixture
    if (is.na(trace_isolate1) || is.na(trace_isolate2)) {
        cat("One of the isolate trace is missed")
        next # Skip if either of the isolate trace is missed
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
    print(i)
    cat((proc.time() - tt)[3], "seconds\n")

}

save(pairs_trace, file = "~/Dropbox/lab/invasion-network/data/temp/CASEU_six_plates_result.Rdata")

if (FALSE) {
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

}

if (FALSE) {

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

