#' Read the CASEU result of random network
#' Batch 1: plate layout T3 C P2
#' Important information before starting: The sample 3, 4, 5, 12, 13, 27, 29, 33, 73, 85, 86, 91, and 94
#' arrived Genewiz in dry tube, so the samples are not included in the batch.
library(tidyverse)
library(CASEU)

read_trace_matrix <- function(abif.file) {
    x <- sangerseqR::read.abif(abif.file) %>% sangerseqR::sangerseq()
    return(x@traceMatrix)
}

folder_directory_RN1 <- "~/Dropbox/lab/emergent-coexistence/data/raw/Sanger/CASEU_RN1/Caseu_RN_1_30-322411127_ab1/"
available_sanger <- list.files(here::here(paste0(folder_directory_RN1))) %>% gsub("-27F.ab1", "", .)
if (FALSE) {
    plates_random <- read_csv("~/Dropbox/lab/emergent-coexistence/data/output/plates_random.csv", col_types = cols())
    plates_RN1 <- plates_random %>%
        filter(PlateLayout == "C", MixPlate == "P2") %>%
        mutate(Well = factor(Well, paste0(rep(LETTERS[1:8], 12), sprintf("%02d", rep(1:12, each = 8))))) %>%
        arrange(Well) %>%
        mutate(Sample = 1:96) %>%
        mutate(Time = "T3") %>%
        select(Sample, everything())
    write_csv(plates_RN1, "~/Dropbox/lab/emergent-coexistence/data/raw/Sanger/CASEU_RN1/plates_RN1.csv")

}
plates_RN1 <- read_csv("~/Dropbox/lab/emergent-coexistence/data/raw/Sanger/CASEU_RN1/plates_RN1.csv", col_types = cols())

# Isolate trace
cat("\nStart reading isolate traces")
trace_isolate <- plates_RN1 %>%
    mutate(FileName = paste0(folder_directory_RN1, Sample, "-27F.ab1")) %>%
    #filter(Isolate1 == Isolate2) %>%
    filter(Isolate1 == Isolate2, MixIsolate == F) %>%
    select(FileName, Sample, Community, Isolate = Isolate1, MixIsolate, Time) %>%
    # Remove unavailable samples
    filter(Sample %in% available_sanger) %>%
    rowwise() %>%
    mutate(Trace = list(read_trace_matrix(FileName))) %>%
    mutate(TraceLength = nrow(Trace))
cat("\nFinish reading", nrow(trace_isolate), "isolate traces")

# Mixture trace
cat("\nStart reading mixture traces")
trace_mixture <- plates_RN1 %>%
    filter(MixIsolate == T, Isolate1 != Isolate2) %>%
    mutate(FileName = paste0(folder_directory_RN1, Sample, "-27F.ab1")) %>%
    select(FileName, Sample, Community, Isolate1, Isolate2, Isolate1Freq, Isolate2Freq) %>%
    # Remove unavailable samples
    filter(Sample %in% available_sanger) %>%
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

    # Skip mixture with missing isolate electropherogram
    if (length(isolate1_index) == 0 | length(isolate2_index) == 0) next

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
    CASEU_RN1 <- trace_mixture %>% select(FileName,Community, Isolate1, Isolate2, Isolate1Freq, Isolate2Freq, Isolate1FreqPredicted)
    write_csv(CASEU_RN1, "~/Dropbox/lab/emergent-coexistence/data/temp/CASEU_RN1.csv")
    cat("\nrow", i, "out of", nrow(trace_mixture), "is finished")
}

CASEU_RN1_trace_isolate <- trace_isolate
CASEU_RN1_trace_mixture <- trace_mixture
save(CASEU_RN1_trace_isolate, CASEU_RN1_trace_mixture, file = "~/Dropbox/lab/emergent-coexistence/data/temp/CASEU_RN1.Rdata")




if (FALSE) {
    ### Parse mixture column
    CASEU_RN1 <- mixture_list %>%
        tidyr::separate(col = Mixture, sep = "_", remove = F, into = c("Community", "Isolate1Freq", "Isolate2Freq", "Isolate1", "Isolate2") ) %>%
        mutate(
            Isolate1 = as.numeric(as.character(Isolate1)),
            Isolate2 = as.numeric(as.character(Isolate2)),
            Isolate1Freq = as.numeric(as.character(Isolate1Freq))/100,
            Isolate2Freq = as.numeric(as.character(Isolate2Freq))/100) %>%
        # Transfer
        mutate(Transfer = "T3")

    #### Available sanger seuqencing result
    CASEU_RN1$DataAvailability <- CASEU_RN1$Community == "RanAss1" & !(CASEU_RN1$Isolate1 %in% 5:6 | CASEU_RN1$Isolate2 %in% 5:6 | CASEU_RN1$Isolate1 == CASEU_RN1$Isolate2)

    ## Fit Sanger sequnece electorpheogram of mixture by using CASEU packages. This may take a few minutes.
    ### Keep raw CASEU outputs

    caseu_prediction <- rep(list(NA), nrow(mixture_list))
    names(caseu_prediction) <- mixture_list$Mixture
    tt <- proc.time()

    for (i in 1:length(caseu_prediction)) {
        #for (i in 1:2) {
        community <- CASEU_RN1$Community[i]
        isolate1 <- CASEU_RN1$Isolate1[i]
        isolate2 <- CASEU_RN1$Isolate2[i]

        trace_mixture <- traces_mixture[[i]]
        trace_isolate1 <- traces_isolate[[match(paste0(community, "_", isolate1), names(traces_isolate))]]
        trace_isolate2 <- traces_isolate[[match(paste0(community, "_", isolate2), names(traces_isolate))]]

        if (isolate1 == isolate2) next
        if (is.na(traces_mixture[[i]])) next
        if (CASEU_RN1$DataAvailability[i]) {
            caseu_prediction[[i]] <- CASEU::fitSangerMixture( # CASEU package function
                mixture = trace_mixture, # Mixture trace
                knots = seq(1500, min(c(nrow(trace_mixture), nrow(trace_isolate1), nrow(trace_isolate2))), by=1500),
                components = list(trace_isolate1, trace_isolate2)
            )
        }

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
        rbindlist(idcol = "Mixture")

    ### Format the CASEU df
    switch_pairwise_column <- function (df, bypair = T) {
        if (any(is.factor(df$Isolate1))) df$Isolate1 <- as.numeric(df$Isolate1); df$Isolate2 <- as.numeric(df$Isolate2)
        if ("Isolate1FreqPredicted" %in% colnames(df)) {
            if (bypair == T) {
                temp_index <- df$Isolate1 > df$Isolate2
                df[temp_index, c("Isolate1", "Isolate2", "Isolate1Freq", "Isolate2Freq", "Isolate1FreqPredicted", "Isolate2FreqPredicted")] <-
                    df[temp_index, c("Isolate2", "Isolate1", "Isolate2Freq", "Isolate1Freq", "Isolate2FreqPredicted", "Isolate1FreqPredicted")]

                df %>% arrange(Isolate1, Isolate2, Isolate1Freq) %>% return()
            } else if (bypair == F) {
                temp_index <- df$Isolate1Freq == 5
                df[temp_index, c("Isolate1", "Isolate2", "Isolate1Freq", "Isolate2Freq", "Isolate1FreqPredicted", "Isolate2FreqPredicted")] <-
                    df[temp_index, c("Isolate2", "Isolate1", "Isolate2Freq", "Isolate1Freq", "Isolate2FreqPredicted", "Isolate1FreqPredicted")]

                df %>% arrange(Isolate1Freq, Isolate1, Isolate2) %>% return()
            }
        } else {

            if (bypair == T) {
                temp_index <- df$Isolate1 > df$Isolate2
                df[temp_index, c("Isolate1", "Isolate2", "Isolate1Freq", "Isolate2Freq")] <-
                    df[temp_index, c("Isolate2", "Isolate1", "Isolate2Freq", "Isolate1Freq")]

                df %>% arrange(Isolate1, Isolate2, Isolate1Freq) %>% return()
            } else if (bypair == F) {
                temp_index <- df$Isolate1Freq == 5
                df[temp_index, c("Isolate1", "Isolate2", "Isolate1Freq", "Isolate2Freq")] <-
                    df[temp_index, c("Isolate2", "Isolate1", "Isolate2Freq", "Isolate1Freq")]

                df %>% arrange(Isolate1Freq, Isolate1, Isolate2) %>% return()
            }
        }
    }
    CASEU_RN1 <- CASEU_RN1 %>%
        left_join(temp) %>%
        select(-Mixture, -DataAvailability, -Sample) %>%
        filter(!is.na(Isolate1FreqPredicted)) %>%
        select(Community, Isolate1, Isolate2, Isolate1Freq, Isolate2Freq, Isolate1FreqPredicted, Isolate2FreqPredicted) %>%
        arrange(Isolate1, Isolate2, Isolate1Freq, Isolate2Freq) %>%
        # Switch isolates so that isolate2 is always greater than isolate1
        switch_pairwise_column()


    # Save raw CASEU output
    save(traces_isolate, file = here::here("data/temp/CASEU_RN_isolates_trace.Rdata"))
    CASEU_RN1_raw_output <- caseu_prediction # Raw CASEU output
    CASEU_RN1_mixture_list <- mixture_list # list of the mixture
    #save(CASEU_RN1_raw_output, CASEU_RN1_mixture_list, file = here::here("data/temp/CASEU_RN1_raw_output.Rdata"))


    fwrite(CASEU_RN1, file = here::here("data/temp/CASEU_RN1.csv"))
}
