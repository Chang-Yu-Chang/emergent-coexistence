 #' This script runs the analysis of CASEU pilot2 Sanger sequences raw data
library(CASEU)
library(sangerseqR)

# Read trace matrices for isolates and mixtures from Sanger sequences
read_trace_matrix <- function(abif.file) {
    x <- sangerseq(read.abif(abif.file))
    return(x@traceMatrix)
}

folder_directory <- "~/Dropbox/lab/emergent-coexistence/data/raw/Sanger/CASEU_pilot2/30-280311837_ab1/" # Replace it with the correct folder directory

genewiz_pilot2 <- fread("~/Dropbox/lab/emergent-coexistence/data/raw/Sanger/CASEU_pilot2/protocol_20190910_Sanger_seq_prep-genewiz_table_CYC.csv")

# Isolate trace
trace_isolate <- data.table(FileName = paste0(folder_directory, genewiz_pilot2$`DNA Name`, "-27F.ab1"))
trace_isolate <- trace_isolate[FileName %like% 'single']
trace_isolate[, temp := gsub('(-27F.ab1)|(-27F_R.ab1)', '', basename(FileName))]
trace_isolate[, c("temp1", "Community", "temp2", "Isolate") := tstrsplit(temp, '_')]
trace_isolate$Trace <- NA

# Choose isolates
for (i in 1:nrow(trace_isolate)) trace_isolate$Trace[i] <- list(read_trace_matrix(trace_isolate$FileName[i]))


# Mixture trace
trace_mixture <- data.table(FileName = paste0(folder_directory, genewiz_pilot2$`DNA Name`, "-27F.ab1"))
trace_mixture <- trace_mixture[!(FileName %like% 'single')]
trace_mixture[, temp := gsub('(-27F.ab1)|(-27F_R.ab1)', '', basename(FileName))]
trace_mixture[, c("temp1", "Community", "Isolate1Freq", "Isolate2Freq", "Isolate") := tstrsplit(temp, '_')]
trace_mixture[, c("Isolate1", "Isolate2") := tstrsplit(Isolate, '')]
trace_mixture$Trace <- NA

for (i in 1:nrow(trace_mixture)) trace_mixture$Trace[i] <- list(read_trace_matrix(trace_mixture$FileName[i]))
trace_mixture$TraceLength  <- unlist(lapply(trace_mixture$Trace, nrow))
trace_mixture[, CASEU := NA]
trace_mixture[, Isolate1FreqPredicted:= NA]
trace_mixture[, R2 := NA]


# Fit mixture Sanger electropherogram  using CASEU packages. This may take a few minutes.
cat("\nStart caseu fitting", nrow(trace_mixture), "pairs")
for (i in 1:nrow(trace_mixture)) {
    isolate1_index <- which(trace_isolate$Community == trace_mixture$Community[i] &
                            trace_isolate$Isolate == trace_mixture$Isolate1[i])
    isolate2_index <- which(trace_isolate$Community == trace_mixture$Community[i] &
                            trace_isolate$Isolate == trace_mixture$Isolate2[i])
    # CASEU output
    trace_mixture$CASEU[i] <- list(CASEU::fitSangerMixture(
        mixture = trace_mixture$Trace[[i]],
        components = list(trace_isolate$Trace[[isolate1_index]], trace_isolate$Trace[[isolate2_index]])
    ))

    # Extract predicted fraction
    trace_mixture$Isolate1FreqPredicted[i] <- trace_mixture$CASEU[[i]]$frac[1]
    trace_mixture$R2[i] <- trace_mixture$CASEU[[i]]$r2

    # Output intermediate
    cat("\nrow", i, "out of", nrow(trace_mixture), "is finished")
}

# Output
isolates_ID_match <- fread('data/temp/isolates_ID_match.csv')

aux <- isolates_ID_match
names(aux) <- c('ExpID_1', 'ID_1', 'Community', 'Isolate1')
aux[, Isolate1 := as.character(Isolate1)]
trace_mixture <- merge(trace_mixture, aux, by=c('Community', 'Isolate1'), all.x=TRUE)

aux <- isolates_ID_match
names(aux) <- c('ExpID_2', 'ID_2', 'Community', 'Isolate2')
aux[, Isolate2 := as.character(Isolate2)]
trace_mixture <- merge(trace_mixture, aux, by=c('Community', 'Isolate2'), all.x=TRUE)

CASEU_pilot2_trace_isolate <- trace_isolate
CASEU_pilot2_trace_mixture <- trace_mixture
save(CASEU_pilot2_trace_isolate, CASEU_pilot2_trace_mixture,
     file = "~/Dropbox/lab/emergent-coexistence/data/temp/CASEU_pilot2.Rdata")



