#' This script runs the analysis of CASEU pilot1 Sanger sequences raw data
#' Genewiz repeat the sanger twice so the one isolate and a couple of pairs have repeated sequence data
library(tidyverse)
library(CASEU)
library(sangerseqR)

# Read trace matrices for isolates and mixtures from Sanger sequences
read_trace_matrix <- function(abif.file) {
    x <- sangerseq(read.abif(abif.file))
    return(x@traceMatrix)
}

# Isolates
trace_isolate <- data.table(FileName = paste0('~/Dropbox/lab/emergent-coexistence/data/raw/Sanger/CASEU_pilot1/Caseu_pilot_30-265554768_ab1/',
                                              LETTERS[1:4],
                                              '-27F.ab1'),
                            Isolate = LETTERS[1:4],
                            Trace=NA)

# Read isolate traces
for (i in 1:4) trace_isolate$Trace[i] <- list(read_trace_matrix(trace_isolate$FileName[i]))
trace_isolate$Trace[1] <- list(read_trace_matrix(paste0("~/Dropbox/lab/emergent-coexistence/data/raw/Sanger/CASEU_pilot1/[Repeats]Caseu_pilot_30-269695307_ab1/A-27F_R.ab1"))) ## Repeat file by Genewiz

#djordje - why repeat file here? 

# Mixtures full list
list_mixture <- fread('~/Dropbox/lab/emergent-coexistence/data/raw/Sanger/CASEU_pilot1/protocol_20190813_Sanger_seq_prep-genewiz_table.csv')
list_mixture <- list_mixture$`DNA Name`
list_mixture <- gsub('  ', '', list_mixture)
list_mixture <- list_mixture[list_mixture %like% '_']
# Amplicon treatment uses 3:99 because the small volume. Replace by 3 such that file names match
list_mixture <- gsub('3', '1', list_mixture)
list_mixture <- paste0('~/Dropbox/lab/emergent-coexistence/data/raw/Sanger/CASEU_pilot1/Caseu_pilot_30-265554768_ab1/', list_mixture, '-27F.ab1')

#djordje - why replace 3 by 1? 

## List of mixture repeated by Genewiz
aux <- list.files("~/Dropbox/lab/emergent-coexistence/data/raw/Sanger/CASEU_pilot1/[Repeats]Caseu_pilot_30-270218908_ab1", full.names = T)
aux <- c(aux, '~/Dropbox/lab/emergent-coexistence/data/raw/Sanger/CASEU_pilot1/[Repeats]Caseu_pilot_30-269695307_ab1/OD_50_50_AC-27F_R.ab1')
list_mixture_repeat <- aux[!(aux %like% "/B")]


# build mixtures dataframe
trace_mixture <- data.table(FileName = c(list_mixture, list_mixture_repeat))
trace_mixture[, Repeat := ifelse(FileName %like% '_R.ab1', T, F)]
trace_mixture[, temp := gsub('(-27F.ab1)|(-27F_R.ab1)', '', basename(FileName))]
trace_mixture[, c("Treatment", "Isolate1Freq", "Isolate2Freq", "Isolate") := tstrsplit(temp, '_')]
trace_mixture[, c("Isolate1", "Isolate2") := tstrsplit(Isolate, '')]

trace_mixture[, Trace := NA]
for (i in 1:nrow(trace_mixture)) trace_mixture$Trace[i] <- list(read_trace_matrix(trace_mixture$FileName[i]))
trace_mixture$TraceLength  <- unlist(lapply(trace_mixture$Trace, nrow))
trace_mixture[, CASEU := NA]
trace_mixture[, Isolate1FreqPredicted:= NA]
trace_mixture[, R2 := NA]

trace_mixture <- trace_mixture[!(Isolate %like% 'A')]

# Fit mixtures Sanger electropherogram using CASEU packages. This may take a few minutes.
for (i in 1:nrow(trace_mixture)) {
    # CASEU output
    trace_mixture$CASEU[i] <- list(CASEU::fitSangerMixture(
        mixture = trace_mixture$Trace[i][[1]],
        components = list(trace_isolate$Trace[which(trace_isolate$Isolate == trace_mixture$Isolate1[i])][[1]],
                          trace_isolate$Trace[which(trace_isolate$Isolate == trace_mixture$Isolate2[i])][[1]])))

    # Extract predicted fraction
    trace_mixture$Isolate1FreqPredicted[i] <- trace_mixture$CASEU[[i]]$frac[1]
    trace_mixture$R2[i] <- trace_mixture$CASEU[[i]]$r2
    cat("\nrow", i, "out of", nrow(trace_mixture), "is finished")
}

# Output
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







isolates_ID_match <- fread('data/temp/isolates_ID_match.csv')
isolates_ID_match <- isolates_ID_match[ID %in% c('440', '287')]

aux <- isolates_ID_match
names(aux) <- paste0(names(aux), '_isol1')
aux$Isolate1 <- c('C', 'D')
trace_mixture <- merge(trace_mixture, aux, by='Isolate1', all.x=TRUE)

aux <- isolates_ID_match
names(aux) <- paste0(names(aux), '_isol2')
aux$Isolate2 <- c('C', 'D')
trace_mixture <- merge(trace_mixture, aux, by='Isolate2', all.x=TRUE)



CASEU_pilot1_trace_isolate <- trace_isolate
CASEU_pilot1_trace_mixture <- trace_mixture
save(CASEU_pilot1_trace_isolate, CASEU_pilot1_trace_mixture, file = "~/Dropbox/lab/emergent-coexistence/data/temp/CASEU_pilot1.Rdata")
