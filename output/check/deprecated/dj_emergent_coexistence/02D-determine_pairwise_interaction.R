


# CASEU results ----
load("~/Dropbox/lab/emergent-coexistence/data/temp/CASEU_pilot2.Rdata")
pilot2 <- CASEU_pilot2_trace_mixture
pilot2 <- pilot2[, .(Community, Isolate1, Isolate2,
                     ExpID_1, ID_1, ExpID_2, ID_2,
                     Isolate1Freq, Isolate2Freq, CASEU,
                     Isolate1FreqPredicted, R2)]

load("~/Dropbox/lab/emergent-coexistence/data/temp/CASEU_pilot3.Rdata")
pilot3 <- CASEU_pilot3_trace_mixture
pilot3 <- pilot3[, .(Community, Isolate1, Isolate2,
                     ExpID_1, ID_1, ExpID_2, ID_2,
                     Isolate1Freq, Isolate2Freq, CASEU,
                     Isolate1FreqPredicted, R2)]

load("~/Dropbox/lab/emergent-coexistence/data/temp/CASEU_pilot4.Rdata")
pilot4 <- CASEU_pilot4_trace_mixture
pilot4 <- pilot4[, .(Community, Isolate1, Isolate2,
                     ExpID_1, ID_1, ExpID_2, ID_2,
                     Isolate1Freq, Isolate2Freq, CASEU,
                     Isolate1FreqPredicted, R2)]

load("~/Dropbox/lab/emergent-coexistence/data/temp/CASEU_sixplates.Rdata")
pilot6p <- CASEU_sixplates_trace_mixture
pilot6p <- pilot6p[, .(Community, Isolate1, Isolate2,
                       ExpID_1, ID_1, ExpID_2, ID_2,
                       Isolate1Freq, Isolate2Freq, CASEU,
                       Isolate1FreqPredicted, R2)]


d <- rbind(pilot2, pilot3, pilot4, pilot6p)
d <- d[!is.na(ExpID_1) & !is.na(ExpID_2)]
d[, Isolate1 := as.integer(Isolate1)]
d[, Isolate2 := as.integer(Isolate2)]

# assign pair number 
pairs_ID <- fread("~/Dropbox/lab/emergent-coexistence/data/temp/pairs_ID.csv")

d1 <- d[Isolate1<Isolate2]
d1 <- merge(d1, pairs_ID, all.x=TRUE)


d2 <- d[Isolate1>Isolate2]
names(d2) <- c('Community',
               'Isolate2', 'Isolate1', 'ExpID_2', 'ID_2', 'ExpID_1', 'ID_1',
               'Isolate2Freq', 'Isolate1Freq', 'CASEU',
               'Isolate1FreqPredicted', 'R2')
d2[, Isolate1FreqPredicted := 1-Isolate1FreqPredicted]
d2 <- merge(d2, pairs_ID, all.x=TRUE)

d <- rbind(d1, d2)
d <- d[R2>0.8]
d[, Isolate1Freq := as.numeric(Isolate1Freq)/100]
d[, Isolate2Freq := as.numeric(Isolate2Freq)/100]

# Determine the interaction types by fitness functions ----
# Table for determining interaction types
interaction_type_three <- data.table(
    FromRare = rep(c(1, -1, 0), each = 9),
    FromMedium = rep(rep(c(1, -1, 0), each = 3), 3),
    FromAbundant = rep(c(1, -1, 0), 9),
    InteractionType = NA,
    InteractionTypeFiner = NA)

interaction_type_two <- data.table(
    FromRare = rep(c(1, -1, 0), each = 3),
    FromMedium = rep(NA, 9),
    FromAbundant = rep(c(1, -1, 0), 3),
    InteractionType = NA,
    InteractionTypeFiner = NA)

interaction_type <- rbind(interaction_type_three, interaction_type_two)

## Assign interaction types to combinations of frequency changes signs
interaction_type$InteractionType[c(1, 14, 28, 32, 10, 13, 31)] <- "exclusion"
interaction_type$InteractionType[c(2, 3, 5, 8, 9, 23, 26, 29, 30, 33, 4, 11, 12, 15, 17, 20, 34, 35)] <- "coexistence"
interaction_type$InteractionType[c(27, 36)] <- "neutrality"


## Assign finer interaction types to combinations of frequency changes signs
interaction_type$InteractionTypeFiner[c(1, 14, 28, 32)] <- "competitive exclusion"
interaction_type$InteractionTypeFiner[c(10, 13, 31)] <- "mutual exclusion"
interaction_type$InteractionTypeFiner[c(2, 3, 5, 8, 9, 23, 26, 29, 30, 33)] <- "stable coexistence"
interaction_type$InteractionTypeFiner[c(4, 11, 12, 15, 17, 20, 34, 35)] <- "frequency-dependent coexistence"
interaction_type$InteractionTypeFiner[c(27,36)] <- "neutrality"

interaction_type <- interaction_type %>%  mutate(FreqFunc = paste(FromRare, FromMedium, FromAbundant, sep = "_"))

getCompetitionOutcome <- function(x)
{
    r <- x[1]
    r <- r[, .(Community, Isolate1, Isolate2, ExpID_1, ID_1, ExpID_2, ID_2, PairID)]
    r$FromRare <- NA
    r$FromMedium <- NA
    r$FromAbundant <- NA

    r$FreqIniRare <- NA
    r$FreqFinRare <- NA
    r$FreqIniMedium <- NA
    r$FreqFinMedium <- NA
    r$FreqIniAbundant <- NA
    r$FreqFinAbundant <- NA
    
    for (i in 1:nrow(x)) {
        if (x$Isolate1Freq[i]==.05) {

            r$FreqIniRare <- x$Isolate1Freq[i]
            r$FreqFinRare <- x$Isolate1FreqPredicted[i]
            
            if (x$Isolate1FreqPredicted[i]>.05) r$FromRare <- 1
            if (x$Isolate1FreqPredicted[i]<.05) r$FromRare <- -1
            if (x$Isolate1FreqPredicted[i]==.05) r$FromRare <- 0
        }

        if (x$Isolate1Freq[i]==.5) {

            r$FreqIniMedium <- x$Isolate1Freq[i]
            r$FreqFinMedium <- x$Isolate1FreqPredicted[i]

            if (x$Isolate1FreqPredicted[i]>.5) r$FromMedium <- 1
            if (x$Isolate1FreqPredicted[i]<.5) r$FromMedium <- -1
            if (x$Isolate1FreqPredicted[i]==.5) r$FromMedium <- 0
        }
        
        if (x$Isolate1Freq[i]==.95) {

            r$FreqIniAbundant <- x$Isolate1Freq[i]
            r$FreqFinAbundant <- x$Isolate1FreqPredicted[i]
            
            if (x$Isolate1FreqPredicted[i]>.95) r$FromAbundant <- 1
            if (x$Isolate1FreqPredicted[i]<.95) r$FromAbundant <- -1
            if (x$Isolate1FreqPredicted[i]==.95) r$FromAbundant <- 0
        }
    }
    return(r)
}

res <- rbindlist(lapply(split(d, by='PairID'), getCompetitionOutcome))
res <- merge(res, interaction_type, all.x=TRUE)        

pairs_dj <- fread('data/temp/pairs_ID.csv')

res[, s1 := ifelse(ID_1>ID_2, ID_2, ID_1)] 
res[, s2 := ifelse(ID_1>ID_2, ID_1, ID_2)] 

res <- merge(pairs_dj, res, by=c('s1', 's2', 'Community'), all.x=TRUE)

