library(data.table)
library(tidyverse)
library(cowplot)
source(here::here("analysis/00-metadata.R"))

pairs_ID <- read_csv(paste0(folder_data, "temp/00c-pairs_ID.csv"), show_col_types = F)
pairs_boots <- read_csv(paste0(folder_data, "temp/07-pairs_boots.csv"), show_col_types = F)
pairs_freq <- read_csv(paste0(folder_data, "temp/25-pairs_freq.csv"), show_col_types = F)

boots_T0 <- pairs_boots %>% filter(Time == "T0") %>% as.data.table
boots_T0 <- boots_T0[, .(Community, Isolate1, Isolate2, Isolate1InitialODFreq, BootstrapID, Isolate1CFUFreq)]
boots_T8 <- pairs_boots %>% filter(Time == "T8") %>% as.data.table
pairs_freq <- pairs_freq %>% select(PairFreqID, Community, Isolate1, Isolate2, Isolate1InitialODFreq, Time, Isolate1CFUFreqMean)

#
pairs_I1 <- dcast(pairs_freq, PairFreqID + Community + Isolate1 + Isolate2 + Isolate1InitialODFreq ~ Time, value.var = "Isolate1CFUFreqMean")
names(pairs_I1)[6:7] <- c('T0_I1_mean', 'T8_I1_mean')
pairs_I1 <- as.data.table(pairs_I1)
pairs_I1[, change := NA]
pairs_I1 <- pairs_I1[!is.na(T8_I1_mean)]

# assign how frequency of isolate 1 changes
for (i in 1:nrow(pairs_I1)) {

    aux0 <- boots_T0[Community==pairs_I1$Community[i] &
                     Isolate1==pairs_I1$Isolate1[i] &
                     Isolate2==pairs_I1$Isolate2[i] &
                     Isolate1InitialODFreq==pairs_I1$Isolate1InitialODFreq[i]]
    aux0 <- aux0$Isolate1CFUFreq

    aux8 <- boots_T8[Community==pairs_I1$Community[i] &
                     Isolate1==pairs_I1$Isolate1[i] &
                     Isolate2==pairs_I1$Isolate2[i] &
                     Isolate1InitialODFreq==pairs_I1$Isolate1InitialODFreq[i]]
    aux8 <- aux8$Isolate1CFUFreq

    sgnf <- wilcox.test(aux0, aux8)$p.value

    if (mean(aux0)<mean(aux8) & sgnf<0.001) {
        pairs_I1$change[i] <- 'up'
    } else if (mean(aux0)>mean(aux8) & sgnf<0.001) {
        pairs_I1$change[i] <- 'dw'
    } else {
        pairs_I1$change[i] <- 'nc'
    }
}

# assign competition outcome
pairs_outcome <- unique(pairs_I1[, 2:4])
pairs_outcome[, outcome := NA]

#read_csv(paste0(folder_data, "temp/26-pairs_freq.csv"), show_col_types = F)
#isol <- fread(paste0(folder_data, "output/pairs_remained.csv"))
isol <- fread(paste0(folder_data, "temp/00c-pairs_ID.csv"))
pairs_outcome <- merge(isol[, 3:5], pairs_outcome, all.x=TRUE)

for (i in 1:nrow(pairs_outcome)) {

    aux <- pairs_I1[Community==pairs_outcome$Community[i] &
                    Isolate1==pairs_outcome$Isolate1[i] &
                    Isolate2==pairs_outcome$Isolate2[i]]
    aux <- aux[order(Isolate1InitialODFreq)]

    # Skip no-colony pairs
    if (nrow(aux) < 3) next

    # bootstraps
    aux0 <- boots_T0[Community==pairs_outcome$Community[i] &
                     Isolate1==pairs_outcome$Isolate1[i] &
                     Isolate2==pairs_outcome$Isolate2[i]]
    aux0 <- dcast(aux0, BootstrapID~Isolate1InitialODFreq, value.var='Isolate1CFUFreq')

    aux8 <- boots_T8[Community==pairs_outcome$Community[i] &
                     Isolate1==pairs_outcome$Isolate1[i] &
                     Isolate2==pairs_outcome$Isolate2[i]]
    aux8 <- dcast(aux8, BootstrapID~Isolate1InitialODFreq, value.var='Isolate1CFUFreq')

    eq <- apply(as.matrix(aux8[,-1]), 1, mean)
    eq.m <- mean(eq)

    # does freq 50 start above, below or indistinguishable from eq?
    starts.at.eq <- wilcox.test(aux0$`50`, eq)$p.value>0.05

    if (all(aux$T8_I1_mean==0) | all(aux$T8_I1_mean==1)) {
        pairs_outcome$outcome[i] <- '1-exclusion' # exctinction of one isolate in all three expts
    } else if (all(aux$change=='up') |
               all(aux$change=='dw')) {
        pairs_outcome$outcome[i] <- '2-exclusion' # one goes towards extinction in all three expts
    } else if (aux$T0_I1_mean[1]<eq.m & aux$change[1]=='up' &
               aux$T0_I1_mean[3]>eq.m & aux$change[3]=='dw' &
               ((!starts.at.eq & aux$T0_I1_mean[2]>eq.m & aux$change[2]=='dw') |
                (!starts.at.eq & aux$T0_I1_mean[2]<eq.m & aux$change[2]=='up') |
                (starts.at.eq & aux$change[2]=='nc')) &
               all(aux$T8_I1_mean>0) & all(aux$T8_I1_mean<1)) { # all go towards equilibrium
        if (starts.at.eq & aux$change[2]=='nc') print ("50 stays at eq")
        pairs_outcome$outcome[i] <- '3-coexistence'
    } else if (all(aux$T8_I1_mean>0) & all(aux$T8_I1_mean<1)) {
        pairs_outcome$outcome[i] <- '4-coexistence' # coexist without freq dep
    } else {
        pairs_outcome$outcome[i] <- '5-inconclusive'
    }
}

pairs_outcome <- pairs_outcome[Community!='C10R2']

# Determine the direction of exclusion
pairs_outcome <- pairs_outcome %>%
    left_join(pairs_ID)

# Determine direction of exclusion
pairs_exclusion <- pairs_outcome %>%
    filter(outcome %in% c("1-exclusion", "2-exclusion")) %>%
    left_join(pairs_ID)
pairs_freq_exclusion <- pairs_freq %>%
    left_join(pairs_ID) %>%
    select(PairID, Community, Isolate1, Isolate2, Isolate1InitialODFreq, Time, Isolate1CFUFreqMean) %>%
    pivot_wider(names_from = Time, values_from = Isolate1CFUFreqMean) %>%
    drop_na(T8) %>%
    filter(PairID %in% pairs_exclusion$PairID) %>%
    mutate(Sign = sign(T8-T0)) %>%
    group_by(PairID) %>%
    summarize(Sign = case_when(all(Sign==1) ~ 1, all(Sign == -1) ~-1))

pairs_ID_winner <- pairs_freq_exclusion$PairID[pairs_freq_exclusion$Sign == 1]
pairs_ID_loser <- pairs_freq_exclusion$PairID[pairs_freq_exclusion$Sign == -1]
pairs_ID_others <- pairs_outcome$PairID[pairs_outcome$outcome %in% c("3-coexistence", "4-coexistence", "5-inconclusive")]

pairs_outcome <- pairs_outcome %>%
    mutate(From = case_when( # For network
        PairID %in% pairs_ID_winner ~ Isolate1,
        PairID %in% pairs_ID_loser ~ Isolate2,
        PairID %in% pairs_ID_others ~ Isolate1,
    )) %>%
    mutate(To = case_when( # For network
        PairID %in% pairs_ID_winner ~ Isolate2,
        PairID %in% pairs_ID_loser ~ Isolate1,
        PairID %in% pairs_ID_others ~ Isolate2,
    )) %>% # For plotting the frequency
    mutate(Isolate1IsLoser = case_when(
        PairID %in% pairs_ID_winner ~ F,
        PairID %in% pairs_ID_loser ~ T,
        PairID %in% pairs_ID_others ~ NA,
    ))


write_csv(pairs_outcome, paste0(folder_data, "temp/26-pairs_outcome.csv"))
