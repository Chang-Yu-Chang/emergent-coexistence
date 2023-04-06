
# 
pairs <- fread("pairs.csv")
pairs_freq <- fread("93-pairs_freq.csv")

boots_T0 <- fread("06-pairs_T0_boots.csv")
boots_T0 <- boots_T0[, .(Community, Isolate1, Isolate2,
                         Isolate1InitialODFreq, Isolate2InitialODFreq, BootstrapID, Isolate1CFUFreq)]

boots_T8 <- fread("06-pairs_T8_boots.csv")



# 
pairs_I1 <- dcast(pairs_freq, PairFreqID + Community + Isolate1 + Isolate2 + Isolate1InitialODFreq ~ Time,
                  value.var = "Isolate1CFUFreqMean")
names(pairs_I1)[6:7] <- c('T0_I1_mean', 'T8_I1_mean')

pairs_sd <- dcast(pairs_freq, PairFreqID + Community + Isolate1 + Isolate2 + Isolate1InitialODFreq ~ Time,
                  value.var = "Isolate1CFUFreqSd")
names(pairs_sd)[6:7] <- c('T0_I1_SD', 'T8_I1_SD')

pairs_I1 <- merge(pairs_I1, pairs_sd)

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

isol <- fread('pairs_remained.csv')
pairs_outcome <- merge(isol[, 2:4], pairs_outcome, all.x=TRUE)

for (i in 1:nrow(pairs_outcome)) {

    aux <- pairs_I1[Community==pairs_outcome$Community[i] &
                    Isolate1==pairs_outcome$Isolate1[i] &
                    Isolate2==pairs_outcome$Isolate2[i]]
    aux <- aux[order(Isolate1InitialODFreq)]    
    
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
fwrite(pairs_outcome, 'pairs_outcome.csv')

