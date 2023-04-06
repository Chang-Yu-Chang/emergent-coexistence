library(tidyverse)
library('ggstats')

# ESVs
isol <- fread('isolates_remained.csv')
isol <- isol[,c('Community', 'Isolate', 'ESV', 'BasePairMismatch', 'ID')]

# 
communities_abundance <- fread("Emergent_Comunity_Data.csv")
communities_abundance[, Community := factor(paste0("C", Inoculum, "R", Replicate),
                                            paste0("C", rep(1:12, each = 8), "R", rep(1:8, 12)))]

communities_abundance <- communities_abundance[, .(ESV_ID, ESV, Relative_Abundance,
                                                   Carbon_Source, Community, Transfer)]
communities_abundance <- communities_abundance[Carbon_Source == "Glucose" | Transfer %in% 9:12]

# average 4 last transfers
eq_abundance <- communities_abundance[, .(eq_abundance = mean(Relative_Abundance, na.rm=TRUE)),
                                          by=c('ESV_ID', 'ESV', 'Community')]

# merge by ESV sequence
isol <- merge(isol, eq_abundance, by=c('ESV', 'Community'), all.x=TRUE)
isol[, ESV:=NULL]

# determine coexistence score for each isolate
outcomes <- fread('pairs_outcome.csv')
outcomes <- outcomes[outcome %like% 'coexistence' | outcome %like% 'exclusion']
outcomes[, outc := ifelse(outcome %like% 'coexistence', 'coexistence', 'exclusion')]

# compute with averaged ranks 
isol[, rank_ab := rank(-eq_abundance), by='Community']

isol1 <- isol
names(isol1)[c(2, 7)] <- c('Isolate1', 'rank_ab_1')
isol1 <- isol1[, c(1, 2,7)]

isol2 <- isol
names(isol2)[c(2, 7)] <- c('Isolate2', 'rank_ab_2')
isol2 <- isol2[, c(1,2,7)]

res <- merge(outcomes, isol1, by=c('Community', 'Isolate1'), all.x=TRUE)
res <- merge(res, isol2, by=c('Community', 'Isolate2'), all.x=TRUE)

res[, involves_top := rank_ab_1<2 | rank_ab_2<2]


# bootstrap because 1 ESV sometimes maps to multiple isolates 
r <- matrix(NA, ncol=3, nrow=1000)
# sample one isolate/ESV
for (i in 1:1000) { 

    isol[, rank_ab := rank(-eq_abundance, ties.method='random'), by='Community']

    isol1 <- isol
    names(isol1)[c(2, 7)] <- c('Isolate1', 'rank_ab_1')
    isol1 <- isol1[, c(1, 2,7)]
    
    isol2 <- isol
    names(isol2)[c(2, 7)] <- c('Isolate2', 'rank_ab_2')
    isol2 <- isol2[, c(1,2,7)]

    res <- merge(outcomes, isol1, by=c('Community', 'Isolate1'), all.x=TRUE)
    res <- merge(res, isol2, by=c('Community', 'Isolate2'), all.x=TRUE)

    res[, involves_top := rank_ab_1<2 | rank_ab_2<2]   

    res <- res[,table(involves_top, outc)]
    
    r[i,1] <- res[2,1]/sum(res[2,])
    r[i,2] <- res[1,1]/sum(res[1,])
    r[i, 3] <- fisher.test(res)$p.value
}

colnames(r) <- c('Rank = 1', 'Rank > 1', 'p')
res <- as.data.table(r)
res[, id:=1:1000]
res <- melt(res[,-3], id.vars='id')
res <- res[, .(m = mean(value), sd = sd(value)), by='variable']

p <- ggplot(res, aes(value, fill=variable)) +
  geom_histogram(bins=15, color='black') +
  theme_classic() + xlim(0, .4) +
  labs(x = "% Coexistence among pairs", y = 'Count (# bootstrap samples)') +
  scale_fill_discrete(name = 'Pairs containing:')

ggsave(here::here('Fig_S15.png'), p, width = 6, height = 4)

