pairs <- fread('pairs_remained.csv')
pairs <- pairs[,c(1:4, 9, 11, 12, 17, 19, 20, 23, 24, 25)]

outc <- fread('pairs_outcome.csv')

pairs <- merge(outc, pairs, all.x=TRUE)
pairs[, out := ifelse(outcome %like% 'coexistence', 'coexistence',
               ifelse(outcome %like% 'excl', 'exclusion', 'inconclusive'))]


# Fig S17
png('Fig_S17.png', 500, 400)
ggplot(pairs[out!='inconclusive'], aes(x = Mismatch, fill=out)) +
  geom_histogram(alpha=.4, position="identity", bins=20) + 
  theme_classic() + labs(x='Mismatch (bp)', y=' Number of pairs') +
  scale_fill_manual(values=c('blue', 'red'), name = "Outcome")
dev.off()

with(pairs[out!='inconclusive'], wilcox.test(Mismatch~out))

