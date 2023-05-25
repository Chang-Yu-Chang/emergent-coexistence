library(tidyverse)
library(ggstats)
library(data.table)
source(here::here("processing_scripts/00-metadata.R"))

# ESVs
#pairs <- fread('pairs_remained.csv')
pairs <- fread(paste0(folder_data, "output/pairs_remained.csv"))
pairs <- pairs[,c(1:4, 9, 11, 12, 17, 19, 20, 23, 24, 25)]

#outc <- fread('pairs_outcome.csv')
outc <- fread(paste0(folder_data, 'temp/26-pairs_outcome.csv'))

pairs <- merge(outc, pairs, all.x=TRUE)
pairs[, out := ifelse(outcome %like% 'coexistence', 'coexistence',
               ifelse(outcome %like% 'excl', 'exclusion', 'inconclusive'))]


# Fig S17
p <- ggplot(pairs[out!='inconclusive'], aes(x = Mismatch, fill=out)) +
  geom_histogram(alpha=.4, position="identity", bins=20) +
  theme_classic() + labs(x='Mismatch (bp)', y=' Number of pairs') +
  scale_fill_manual(values=c('blue', 'red'), name = "Outcome")

ggsave(here::here('plots/FigS17.png'), p, width = 6, height = 4)

with(pairs[out!='inconclusive'], wilcox.test(Mismatch~out))

