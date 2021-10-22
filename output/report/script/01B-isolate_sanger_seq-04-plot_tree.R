library(tidyverse)
library(ggtree)

load(here::here("data/temp/isolates_sanger_seq.Rdata"))
isolates_seq

p <- ggtree(tree2)

ggsave(here::here("output/report/figure/01B-tree.png"), p)

