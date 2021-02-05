#' Compare pairwise phylogenetic identities to the outcomes of pairwise competition.
#' The pairwise phylogenetic identities includes:
#' 1. RDP taxonomy
#' 2. 16S
#' 3. Branch distances on a tree
library(tidyverse)
library(data.table)

pairs <- fread(here::here("data/output/pairs.csv"))
interaction_type <- c("exclusion", "coexistence", "neutrality", "mutual exclusion", "frequency-dependent\ncoexistence")
myColor = c("#DB7469", "#557BAA", "#8650C4", "red", "blue")
names(myColor) <- interaction_type


# Plot the pairs count----
## RDP Fermenter vs. fraction of coexistence ----
### Fraction
pairs_interaction_fermenter <- pairs %>%
  filter(PairFermenter != "") %>%
  filter(InteractionType != "", !is.na(PairFermenter)) %>%
  filter(InteractionType %in% c("coexistence", "exclusion")) %>%
  group_by(PairFermenter, InteractionType) %>%
  summarize(Count = n()) %>%
  spread(InteractionType, Count) %>%
  mutate(FracCoext = coexistence / (coexistence + exclusion))


### Plot
my_factor_fermenter <- apply(pairs_interaction_fermenter[,c("coexistence", "exclusion")], 1, sum) %>% max
p_pairs_interaction_fermenter <- pairs %>%
  filter(PairFermenter != "") %>%
  filter(InteractionType != "", !is.na(PairFermenter)) %>%
  filter(InteractionType %in% c("coexistence", "exclusion")) %>%
  ggplot() +
  geom_bar(aes(x = PairFermenter, fill = InteractionType, group = InteractionType), color = "black", position = "stack") +
  geom_line(data = pairs_interaction_fermenter, aes(x = PairFermenter, y = FracCoext * my_factor_fermenter, group = 1)) +
  geom_point(data = pairs_interaction_fermenter, aes(x = PairFermenter, y = FracCoext * my_factor_fermenter)) +
  scale_x_discrete(labels = c("Fermenter\nFermenter", "Fermenter\nNon-fermenter", "Non-fermenter\nNon-fermenter")) +
  scale_y_continuous(sec.axis = sec_axis(trans = ~./my_factor_fermenter, name = "fraction of coexistence"), limits = c(0, my_factor_fermenter)) +
  scale_fill_manual(values = myColor) +
  theme_bw() +
  theme(legend.position = "top", legend.title = element_blank()) +
  labs(x = "")

## RDP Family vs. fraction of coexistence ----
### Fraction
pairs_interaction_family <- pairs %>%
  filter(PairFamily != "") %>%
  filter(InteractionType != "", !is.na(PairFamily)) %>%
  filter(InteractionType %in% c("coexistence", "exclusion")) %>%
  group_by(PairFamily, InteractionType) %>%
  summarize(Count = n()) %>%
  spread(InteractionType, Count) %>%
  mutate(FracCoext = coexistence / (coexistence + exclusion))


### Plot
my_factor_family <- apply(pairs_interaction_family[,c("coexistence", "exclusion")], 1, sum) %>% max
p_pairs_interaction_family <- pairs %>%
  filter(PairFamily != "") %>%
  filter(InteractionType != "", !is.na(PairFamily)) %>%
  filter(InteractionType %in% c("coexistence", "exclusion")) %>%
  ggplot() +
  geom_bar(aes(x = PairFamily, fill = InteractionType, group = InteractionType), color = "black", position= "stack") +
  geom_line(data = pairs_interaction_family, aes(x = PairFamily, y = FracCoext * my_factor_family, group = 1)) +
  geom_point(data = pairs_interaction_family, aes(x = PairFamily, y = FracCoext * my_factor_family)) +
  scale_x_discrete(labels = c("Enterobacteriaceae\nPseudomonadaceae", "Enterobacteriaceae\nPseudomonadaceae", "Pseudomonadaceae\nPseudomonadaceae")) +
  scale_y_continuous(sec.axis = sec_axis(trans = ~./my_factor_family, name = "fraction of coexistence"), limits = c(0, my_factor_family)) +
  scale_fill_manual(values = myColor) +
  theme_bw() +
  theme(legend.position = "top", legend.title = element_blank()) +
  labs(x = "")


## 16S sequence difference vs. fraction of coexistence----
### Cut the continuous sequence difference into categories
pairs_interaction_16S_cut <-  pairs %>%
  filter(!is.na(SeqDifference)) %>%
  filter(InteractionType %in% c("coexistence", "exclusion")) %>%
  mutate(InteractionType = factor(InteractionType)) %>%
  #  mutate(SeqDifference = cut(SeqDifference, breaks = c(seq(0, 30, 10), seq(60, 120, 10), 150, Inf))) %>%
  mutate(SeqDifference = cut(SeqDifference, breaks = seq(0, 160, 10), include.lowest = T)) %>%
  group_by(SeqDifference, InteractionType) %>%
  summarize(Count = n()) %>%
  spread(InteractionType, Count, fill = 0) %>%
  mutate(FracCoext = coexistence / (coexistence + exclusion)) %>%
  mutate(Count = sum(coexistence, exclusion))

### Plot
p_pairs_interaction_16S_cut <- pairs_interaction_16S_cut %>%
  gather("yAxisData", "yAxisValue", 4:5) %>%
  ggplot(aes(x = SeqDifference, y = yAxisValue)) +
  geom_col(col = "black", fill = "white") +
#  scale_y_continuous(limits = c(0,1)) +
  facet_grid(yAxisData~., scales = "free_y") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90)) +
  labs(x = "pairwise differernces in 16S sequences (bp)", y = "fraction of coexistence")

if (FALSE) {


## Tree distance vs. fraction of coexistence ----
### Cut the continuous sequence difference into categories
pairs_interaction_tree_cut <-
  pairs %>%
  filter(!is.na(PairTreeDistance)) %>%
  filter(InteractionType %in% c("coexistence", "exclusion")) %>%
  mutate(InteractionType = factor(InteractionType)) %>%
  mutate(PairTreeDistance = cut(PairTreeDistance, breaks = seq(-0.1, 0.5, .05), include.lowest = T)) %>%
  group_by(PairTreeDistance, InteractionType) %>%
  summarize(Count = n()) %>%
  spread(InteractionType, Count, fill = 0) %>%
  mutate(FracCoext = coexistence / (coexistence + exclusion)) %>%
  mutate(Count = sum(coexistence, exclusion))

### Plot
p_pairs_interaction_tree_cut <- pairs_interaction_tree_cut %>%
  gather("yAxisData", "yAxisValue", 4:5) %>%
  ggplot(aes(x = PairTreeDistance, y = yAxisValue)) +
  geom_col(col = "black", fill = "white") +
  #  scale_y_continuous(limits = c(0,1)) +
  facet_grid(yAxisData~., scales = "free_y") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90)) +
  labs(x = "pairwise branch distance", y = "fraction of coexistence")

p_pairs_interaction_tree_cut
}

