#' Relative abundances within pairs at T8
#' `data/temp/pairs_freq_taxonomy` is from the R script `output/report/script/02-pairs-08-pairs_interaction_taxonomy.R`
library(tidyverse)
library(data.table)
pairs_melted <- fread(here::here("data/output/pairs_melted.csv"))

# Switch fermenter or family ----
## Fermenter and non-fermenter pairs
#' Switch E and R so that in fermenter and non-fermenter pairs, isolate1 is
#' always fermenter and isolate 2 is always non-fermenter
pairs_freq_taxonomy_fermenter <-
  pairs_melted %>%
  mutate(Isolate2MeasuredFreq = 1 - Isolate1MeasuredFreq) %>%
  as.data.table()

temp_index <- pairs_freq_taxonomy_fermenter$Fermenter1 == F & pairs_freq_taxonomy_fermenter$Fermenter2 == T

pairs_freq_taxonomy_fermenter[temp_index, c("Isolate1", "Isolate2", "Isolate1InitialODFreq", "Isolate2InitialODFreq",
  "Isolate1MeasuredFreq", "Isolate2MeasuredFreq", "Fermenter1", "Fermenter2", "Family1", "Family2", "Genus1", "Genus2")] <-
  pairs_freq_taxonomy_fermenter[temp_index, c("Isolate2", "Isolate1", "Isolate2InitialODFreq", "Isolate1InitialODFreq",
    "Isolate2MeasuredFreq", "Isolate1MeasuredFreq", "Fermenter2", "Fermenter1", "Family2", "Family1", "Genus2", "Genus1")]


## Enterobacteriaceae and Pseudomonadaceae pairs
#' Switch E and R so that in EP pairs, isolate1 is always Enterobacteriaceae and isolate 2 is always Pseudomonadaceae
pairs_freq_taxonomy_family <-
  pairs_melted %>%
  mutate(Isolate2MeasuredFreq = 1 - Isolate1MeasuredFreq) %>%
  as.data.table()

temp_index <- pairs_freq_taxonomy_family$Family1 == "Pseudomonadaceae" & pairs_freq_taxonomy_family$Family2 == "Enterobacteriaceae"
pairs_freq_taxonomy_family[temp_index, c("Isolate1", "Isolate2", "Isolate1InitialODFreq", "Isolate2InitialODFreq",
  "Isolate1MeasuredFreq", "Isolate2MeasuredFreq", "Fermenter1", "Fermenter2", "Family1", "Family2", "Genus1", "Genus2")] <-
  pairs_freq_taxonomy_family[temp_index, c("Isolate2", "Isolate1", "Isolate2InitialODFreq", "Isolate1InitialODFreq",
    "Isolate2MeasuredFreq", "Isolate1MeasuredFreq", "Fermenter2", "Fermenter1", "Family2", "Family1", "Genus2", "Genus1")]


#
temp_df <- pairs %>% as_tibble() %>%
  filter(InteractionType == "coexistence") %>%
  select(Community, Isolate1, Isolate2, InteractionType) %>%
  unite(col = "temp", Community, Isolate1, Isolate2)

pairs_freq_taxonomy_fermenter %>% as_tibble() %>%
  filter(Time == "T8") %>%
  # Filter only the coexitence pairs
  unite(col = "temp", Community, Isolate1, Isolate2) %>%
  right_join(temp_df, by = "temp") %>%
  separate(col = temp, into = c("Community", "Isolate1", "Isolate2")) %>%
  #
  filter(!Isolate1MeasuredFreq %in% c(0,1)) %>%
  ggplot() +
  geom_histogram(aes(x = Isolate1MeasuredFreq), col = "black", fill = "white") +
  facet_grid(PairFermenter~., scale = "free") +
  theme_bw() +
  labs("frequencies of isolate 1 at T8")



save(pairs_freq_taxonomy_fermenter, pairs_freq_taxonomy_family, file = here::here("data/temp/pairs_hist_plots.Rdata"))



