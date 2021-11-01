#' Pick random isolates for experiment

library(tidyverse)
library(data.table)
isolates <- read_csv(here::here("data/output/isolates.csv"))
jean_isolates_RDP <- read_csv(here::here("data/temp/jean_isolates_RDP.csv"))

# Read alignment data
# pairs_sanger_alignment <- read_csv(here::here("data/temp/pairs_sanger_alignment.csv"))
# jean_pairs_sanger_alignment <- read_csv(here::here("data/temp/jean_pairs_sanger_alignment.csv"))

# Randomly pick isolates ----
# 16 isolates from self-assembled community isolates. 8 isolates in each of the new assembled communities
# Randomly pick 8 from 13 communities, and pick one isolate in each picked community
## First community
temp_df <- data.frame(ExpID = c("1.4.A.1", "1.6.A.5", "2.6.A.3", "2.8.A.1", "8.4.A.2", "10.2.A.1", "11.1.B.1", "11.2.B.2"), stringsAsFactors = F)

## Repeat it for the second community
temp_df2 <- data.frame(ExpID = c("1.6.A.6", "1.7.A.1.2", "2.8.A.3", "4.1.A.2", "8.4.A.1", "10.2.A.2", "11.1.B.3", "1.2.A.5"), stringsAsFactors = F)

# Randomly pick another 16 isolates from Jean's isolates
jean_isolates_RDP <- jean_isolates_RDP %>%
  filter(Genus %in% c("Enterobacter", "Pseudomonas", "Pantoea", "Citrobacter", "Raoultella", "Kluyvera")) %>%
  mutate(ExpID = ordered(ExpID, paste0("JVN", 1:100)))

## Random community 1
temp_df3 <-  jean_isolates_RDP %>%
  filter(ExpID %in% paste0("JVN", c(1,2,7,8,21,25,30, 46))) %>% # 30 and 46 are Pseudomonas
  arrange(ExpID) %>% as_tibble()

## Random community 2
temp_df4 <- jean_isolates_RDP %>%
  filter(ExpID %in% paste0("JVN", c(3, 11, 12, 14, 23, 28, 35, 41))) %>% # 14, 23, 35 are Pseudomonas
  arrange(ExpID) %>% as_tibble()


# Rbind the new isolates list ----
## Each community assembly
isolates_assembly_across <-
  temp_df %>%
  bind_rows(temp_df2) %>%
  left_join(isolates, by = c("ExpID")) %>%
  mutate(Assembly = "across-community",
    Community = rep(paste0("AcrAss", 1:2), each = 8)) %>%
  mutate(Isolate = rep(1:8, 2)) %>%
  select(Assembly, ExpID, Community, Isolate, Family, Genus)

isolates_assembly_random <-
  temp_df3 %>%
  bind_rows(temp_df4) %>%
  mutate(Assembly = rep("random-assembly", 16),
    Community = rep(paste0("RanAss", 1:2), each = 8),
    Isolate = rep(1:8, 2)) %>%
    select(Assembly, ExpID, Community, Isolate, Family, Genus)

#
isolates_random <- isolates_assembly_across %>%
  bind_rows(isolates_assembly_random) %>%
  replace_na(replace = list(Isolate = ""))

#
isolates_random$Fermenter <- ifelse(isolates_random$Family %in% c("Pseudomonadaceae", "Moraxellaceae", "Xanthomonadaceae", "Alcaligenaceae", "Comamonadaceae"), F, ifelse(isolates_random$Family %in% c("Enterobacteriaceae", "Aeromonadaceae"), T, NA))

write_csv(isolates_random, file = here::here("data/output/isolates_random.csv"))


if (FALSE) {


# Check the pairwise alignment
pairs_assembly <-
isolates_assembly  %>%
  group_by(AssemblyCommunity) %>%
  group_split() %>%
  lapply(function(x){
    x %>%
      pull(ExpID) %>%
      combn(2) %>% t() %>%
      as.data.frame() %>% setNames(c("ExpID1", "ExpID2")) %>%
      mutate(AssemblyCommunity = unique(x$AssemblyCommunity))
  }) %>%
  rbindlist()

pairs_assembly %>%
  left_join(bind_rows(pairs_sanger_alignment, jean_pairs_sanger_alignment)) %>%
  filter(AlignmentType == "local") %>%
  # At least 2 base pair difference
  filter(BasePairMismatch <= 1) %>%
  as_tibble()

}
