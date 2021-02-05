#'  This script takes the determined outcomes of pairwise competition by chromogenic plate and CASEU pilot2
#'  and returns a list of ambiguous pairs that have to be sequenced at CASEU pilot3.


# Pairs clarified from chromogenic agar plates ----
pairs_chromogenic <-
# fread(root$find_file("output/protocol/tab_fig/protocol_20190912_plate_ambiguous_pairs-plating_pairs.csv")) %>%
  fread(root$find_file("data/raw/pairwise_competition/result_pairwise_competition_chromogenic.csv")) %>%
  filter(!is.na(ColonyCount), ColonyCount != 0) %>%
  distinct() %>%
  select(Community, Isolate1, Isolate2, Isolate1Freq, Isolate2Freq)

temp_index <- pairs_chromogenic$Isolate1 > pairs_chromogenic$Isolate2
pairs_chromogenic[temp_index, c("Isolate1", "Isolate2", "Isolate1Freq", "Isolate2Freq")] <- pairs_chromogenic[temp_index, c("Isolate2", "Isolate1", "Isolate2Freq", "Isolate1Freq")]

# Pairs clarified from CASEU pilot2 Sanger seqeuncing ----
pairs_CASEU_pilot2 <- fread(root$find_file("data/temp/CASEU_pilot2.csv")) %>%
  filter(Treatment == "C11R1") %>%
  mutate(Community = Treatment, Isolate1 = as.numeric(Isolate1), Isolate2 = as.numeric(Isolate2),
    Isolate1Freq = Isolate1Freq * 100, Isolate2Freq = Isolate2Freq * 100) %>%
  select(Community, Isolate1, Isolate2, Isolate1Freq, Isolate2Freq, Isolate1FreqPredicted, Isolate2FreqPredicted)

temp_index <- pairs_CASEU_pilot2$Isolate1 > pairs_CASEU_pilot2$Isolate2
pairs_CASEU_pilot2[temp_index, c("Isolate1", "Isolate2", "Isolate1Freq", "Isolate2Freq")] <- pairs_CASEU_pilot2[temp_index, c("Isolate2", "Isolate1", "Isolate2Freq", "Isolate1Freq")]


# The rest ambiguous pairs for CASEU pilot3 sequencing ----
# pairs_ambiguous_sequencing <- pairs_ambiguous %>%
#   anti_join(pairs_chromogenic, by = c("Community", "Isolate1", "Isolate2", "Isolate1Freq", "Isolate2Freq")) %>%
#   anti_join(pairs_CASEU_pilot2, by = c("Community", "Isolate1", "Isolate2", "Isolate1Freq", "Isolate2Freq")) %>%
#   arrange(Community, Isolate1, Isolate2)


"
Check manual on TSA plate and make a manual seuqenicng table for CASEU pilot3
Only sequence a pairs with all of its three initial frequencies
"

temp_df <- data.frame(
  Community = c(rep("C11R2", 5), rep("C11R5", 2), rep("C1R4", 3), rep("C1R6", 2), "C4R1", "C7R1"),
  Isolate1 = c(3,  3,  7,  6, 2, 1,2,2,2,4,2,1,1,3),
  Isolate2 = c(7, 12, 12, 11, 8, 3,4,4,5,5,5,2,3,4)
)

pairs_ambiguous_sequencing <- rbind(temp_df, temp_df, temp_df) %>%
  mutate(Isolate1Freq = rep(c(5,50,95), each = nrow(temp_df)),
    Isolate2Freq = rep(c(95,50,5), each = nrow(temp_df)))




# Match ambiguous pairs to well positon on plates ----
## Read DW96 plate layout
plate_df <- fread(root$find_file("data/temp/plate_df.csv"))

## Switch the isolate1 and isolate2 since the P1 is 50:50 and
## rows and columns in P2 P3 are for 95 and 5 respectively
pairs_ambiguous_sequencing[pairs_ambiguous_sequencing$Isolate1Freq == 5, c("Isolate1", "Isolate2", "Isolate1Freq", "Isolate2Freq")] <- pairs_ambiguous_sequencing[pairs_ambiguous_sequencing$Isolate1Freq == 5, c("Isolate2", "Isolate1", "Isolate2Freq", "Isolate1Freq")]

## Mutate the isolate variables
pairs_ambiguous_sequencing <- pairs_ambiguous_sequencing %>% mutate(Isolate1 = factor(Isolate1), Isolate2 = factor(Isolate2))

## Join plates
pairs_ambiguous_sequencing_on_DW96 <- pairs_ambiguous_sequencing %>%
  left_join(plate_df, by = c("Community", "Isolate1", "Isolate2", "Isolate1Freq", "Isolate2Freq")) %>%
  mutate(Plate = ifelse((Isolate1Freq != 50 & PlateLayout != "C11R1") , "P2", "P1"),
    Community = ordered(Community, community_name)) %>%
  distinct(Community, Isolate1, Isolate2, Isolate1Freq, Isolate2Freq, .keep_all = T) %>%
  select(-MixIsolate)
