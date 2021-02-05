#' Read and match isolates' 16S to pairs
library(tidyverse)
library(data.table)
# Read data ----
## Read isolates' RDP and ID match
isolates_RDP <- fread(here::here("data/temp/isolates_RDP.csv"))
isolates_ID_match <- fread(here::here("data/temp/isolates_ID_match.csv"))

## Read pairs
pairs_ID <- fread(here::here("data/temp/pairs_ID.csv"))
#pairs_freq_func <- fread(here::here("data/temp/pairs_freq_func.csv"))

## Match isolates' information
isolates_RDP_ID <- isolates_ID_match %>%
  left_join(isolates_RDP, by = c("ID")) %>%
  filter(Community %in% communities_name) %>%
  select(ExpID, ID, Community, Isolate, Family, Genus, GenusScore, Fermenter, Sequence) %>%
  as_tibble()


# Match the isolates taxonomic information to pairs ----

## R function for computing pairwise alignment difference of two sequences
sequence1 <- isolates_RDP_ID$Sequence[1]
sequence2 <- isolates_RDP_ID$Sequence[2]
get_aligment_difference <- function(sequence1, sequence2, print = FALSE, type = "local") {
  algn <- Biostrings::pairwiseAlignment(sequence1, sequence2, type = "local")

  temp_df <- tibble(
    Sequence1 = as.character(Biostrings::pattern(algn)) %>% strsplit('') %>% unlist(),
    Sequence2 = as.character(IRanges::subject(algn))  %>% strsplit('') %>% unlist()
  )

  bp_diff <- temp_df %>%
    filter(Sequence1 != Sequence2, Sequence1 != "N", Sequence2 != "N") %>%
    nrow()

  if (print == TRUE) return(algn) else return(c(bp_diff, nchar(as.character(Biostrings::pattern(algn)))))
}

get_aligment_difference(isolates_RDP_ID$Sequence[1], isolates_RDP_ID$Sequence[2], print = F)


## Loop through 183 pairs and compute 16S base pair difference. Skip C10R2 isolate 3 that has no sequence
pairs_ID$SeqDifference <- NA
pairs_ID$SeqLength <- NA

for (i in 1:nrow(pairs_ID)) {
  seq1 <- isolates_RDP_ID$Sequence[which(isolates_RDP_ID$Community == pairs_ID$Community[i] & isolates_RDP_ID$Isolate == pairs_ID$Isolate1[i])]
  seq2 <- isolates_RDP_ID$Sequence[which(isolates_RDP_ID$Community == pairs_ID$Community[i] & isolates_RDP_ID$Isolate == pairs_ID$Isolate2[i])]

  if (!is.na(seq1) & !is.na(seq2)) {
    align_diff <- get_aligment_difference(seq1, seq2)
    pairs_ID$SeqDifference[i] <- align_diff[1]
    pairs_ID$SeqLength[i] <- align_diff[2]
  }
  cat(paste0(i, " "))
}

pairs_16S <- pairs_ID

#
fwrite(pairs_16S, file = here::here("data/temp/pairs_16S.csv"))
