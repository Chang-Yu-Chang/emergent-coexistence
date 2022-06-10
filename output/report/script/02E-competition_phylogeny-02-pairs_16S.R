#' Read and match isolates' 16S to pairs
library(tidyverse)

isolates_RDP <- read_csv("~/Dropbox/lab/emergent-coexistence/data/temp/isolates_RDP.csv", col_types = cols())
isolates_ID_match <- read_csv("~/Dropbox/lab/emergent-coexistence/data/temp/isolates_ID_match.csv", col_types = cols())
pairs_ID <- read_csv("~/Dropbox/lab/emergent-coexistence/data/temp/pairs_ID.csv", col_types = cols())
communities <- read_csv("~/Dropbox/lab/emergent-coexistence/data/output/communities.csv", col_types = cols())


# Match isolates' information
isolates_RDP_ID <- isolates_ID_match %>%
  left_join(isolates_RDP) %>%
  filter(Community %in% communities$Community) %>%
  select(ExpID, ID, Community, Isolate, Family, Genus, GenusScore, Fermenter, Sequence) %>%
  as_tibble()


# Match the isolates taxonomic information to pairs
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
pairs_16S <- pairs_ID %>%
    left_join(rename_with(isolates_ID_match, ~ paste0(., "1"), !contains("Community"))) %>%
    left_join(rename_with(isolates_ID_match, ~ paste0(., "2"), !contains("Community"))) %>%
    left_join(rename_with(isolates_RDP, ~ paste0(., "1"), !contains("Community"))) %>%
    left_join(rename_with(isolates_RDP, ~ paste0(., "2"), !contains("Community"))) %>%
    rowwise() %>%
    mutate(temp = get_aligment_difference(Sequence1, Sequence2) %>% list) %>%
    mutate(SequenceDifference = temp[1], SequenceLength = temp[2]) %>%
    select(PairID, Community, Isolate1, Isolate2, SequenceDifference, SequenceLength)

#
write_csv(pairs_16S, "~/Dropbox/lab/emergent-coexistence/data/temp/pairs_16S.csv")
















