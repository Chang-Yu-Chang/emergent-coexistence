#' Match isolates' 16S sequences to community ESV
library(tidyverse)
library(data.table)


# Read data ----
communities_abundance_syl <- read_csv(here::here("data/temp/communities_abundance_syl.csv")) # Data from Sylvie
isolates_ID_match <- read_csv(here::here("data/temp/isolates_ID_match.csv"))
isolates_RDP <- read_csv(here::here("data/temp/isolates_RDP.csv"))
isolates_ID_RDP <- isolates_ID_match %>% filter(!grepl("Ass", Community)) %>% left_join(isolates_RDP) %>% select(ExpID, ID, Community, Isolate, Family, Genus, Sequence) %>% as_tibble()

if (FALSE) {
# C2R4 isolates, hardcoded
temp <- list.files(here::here("data/raw/Sanger/Sanger_ES_C2R4_SE"), "merged.fasta", full.names = T)
temp_df <- tibble(ExpID = c("E", "P", "R", "C"), ID = c(160, 162, 165, 169), Community = "C2R4",
  Family = c("Enterobacteriaceae", "Pseudomonadaceae", "Pseudomonadaceae", "Enterobacteriaceae"),
  Genus = c("Enterobacter", "Pseudomonas", "Raoultella", "Citrobacter"), Sequence =
    c(seqinr::read.fasta(temp[1], as.string = T, seqonly = T, seqtype = "DNA")[[1]],
      seqinr::read.fasta(temp[2], as.string = T, seqonly = T, seqtype = "DNA")[[1]],
      seqinr::read.fasta(temp[3], as.string = T, seqonly = T, seqtype = "DNA")[[1]],
      seqinr::read.fasta(temp[4], as.string = T, seqonly = T, seqtype = "DNA")[[1]]))

isolates_Sanger_ID <- fread(here::here("output/report/data/01E-isolates_SangerID_SE.csv")) %>% transmute(ID = SangerID); isolates_Sanger_ID$ID[isolates_Sanger_ID$ID == 523] <- 162
isolates_RDP <- isolates_RDP %>% bind_rows(temp_df)


  isolates_ID_match %>%
  left_join(isolates_RDP) %>%
  as_tibble()

}

#' To use community ESVs, here are some criteria for sequence alignment:
#' 1. Only use ESVs from transfer 12
#' 2. Use the communities assembled in glucose media
#' 3. Remove duplicates
if (FALSE) {
communities_abundance_align <- communities_abundance %>%
  filter(!grepl("AA", SampleID)) %>%
  filter(Transfer == 12) %>%
  mutate(Community = ordered(Community, levels = c(paste0("C", 1:12, "Rpool"), paste0("C", rep(1:12, each = 8), "R", rep(1:8, 12))))) %>%
  group_by(SampleID, Community) %>%
  arrange(Community, SampleID)
}

communities_abundance_syl_align <- communities_abundance_syl %>%
  filter(CarbonSource == "Glucose", Transfer == 12)


# Match isolate Sanger seq to community ESV ----
## R function for pairwise alignment and counting bp gap and mismatch
alignment_bp_count <- function(seq1, seq2, type = "local"){
  # Make pairwise alignment
  algn <- Biostrings::pairwiseAlignment(seq1, seq2, type = type)

  # Make a df
  temp_df <- tibble(
    s1 = strsplit(as.character(Biostrings::pattern(algn)), '')[[1]],
    s2 = strsplit(as.character(Biostrings::subject(algn)), '')[[1]]
  )

  # Find difference; ignore N (which means any nucleotide)
  temp_df <- temp_df %>% mutate(Diff = ((s1 != s2) & s1 != 'N' & s2 != 'N'))

  # Find gaps
  ngap <- sum(temp_df$s1 == "-" | temp_df$s2 == "-")

  # Find mismatch
  nmis <- sum((temp_df$s1 != temp_df$s2) & temp_df$s1 != "-" & temp_df$s2 != "-" & temp_df$s1 != "N" & temp_df$s2 != "N")

  # Length of consensus
  nlength <- nrow(temp_df)

  #
  c(ConsensusLength = nlength, BasePairGap = ngap, BasePairMismatch = nmis, AlignmentScore = algn@score) %>%
    return()
}
## Example code for match one pair of sequences
# seq1 <- isolates_ID_RDP$Sequence[1]
# seq2 <- isolates_ID_RDP$Sequence[3]
# alignment_bp_count(seq1, seq2) # R function from `sequence_pairwise_alignment.R`
# printPairwiseAlignment(Biostrings::pairwiseAlignment(seq1, seq2))

## R function for align isolates' Sanger sequences to ESVs within communities, based on different alignment methods
align_sequence <- function(communities_abundance_align, isolates_ID_RDP, alignment_type = "global") {
  # Community names. All communties with isolates
  comm_name <- unique(isolates_ID_RDP$Community)

  # Empty list for alignment
  sequences_alignment_list_func <- rep(list(NA), length(comm_name)) %>% setNames(comm_name)

  # Loop through comunities
  for (i in 1:length(comm_name)) {
    # Subset for community ESV
    temp_comm <- communities_abundance_align %>%
      filter(Community == comm_name[i]) %>%
      select(SampleID, Community, Family, RelativeAbundance, CommunityESVID, ESV)

    # Isolate sequences
    temp_iso <- isolates_ID_RDP %>%
      filter(Community == comm_name[i]) %>%
      # Change variable name used in isolates to avoid duplicate variable name
      mutate(IsolateID = ID, IsolateGenus = Genus, IsolateSangerSequence = Sequence) %>%
      select(Community, ExpID, IsolateID, IsolateGenus, IsolateSangerSequence)

    # Merge two dfs: it will auto create df with nrow(temp_df) = nrow(temp_comm) * nrow(temp_iso)
    temp_df <- temp_comm %>% left_join(temp_iso, by = "Community")
    cat("\n", comm_name[i], "\n number of ESVs = ", nrow(temp_comm), "\n number of isolates = ", nrow(temp_iso), "\n pairs = ", nrow(temp_df), "\n")

    # Add five new columns to temp_df
    temp_df$AlignmentType <- "NA"
    temp_df[,c("ConsensusLength", "BasePairGap", "BasePairMismatch", "AlignmentScore")] <- NaN

    # Pairwise alignment
    for (j in 1:nrow(temp_df)) {
      # Skip isolates that don't have Sanger sequences
      if (temp_df$IsolateSangerSequence[j] == "" | is.na(temp_df$IsolateSangerSequence[j])) next

      # Alignment
      result <- alignment_bp_count(temp_df$ESV[j], temp_df$IsolateSangerSequence[j], type = alignment_type)
      temp_df[j,c("AlignmentType","ConsensusLength", "BasePairGap", "BasePairMismatch", "AlignmentScore")] <-
        tibble(AlignmentType = alignment_type,
               ConsensusLength = result[1],
               BasePairGap = result[2],
               BasePairMismatch = result[3],
               AlignmentScore = result[4])

      #
      if ((j %% 10) == 0) cat(j, " ")
    }

    # Assignment
    sequences_alignment_list_func[[i]] <- temp_df
    cat("\n")
  }

  # Merge list
  sequences_alignment <- rbindlist(sequences_alignment_list_func[!is.na(sequences_alignment_list_func)])  %>% as_tibble

  return(sequences_alignment)
}

## Align Sanger to ESVs
alignment_type <- c("global", "local", "overlap", "global-local", "local-global")
sequences_alignment_list <- rep(list(NA), length(alignment_type))
for (i in 1:length(alignment_type)) {
  sequences_alignment_list[[i]] <- align_sequence(communities_abundance_syl_align, isolates_ID_RDP, alignment_type = alignment_type[i])
  print(alignment_type[i])
}
## Merge list
sequences_alignment_syl <- rbindlist(sequences_alignment_list) %>% as_tibble

fwrite(sequences_alignment_syl, here::here("data/temp/sequences_alignment_syl.csv"))
















