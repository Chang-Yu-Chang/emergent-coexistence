#' This script aligns the isolate Sanger and community ESV
#'
#' The alignment is performed within the 12 communities

library(tidyverse)
library(cowplot)
source(here::here("processing_scripts/00-metadata.R"))

factorize_communities <- function (x) x %>% mutate(Community = factor(Community, paste0("C", rep(1:12, each = 8), "R", rep(1:8, 12))))
communities <- read_csv(paste0(folder_data, "temp/00c-communities.csv"), col_types = cols())

# 1. Read community ESV sequence ----
communities_abundance <- read_csv(paste0(folder_data, "temp/14-communities_abundance.csv"), show_col_types = F) %>% factorize_communities %>%
    # Tidy up the variable names
    rename(CommunityESVID = ESV_ID, SampleID = Sample_ID, RelativeAbundance = Relative_Abundance, CarbonSource = Carbon_Source, ESVFamily = Family, ESVGenus = Genus)
communities_abundance_12comm <- communities_abundance %>% filter(Community %in% communities$Community)
communities_abundance_T12 <- communities_abundance_12comm %>% filter(Transfer == 12)

nrow(communities_abundance_T12) # 112 ESVs across the 13 communities
communities_abundance_T12 %>% filter(Community != "C10R2") %>% nrow() # 102 ESVs after removing C10R2

# Number of ESVs within communities
communities_abundance_T12 %>%
    group_by(Community) %>%
    count(name = "Count") %>%
    pull(Count) %>% range()

# 2. Read isolate Sanger sequences ----
isolates_RDP <- read_csv(paste0(folder_data, "temp/12-isolates_RDP.csv"), col_types = cols()) %>%
    select(ExpID, ID, Community, Isolate, Family, Genus, Sequence)

# 3. Alignment ----
# Align isolate 16S sequences to community ESV
count_alignment_bp <- function(seq1, seq2, type = "local"){
    #' R function for pairwise alignment and counting bp gap and mismatch
    #' Modified from Djordje's code

    # Make pairwise alignment
    algn <- Biostrings::pairwiseAlignment(seq1, seq2, type = type)

    # Table of two aligned sequences
    algn_t <- tibble(
        s1 = strsplit(as.character(algn@pattern), '')[[1]],
        s2 = strsplit(as.character(algn@subject), '')[[1]]
    ) %>%
        # Find difference; ignore N (which means any nucleotide)
        mutate(Diff = ((s1 != s2) & s1 != 'N' & s2 != 'N'))

    # Find gaps
    ngap <- sum(algn_t$s1 == "-" | algn_t$s2 == "-")

    # Find mismatch
    nmis <- sum((algn_t$s1 != algn_t$s2) & algn_t$s1 != "-" & algn_t$s2 != "-" & algn_t$s1 != "N" & algn_t$s2 != "N")

    # Length of consensus
    nlength <- nrow(algn_t)

    #
    return(c(
        ConsensusLength = nlength,
        BasePairGap = ngap,
        BasePairMismatch = nmis,
        AlignmentScore = algn@score
    ))
}

## Example code for match one pair of sequences
# seq1 <- isolates_RDP$Sequence[1]
# seq2 <- isolates_RDP$Sequence[3]
# count_alignment_bp(seq1, seq2) # R function from `sequence_pairwise_alignment.R`

sequences_alignment_list <- rep(list(NA), nrow(communities))
for (i in 1:nrow(communities)) {
    comm_ESV <- communities_abundance_T12 %>% filter(Community == communities$Community[i])
    comm_Sanger <- isolates_RDP %>% filter(Community == communities$Community[i])

    # Merge two dfs: it will auto create df with nrow(comm_seqs) = nrow(comm_ESV) * nrow(comm_Sanger)
    comm_seqs <- full_join(comm_ESV, comm_Sanger, by = join_by(Community), multiple = "all")
    cat("\n", communities$Community[i], "\n number of ESVs = ", nrow(comm_ESV), "\n number of isolates_RDP = ", nrow(comm_Sanger), "\n number of matches = ", nrow(comm_seqs), "\n")
    comm_seqs[,c("AlignmentType", "ConsensusLength", "BasePairGap", "BasePairMismatch", "AlignmentScore")] <- NA

    # Pairwise alignment
    for (j in 1:nrow(comm_seqs)) {
        # Skip isolates_RDP that don't have Sanger sequences
        if (comm_seqs$Sequence[j] == "" | is.na(comm_seqs$Sequence[j])) next

        # Alignment
        result <- count_alignment_bp(comm_seqs$ESV[j], comm_seqs$Sequence[j], type = "local")
        comm_seqs[j,c("AlignmentType","ConsensusLength", "BasePairGap", "BasePairMismatch", "AlignmentScore")] <-
            tibble(AlignmentType = "local",
                   ConsensusLength = result[1],
                   BasePairGap = result[2],
                   BasePairMismatch = result[3],
                   AlignmentScore = result[4])

        #
        if ((j %% 10) == 0) cat(j, " ")
    }

    sequences_alignment_list[[i]] <- comm_seqs

}
sequences_alignment <- bind_rows(sequences_alignment_list)
write_csv(sequences_alignment, paste0(folder_data, "temp/16-sequences_alignment.csv"))


# 4. For each isolate, save the matched ESV abundnance ----
sequences_alignment <- read_csv(paste0(folder_data, "temp/16-sequences_alignment.csv"), show_col_types = F)
# Apply filter and find match
algn_Sanger_ESV <- sequences_alignment %>%
    filter(ConsensusLength >= 200, BasePairMismatch <= 4) %>% # 133 alignments
    # For each Sanger, find the ESV which it has least base pair mismatch with
    group_by(ExpID) %>%
    filter(BasePairMismatch == min(BasePairMismatch)) %>%
    # When a Sanger has two ESV that both has zero, pick the first one
    slice(1) %>%
    ungroup()

nrow(algn_Sanger_ESV) # The number of data points should be 62
table(algn_Sanger_ESV$BasePairMismatch)
# 0  1  2  3  4
# 40 11  6  2  3

isolates_abundance <- algn_Sanger_ESV %>%
    select(Community, RelativeAbundance, CommunityESVID, ESV, ESVFamily, ESVGenus, ExpID, ID, Isolate, Family, Genus, Sequence, AlignmentType, ConsensusLength, BasePairGap, BasePairMismatch, AlignmentScore)
write_csv(isolates_abundance, paste0(folder_data, "temp/16-isolates_abundance.csv"))





# 5. Check if the number of alignments is correct ----
n_Sanger_comm <- isolates_RDP %>% drop_na(Sequence) %>% group_by(Community) %>% count(name = "n_Sanger") # Number of sanger sequences per community
n_ESV_comm <- communities_abundance_T12 %>% drop_na(ESV) %>% group_by(Community) %>% count(name = "n_ESV") # Number of ESVs per community
n_align_comm <- n_Sanger_comm %>%
    left_join(n_ESV_comm) %>%
    mutate(n_algn = n_Sanger * n_ESV)

sum(n_align_comm$n_algn) # expected number of alignments is 599
nrow(sequences_alignment) # actual number of alignments is 599

# Check the number of alignments pass the first criteria
sequences_alignment %>%
    filter(ConsensusLength >= 200, BasePairMismatch <= 4) %>%
    nrow # 129 alignments pass


# The isolates_RDP where its alignments are all dropped because of low-alignment quality
sequences_alignment %>%
    filter(ConsensusLength >= 200, BasePairMismatch <= 4) %>%
    distinct(Community, ExpID) %>%
    mutate(PassFirstFilter = T) %>%
    right_join(isolates_RDP) %>%
    filter(is.na(PassFirstFilter)) # three isolates haave low alignment quality
