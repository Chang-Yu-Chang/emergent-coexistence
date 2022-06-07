#' Match isolates' 16S sequences to community ESV
library(tidyverse)

# ESV abundance in community
communities_abundance <- read_csv("~/Dropbox/lab/emergent-coexistence/data/temp/communities_abundance.csv", col_types = cols()) # Data curated by Sylvie
isolates_ID_match <- read_csv("~/Dropbox/lab/emergent-coexistence/data/temp/isolates_ID_match.csv", col_types = cols())
isolates_RDP <- read_csv("~/Dropbox/lab/emergent-coexistence/data/temp/isolates_RDP.csv", col_types = cols())
isolates_ID_RDP <- isolates_ID_match %>% left_join(isolates_RDP) %>% select(ExpID, ID, Community, Isolate, Family, Genus, Sequence)

# Match isolate Sanger seq to community ESV ----
## R function for pairwise alignment and counting bp gap and mismatch
alignment_bp_count <- function(seq1, seq2, type = "local"){
    # Make pairwise alignment
    algn <- Biostrings::pairwiseAlignment(seq1, seq2, type = type)

    # Make a df
    temp <- tibble(
        s1 = strsplit(as.character(Biostrings::pattern(algn)), '')[[1]],
        s2 = strsplit(as.character(Biostrings::subject(algn)), '')[[1]]
    )

    # Find difference; ignore N (which means any nucleotide)
    temp <- temp %>% mutate(Diff = ((s1 != s2) & s1 != 'N' & s2 != 'N'))

    # Find gaps
    ngap <- sum(temp$s1 == "-" | temp$s2 == "-")

    # Find mismatch
    nmis <- sum((temp$s1 != temp$s2) & temp$s1 != "-" & temp$s2 != "-" & temp$s1 != "N" & temp$s2 != "N")

    # Length of consensus
    nlength <- nrow(temp)

    #
    c(ConsensusLength = nlength, BasePairGap = ngap, BasePairMismatch = nmis, AlignmentScore = algn@score) %>%
        return()
}
## Example code for match one pair of sequences
# seq1 <- isolates_ID_RDP$Sequence[1]
# seq2 <- isolates_ID_RDP$Sequence[3]
# alignment_bp_count(seq1, seq2) # R function from `sequence_pairwise_alignment.R`

## R function for align isolates' Sanger sequences to ESVs within communities, based on different alignment methods
align_sequence <- function(communities_abundance, isolates_ID_RDP, alignment_type = "global") {
    # Community names. All communities with isolates
    comm_name <- unique(isolates_ID_RDP$Community)

    # Empty list for alignment
    sequences_alignment_list_func <- rep(list(NA), length(comm_name)) %>% setNames(comm_name)

    # Loop through communities
    for (i in 1:length(comm_name)) {

        # Subset for community ESV
        temp_comm <- communities_abundance %>%
            filter(Community == comm_name[i])

        # Isolate sequences
        temp_iso <- isolates_ID_RDP %>%
            filter(Community == comm_name[i])

        # Merge two dfs: it will auto create df with nrow(temp) = nrow(temp_comm) * nrow(temp_iso)
        temp <- temp_comm %>% left_join(temp_iso, by = "Community")
        cat("\n", comm_name[i], "\n number of ESVs = ", nrow(temp_comm), "\n number of isolates = ", nrow(temp_iso), "\n pairs = ", nrow(temp), "\n")

        # Add five new columns to temp
        temp$AlignmentType <- "NA"
        temp[,c("ConsensusLength", "BasePairGap", "BasePairMismatch", "AlignmentScore")] <- NaN

        # Pairwise alignment
        for (j in 1:nrow(temp)) {
            # Skip isolates that don't have Sanger sequences
            if (temp$Sequence[j] == "" | is.na(temp$Sequence[j])) next

            # Alignment
            result <- alignment_bp_count(temp$ESV[j], temp$Sequence[j], type = alignment_type)
            temp[j,c("AlignmentType","ConsensusLength", "BasePairGap", "BasePairMismatch", "AlignmentScore")] <-
                tibble(AlignmentType = alignment_type,
                       ConsensusLength = result[1],
                       BasePairGap = result[2],
                       BasePairMismatch = result[3],
                       AlignmentScore = result[4])

            #
            if ((j %% 10) == 0) cat(j, " ")
        }

        # Assignment
        sequences_alignment_list_func[[i]] <- temp
        cat("\n")
    }

    # Merge list
    sequences_alignment <- bind_rows(sequences_alignment_list_func[!is.na(sequences_alignment_list_func)])  %>% as_tibble

    return(sequences_alignment)
}

## Align Sanger to ESVs
alignment_type <- c("global", "local", "overlap", "global-local", "local-global")
sequences_alignment_list <- rep(list(NA), length(alignment_type))
for (i in 1:length(alignment_type)) {
    sequences_alignment_list[[i]] <- align_sequence(communities_abundance, isolates_ID_RDP, alignment_type = alignment_type[i])
    print(alignment_type[i])
}
## Merge list
sequences_alignment <- bind_rows(sequences_alignment_list) %>% as_tibble
write_csv(sequences_alignment, "~/Dropbox/lab/emergent-coexistence/data/temp/sequences_alignment.csv")
















