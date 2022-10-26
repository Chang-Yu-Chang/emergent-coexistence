#' This script reads community OTU table from the curated community files from Sylvie
#'
#' 1.

library(tidyverse)
source(here::here("analysis/00-metadata.R"))

communities <- read_csv(paste0(folder_data, "temp/00c-communities.csv"), col_types = cols())

# 1. Read the community ESV abundance ----
#' Note that this csv only has transfer 12
communities_abundance <- read_csv(paste0(folder_data, "raw/community_ESV/Emergent_Simplicity_Equilibrium_Data.csv"), col_types = cols()) %>%
    # Tidy up the variable names
    rename(CommunityESVID = ESV_ID, SampleID = Sample_ID, RelativeAbundance = Relative_Abundance, CarbonSource = Carbon_Source,
           ESVFamily = Family, ESVGenus = Genus) %>%
    # The inoculum, replicate, and sample ID in this df is zero-based, so one has to convert it into a one-based system
    mutate(Community = paste0("C", Inoculum+1, "R", Replicate+1)) %>%
    filter(Community %in% communities$Community) %>%
    mutate(Community = factor(Community, levels = communities$Community)) %>%
    filter(CarbonSource == "Glucose") %>%
    arrange(Community, CarbonSource, desc(RelativeAbundance)) %>%
    select(SampleID, Community, RelativeAbundance, CommunityESVID, ESV, ESVFamily, ESVGenus)

write_csv(communities_abundance, paste0(folder_data, "temp/13-communities_abundance.csv"))


# 2. Align isolate 16S sequences to community ESV ----
# ESV abundance in community
isolates_ID <- read_csv(paste0(folder_data, "temp/00c-isolates_ID.csv"), col_types = cols())
isolates_RDP <- read_csv(paste0(folder_data, "temp/12-isolates_RDP.csv"), col_types = cols())
isolates_ID_RDP <- isolates_ID %>% left_join(isolates_RDP) %>% select(ExpID, ID, Community, Isolate, Family, Genus, Sequence)

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
write_csv(sequences_alignment, paste0(folder_data, "temp/13-sequences_alignment.csv"))


# 3. Pick the best-matched ESV-isolate pairs and plot the relative abundances ----

# Find the match with highest alignment score; 68 rows
sequences_abundance_list <- rep(list(NA), 3)
allow_mismatch <- c(0:2, Inf)

for (i in 1:4) {
    sequences_abundance_list[[i]] <- sequences_alignment %>%
        filter(AlignmentType == "local") %>%
        # Filter for BasePairMatch
        filter(BasePairMismatch <= allow_mismatch[i]) %>%
        # For each Sanger, find the Sanger-ESV match with highest alignment score
        group_by(AlignmentType, ExpID) %>%
        arrange(desc(AlignmentScore)) %>%
        slice(1) %>%
        ungroup() %>%
        # For the scenario that two (or more) Sangers match to the same ESV, choose the match with the higher alignment score
        group_by(AlignmentType, Community, CommunityESVID) %>%
        arrange(desc(AlignmentScore)) %>%
        slice(1) %>%
        ungroup() %>%
        arrange(Community) %>%
        # Specify mismatch allowed
        mutate(AllowMismatch = allow_mismatch[i])
}

isolates_abundance <- bind_rows(sequences_abundance_list) %>%
    mutate(AlignmentType = factor(AlignmentType, levels = c("global", "local", "overlap", "global-local", "local-global"))) %>%
    arrange(AlignmentType, AllowMismatch, Community) %>%
    filter(AlignmentType == "local", AllowMismatch == "Inf") %>%
    # Select only necessary variables
    select(Community, ExpID, RelativeAbundance, CommunityESVID) %>%
    # Match it to isolate
    right_join(isolates_ID, by = c("Community", "ExpID")) %>%
    select(Community, Isolate, ExpID, RelativeAbundance) %>%
    group_by(Community) %>%
    mutate(RankRelativeAbundance = rank(-RelativeAbundance))


write_csv(isolates_abundance, paste0(folder_data, "temp/13-isolates_abundance.csv"))




















