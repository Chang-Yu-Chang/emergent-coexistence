#' This script reads community OTU table from the curated community files from Sylvie

library(tidyverse)
source(here::here("analysis/00-metadata.R"))

communities <- read_csv(paste0(folder_data, "temp/00c-communities.csv"), col_types = cols())

# 1. Read the community ESV abundance ----
communities_abundance <- read_csv(paste0(folder_data, "raw/community_ESV/Emergent_Comunity_Data.csv"), show_col_types = F) %>%
    filter(Transfer == 12) %>%
    mutate(Community = paste0("C", Inoculum, "R", Replicate)) %>%
    # Tidy up the variable names
    rename(CommunityESVID = ESV_ID, SampleID = Sample_ID, RelativeAbundance = Relative_Abundance,
           CarbonSource = Carbon_Source, ESVFamily = Family, ESVGenus = Genus) %>%
    filter(Community %in% communities$Community) %>%
    mutate(Community = factor(Community, levels = communities$Community)) %>%
    filter(CarbonSource == "Glucose") %>%
    arrange(Community, CarbonSource, desc(RelativeAbundance)) %>%
    select(SampleID, Community, RelativeAbundance, CommunityESVID, ESV, ESVFamily, ESVGenus)

write_csv(communities_abundance, paste0(folder_data, "temp/14-communities_abundance.csv"))

# 2. Align isolate 16S sequences to community ESV ----
# ESV abundance in community
isolates_RDP <- read_csv(paste0(folder_data, "temp/12-isolates_RDP.csv"), col_types = cols()) %>%
    select(ExpID, ID, Community, Isolate, Family, Genus, Sequence)

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
# seq1 <- isolates_RDP$Sequence[1]
# seq2 <- isolates_RDP$Sequence[3]
# alignment_bp_count(seq1, seq2) # R function from `sequence_pairwise_alignment.R`

## R function for align isolates' Sanger sequences to ESVs within communities, based on different alignment methods
align_sequence <- function(communities_abundance, isolates_RDP, alignment_type = "local") {
    # Community names. All communities with isolates
    comm_name <- unique(isolates_RDP$Community)

    # Empty list for alignment
    sequences_alignment_list_func <- rep(list(NA), length(comm_name)) %>% setNames(comm_name)

    # Loop through communities
    for (i in 1:length(comm_name)) {

        # Subset for community ESV
        temp_comm <- communities_abundance %>%
            filter(Community == comm_name[i])

        # Isolate sequences
        temp_iso <- isolates_RDP %>%
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
#alignment_type <- c("global", "local", "overlap", "global-local", "local-global")
alignment_type <- c("local")
sequences_alignment_list <- rep(list(NA), length(alignment_type))
for (i in 1:length(alignment_type)) {
    sequences_alignment_list[[i]] <- align_sequence(communities_abundance, isolates_RDP, alignment_type = alignment_type[i])
    print(alignment_type[i])
}
## Merge list
sequences_alignment <- bind_rows(sequences_alignment_list) %>% as_tibble
write_csv(sequences_alignment, paste0(folder_data, "temp/14-sequences_alignment.csv"))


# 3. Pick the best-matched ESV-isolate pairs with up to 0-2 mismatch ----
isolates_ID <- read_csv(paste0(folder_data, "temp/00c-isolates_ID.csv"), col_types = cols())
sequences_alignment <- read_csv(paste0(folder_data, "temp/14-sequences_alignment.csv"), show_col_types = F)

# Find the match with highest alignment score; 68 rows
sequences_abundance <- sequences_alignment %>%
    filter(AlignmentType == "local") %>%
    # Filter for BasePairMatch
    filter(BasePairMismatch <= "Inf") %>%
    # For each Sanger, find the Sanger-ESV match with highest alignment score
    group_by(AlignmentType, ExpID) %>%
    arrange(desc(AlignmentScore)) %>%
    slice(1) %>%
    ungroup() %>%
    # When two (or more) Sangers match the same ESV, choose the Sanger with the higher alignment score
    group_by(AlignmentType, Community, CommunityESVID) %>%
    arrange(desc(AlignmentScore)) %>%
    slice(1) %>%
    ungroup() %>%
    arrange(Community) %>%
    # Specify mismatch allowed
    mutate(AllowMismatch = allow_mismatch[i])

isolates_abundance <- sequences_abundance %>%
    #mutate(AlignmentType = factor(AlignmentType, levels = c("global", "local", "overlap", "global-local", "local-global"))) %>%
    arrange(AlignmentType, AllowMismatch, Community) %>%
    filter(AlignmentType == "local", AllowMismatch == "Inf") %>%
    select(Community, ExpID, everything()) %>%
    # Match it to isolate
    left_join(isolates_ID, by = c("Community", "ExpID", "ID", "Isolate")) %>%
    select(Community, Isolate, ExpID, Family, Genus, CommunityESVID, RelativeAbundance, everything()) %>%
    #arrange(Community, Isolate) %>%
    group_by(Community) %>%
    mutate(RankRelativeAbundance = rank(-RelativeAbundance)) %>%
    ungroup()

write_csv(isolates_abundance, paste0(folder_data, "temp/14-isolates_abundance.csv"))

# # 4. Keep all Sangers that match any ESV with up to 0-2 mismatch ----
# sequences_alignment <- read_csv(paste0(folder_data, "temp/14-sequences_alignment.csv"), show_col_types = F)
#
# # Find the match with highest alignment score; 68 rows
# sequences_abundance_list <- rep(list(NA), 3)
# allow_mismatch <- c(0:2, Inf)
#
# for (i in 1:4) {
#     sequences_abundance_list[[i]] <- sequences_alignment %>%
#         filter(AlignmentType == "local") %>%
#         # Filter for BasePairMatch
#         filter(BasePairMismatch <= allow_mismatch[i]) %>%
#         # For each Sanger, find the Sanger-ESV match with highest alignment score
#         group_by(AlignmentType, ExpID) %>%
#         arrange(desc(AlignmentScore)) %>%
#         slice(1) %>%
#         ungroup() %>%
#         arrange(Community) %>%
#         # Specify mismatch allowed
#         mutate(AllowMismatch = allow_mismatch[i])
# }
#
# isolates_abundance_loose <- bind_rows(sequences_abundance_list) %>%
#     mutate(AlignmentType = factor(AlignmentType, levels = c("global", "local", "overlap", "global-local", "local-global"))) %>%
#     arrange(AlignmentType, AllowMismatch, Community) %>%
#     filter(AlignmentType == "local", AllowMismatch == "Inf") %>%
#     # Select only necessary variables
#     select(Community, ExpID, RelativeAbundance, CommunityESVID) %>%
#     # Match it to isolate
#     right_join(isolates_ID, by = c("Community", "ExpID")) %>%
#     select(Community, Isolate, ExpID, RelativeAbundance)
#
# write_csv(isolates_abundance_loose, paste0(folder_data, "temp/14-isolates_abundance_loose.csv"))


# 5. Plot ESV-Sanger alignment ----
isolates <- read_csv(paste0(folder_data, "output/isolates.csv"), show_col_types = F)
pairs <- read_csv(paste0(folder_data, "output/pairs.csv"), show_col_types = F)
communities <- read_csv(paste0(folder_data, "temp/00c-communities.csv"), show_col_types = F)

sequences_alignment <- read_csv(paste0(folder_data, "temp/14-sequences_alignment.csv"), show_col_types = F)
isolates_abundance <- read_csv(paste0(folder_data, "temp/14-isolates_abundance.csv"), show_col_types = F) %>%
    left_join(communities) %>%
    mutate(Community = factor(Community, communities$Community))
isolate_ESV_1on1_match <- read_csv(paste0(folder_data, "temp/14-isolate_ESV_1on1_match.csv"), show_col_types = F) %>%
    left_join(communities) %>%
    mutate(Community = factor(Community, communities$Community))

p1 <- sequences_alignment %>%
    left_join(communities) %>%
    mutate(Community = factor(Community, communities$Community)) %>%
    filter(AlignmentType == "local") %>%
    drop_na() %>%
    ggplot() +
    geom_tile(aes(x = ExpID, y = CommunityESVID, fill = BasePairMismatch)) +
    # Each isolate has one ESV match
    geom_tile(data = isolate_ESV_1on1_match, aes(x = ExpID, y = CommunityESVID), fill = NA, color = grey(0.1), linewidth = .5) +
    # Best match
    geom_tile(data = isolates_abundance, aes(x = ExpID, y = CommunityESVID), fill = NA, color = "red", linewidth = .5) +
    facet_wrap(.~CommunityLabel, scales = "free", ncol = 4) +
    scale_fill_gradient(low = grey(0.9), high = "darkblue") +
    scale_x_discrete(expand = c(0,0)) +
    scale_y_discrete(expand = c(0,0)) +
    theme_classic() +
    theme(
        panel.border = element_rect(fill = NA, color = 1, linewidth = 0.5),
        axis.text.x = element_text(angle = 90, vjust = .5, hjust = 1, size = 7),
        axis.text.y = element_text(size = 6),
        strip.background = element_blank(),
        legend.position = c(.6, .05),
        plot.title = element_text(hjust = 0.5)
    ) +
    guides(fill = guide_legend(nrow = 1, title = "bp difference")) +
    labs(x = "", y = "ESV") +
    ggtitle("community")













