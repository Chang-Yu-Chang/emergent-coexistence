#' This script reads community OTU table from the curated ESV table files from Jean's ISME paper
#' and align the ESV to Sanger sequence

library(tidyverse)
library(cowplot)
source(here::here("analysis/00-metadata.R"))


# 0.1. Read the community ESV abundance ----
communities <- read_csv(paste0(folder_data, "temp/00c-communities.csv"), col_types = cols())
communities_abundance <- read_csv(paste0(folder_data, "raw/community_ESV/Emergent_Comunity_Data.csv"), show_col_types = F) %>%
    mutate(Community = paste0("C", Inoculum, "R", Replicate)) %>%
    # Tidy up the variable names
    rename(CommunityESVID = ESV_ID, SampleID = Sample_ID, RelativeAbundance = Relative_Abundance,
           CarbonSource = Carbon_Source, ESVFamily = Family, ESVGenus = Genus) %>%
    filter(Community %in% communities$Community) %>%
    mutate(Community = factor(Community, communities$Community)) %>%
    filter(CarbonSource == "Glucose") %>%
    arrange(Community, CarbonSource, desc(RelativeAbundance)) %>%
    select(SampleID, Community, Transfer, RelativeAbundance, CommunityESVID, ESV, ESVFamily, ESVGenus)

communities_abundance_T12 <- communities_abundance %>%
    filter(Transfer == 12)

write_csv(communities_abundance, paste0(folder_data, "temp/32-communities_abundance.csv"))
write_csv(communities_abundance_T12, paste0(folder_data, "temp/32-communities_abundance_T12.csv"))

# 0.2. Read isolate Sanger sequences ----
isolates_RDP <- read_csv(paste0(folder_data, "temp/12-isolates_RDP.csv"), col_types = cols()) %>%
    select(ExpID, ID, Community, Isolate, Family, Genus, Sequence)

# 1. Align isolate 16S sequences to community ESV ----
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
    cat("\n", communities$Community[i], "\n number of ESVs = ", nrow(comm_ESV), "\n number of isolates = ", nrow(comm_Sanger), "\n number of matches = ", nrow(comm_seqs), "\n")
    comm_seqs[,c("AlignmentType", "ConsensusLength", "BasePairGap", "BasePairMismatch", "AlignmentScore")] <- NA

    # Pairwise alignment
    for (j in 1:nrow(comm_seqs)) {
        # Skip isolates that don't have Sanger sequences
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
write_csv(sequences_alignment, paste0(folder_data, "temp/32-sequences_alignment.csv"))


# 2. Checking ----
communities <- read_csv(paste0(folder_data, "temp/00c-communities.csv"), show_col_types = F)
isolates <- read_csv(paste0(folder_data, "output/isolates.csv"), show_col_types = F)
sequences_alignment <- read_csv(paste0(folder_data, "temp/32-sequences_alignment.csv"), show_col_types = F)
communities_abundance_T12 <- read_csv(paste0(folder_data, "temp/32-communities_abundance_T12.csv"), show_col_types = F)

# 2.1 Check if the number of alignment is correct ----
n_Sanger_comm <- isolates %>% drop_na(Sequence) %>% group_by(Community) %>% count(name = "n_Sanger") # Number of sanger sequences per community
n_ESV_comm <- communities_abundance_T12 %>% drop_na(ESV) %>% group_by(Community) %>% count(name = "n_ESV") # Number of ESVs per community
n_align_comm <- n_Sanger_comm %>%
    left_join(n_ESV_comm) %>%
    mutate(n_algn = n_Sanger * n_ESV)

sum(n_align_comm$n_algn) # expected number of alignments
nrow(sequences_alignment) # actual number of alignments

# Check the bp distribution
sequences_alignment %>%
    ggplot() +
    geom_histogram(aes(x = BasePairMismatch), binwidth = 1, color = "black", fill = "white") +
    scale_x_continuous(breaks = seq(0, 50, 5)) +
    theme_classic() +
    theme() +
    guides() +
    labs()

# Check the consensus length distribution
all(str_count(sequences_alignment$ESV) == 233) # All ESVs have 233 bp
table(sequences_alignment$ConsensusLength) # Some alignment has very short consensus
sort(unique(sequences_alignment$ConsensusLength))
sequences_alignment %>%
    ggplot() +
    geom_histogram(aes(x = ConsensusLength), binwidth = 1, color = "black", fill = "white") +
    #geom_point(aes(x = ConsensusLength, y = str_count(ESV)), size = 2, shape = 21) +
    #scale_x_continuous(breaks = seq(0, 50, 5)) +
    theme_classic() +
    theme() +
    guides() +
    labs()

# Check the alignment score distribution
sequences_alignment %>%
    ggplot() +
    geom_histogram(aes(x = AlignmentScore), binwidth = 5, color = "black", fill = "white") +
    theme_classic() +
    theme() +
    guides() +
    labs()


# 2.2 For each isolate Sanger, find the ESVs that it has the least base pair mismatch ----
algn_Sanger_ESV <- sequences_alignment %>%
    filter(ConsensusLength > 160) %>%
    # For each Sanger, find the ESV which it has least base pair mismatch with
    group_by(ExpID) %>%
    filter(BasePairMismatch == min(BasePairMismatch)) %>%
    # When a Sanger has two ESV that both has zero, pick the first one
    slice(1) %>%
    ungroup()

nrow(algn_Sanger_ESV) # The number of data points should be 68 because of 68 Sanger sequences

algn_Sanger_ESV %>% # Number of mismatches for each isolate Sanger
    group_by(BasePairMismatch) %>%
    count()
# algn_Sanger_ESV %>%
#     filter(BasePairMismatch > 4) %>%
#     select(Community, ExpID, Genus, CommunityESVID, BasePairMismatch)


# Heatmap, ESV and sanger alignment
isolates_algn_bp_mismatch <- isolates %>%
    select(Community, ExpID) %>%
    left_join(algn_Sanger_ESV) %>%
    mutate(Community = factor(Community, communities$Community)) %>%
    mutate(ExpIDGenus = paste0(ExpID, "-", Genus))

p <- sequences_alignment %>%
    mutate(Community = factor(Community, communities$Community)) %>%
    mutate(ExpIDGenus = paste0(ExpID, "-", Genus)) %>%
    ggplot() +
    geom_tile(aes(x = ExpIDGenus, y = CommunityESVID, fill = BasePairMismatch)) +
    # Each Sanger has one ESV match
    geom_tile(data = isolates_algn_bp_mismatch, aes(x = ExpIDGenus, y = CommunityESVID), fill = NA, color = "gold", linewidth = 1) +
    facet_wrap(.~Community, scales = "free", ncol = 4) +
    scale_fill_gradient(low = RColorBrewer::brewer.pal(n = 9, "Blues")[9], high = grey(0.9)) +
    scale_x_discrete(expand = c(0,0)) +
    scale_y_discrete(expand = c(0,0)) +
    theme_classic() +
    theme(
        panel.border = element_rect(fill = NA, color = 1, linewidth = 0.5),
        axis.text.x = element_text(angle = 30, vjust = 1.1, hjust = 1, size = 7),
        axis.text.y = element_text(size = 7),
        legend.position = c(.6, .05)
    ) +
    guides(fill = guide_legend(nrow = 1)) +
    labs()

ggsave(paste0(folder_data, "temp/32-01-ESV_Sanger_bp_mismatch.png"), p, width = 12, height = 12)


# Heatmap, ESV and sanger alignment, binned by mismatch
p <- sequences_alignment %>%
    mutate(Community = factor(Community, communities$Community)) %>%
    mutate(ExpIDGenus = paste0(ExpID, "-", Genus)) %>%
    mutate(BasePairMismatchCategory = case_when(
        BasePairMismatch == 0 ~ "0",
        BasePairMismatch == 1 ~ "1",
        BasePairMismatch == 2 ~ "2",
        BasePairMismatch == 3 ~ "3",
        BasePairMismatch == 4 ~ "4",
        # BasePairMismatch >= 4 & BasePairMismatch <= 10 ~ "4-10",
        # BasePairMismatch >= 11 & BasePairMismatch <= 29 ~ "11-29",
        # BasePairMismatch >= 30 ~ ">=30"
        BasePairMismatch >= 5 ~ ">=5"
    #) %>% ordered(c(0:3, "4-10", "11-29", ">=30"))
    ) %>% ordered(c(0:4, ">=5"))
    ) %>%
    ggplot() +
    geom_tile(aes(x = ExpIDGenus, y = CommunityESVID, fill = BasePairMismatchCategory)) +
    # Each Sanger has one ESV match
    geom_tile(data = isolates_algn_bp_mismatch, aes(x = ExpIDGenus, y = CommunityESVID), fill = NA, color = "gold", linewidth = 1) +
    facet_wrap(.~Community, scales = "free", ncol = 4) +
    scale_fill_manual(values = rev(RColorBrewer::brewer.pal(n = 6, "Blues"))) +
    scale_x_discrete(expand = c(0,0)) +
    scale_y_discrete(expand = c(0,0)) +
    theme_classic() +
    theme(
        panel.border = element_rect(fill = NA, color = 1, linewidth = 0.5),
        axis.text.x = element_text(angle = 30, vjust = 1.1, hjust = 1, size = 7),
        axis.text.y = element_text(size = 7),
        legend.position = c(.6, .05)
    ) +
    guides(fill = guide_legend(nrow = 1)) +
    labs()

ggsave(paste0(folder_data, "temp/32-02-ESV_Sanger_bp_mismatch_binned.png"), p, width = 12, height = 12)

#
p <- sequences_alignment %>%
    mutate(Community = factor(Community, communities$Community)) %>%
    mutate(ExpIDGenus = paste0(ExpID, "-", Genus)) %>%
    ggplot() +
    geom_tile(aes(x = ExpIDGenus, y = CommunityESVID, fill = ConsensusLength)) +
    # Each Sanger has one ESV match
    geom_tile(data = isolates_algn_bp_mismatch, aes(x = ExpIDGenus, y = CommunityESVID), fill = NA, color = "gold", linewidth = 1) +
    facet_wrap(.~Community, scales = "free", ncol = 4) +
    scale_fill_gradient(low = RColorBrewer::brewer.pal(n = 9, "Blues")[9], high = grey(0.9)) +
    scale_x_discrete(expand = c(0,0)) +
    scale_y_discrete(expand = c(0,0)) +
    theme_classic() +
    theme(
        panel.border = element_rect(fill = NA, color = 1, linewidth = 0.5),
        axis.text.x = element_text(angle = 30, vjust = 1.1, hjust = 1, size = 7),
        axis.text.y = element_text(size = 7),
        legend.position = c(.6, .05)
    ) +
    guides(fill = guide_legend(nrow = 1)) +
    labs()

ggsave(paste0(folder_data, "temp/32-03-ESV_Sanger_bp_concensus_length.png"), p, width = 12, height = 12)


# 2.3 For each isolate Sanger, find the ESVs that it has the highest alignment score ----
algn_Sanger_ESV <- sequences_alignment %>%
    # For each Sanger, find the ESV which it has least base pair mismatch with
    group_by(ExpID) %>%
    filter(AlignmentScore == max(AlignmentScore)) %>%
    # When a Sanger has two ESV that both has zero, pick the first one
    slice(1) %>%
    ungroup()

nrow(algn_Sanger_ESV) # The number of data points should be 66 because of 66 Sanger sequences

# Heatmap, ESV and sanger alignment
isolates_algn_score <- isolates %>%
    select(Community, ExpID) %>%
    left_join(algn_Sanger_ESV) %>%
    mutate(Community = factor(Community, communities$Community)) %>%
    mutate(ExpIDGenus = paste0(ExpID, "-", Genus))

p <- sequences_alignment %>%
    mutate(Community = factor(Community, communities$Community)) %>%
    mutate(ExpIDGenus = paste0(ExpID, "-", Genus)) %>%
    ggplot() +
    geom_tile(aes(x = ExpIDGenus, y = CommunityESVID, fill = AlignmentScore)) +
    # Each Sanger has one ESV match
    geom_tile(data = isolates_algn_score, aes(x = ExpIDGenus, y = CommunityESVID), fill = NA, color = "gold", linewidth = 1) +
    facet_wrap(.~Community, scales = "free", ncol = 4) +
    scale_fill_gradient(low = RColorBrewer::brewer.pal(n = 9, "Blues")[9], high = grey(0.9)) +
    scale_x_discrete(expand = c(0,0)) +
    scale_y_discrete(expand = c(0,0)) +
    theme_classic() +
    theme(
        panel.border = element_rect(fill = NA, color = 1, linewidth = 0.5),
        axis.text.x = element_text(angle = 30, vjust = 1.1, hjust = 1, size = 7),
        axis.text.y = element_text(size = 7),
        legend.position = c(.6, .05)
    ) +
    guides(fill = guide_legend(nrow = 1)) +
    labs()

ggsave(paste0(folder_data, "temp/32-04-ESV_Sanger_score.png"), p, width = 12, height = 12)

# 2.4 For each ESV, find the Sanger with least bp mismatch ----
algn_Sanger_ESV1 <- algn_Sanger_ESV %>%
    group_by(Community, CommunityESVID) %>%
    arrange(BasePairMismatch) %>%
    slice(1)

isolates_algn_bp_mismatch1 <- isolates %>%
    left_join(algn_Sanger_ESV1) %>%
    mutate(Community = factor(Community, communities$Community)) %>%
    mutate(ExpIDGenus = paste0(ExpID, "-", Genus))

# Heatmap, ESV and sanger alignment, binned by mismatch
p <- sequences_alignment %>%
    mutate(Community = factor(Community, communities$Community)) %>%
    mutate(ExpIDGenus = paste0(ExpID, "-", Genus)) %>%
    mutate(BasePairMismatchCategory = case_when(
        BasePairMismatch == 0 ~ "0",
        BasePairMismatch == 1 ~ "1",
        BasePairMismatch == 2 ~ "2",
        BasePairMismatch == 3 ~ "3",
        BasePairMismatch == 4 ~ "4",
        # BasePairMismatch >= 4 & BasePairMismatch <= 10 ~ "4-10",
        # BasePairMismatch >= 11 & BasePairMismatch <= 29 ~ "11-29",
        # BasePairMismatch >= 30 ~ ">=30"
        BasePairMismatch >= 5 ~ ">=5"
        #) %>% ordered(c(0:3, "4-10", "11-29", ">=30"))
    ) %>% ordered(c(0:4, ">=5"))
    ) %>%
    ggplot() +
    geom_tile(aes(x = ExpIDGenus, y = CommunityESVID, fill = BasePairMismatchCategory)) +
    # Each Sanger has one ESV match
    geom_tile(data = isolates_algn_bp_mismatch, aes(x = ExpIDGenus, y = CommunityESVID), fill = NA, color = "gold", linewidth = 1) +
    # Each ESV has one Sanger
    geom_tile(data = isolates_algn_bp_mismatch1, aes(x = ExpIDGenus, y = CommunityESVID), fill = NA, color = "red", linewidth = 1) +
    facet_wrap(.~Community, scales = "free", ncol = 4) +
    scale_fill_manual(values = rev(RColorBrewer::brewer.pal(n = 6, "Blues"))) +
    scale_x_discrete(expand = c(0,0)) +
    scale_y_discrete(expand = c(0,0)) +
    theme_classic() +
    theme(
        panel.border = element_rect(fill = NA, color = 1, linewidth = 0.5),
        axis.text.x = element_text(angle = 30, vjust = 1.1, hjust = 1, size = 7),
        axis.text.y = element_text(size = 7),
        legend.position = c(.6, .05)
    ) +
    guides(fill = guide_legend(nrow = 1)) +
    labs()

ggsave(paste0(folder_data, "temp/32-05-ESV_Sanger_bp_mismatch_binned.png"), p, width = 12, height = 12)

# 2.5 plot the isolate abundance ----
# Check if the number is correct. Facet by community
algn_Sanger_ESV1 %>%
    mutate(Community = factor(Community, communities$Community)) %>%
    group_by(Community) %>% count()

p <- algn_Sanger_ESV1 %>%
    mutate(Community = factor(Community, communities$Community)) %>%
    ggplot() +
    geom_col(aes(x = CommunityESVID, y = RelativeAbundance, fill = ESVGenus)) +
    scale_fill_manual(values = c(RColorBrewer::brewer.pal(8, "Set1"), RColorBrewer::brewer.pal(8, "Set2"))) +
    facet_wrap(~Community, ncol = 3, scales = "free_x") +
    theme_classic() +
    theme(
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)
    ) +
    guides() +
    labs()

ggsave(paste0(folder_data, "temp/32-06-matched_ESV_abundance_comm.png"), p, width = 12, height = 12)

# Total abundance by community
p1 <- algn_Sanger_ESV1 %>%
    mutate(Community = factor(Community, communities$Community)) %>%
    ggplot() +
    geom_col(aes(x = Community, y = RelativeAbundance, fill = ESVGenus, group = CommunityESVID), color = "black") +
    scale_fill_manual(values = c(RColorBrewer::brewer.pal(8, "Set1"), RColorBrewer::brewer.pal(8, "Set2"))) +
    #facet_wrap(~Community, ncol = 3, scales = "free_x") +
    theme_classic() +
    theme(
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        panel.border = element_rect(color = 1, fill = NA, linewidth = 0.5)
    ) +
    guides() +
    labs() +
    ggtitle("ESV genus")

p2 <- algn_Sanger_ESV1 %>%
    mutate(Community = factor(Community, communities$Community)) %>%
    ggplot() +
    geom_col(aes(x = Community, y = RelativeAbundance, fill = Genus, group = CommunityESVID), color = "black") +
    scale_fill_manual(values = c(RColorBrewer::brewer.pal(8, "Set1"), RColorBrewer::brewer.pal(8, "Set2"))) +
    #facet_wrap(~Community, ncol = 3, scales = "free_x") +
    theme_classic() +
    theme(
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        panel.border = element_rect(color = 1, fill = NA, linewidth = 0.5)
    ) +
    guides() +
    labs() +
    ggtitle("Sanger genus")
p <- plot_grid(p1, p2, ncol = 1, align = "v", axis = "lr")
ggsave(paste0(folder_data, "temp/32-07-matched_ESV_abundance_color.png"), p, width = 6, height = 8)


# 2.6 Plot by family ----
p <- algn_Sanger_ESV1 %>%
    mutate(Community = factor(Community, communities$Community)) %>%
    group_by(Community) %>%
    arrange(Community, desc(RelativeAbundance)) %>%
    mutate(CommunityESVID = factor(CommunityESVID)) %>%
    #mutate(Family = factor(CommunityESVID)) %>%
    ggplot() +
    geom_col(aes(x = Community, y = RelativeAbundance, fill = Family, group = CommunityESVID), color = "black") +
    scale_fill_manual(values = c(RColorBrewer::brewer.pal(8, "Set1"), RColorBrewer::brewer.pal(8, "Set2"))) +
    #facet_wrap(~Community, ncol = 3, scales = "free_x") +
    theme_classic() +
    theme(
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        panel.border = element_rect(color = 1, fill = NA, linewidth = 0.5)
    ) +
    guides() +
    labs()

ggsave(paste0(folder_data, "temp/32-08-matched_ESV_abundance_family.png"), p, width = 6, height = 4)

# 2.7 Check if the total number of ESV adds up ----
n_ESV <- communities_abundance_T12 %>%
    mutate(Community = factor(Community, communities$Community)) %>%
    group_by(Community) %>% count(name = "n_ESV")
n_ESV_matched <- algn_Sanger_ESV1 %>%
    mutate(Community = factor(Community, communities$Community)) %>%
    group_by(Community) %>% count(name = "n_ESV_matched")
n_ESV %>%
    left_join(n_ESV_matched) %>%
    ungroup() %>%
    summarise(n_ESV = sum(n_ESV), n_ESV_matched = sum(n_ESV_matched)) # total 112 ESVs, 47 matched to Sanger

algn_Sanger_ESV1 %>%
    group_by(Community) %>%
    summarize(TotalAbundance = sum(RelativeAbundance)) %>%
    summarize(MeanTotalAbundance = mean(TotalAbundance)) # Matched ESVs represent 89.5% of abundance


algn_ESV <- algn_Sanger_ESV %>%
    # When two (or more) Sangers match the same ESV, choose the Sanger with the lowest mismatch
    group_by(Community, CommunityESVID) %>%
    arrange(BasePairMismatch) %>%
    slice(1) %>%
    ungroup() %>%
    arrange(Community)

# Check the alignment of the two strains
#sequences_alignment %>%
algn_Sanger_ESV %>%
    filter(ExpID %in% c("10.2.C.4", "2.6.A.5"))



# Store the isolate abundance data ----

# Store csv where each Sanger has one ESV
isolates_abundance_all_sanger <- algn_Sanger_ESV
write_csv(isolates_abundance_all_sanger, paste0(folder_data, "temp/32-isolates_abundance_all_sanger.csv"))

# Store csv where one Sanger match one ESV
isolates_abundance <- algn_Sanger_ESV1
write_csv(isolates_abundance, paste0(folder_data, "temp/32-isolates_abundance.csv"))


# isolates_abundance <- sequences_abundance %>%
#     arrange(AlignmentType, AllowMismatch, Community) %>%
#     filter(AlignmentType == "local", AllowMismatch == "Inf") %>%
#     select(Community, ExpID, everything()) %>%
#     # Match it to isolate
#     left_join(isolates_ID, by = c("Community", "ExpID", "ID", "Isolate")) %>%
#     select(Community, Isolate, ExpID, Family, Genus, CommunityESVID, RelativeAbundance, everything()) %>%
#     #arrange(Community, Isolate) %>%
#     group_by(Community) %>%
#     mutate(RankRelativeAbundance = rank(-RelativeAbundance)) %>%
#     ungroup()
#write_csv(isolates_abundance, paste0(folder_data, "temp/32-isolates_abundance.csv"))













