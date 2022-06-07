#' Read isolate 16S sequences and assign taxonomy using RDP

library(tidyverse)
library(rRDPData)
library(rRDP)

isolates_ID_match <- read_csv("~/Dropbox/lab/emergent-coexistence/data/temp/isolates_ID_match.csv", col_types = cols())

# Read 16S sequence
isolates_16S <- read_csv("~/Dropbox/lab/emergent-coexistence/data/temp/isolates_16S.csv", col_types = cols()) %>%
    right_join(isolates_ID_match) %>%
    mutate(Isolate = as.character(Isolate)) %>%
    select(Assembly, ExpID, ID, Community, Isolate, Sequence)


# Two isolates using caseu sanger sequences: C10R2 isolate 3 (10.2.C.4) and C2R6 isolate 4 (2.6.A.5)
folder_directory <- "~/Dropbox/lab/emergent-coexistence/data/raw/Sanger/CASEU_six_plates/all_ab1/"
read_sequence <- function(Seq_ID) {
    paste0(folder_directory, Seq_ID, "-27F.ab1") %>%
        sangerseqR::read.abif() %>%
        sangerseqR::sangerseq() %>%
        sangerseqR::primarySeq(string = T)
}

temp <- tibble(Community = c("C2R6", "C10R2"), Isolate = c("4","3"), ExpID = c("2.6.A.5", "10.2.C.4"),
               Seq_ID = c("B2_T7_444_P2_28", "B2_T7_933_P2_95")) %>%
    rowwise() %>%
    mutate(Sequence = read_sequence(Seq_ID)) %>%
    select(-Seq_ID)

isolates_16S$Sequence[isolates_16S$ExpID == "2.6.A.5"] <- temp$Sequence[temp$ExpID == "2.6.A.5"]
isolates_16S$Sequence[isolates_16S$ExpID == "10.2.C.4"] <- temp$Sequence[temp$ExpID == "10.2.C.4"]


# Make DNA string set object
isolates_seq_set <- Biostrings::DNAStringSet(isolates_16S$Sequence,use.names = T)
names(isolates_seq_set) <- isolates_16S$ExpID

# Use rdp for classification (this needs package rRDPData)
pred <- predict(rdp(), isolates_seq_set, confidence = 0)
conf_score <- attr(pred, "confidence") %>% as.data.frame()
colnames(pred) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus")
colnames(conf_score) <- paste0(c("Domain", "Phylum", "Class", "Order", "Family", "Genus"), "Score")
pred$ExpID <- rownames(pred)
conf_score$ExpID <- rownames(pred)

# Join the predicted taxonomy and isolates information
isolates_RDP <- left_join(isolates_16S, pred, by = "ExpID") %>% left_join(conf_score, "ExpID")


#' Fermenter or non-fermenter based on literatures
#' - Fermenter: Enterobacteriaceae, Aeromonadaceae
#' - Non-fermenter: Pseudomonadaceae, Moraxellaceae, Xanthomonadaceae, Alcaligenaceae, Comamonadaceae
isolates_RDP$Fermenter <- ifelse(isolates_RDP$Family %in% c("Pseudomonadaceae", "Moraxellaceae", "Xanthomonadaceae", "Alcaligenaceae", "Comamonadaceae"), F, ifelse(isolates_RDP$Family %in% c("Enterobacteriaceae", "Aeromonadaceae", "Bacillaceae 1"), T, NA))

#' Gram-positive or gram-negative
#' It takes additional step to lyse the cell wall of Fram-positive strains in DNA extraction.
isolates_RDP$GramPositive <- ifelse(isolates_RDP$Family %in% c("Pseudomonadaceae", "Moraxellaceae", "Xanthomonadaceae", "Alcaligenaceae", "Comamonadaceae", "Enterobacteriaceae", "Aeromonadaceae", "Bacillaceae 1"), F, ifelse(isolates_RDP$Family %in% c(""), T, NA))

# Remove contamination
isolates_RDP <- isolates_RDP %>%
    select(ExpID, Fermenter, GramPositive, Family, Genus, GenusScore, Sequence) %>%
    filter(!Genus == "Staphylococcus")

# Remove unassgined families
isolates_RDP <- isolates_RDP %>%
    filter(!is.na(Fermenter), !is.na(GramPositive))

write_csv(isolates_RDP, "~/Dropbox/lab/emergent-coexistence/data/temp/isolates_RDP.csv")

