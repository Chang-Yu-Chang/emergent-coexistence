#' This scripts reads isolate 16S sequences and assigns taxonomy using RDP

if (FALSE) {
    BiocManager::install("rRDPData")
    BiocManager::install("rRDP")
}

library(tidyverse)
library(rRDPData)
library(rRDP)
library(Biostrings)
source(here::here("processing_scripts/00-metadata.R"))

isolates_ID <- read_csv(paste0(folder_data, "temp/00c-isolates_ID.csv"), col_types = cols())

# Read 16S sequence
isolates_16S <- read_csv(paste0(folder_data, "temp/11-isolates_16S.csv"), col_types = cols()) %>%
    right_join(isolates_ID) %>%
    mutate(Isolate = as.character(Isolate)) %>%
    select(ExpID, ID, Community, Isolate, Sequence) %>%
    filter(!is.na(Sequence))

# Make DNA string set object
isolates_seq_set <- DNAStringSet(isolates_16S$Sequence, use.names = T)
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
#' It takes additional step to lyse the cell wall of gram-positive strains in DNA extraction.
isolates_RDP$GramPositive <- ifelse(isolates_RDP$Family %in% c("Pseudomonadaceae", "Moraxellaceae", "Xanthomonadaceae", "Alcaligenaceae", "Comamonadaceae", "Enterobacteriaceae", "Aeromonadaceae", "Bacillaceae 1"), F, ifelse(isolates_RDP$Family %in% c(""), T, NA))

write_csv(isolates_RDP, paste0(folder_data, "temp/12-isolates_RDP.csv"))













