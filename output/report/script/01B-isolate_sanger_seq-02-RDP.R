#' Read isolate 16S sequences and assign taxonomy using RDP

library(tidyverse)
library(data.table)
library(rRDPData)
library(rRDP)

isolates_ID_match <- fread(here::here("data/temp/", "isolates_ID_match.csv"))
isolates_16S <- fread(here::here("data/temp/", "isolates_16S.csv"))

# Djordje's isolates----
# Use the default classifier
isolate_seq_set <- Biostrings::DNAStringSet(isolates_16S$Sequence,use.names = T)
names(isolate_seq_set) <- isolates_16S$ID # Rename the sequence

# Use rdp for classification (this needs package rRDPData)
pred <- predict(rdp(), isolate_seq_set, confidence = 0) #  Predict 16s
conf_score <- attr(pred, "confidence") %>% as.data.frame()  # Confidence score
colnames(pred ) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus")
colnames(conf_score) <- paste0(c("Domain", "Phylum", "Class", "Order", "Family", "Genus"), "Score")
pred$ID <- rownames(pred) %>% as.numeric()
conf_score$ID <- rownames(pred) %>% as.numeric()

# Join the predicted taxonomy and isolates information
isolates_RDP <- left_join(isolates_16S, pred, "ID") %>% left_join(conf_score, "ID")

#' Fermenter or non-fermenter based on literatures
#' - Fermenter: Enterobacteriaceae, Aeromonadaceae
#' - Non-fermenter: Pseudomonadaceae, Moraxellaceae, Xanthomonadaceae, Alcaligenaceae, Comamonadaceae
isolates_RDP$Fermenter <- ifelse(isolates_RDP$Family %in% c("Pseudomonadaceae", "Moraxellaceae", "Xanthomonadaceae", "Alcaligenaceae", "Comamonadaceae"), F, ifelse(isolates_RDP$Family %in% c("Enterobacteriaceae", "Aeromonadaceae"), T, NA))

#' Gram-positive or gram-negative
#' It takes additional step to lyse the cell wall of Fram-positive strains in DNA extraction.
isolates_RDP$GramPositive <- ifelse(isolates_RDP$Family %in% c("Pseudomonadaceae", "Moraxellaceae", "Xanthomonadaceae", "Alcaligenaceae", "Comamonadaceae", "Enterobacteriaceae", "Aeromonadaceae"), F, ifelse(isolates_RDP$Family %in% c(""), T, NA))

# Remove contamination
isolates_RDP <- isolates_RDP %>%
    select(ID, Fermenter, GramPositive, Family, Genus, GenusScore, Sequence) %>%
    filter(!Genus == "Staphylococcus")

# Save 198 isolates with taxonomic information in a csv file `data/temp/isolates_RDP.csv`
fwrite(isolates_RDP, here::here("data/temp/isolates_RDP.csv"))



# Jean's isolates ----
ID_jean <- fread(here::here("data/raw/Sanger/jean_isolate_ids.csv")) %>%
    filter(!is.na(V2)) %>%
    select(V1, V2) %>%
    setNames(c("ExpID", "ID"))
# Use the default classifier
jean_isolates_16S <- ID_jean %>% left_join(isolates_16S) %>% filter(!is.na(Sequence))
isolate_seq_set <- Biostrings::DNAStringSet(jean_isolates_16S$Sequence)
names(isolate_seq_set) <- jean_isolates_16S$ID # Rename the sequence

# Use rdp for classification (this needs package rRDPData)
pred <- predict(rdp(), isolate_seq_set, confidence = 0) #  Predict 16s
conf_score <- attr(pred, "confidence") %>% as.data.frame()  # Confidence score
colnames(pred ) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus")
colnames(conf_score) <- paste0(c("Domain", "Phylum", "Class", "Order", "Family", "Genus"), "Score")
pred$ID <- rownames(pred) %>% as.numeric()
conf_score$ID <- rownames(pred) %>% as.numeric()

# Join the predicted taxonomy and isolates information
jean_isolates_RDP <- left_join(jean_isolates_16S, pred, "ID") %>% left_join(conf_score, "ID")
#' - Fermenter: Enterobacteriaceae, Aeromonadaceae
#' - Non-fermenter: Pseudomonadaceae, Moraxellaceae, Xanthomonadaceae, Alcaligenaceae, Comamonadaceae
jean_isolates_RDP$Fermenter <- ifelse(jean_isolates_RDP$Family %in% c("Pseudomonadaceae", "Moraxellaceae", "Xanthomonadaceae", "Alcaligenaceae", "Comamonadaceae"), F, ifelse(isolates_RDP$Family %in% c("Enterobacteriaceae", "Aeromonadaceae"), T, NA))

#' Gram-positive or gram-negative
#' It takes additional step to lyse the cell wall of Fram-positive strains in DNA extraction.
jean_isolates_RDP$GramPositive <- ifelse(jean_isolates_RDP$Family %in% c("Pseudomonadaceae", "Moraxellaceae", "Xanthomonadaceae", "Alcaligenaceae", "Comamonadaceae", "Enterobacteriaceae", "Aeromonadaceae"), F, ifelse(isolates_RDP$Family %in% c(""), T, NA))

# Remove contamination
jean_isolates_RDP <- jean_isolates_RDP %>%
    select(ExpID, ID, Fermenter, GramPositive, Family, Genus, GenusScore, Sequence) %>%
    filter(!Genus == "Staphylococcus")

fwrite(jean_isolates_RDP, here::here("data/temp/jean_isolates_RDP.csv"))





