#' Make DNA string sets for isolate 16S Sangers, multiple alignment, and make basic tree

library(tidyverse)
library(msa)
library(Biostrings)

# Make DNA string sets
isolates_RDP <- read_csv("~/Dropbox/lab/emergent-coexistence/data/temp/isolates_RDP.csv", col_types = cols())
isolates_seq <- DNAStringSet(isolates_RDP$Sequence)
names(isolates_seq) <- isolates_RDP$ExpID

# Multiple alignment on isolates across all communities. It takes a few seconds.
aln <- msa(isolates_seq, method = "ClustalW")

# Make basic tree
aln2 <- msa::msaConvert(aln, type = "seqinr::alignment")
d <- seqinr::dist.alignment(aln2, "identity")
tree <- ape::nj(d) # Neighbor-joining tree estimation

# Append trait data
isolates_ID_match <- read_csv("~/Dropbox/lab/emergent-coexistence/data/temp/isolates_ID_match.csv", col_types = cols())
tree2 <- as_tibble(tree) %>%
    mutate(ExpID = label) %>%
    full_join(isolates_RDP, by = "ExpID") %>%
    full_join(isolates_ID_match, by = "ExpID") %>%
    tidytree::as.treedata()

# Save
save(aln, aln2, tree, tree2, isolates_seq, file = "~/Dropbox/lab/emergent-coexistence/data/temp/isolates_sanger_seq.Rdata")




















