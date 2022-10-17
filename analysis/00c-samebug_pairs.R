#' This script takes the outcome from Jean's script to identify pairs of isolates
#' with 0 mismatch and remove them from the list of pairs

library(tidyverse)

folder_main <- "~/Dropbox/lab/emergent-coexistence/plate_scan_pipeline/"
folder_script <- "~/Desktop/lab/emergent-coexistence/analysis/"

# 1. Create a folder such that Jean's script can work ----
## Sanger ID of isolates. There are 66 of them
ids <- read_csv("~/Dropbox/lab/emergent-coexistence/data/raw/pairwise_competition/isolates1.csv", show_col_types = F) %>%
    select(ExpID, ID, Community, Isolate) %>%
    pull(ID) %>%
    sort()
## To move ab1 of the isolates I used in the competition to one folder
## Go to ~/Dropbox/lab/emergent-coexistence/data/raw/Sanger and copy-paste the following commands to the commandline
for (i in 1:length(ids)) {
    if (ids[i] %in% 253:369) {
        cat("\ncp sanger_seq_16S_20180702/", ids[i], "-16S-rRNA-SeqF.ab1 sanger_seq_16S_communities/", ids[i], "-16S-rRNA-SeqF.ab1", sep = "")
        cat("\ncp sanger_seq_16S_20180702/", ids[i], "-16S-rRNA-SeqR.ab1 sanger_seq_16S_communities/", ids[i], "-16S-rRNA-SeqR.ab1", sep = "")
    }
    if (ids[i] %in% 370:462) {
        cat("\ncp sanger_seq_16S_20180726/", ids[i], "-16S-rRNA-SeqF.ab1 sanger_seq_16S_communities/", ids[i], "-16S-rRNA-SeqF.ab1", sep = "")
        cat("\ncp sanger_seq_16S_20180726/", ids[i], "-16S-rRNA-SeqR.ab1 sanger_seq_16S_communities/", ids[i], "-16S-rRNA-SeqR.ab1", sep = "")
    }
}

# 2. Run jean's python script to obtain the mismatch matrix ----
#' I modified the script so that it only takes the folder I created in the previous step
#' in commandline where 00-16s_mismatch.py is stored, execture `python 00-16s_mismatch.py`
#' the output matrix is saved as ~/Dropbox/lab/emergent-coexistence/data/raw/Sanger/Temp/mismatch_matrix_communities.csv
#' where the column and row names are the Sanger ID


# 3. Match the mismatch matrix to my internal ID for isoaltes
mismatch_matrix <- read_csv("~/Dropbox/lab/emergent-coexistence/data/raw/Sanger/Temp/mismatch_matrix_communities.csv", show_col_types = F) %>%
    pivot_longer(cols = -`...1`) %>%
    rename(ID_row = `...1`, ID_col = name, Mismatch = value) %>%
    mutate(ID_col = as.numeric(ID_col))

isolates_ID_match <- read_csv("~/Dropbox/lab/emergent-coexistence/data/raw/pairwise_competition/isolates1.csv", col_types = cols()) %>%
    mutate(Assembly = "self_assembly") %>%
    select(ID, Community, Isolate)

pairs_mismatch <- mismatch_matrix %>%
    filter(ID_row != ID_col) %>%
    # Find within-community pairs
    left_join(select(isolates_ID_match, ID_row = ID, Community1 = Community), by = "ID_row") %>%
    left_join(select(isolates_ID_match, ID_col = ID, Community2 = Community), by = "ID_col") %>%
    filter(Community1 == Community2) %>%
    select(Community = Community1, ID_row, ID_col, Mismatch) %>%
    arrange(Community) %>%
    # Reomve duplicate pairs
    left_join(select(isolates_ID_match, ID_row = ID, Isolate_row = Isolate), by = "ID_row") %>%
    left_join(select(isolates_ID_match, ID_col = ID, Isolate_col = Isolate), by = "ID_col") %>%
    filter(Isolate_row < Isolate_col) %>%
    rename_with(~str_replace(.,"_row", "1"), ends_with("_row")) %>%
    rename_with(~str_replace(.,"_col", "2"), ends_with("_col")) %>%
    select(Community, starts_with("Isolate"), starts_with("ID"), Mismatch) %>%
    mutate(Community = factor(Community, paste0("C", rep(1:12, each = 8), "R", 1:8))) %>%
    arrange(Community, Isolate1, Isolate2)

write_csv(pairs_mismatch, paste0(folder_main, "meta/00c-pairs_mismatch.csv"))

