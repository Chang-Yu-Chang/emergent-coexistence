#' This script takes the outcome from Jean's script to identify pairs of isolates
#' with 0 mismatch and remove them from the list of pairs

library(tidyverse)
source(here::here("analysis/00-metadata.R"))

# 1. Run jean's python script to obtain the mismatch matrix ----
#' I modified the script so that it only takes the folder I created in the previous step
#' in commandline where 00-16s_mismatch.py is stored, execute `python 21-pairwise_16s_mismatch.py PATH_TO_DATA/`
#' the output matrix is saved as `PATH_TO_DATA/temp/21-mismatch_matrix_communities.csv`
#' where the column and row names are the Sanger ID


# 2. Match the mismatch matrix to my internal ID for isoaltes
mismatch_matrix <- read_csv(paste0(folder_data, "temp/21-mismatch_matrix_communities.csv"), show_col_types = F) %>%
    pivot_longer(cols = -`...1`) %>%
    rename(ID_row = `...1`, ID_col = name, Mismatch = value) %>%
    mutate(ID_row = as.character(ID_row)) %>%
    mutate(ID_col = as.character(ID_col))

isolates_ID <- read_csv(paste0(folder_data, "temp/00c-isolates_ID.csv"), show_col_types = F) %>%
    select(ID, Community, Isolate)

pairs_mismatch <- mismatch_matrix %>%
    filter(ID_row != ID_col) %>%
    # Find within-community pairs
    left_join(select(isolates_ID, ID_row = ID, Community1 = Community), by = "ID_row") %>%
    left_join(select(isolates_ID, ID_col = ID, Community2 = Community), by = "ID_col") %>%
    filter(Community1 == Community2) %>%
    select(Community = Community1, ID_row, ID_col, Mismatch) %>%
    arrange(Community) %>%
    # Reomve duplicate pairs
    left_join(select(isolates_ID, ID_row = ID, Isolate_row = Isolate), by = "ID_row") %>%
    left_join(select(isolates_ID, ID_col = ID, Isolate_col = Isolate), by = "ID_col") %>%
    filter(Isolate_row < Isolate_col) %>%
    rename_with(~str_replace(.,"_row", "1"), ends_with("_row")) %>%
    rename_with(~str_replace(.,"_col", "2"), ends_with("_col")) %>%
    select(Community, starts_with("Isolate"), starts_with("ID"), Mismatch) %>%
    mutate(Community = factor(Community, paste0("C", rep(1:12, each = 8), "R", 1:8))) %>%
    arrange(Community, Isolate1, Isolate2)

write_csv(pairs_mismatch, paste0(folder_data, "temp/22-pairs_mismatch.csv"))
