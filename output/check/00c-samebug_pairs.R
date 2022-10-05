#' This script takes the outcome from Jean's script to identify pairs of isolates
#' with 0 mismatch and remove them from the list of pairs

library(tidyverse)

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


# pairs_samebug <- pairs_mismatch %>% filter(Mismatch == 0)
# pairs_not_samebug <- pairs_mismatch %>% filter(Mismatch == 0)

write_csv(pairs_mismatch, "pairs_mismatch.csv")


#
#
# #
# isolates_ID_match <- read_csv("temp/isolates_ID_match.csv", show_col_types = F)
# pairs_ID <- read_csv("temp/pairs_ID.csv", show_col_types = F)

# Isolates from self-assembled communities



# # construct pairs
# pairs_ID <- rbindlist(lapply(split(isolates_ID_match, by='Community'),
#                              function(x) cbind(as.data.frame(t(combn(x$ID, 2))),
#                                                Community=x$Community[1])))
# pairs_ID [, s1 := ifelse(V1<V2, V1, V2)]
# pairs_ID [, s2 := ifelse(V1<V2, V2, V1)]
# pairs_ID  <- pairs_ID [, .(s1, s2, Community)]
# pairs_ID <- unique(pairs_ID)


# # identify identical isolates within community
# mismatches <- fread('~/Dropbox/lab/emergent-coexistence/data/raw/Sanger/mismatch_matrix_raw.csv', header=TRUE)
# mismatches <- melt(mismatches, id.vars='V1')
# mismatches[, variable := as.numeric(as.character(variable))]
# mismatches[, s1 := ifelse(V1<variable, V1, variable)]
# mismatches[, s2 := ifelse(V1<variable, variable, V1)]
# mismatches <- mismatches[, .(s1, s2, value)]
# mismatches <- mismatches[s1<s2]
# names(mismatches)[3] <- 'mismatch'
# mismatches <- unique(mismatches)
#
# pairs_ID <- merge(pairs_ID, mismatches, all.x=TRUE)
#
# # Manually identified duplicate isolates to remove from dataset
# trmi <- c(462, 355, 356, 461, 452, 446, 305, 435, 444, 348, 460, 454)
#
# isolates_ID_match <- isolates_ID_match[ID %nin% trmi]
# names(isolates_ID_match)[4] <- 'CYC_isolate'
#
# pairs_ID  <- pairs_ID[s1 %nin% trmi & s2 %nin% trmi]
#
#
# #
# fwrite(isolates_ID_match, "data/temp/isolates_ID_match.csv")
# fwrite(pairs_ID, "data/temp/pairs_ID.csv")
# # write_csv(communities, "~/Dropbox/lab/emergent-coexistence/data/output/communities.csv")
#




# Read the pair ID and isolate ID
pairs_ID_isolate <- pairs_ID %>%
    mutate(Community = factor(Community, paste0("C", rep(1:12, each = 8), "R", rep(1:8, 12)))) %>%
    rename(ID1 = s1, ID2 = s2) %>%
    left_join(select(isolates_ID_match, ID1 = ID, Isolate1 = CYC_isolate)) %>%
    left_join(select(isolates_ID_match, ID2 = ID, Isolate2 = CYC_isolate)) %>%
    arrange(Community, Isolate1, Isolate2)

# Switch the isolate's ID columns such that isolate1 always has smaller number than isolate2
temp_index <- pairs_ID_isolate$Isolate1 > pairs_ID_isolate$Isolate2

temp <- pairs_ID_isolate$ID1[temp_index]
pairs_ID_isolate$ID1[temp_index] <- pairs_ID_isolate$ID2[temp_index]
pairs_ID_isolate$ID2[temp_index] <- temp

temp <- pairs_ID_isolate$Isolate1[temp_index]
pairs_ID_isolate$Isolate1[temp_index] <- pairs_ID_isolate$Isolate2[temp_index]
pairs_ID_isolate$Isolate2[temp_index] <- temp

pairs_ID_isolate <- pairs_ID_isolate %>% arrange(Community, Isolate1, Isolate2)

## Dont Uncomment this trunk as the csv is later manually edited
# pairs_ID_isolate %>%
#     mutate(Pair = 1:n()) %>%
#     select(Pair, everything()) %>%
#     write_csv("output/check/temp/pairs_ID_isolate.csv")
#pairs_ID_isolate %>% filter(Community == "C11R5")


