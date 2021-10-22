#' Match pairwise Sanger sequences
library(tidyverse)

# Read data
isolates_RDP <- read_csv(here::here("data/temp/isolates_RDP.csv")) %>% as_tibble()
#jean_isolates_RDP <- fread(here::here("data/temp/jean_isolates_RDP.csv")) %>% as_tibble()

# Align self-assembled communities isolates
## Call alignment function
source(here::here("misc/sequence_pairwise_alignment.R"))

## Create an empty dataframe
alignment_type <- c("global", "local", "overlap", "global-local", "local-global")

pairs_sanger_alignment <-
  t(combn(isolates_RDP$ExpID, 2)) %>%
  as.data.frame() %>%
  setNames(c("ExpID1", "ExpID2")) %>%
  mutate(ConsensusLength = NA, BasePairGap = NA, BasePairMismatch = NA, AlignmentScore = NA) %>%
  as_tibble()

pairs_sanger_alignment_list <- rep(list(pairs_sanger_alignment), length(alignment_type))

## Alignment
tt <- proc.time()
for (i in 1:length(alignment_type)) {
  for (j in 1:nrow(pairs_sanger_alignment)) {
    pairs_sanger_alignment_list[[i]][j, c("ConsensusLength", "BasePairGap", "BasePairMismatch", "AlignmentScore")] <-
      alignment_bp_count(
        seq1 = isolates_RDP$Sequence[match(unlist(pairs_sanger_alignment_list[[i]][j, "ExpID1"]), isolates_RDP$ExpID)],
        seq2 = isolates_RDP$Sequence[match(unlist(pairs_sanger_alignment_list[[i]][j, "ExpID2"]), isolates_RDP$ExpID)],
        type = alignment_type[i])
    if ((j %% 100) == 0)  cat("\n", j)

  }
  cat((proc.time() - tt)[1])
  cat("\n\n", alignment_type[i])
}

proc.time() - tt

pairs_sanger_alignment <-
  pairs_sanger_alignment_list %>%
  setNames(alignment_type) %>%
  rbindlist(idcol = "AlignmentType")




fwrite(pairs_sanger_alignment, root$find_file("data/temp/pairs_sanger_alignment.csv"))
fwrite(jean_pairs_sanger_alignment, root$find_file("data/temp/jean_pairs_sanger_alignment.csv"))




if (FALSE) {
    # Align random isolates from Jean
    jean_isolates_RDP

    jean_pairs_sanger_alignment <-
        t(combn(jean_isolates_RDP$ExpID, 2)) %>%
        as.data.frame() %>%
        setNames(c("ExpID1", "ExpID2")) %>%
        mutate(ConsensusLength = NA, BasePairGap = NA, BasePairMismatch = NA, AlignmentScore = NA) %>%
        as_tibble()

    jean_pairs_sanger_alignment_list <- rep(list(jean_pairs_sanger_alignment), length(alignment_type))

    ## Alignment
    tt <- proc.time()
    for (i in 1:length(alignment_type)) {
        for (j in 1:nrow(jean_pairs_sanger_alignment)) {
            jean_pairs_sanger_alignment_list[[i]][j, c("ConsensusLength", "BasePairGap", "BasePairMismatch", "AlignmentScore")] <-
                alignment_bp_count(
                    seq1 = jean_isolates_RDP$Sequence[match(unlist(jean_pairs_sanger_alignment_list[[i]][j, "ExpID1"]), jean_isolates_RDP$ExpID)],
                    seq2 = jean_isolates_RDP$Sequence[match(unlist(jean_pairs_sanger_alignment_list[[i]][j, "ExpID2"]), jean_isolates_RDP$ExpID)],
                    type = alignment_type[i])
            if ((j %% 100) == 0)  cat("\n", j)

        }
        cat((proc.time() - tt)[1])
        cat("\n\n", alignment_type[i])
    }

    proc.time() - tt

    jean_pairs_sanger_alignment <-
        jean_pairs_sanger_alignment_list %>%
        setNames(alignment_type) %>%
        rbindlist(idcol = "AlignmentType")

    #
}
