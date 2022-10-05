#' Generate isolate and community metadata files
#' Two sets of isolates used: self-assembled community isolates and Jean's natural isolates
library(tidyverse)
library(data.table)

# Isolates from self-assembled communities
isolates_ID_match <- read_csv("~/Dropbox/lab/emergent-coexistence/data/raw/pairwise_competition/isolates1.csv", col_types = cols()) %>%
#isolates_ID_match <- read_csv("~/Dropbox/lab/emergent-coexistence/data/raw/pairwise_competition/isolates1_cp.csv", col_types = cols()) %>%
  mutate(Assembly = "self_assembly") %>%
  select(ExpID, ID, Community, Isolate) %>% as.data.table()


# construct pairs
pairs_ID <- rbindlist(lapply(split(isolates_ID_match, by='Community'),
                             function(x) cbind(as.data.frame(t(combn(x$ID, 2))),
                                               Community=x$Community[1])))
pairs_ID [, s1 := ifelse(V1<V2, V1, V2)]
pairs_ID [, s2 := ifelse(V1<V2, V2, V1)]
pairs_ID  <- pairs_ID [, .(s1, s2, Community)]
pairs_ID <- unique(pairs_ID)


# identify identical isolates within community

mismatches <- fread('~/Dropbox/lab/emergent-coexistence/data/raw/Sanger/mismatch_matrix_raw.csv', header=TRUE)
mismatches <- melt(mismatches, id.vars='V1')
mismatches[, variable := as.numeric(as.character(variable))]
mismatches[, s1 := ifelse(V1<variable, V1, variable)]
mismatches[, s2 := ifelse(V1<variable, variable, V1)]
mismatches <- mismatches[, .(s1, s2, value)]
mismatches <- mismatches[s1<s2]
names(mismatches)[3] <- 'mismatch'
mismatches <- unique(mismatches)

pairs_ID <- merge(pairs_ID, mismatches, all.x=TRUE)

# Manually identified duplicate isolates to remove from dataset
trmi <- c(462, 355, 356, 461, 452, 446, 305, 435, 444, 348, 460, 454)

isolates_ID_match <- isolates_ID_match[ID %nin% trmi]
names(isolates_ID_match)[4] <- 'CYC_isolate'

pairs_ID  <- pairs_ID[s1 %nin% trmi & s2 %nin% trmi]


#
fwrite(isolates_ID_match, "data/temp/isolates_ID_match.csv")
fwrite(pairs_ID, "data/temp/pairs_ID.csv")
# write_csv(communities, "~/Dropbox/lab/emergent-coexistence/data/output/communities.csv")
