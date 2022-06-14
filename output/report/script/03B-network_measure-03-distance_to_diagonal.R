#' Count the number of pairwise coexistence as a function of distance to diagonal (|i-j|)
library(tidyverse)
library(tidygraph)
library(igraph)
source(here::here("plotting_scripts/misc.R"))
communities <- read_csv("~/Dropbox/lab/emergent-coexistence/data/output/communities.csv", col_types = cols())
load("~/Dropbox/lab/emergent-coexistence/data/output/communities_network.Rdata")
load("~/Dropbox/lab/emergent-coexistence/data/output/communities_network_randomized.Rdata")


# Observation
communities_diag <- communities_network %>%
    filter(str_detect(Community, "C\\d")) %>%
    rowwise() %>%
    mutate(Diagonal = count_diag_coexistence(Network) %>% mutate(Community) %>% list()) %>%
    pull(Diagonal) %>%
    bind_rows()

# Permutation
communities_diag_randomized_list <- rep(list(NA), nrow(communities_network))
names(communities_diag_randomized_list) <- communities_network$Community
tt <- proc.time()
for (i in 1:nrow(communities_network)) {
    temp_tt <- proc.time()
    communities_diag_randomized_list[[i]] <- lapply(communities_network_randomized$NetworkRandomized[[i]], count_diag_coexistence) %>%
        bind_rows(.id = "Replicate")
    # Print
    cat("\n\n", communities$Community[i])
    cat("\n", (proc.time() - temp_tt)[3], "seconds")
    if (i == nrow(communities_network)) cat("\n\n total time:", (proc.time() - tt)[3], "seconds")
}
communities_diag_randomized <- bind_rows(communities_diag_randomized_list, .id = "Community")

#
write_csv(communities_diag, "~/Dropbox/lab/emergent-coexistence/data/output/communities_diag.csv")
write_csv(communities_diag_randomized, "~/Dropbox/lab/emergent-coexistence/data/output/communities_diag_randomized.csv")










