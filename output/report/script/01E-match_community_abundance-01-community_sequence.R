#' Read and clean up community OTU table
#' There are two version of community data:
#' 1. Raw data from Nanxi (deprecated)
#' 2. Curated community files from Sylvie
#'
#' This script first cleans up the OTU table from Nanxi, which includes the step
#' 1) Read and melt the dataframe and 2) Rename the variables
#' Then work on the curated community abundance data from Sylvie
library(tidyverse)
library(data.table)

# Curated community abundances ----
# Read data
communities_abundance_syl <- fread(here::here("data/raw/community_ESV/Emergent_Simplicity_Equilibrium_Data.csv"))

# Tidy up the variable names
#' The inoculum, replicate, and sample ID in this df is zero-based,
#' so one has to convert it into a one-based system
names(communities_abundance_syl) <- c("CommunityESVID", "ESV", "Order", "Family", "Genus", "SampleID", "RelativeAbundance", "CarbonSource", "Inoculum", "Replicate", "Transfer")
communities_abundance_syl[,"Community"] <-
  paste0("C", communities_abundance_syl$Inoculum + 1,
    "R", communities_abundance_syl$Replicate + 1) # The number here is zero-based
communities_abundance_syl <- communities_abundance_syl %>%
  dplyr::select(SampleID, CarbonSource, Community, Transfer, RelativeAbundance, CommunityESVID, ESV, Order, Family, Genus) %>%
  mutate(Community = ordered(Community, levels = paste0("C", rep(1:12, each = 8), "R", rep(1:8, 12)))) %>%
  arrange(Community, CarbonSource, desc(RelativeAbundance)) %>%
  as_tibble()

fwrite(communities_abundance_syl, here::here("data/temp/communities_abundance_syl.csv"))


if (FALSE) {
    # Raw OTU table from Nanxi ----
    # Read raw data from Rdata file
    load(here::here("data/raw/community_ESV/Glucose_all_from_mergeddataall_20180727.Rdata"))
    d <- data.table(Glucose_all); rm(Glucose_all)
    dim(d) # 474 communities (or replicate) X 4328 (4340 substracts 12)
    colnames(d[,1:12])
    d <- gather(d, "Sequence", "Abundance", -c(1:12)) %>% filter(Abundance != 0) # Melt data.frame and remove 0
    dim(d) # 9545 sequences; now each row represents one sequence in a community.
    d$Rep[is.nan(d$Rep)] <- "pool" # Species pool

    # Specify the variable names, such as community (regional pool), replicate, and transfer.
    communities_abundance <- d %>%
        mutate(
            SampleID = Experiment,
            Community = paste0("C", as.numeric(gsub("[a-zA-Z]", "", Comm1)), "R", Rep),
            Community = ordered(Community, levels = c(paste0("C", 1:12, "Rpool"), paste0("C", rep(1:12, each = 8), "R", rep(1:8, 12)))),
            ESV = Sequence
        ) %>%
        arrange(Community, Transfer, desc(Abundance)) %>%
        # Sequence ID specific to community-sequence
        mutate(CommunityESVID = 1:n()) %>%
        select(SampleID, Community, Transfer, Abundance, CommunityESVID, ESV) %>%
        as_tibble()

    # Calculate relative abundances
    communities_abundance <-
        communities_abundance %>%
        group_by(SampleID, Community, Transfer) %>%
        mutate(RelativeAbundance = round(Abundance / sum(Abundance), 4)) %>%
        #  select(Community, Abundance, RelativeAbundance) %>%
        #  filter(RelativeAbundance >= 0.001) %>%
        ungroup() %>%
        select(SampleID, Community, Transfer, Abundance, RelativeAbundance, CommunityESVID, ESV) %>%
        as_tibble()

    # Community richness
    communities_richness <-
        communities_abundance %>%
        #  filter(Community %in% communities_name) %>%
        mutate(Community = ordered(Community, levels = communities_name),
               Transfer = ordered(Transfer, levels = 1:12)) %>%
        group_by(Community, Transfer) %>%
        summarize(Richness = n()) %>%
        arrange(Community, Transfer)

    #
    fwrite(communities_richness, here::here("data/temp/communities_richness.csv"))
    fwrite(communities_abundance, here::here("data/temp/communities_abundance.csv"))
}
