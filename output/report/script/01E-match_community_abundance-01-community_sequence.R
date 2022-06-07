#' Read and clean up community OTU table from the curated community files from Sylvie
library(tidyverse)

communities <- read_csv("~/Dropbox/lab/emergent-coexistence/data/output/communities.csv", col_types = cols())

# This csv only has transfer 12
communities_abundance <- read_csv("~/Dropbox/lab/emergent-coexistence/data/raw/community_ESV/Emergent_Simplicity_Equilibrium_Data.csv", col_types = cols()) %>%
    # Tidy up the variable names
    rename(CommunityESVID = ESV_ID, SampleID = Sample_ID, RelativeAbundance = Relative_Abundance, CarbonSource = Carbon_Source,
           ESVFamily = Family, ESVGenus = Genus) %>%
    # The inoculum, replicate, and sample ID in this df is zero-based, so one has to convert it into a one-based system
    mutate(Community = paste0("C", Inoculum+1, "R", Replicate+1)) %>%
    filter(Community %in% communities$Community) %>%
    mutate(Community = factor(Community, levels = communities$Community)) %>%
    filter(CarbonSource == "Glucose") %>%
    arrange(Community, CarbonSource, desc(RelativeAbundance)) %>%
    select(SampleID, Community, RelativeAbundance, CommunityESVID, ESV, ESVFamily, ESVGenus)


write_csv(communities_abundance, "~/Dropbox/lab/emergent-coexistence/data/temp/communities_abundance.csv")

