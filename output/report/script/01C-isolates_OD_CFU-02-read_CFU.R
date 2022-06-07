#' Read isolate OD
library(tidyverse)

isolates_OD <- read_csv("~/Dropbox/lab/emergent-coexistence/data/temp/isolates_OD.csv", col_types = cols())

## Isolate colony count of monoculture at T8
isolates_CFU_T8 <- read_csv("~/Dropbox/lab/emergent-coexistence/data/raw/pairwise_competition/colony_count_dilution_factor.csv", col_types = cols()) %>%
    filter(Isolate1Freq == "monoculture") %>% # Monoculture
    filter(!(Community == "C11R2"& Isolate1 == 13)) %>% # Remove staph contamination
    filter(!(Batch == "B2" & Community == "C11R1")) %>% # Remove Community C11R1 in batch B2
    filter(FileName != "D_T8_C1R4_3_-5") %>% # Manually remove the isolates that have two dilution factors
    select(Community, Isolate = Isolate1, ColonyCount, DilutionFactor)

#' Note that there are two isolates that don't have colony at plating dilution
#' factors, so I set the colony count in these plates to 0
isolates_CFU_T8$ColonyCount[isolates_CFU_T8$ColonyCount <= 0] <- 0

# Match CFU and OD
isolates_OD_CFU <- isolates_OD %>%
    mutate(OD620 = AbsMean) %>%
    filter(Transfer == 8, Wavelength == 620) %>%
    select(Community, Isolate, OD620) %>%
    left_join(isolates_CFU_T8, by = c("Community", "Isolate"))
write_csv(isolates_OD_CFU, "~/Dropbox/lab/emergent-coexistence/data/temp/isolates_OD_CFU.csv")

