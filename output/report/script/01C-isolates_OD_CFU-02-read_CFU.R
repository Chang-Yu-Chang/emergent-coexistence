#' Read isolate OD
library(tidyverse)
library(data.table)

isolates_OD <- fread(here::here("data/temp/isolates_OD.csv"))

## Isolate colony count of monoculture at T8
isolates_CFU_T8 <- fread(here::here("data/raw/pairwise_competition/colony_count_dilution_factor.csv")) %>%
  filter(Isolate1Freq == "monoculture") %>% # Monoculture
  filter(!(Community == "C11R2"& Isolate1 == 13)) %>% # Remove contamination
  filter(!(Batch == "B2" & Community == "C11R1")) %>% # Remove Community C11R1 in batch B2
  filter(FileName != "D_T8_C1R4_3_-5") %>% # Manual remove the isolates that have two dilution factors
  #mutate(Isolate = ordered(Isolate1, levels = 1:12)) %>%
  select(Community, Isolate = Isolate1, ColonyCount, DilutionFactor)

"
There are two isolates that don't have colony at plating dilution factors.
Now I set the colony count in these plates to 0
"

isolates_CFU_T8$ColonyCount[isolates_CFU_T8$ColonyCount <= 0] <- 0

# Match CFU and OD ----
isolates_OD_CFU <- isolates_OD %>%
  mutate(OD620 = AbsMean) %>%
  select(Community, Isolate, OD620) %>%
  left_join(isolates_CFU_T8, by = c("Community", "Isolate"))
fwrite(isolates_OD_CFU, file = here::here("data/temp/isolates_OD_CFU.csv"))
