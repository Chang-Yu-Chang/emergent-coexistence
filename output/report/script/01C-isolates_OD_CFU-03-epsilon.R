#' OD and CFU conversion
library(tidyverse)
isolates_OD_CFU <- read_csv("~/Dropbox/lab/emergent-coexistence/data/temp/isolates_OD_CFU.csv", col_types = cols())

# OD-CFU converion coefficient epsilon ----
# From T8 monoculture, calculate the conversion coefficient $\epsilon$ from OD to CFU.
#' $\epsilon_A = \frac{OD_A DF_A v_A}{CFU_A}$, where
#' D is the dilution factor for plating (for instance, $10^{5}$),
#' CFU is the colony count under given dilution factor,
#' OD is the optical density measured at transfer 8 after 48 hours,
#' $v$ is the plating volume (20 uL)

isolates_epsilon <- isolates_OD_CFU %>%
  mutate(Epsilon = OD620 * 10^(-DilutionFactor) * 20 /ColonyCount) %>%
  select(Community, Isolate, ColonyCount, DilutionFactor, OD620, Epsilon)

write_csv(isolates_epsilon, "~/Dropbox/lab/emergent-coexistence/data/temp/isolates_epsilon.csv")
