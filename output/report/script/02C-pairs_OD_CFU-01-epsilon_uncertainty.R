#' This script calculates the uncertainty in epsilon from the uncertainties of OD, CFU, and V

library(tidyverse)

# Read OD and CFU at T8
isolates_OD_CFU <- read_csv("~/Dropbox/lab/emergent-coexistence/data/temp/isolates_OD_CFU.csv", col_types = cols())

# Epsilon uncertainty
isolates_epsilon_uncertainty <-
  isolates_OD_CFU %>%
  mutate( # Variables
    CFU = ColonyCount, OD = OD620, V=20, # Plating volume
    V1 = 10, V2 = 90, # Volume used in serial dilution
    ErrorCFU = sqrt(CFU), ErrorOD = 0.001, ErrorV = 0.4,
    ErrorV1 = 0.4, ErrorV2 = 2
  ) %>%
  mutate( # Variable uncertainties
    n = DilutionFactor, # 4 or 5. shortened to n for conveience
    PartialV1 = n * V1^(n-1) * V2 / (V1+V2)^(n+1), #-n*(V1+V2)^(n-1)*V2/(V1^(n+1)),
    PartialV2 = -n * V1^(n-1) / (V1+V2)^(n+1), #n*(V1+V2)^(n-1)/(V1^n),
    DF = (V1/(V1+V2))^n,
    ErrorDF = sqrt((PartialV1*ErrorV1)^2 + (PartialV2*ErrorV2)^2)
  ) %>%
  select(-n, -V1, -V2, -ErrorV1, -ErrorV2, -PartialV1, -PartialV2) %>% # Remove transient variables
  mutate( # Partial derivatives
    PartialCFU = 1/(OD*DF*V),
    PartialOD = -CFU/(OD^2*DF*V),
    PartialDF = -CFU/(OD*DF^2*V),
    PartialV = -CFU/(OD*CFU*V^2)
  ) %>%
  mutate(
    Epsilon = CFU/(OD*DF*V),
    ErrorEpsilon = sqrt((PartialDF*ErrorDF)^2 + (PartialCFU*ErrorCFU)^2 + (PartialOD*ErrorOD)^2 + (PartialV*ErrorV)^2)
  )


# Set epsilon to NA if OD <= 0 or CFU == 0
isolates_epsilon_uncertainty[which(isolates_epsilon_uncertainty$OD <= 0 | isolates_epsilon_uncertainty$CFU == 0),
                             c("Epsilon", "ErrorEpsilon")] <- NA

# Set epsilon to NA if Colony counts <= 10
isolates_epsilon_uncertainty[isolates_epsilon_uncertainty$ColonyCount <= 10, c("Epsilon", "ErrorEpsilon")] <- NA

#
write_csv(isolates_epsilon_uncertainty, "~/Dropbox/lab/emergent-coexistence/data/temp/isolates_epsilon_uncertainty.csv")

