#' This script calculates the uncertainty in epslion from the uncertainties of OD, CFU, and V

# Read OD and CFU at T8
isolates_OD_CFU <- fread(here::here("data/temp/isolates_OD_CFU.csv"))

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
fwrite(isolates_epsilon_uncertainty, here::here("data/temp/isolates_epsilon_uncertainty.csv"))


if (FALSE) {





  #----
  isolates_epsilon_uncertainty <-
    isolates_OD_CFU %>%
    mutate(
      CFU = ColonyCount, OD = OD620, V=20, # Plating volume
      V1 = 10, V2 = 90, # Volume used in serial dilution
      ErrorCFU = sqrt(CFU), ErrorOD = 0.001, ErrorV = 0.4,
      ErrorV1 = 0.4, ErrorV2 = 2
    ) %>%
    mutate(
      n = DilutionFactor, # 4 or 5. shortened to n for conveience
      PartialV1 = n * V1^(n-1) * V2 / (V1+V2)^(n+1), #-n*(V1+V2)^(n-1)*V2/(V1^(n+1)),
      PartialV2 = -n * V1^(n-1) / (V1+V2)^(n+1), #n*(V1+V2)^(n-1)/(V1^n),
      DF = (V1/(V1+V2))^n,
      ErrorDF = sqrt((PartialV1*ErrorV1)^2 + (PartialV2*ErrorV2)^2)
    ) %>%
    select(-n, -V1, -V2, -ErrorV1, -ErrorV2, -PartialV1, -PartialV2) %>% # Remove transient variables
    mutate(
      PartialDF = CFU/(OD*V),
      PartialCFU = DF/(OD*V),
      PartialOD = (DF*CFU/(OD^2*V)),
      PartialV=(DF*CFU)/(OD*V^2)
    ) %>%
    mutate(
      Epsilon = OD * DF * V / CFU,
      ErrorEpsilon = sqrt((PartialDF*ErrorDF)^2 + (PartialCFU*ErrorCFU)^2 + (PartialOD*ErrorOD)^2 + (PartialV*ErrorV)^2)
    )


  # Apply the calculation to the real data
  # Uncertainty in cell frequency at T0 ----
  #The cellular frequency of isolate A at T0 has the following form:
  #$freq(A)=\frac{\frac{f_1CFU_A}{DF_A v_A}}{\frac{f_1CFU_A}{DF_A v_A} + \frac{f_2CFU_B}{DF_B v_B}}=\frac{1}{1+\frac{f_2 CFU_B DF_A v_A}{f_1 CFU_A DF_B v_B}}$, where the error in this frequency has the form $\delta freq(A)= \sqrt{(\frac{\partial freq(A)}{\partial CFU_A}\delta CFU_A)^2 + (\frac{\partial freq(A)}{\partial CFU_B}\delta CFU_B)^2 + (\frac{\partial freq(A)}{\partial DF_A}\delta DF_A)^2 + (\frac{\partial freq(A)}{\partial DF_B}\delta DF_B)^2 + (\frac{\partial freq(A)}{\partial v_A}\delta v_A)^2 + (\frac{\partial freq(A)}{\partial v_B}\delta v_B)^2}$


  # T0 cell frequency and errors


  temp <- colonies_T0_mono %>%
    select(Community, Isolate, ColonyCount, DilutionFactor) %>%
    arrange(Community) %>%
    right_join(isolates_OD_620, by = c("Community", "Isolate")) %>%
    mutate(OD = 0.1, V = 20, CFU = ColonyCount,
      ErrorOD = 0.001, ErrorV = 0.4) %>%
    # Uncertainty in dilution factor
    mutate(ErrorCFU = sqrt(CFU),
      n = DilutionFactor,  V1 = 10, V2 = 90, ErrorV1 = 0.4, ErrorV2 = 2,
      PartialV1 = -n*(V1+V2)^(n-1)*V2/(V1^(n+1)),
      PartialV2 = n*(V1+V2)^(n-1)/(V1^n),
      DF = ((V1+V2)/V1)^n,
      ErrorDF = sqrt((PartialV1*ErrorV1)^2 + (PartialV2*ErrorV2)^2)) %>%
    select(-n, -V1, -V2, -ErrorV1, -ErrorV2, -PartialV1, -PartialV2) %>%
    select(Community, Isolate, ColonyCount, ErrorCFU, DF, ErrorDF, V, ErrorV)
  temp <- mutate(temp, Isolate1 = Isolate, Isolate2 = Isolate) %>% select(-Isolate)

  # Pairs
  temp_pairs <- left_join(pairs_competition, select(temp, -Isolate2), by = c("Community", "Isolate1")) %>%
    left_join(select(temp, -Isolate1), by = c("Community", "Isolate2"))

  pairs_cell_T0 <-
    mutate(temp_pairs,
      a = (Isolate2Freq*ColonyCount.y*DF.x*V.x) / (Isolate1Freq*ColonyCount.x*DF.y*V.y),
      PartialCFU_A = 1/(1+a)^2 * a / ColonyCount.x,
      PartialCFU_B = -1/(1+a)^2 * a / ColonyCount.y,
      PartialDF_A = -1/(1+a)^2 * a / DF.x,
      PartialDF_B = 1/(1+a)^2 * a / DF.y,
      PartialV_A = -1/(1+a)^2 * a / V.x,
      PartialV_B = 1/(1+a)^2 * a / V.y,
      CellFreq1T0 = 1 / (1 + a),
      ErrorCellFreq1T0 = sqrt((PartialCFU_A*ErrorCFU.x)^2 + (PartialCFU_B*ErrorCFU.y)^2 +
          (PartialDF_A*ErrorDF.x)^2 + (PartialDF_B*ErrorDF.y)^2 +
          (PartialV_A*ErrorV.x)^2 + (PartialV_B*ErrorV.y)^2)
    ) %>%
    select(-starts_with("Partial"))






  # Uncertainty in cell frequency at T8 ----
  # The cell frequency at T8 from CFU is $freq(A)=\frac{CFU_A}{CFU_A+CFU_B}$. The error in this frequency has the following form $\delta freq(A)= \sqrt{(\frac{\partial freq(A)}{\partial CFU_A}\delta CFU_A)^2 + (\frac{\partial freq(A)}{\partial CFU_B}\delta CFU_B)^2}$

  pairs_cell_T8 <- temp_pairs %>%
    select(Community, Isolate1, Isolate2, Isolate1Freq, Isolate2Freq, ColonyCount, ColonyCount1, ColonyCount2) %>%
    mutate(ErrorCFU_A = sqrt(ColonyCount1), ErrorCFU_B = sqrt(ColonyCount2)) %>%
    mutate(PartialCFU_A = ColonyCount2/(ColonyCount1+ColonyCount2)^2,
      PartialCFU_B = -ColonyCount1/(ColonyCount1+ColonyCount2)^2,
      CellFreq1T8 = ColonyCount1 / (ColonyCount1+ColonyCount2),
      ErrorCellFreq1T8 = sqrt((PartialCFU_A*ErrorCFU_A)^2 + (PartialCFU_B*ErrorCFU_B)^2)) %>%
    select(-starts_with("Partial"))

  # Match T0 and T8
  pairs_cell <-
    left_join(pairs_cell_T0, pairs_cell_T8, by = c("Community", "Isolate1", "Isolate2", "Isolate1Freq", "Isolate2Freq")) %>%
    left_join(mutate(pairs, Isolate1 = factor(Isolate1), Isolate2 = factor(Isolate2))) %>%
    select(Community, Isolate1, Isolate2, Isolate1Freq, Isolate2Freq, ColonyCount, ColonyCount1, ColonyCount2,
      CellFreq1T0, CellFreq1T8, ErrorCellFreq1T0, ErrorCellFreq1T8)

  pairs_cell_gathered <- pairs_cell %>%
    gather("Time", "ErrorCellFreq1", 8:9) %>%
    gather("Time2", "CellFreq1", 6:7) %>%
    # Remove mismatched rows that match T0 mean to T8 error
    filter(sub("\\w+T", "", Time) == sub("\\w+T", "", Time2)) %>%
    mutate(Isolate1Freq = factor(Isolate1Freq), Isolate2Freq = factor(Isolate2Freq)) %>%
    mutate(Time = ifelse(Time == "ErrorCellFreq1T0", "T0", "T8")) %>%
    select(-Time2) %>%
    distinct()
}
