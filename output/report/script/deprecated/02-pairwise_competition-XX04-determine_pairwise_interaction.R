#' This script determine the pairwise interactions


if (FALSE) {

  "

1. If the number of colony is greater than 100
- Yes, go to 2.
- No. Not valid.

2. Are the final frequencies of 95:5 and 5:95 both 0?
  - Yes, we can rule out coexistence. The interaction is either mutual exclusion or competitive exclusion.
- No, go to 3.
3. Do they converge to the same frequency? Test the same number of
- Yes, go to 3.1
- No, go to 4.
3.1 Is the convergence of frequencies neutral?
  - No
  "

  # Boostrapping between frequencies at T0 and T8 ----

  #To test whether the increment/decrement in the frequency is significant, here I compare the frequency of isolate A at T0 and T8 by bootstrapping the two frequencies. The frequency values are drawn from normal distribtution which has the mean and standard deviation (error). For each bootstrap, one frequency value from T0 and one frequency value from T8 is drawn accordingly from their distribution and then be compared. The sign of difference of these two values are then recorded for this bootstrap. Then I repeated the simulations for 100000 times.
  b = 100000
  p <- 0.05 # Significant value

  pairs_cell$DifferenceT8T0 <- NA
  pairs_cell$DifferenceT8T0pvalue <- NA

  for (i in 1:nrow(pairs_cell)) {
    df_temp <- data.frame(T0 = rnorm(b, mean = pairs_cell$CellFreq1T0[i], sd = pairs_cell$ErrorCellFreq1T0[i]),
      T8 = rnorm(b, mean = pairs_cell$CellFreq1T8[i], sd = pairs_cell$ErrorCellFreq1T8[i]))
    # p value
    temp <- sum((df_temp$T8 - df_temp$T0) > 0) / b
    pairs_cell$DifferenceT8T0pvalue[i] <- min(temp, 1-temp)

    # Growing or decreasing
    # 1 means that the difference is growing from T0 to T8; 0 for non-significance; -1 for decreasing
    pairs_cell$DifferenceT8T0[i] <- ifelse(temp > 0.5, 1, -1)
    #cat(i); cat(" ")
  }

  # 0 for non-significance
  pairs_cell$DifferenceT8T0[pairs_cell$DifferenceT8T0pvalue > p] <- 0

  # Pairs that have the T8 frequencies from three mixing ratios.
  pairs_cell_T8_mean <-
    pairs_cell %>% select(Community, Isolate1, Isolate2, Isolate1Freq, CellFreq1T8) %>%
    spread(Isolate1Freq, CellFreq1T8) %>%
    filter(!is.na(`5`), !is.na(`50`), !is.na(`95`))

  pairs_cell_T8_sd <-
    pairs_cell %>% select(Community, Isolate1, Isolate2, Isolate1Freq, ErrorCellFreq1T8) %>%
    spread(Isolate1Freq, ErrorCellFreq1T8) %>%
    filter(!is.na(`5`), !is.na(`50`), !is.na(`95`))



  # Interaction classification ----
  ## Step1: Is the colony number N>10?
  pairs_cell %>%
    filter(ColonyCount >= 10)  %>%
    pairs_df_count


  ## Step2: Are the final frequencies both 0?
  pairs_cell_T8_mean %>%
    filter(`5` %in% c(0,1), `95` %in% c(0,1))

  ## Step3: Do they converge to the same frequency?
  pairs_cell_T8_mean_coext <- pairs_cell_T8_mean %>% filter(!(`5` %in% c(0,1) & `95` %in% c(0,1)))
  pairs_cell_T8_sd_coext <- pairs_cell_T8_sd %>% filter(!(`5` %in% c(0,1) & `95` %in% c(0,1)))

  ### Boostrapping between T8 frequencies
  b = 100000
  p <- 0.05 # Significant value

  pairs_cell_T8_mean_coext$DifferenceF95F5 <- NA
  pairs_cell_T8_mean_coext$DifferenceF95F5pvalue <- NA

  for (i in 1:nrow(pairs_cell_T8_mean_coext)) {
    df_temp <- data.frame(F5 = rnorm(b, mean = pairs_cell_T8_mean_coext$`5`[i], sd = pairs_cell_T8_sd_coext$`95`[i]), F95 = rnorm(b, mean = pairs_cell_T8_mean_coext$`95`[i], sd = pairs_cell_T8_sd_coext$`95`[i]))
    # p value
    temp <- sum((df_temp$F95 - df_temp$F5) > 0) / b
    pairs_cell_T8_mean_coext$DifferenceF95F5pvalue[i] <- min(temp, 1-temp)

    # Growing or decreasing
    # 1 means that the difference is growing from T0 to T8; 0 for non-significance; -1 for decreasing
    pairs_cell_T8_mean_coext$DifferenceF95F5[i] <- ifelse(temp > 0.5, 1, -1)
    #cat(i); cat(" ")
  }

  # 0 for non-significance
  pairs_cell_T8_mean_coext$DifferenceF95F5[pairs_cell_T8_mean_coext$DifferenceF95F5pvalue > p] <- 0

  #Number of pairs that have the converged frequencies at T8.

  filter(pairs_cell_T8_mean_coext, DifferenceF95F5 == 0) %>%
    nrow()

  ### Step3.1: Is the convergence neutral?
  freqx = 0.95
  freqy = 0.05

  df <- data.frame(t=1:8, x=NA, y=NA)
  df$x[1] <- freqx; df$y[1] <- freqy

  for (t in 1:7) {
    err_x <- rnorm(1,0,0.02)
    err_y <- rnorm(1,0,0.02)
    df$x[t+1] <- df$x[t] * (1+err_x) / (df$x[t] * (1+err_x) + df$y[t] * (1+err_y))
    df$y[t+1] <- df$y[t] * (1+err_y) / (df$x[t] * (1+err_x) + df$y[t] * (1+err_y))
  }

  df



}

# Classification by the invasive fitness and final frequency ----

#Concanenate the fitness function to identify interaction types.
pairs_interaction <- pairs_cell %>%
  #  mutate(DiffernceT8T0pvalue)
  select(Community, Isolate1, Isolate2, Isolate1Freq, DifferenceT8T0) %>%
  #  group_by(Community, Isolate1, Isolate2) %>%
  filter(!is.na(DifferenceT8T0)) %>%
  spread(Isolate1Freq, DifferenceT8T0) %>%
  filter(!is.na(`5`), !is.na(`50`), !is.na(`95`)) %>%
  mutate(Interaction = paste(`5`, `50`, `95`, sep = ","))

### Competitive exclusion
#- **(-1,-1,-1), (1,1,1)**
filter(pairs_interaction, Interaction %in% c("-1,-1,-1", "1,1,1")) %>% nrow()

### Stable coexistence
#- **(1,1,-1), (1,-1,-1), (1,0,-1)**
filter(pairs_interaction, Interaction %in% c("1,1,-1", "1,-1,-1", "1,0,-1")) %>% nrow()

# Check pairs whose frequency changes meet the stable coexistence criteria but have 0 or 1 at final frequency.
# Only 15 pairs. There are two pairs (C11R1 1-6 and C11R2 8-10) have to be checked since they have at least one mixing ratio going to exclusion.
pairs_interaction_check <- filter(pairs_interaction, Interaction %in% c("1,1,-1", "1,-1,-1", "1,0,-1"))
for (i in 1:nrow(pairs_interaction_check)) {
  pairs_interaction_loop <- filter(pairs_cell, Community == pairs_interaction_check$Community[i],
    Isolate1 == pairs_interaction_check$Isolate1[i], Isolate2 == pairs_interaction_check$Isolate2[i])
  if (any(pairs_interaction_loop$CellFreq1T8 == 0) | any(pairs_interaction_loop$CellFreq1T8 == 1)) print(pairs_interaction_loop)
}


#- **(1,1,0), (0,1,-1)**
# There are 7 pairs converging to a coexistence state. The other 6 pairs have one mixing ratio going to exclusion (final frequency = 1 or 0).
pairs_interaction_check <- filter(pairs_interaction, Interaction %in% c("1,1,0", "0,-1,-1"))
for (i in 1:nrow(pairs_interaction_check)) {
  pairs_interaction_loop <- filter(pairs_cell, Community == pairs_interaction_check$Community[i],
    Isolate1 == pairs_interaction_check$Isolate1[i], Isolate2 == pairs_interaction_check$Isolate2[i])
  if (any(pairs_interaction_loop$CellFreq1T8 == 0) | any(pairs_interaction_loop$CellFreq1T8 == 1)) print(pairs_interaction_loop)
}
nrow(pairs_interaction_check)



## - **(1,0,0), (0,0,-1)**
# 1,0,0
#Only one pair (C2R6 2-4). It is probably stable coexistence but has to check the frequency difference at T8.
filter(pairs_cell, Community == "C2R6", Isolate1 == 2, Isolate2 ==4)


### Neutral
#- **(0,0,0)**

filter(pairs_interaction, Interaction %in% c("0,0,0")) %>% nrow()


### Bistability (two exclusions)
#- **(-1,1,1)**
# C2R6 1-4 and 3-4 pairs both diverge to two states that are competitive exclusion. Bistability.
filter(pairs_interaction, Interaction %in% c("-1,-1,1")) %>% nrow()


### Bistability (coexistence and exclusion)
#- **(-1,1,0)**
#C1R4 1-3 and C4R1 1-2 both have two states: coextence and exclusion, suggesting bistability.
filter(pairs_interaction, Interaction %in% c("-1,1,0")) %>% nrow()
#filter(pairs_cell, Community == "C1R4", Isolate1 == 1, Isolate2 == 3)
#filter(pairs_cell, Community == "C4R1", Isolate1 == 1, Isolate2 == 2)


#- **(-1,1,-1)**
# C2R6 1-3 has two states: coextence and exclusion, suggesting bistability.
filter(pairs_interaction, Interaction %in% c("-1,1,-1"))  %>% nrow()
#filter(pairs_cell, Community == "C2R6", Isolate1 == 1, Isolate2 == 3)


#- **(-1,0,-1)**
# Both C11R1 1-5 and C1R7 5-6 diverge to two states that coexistence and exclusion.

filter(pairs_interaction, Interaction %in% c("-1,0,-1")) %>% nrow()
#filter(pairs_cell, Community == "C11R1", Isolate1 == 1, Isolate2 == 5)
#filter(pairs_cell, Community == "C1R7", Isolate1 == 5, Isolate2 == 6)



### Undefined
#- **(1,1,-1), (1,0,-1), (1,-1,-1)**
#Only 2 pairs. There are two pairs (C11R1 1-6 and C11R2 8-10) have to be checked since they have at least one mixing ratio going to exclusion.
pairs_interaction_check <- filter(pairs_interaction, Interaction %in% c("1,1,-1", "1,-1,-1", "1,0,-1"))
for (i in 1:nrow(pairs_interaction_check)) {
  pairs_interaction_loop <- filter(pairs_cell, Community == pairs_interaction_check$Community[i],
    Isolate1 == pairs_interaction_check$Isolate1[i], Isolate2 == pairs_interaction_check$Isolate2[i])
  if (any(pairs_interaction_loop$CellFreq1T8 == 0) | any(pairs_interaction_loop$CellFreq1T8 == 1)) print(pairs_interaction_loop)
}


#- **(1,1,0), (0,-1,-1)**
#Only 6 pairs are undefined. The 6 pairs have one mixing ratio going to exclusion (final frequency = 1 or 0). The other 7 pairs converging to coexistence.
filter(pairs_interaction, Interaction %in% c("1,1,0", "0,-1,-1")) %>% nrow()



#- **(1,0,0), (0,0,-1)**
#Only 2 are undefined. C1R6 2-4 and C11R2 2-12 have cross-over in the frequency change plot, probably need to plate. The last one pair C2R6 2-4 is probably stable coexistence but has to check the frequency at T8.
filter(pairs_interaction, Interaction %in% c("1,0,0", "0,0,-1")) %>% nrow()
#filter(pairs_cell, Community == "C11R2", Isolate1 == 2, Isolate2 ==12)
#filter(pairs_cell, Community == "C1R6", Isolate1 == 3, Isolate2 ==4)
#filter(pairs_cell, Community == "C2R6", Isolate1 == 2, Isolate2 ==4)



#- **(-1,-1,0)**
#C1R2 has one coexistence and one exclusion state, suggesting probably bistability (coexistence and exclusion) but need to check. C11R2 4-11 has three states: two exclusions and one coexistence, probably bistability (two exclusion).
filter(pairs_interaction, Interaction %in% c("-1,-1,0")) %>%nrow()
#filter(pairs_cell, Community == "C1R2", Isolate1 == 1, Isolate2 == 3)
#filter(pairs_cell, Community == "C11R2", Isolate1 == 4, Isolate2 == 11)





# Interactions ----

### Classify the interactions
#Define the interaction types based on the criteria from the previous section.
pairs_cell_interaction <- as.data.table(pairs_interaction)

# Competitive exclusion
pairs_cell_interaction[Interaction %in% c("1,1,1", "-1,-1,-1"), InteractionType := "exclusion"]

# Stable coexistence
pairs_cell_interaction[Interaction %in% c("1,1,-1", "1,0,-1", "1,-1,-1", "1,1,0", "0,-1,-1", "1,0,0", "0,0,-1"),
  InteractionType := "stable coexistence"]

# Bistability (coexistence and exclusion)
pairs_cell_interaction[Interaction %in% c("-1,1,0", "-1,1,-1", "-1,0,-1"), InteractionType := "bistability (coexistence and exclusion)"]

# Bistability (two exclusions)
pairs_cell_interaction[Interaction %in% c("-1,-1,1"), InteractionType := "bistability (two exclusions)"]

# Neutral
pairs_cell_interaction[Interaction %in% c("0,0,0"), InteractionType := "neutral coexistence"]

# Not defined
pairs_cell_interaction[
  # (1,1,0) or (1,-1,-1)
  (Community == "C10R2" & Isolate1 == 2 & Isolate2 == 3) |
    (Community == "C11R1" & Isolate1 == 2 & Isolate2 == 3) |
    (Community == "C11R2" & Isolate1 == 1 & Isolate2 == 7) |
    (Community == "C11R2" & Isolate1 == 1 & Isolate2 == 12) |
    (Community == "C11R2" & Isolate1 == 3 & Isolate2 == 8) |
    (Community == "C11R2" & Isolate1 == 5 & Isolate2 == 8) |
    # (1,0,0) or (0,0,-1)
    (Community == "C1R6" & Isolate1 == 3 & Isolate2 == 4) |
    (Community == "C11R2" & Isolate1 == 2 & Isolate2 == 12) |
    # (-1,-1,0)
    (Community == "C1R2" & Isolate1 == 1 & Isolate2 == 3) |
    (Community == "C11R2" & Isolate1 == 4 & Isolate2 == 11) |
    # (1,X,-1)
    (Community == "C11R1" & Isolate1 == 1 & Isolate2 == 6) |
    (Community == "C11R2" & Isolate1 == 8 & Isolate2 == 10),
  InteractionType := "undefined"]





# ----
pairs_cell_interaction <- mutate(pairs_cell_interaction, InteractionType = factor(InteractionType, levels = c("exclusion", "stable coexistence", "bistability (coexistence and exclusion)", "bistability (two exclusions)", "neutral coexistence", "undefined")))




