#' Determine the pairwise interaction by the result of competition in three inital frequencies.
#' This script has the following content
#' 1 Compute the sign of frequency changes
#' 2. Determine the interaction types by the fitness function
library(tidyverse)

communities <- read_csv("~/Dropbox/lab/emergent-coexistence/data/output/communities.csv", col_types = cols())
pairs_freq <- read_csv("~/Dropbox/lab/emergent-coexistence/data/output/pairs_freq.csv", col_types = cols())

#
pairs_freq_wide <- pairs_freq %>%
    select(Set, everything()) %>%
    select(-RawDataType) %>%
    # If Isolate1MeasuredFreq is available and the error is NA, set the error to 0
    mutate(ErrorIsolate1MeasuredFreq = case_when(
        !is.na(Isolate1MeasuredFreq) & is.na(ErrorIsolate1MeasuredFreq) ~ 0,
        !is.na(Isolate1MeasuredFreq) & !is.na(ErrorIsolate1MeasuredFreq) ~ ErrorIsolate1MeasuredFreq
    )) %>%
    pivot_wider(names_from = Time, values_from = c(Isolate1MeasuredFreq, ErrorIsolate1MeasuredFreq), names_sep = "") %>%
    mutate(DifferenceTiniTend = NA, DifferenceTiniTendpvalue = NA)


# Compute the sign of frequency changes between T0 and T8 ----
## Boostrapping between frequencies at T0 and T8
#' To test whether the increment/decrement in the frequency is significant,
#' here I compare the frequency of isolate A at T0 and T8 by bootstrapping the
#' two frequencies. The frequency values are drawn from normal distribtution
#' which has the mean and standard deviation (error). For each bootstrap,
#'  one frequency value from T0 and one frequency value from T8 is drawn
#'  accordingly from their distribution and then be compared. The sign of
#'  difference of these two values are then recorded for this bootstrap.
#'  Then I repeated the simulations for 10000 times.


b = 1000
p <- 0.05 # Significant value

names(pairs_freq_wide)

for (i in 1:nrow(pairs_freq_wide)) {
    if (!is.na(pairs_freq_wide$Isolate1MeasuredFreqTend[i])) {
        T0_sd <- ifelse(is.na(pairs_freq_wide$ErrorIsolate1MeasuredFreqTini[i]), 0, pairs_freq_wide$ErrorIsolate1MeasuredFreqTini[i])
        T8_sd <- ifelse(is.na(pairs_freq_wide$ErrorIsolate1MeasuredFreqTend[i]), 0, pairs_freq_wide$ErrorIsolate1MeasuredFreqTend[i])

        # Draw from normal distribution
        df_temp <- data.frame(
            T0 = rnorm(b, mean = pairs_freq_wide$Isolate1MeasuredFreqTini[i], sd = T0_sd),
            T8 = rnorm(b, mean = pairs_freq_wide$Isolate1MeasuredFreqTend[i], sd = T8_sd))

        # p value
        temp <- sum((df_temp$T8 - df_temp$T0) > 0) / b
        pairs_freq_wide$DifferenceTiniTendpvalue[i] <- min(temp, 1-temp)

        # Growing or decreasing
        # 1 means that the difference is growing from T0 to T8; 0 for non-significance; -1 for decreasing
        pairs_freq_wide$DifferenceTiniTend[i] <- ifelse(temp > 0.5, 1, -1)
    }
}

# 0 for non-significance
pairs_freq_wide$DifferenceTiniTend[pairs_freq_wide$DifferenceTiniTendpvalue > p] <- 0

# Determine the direction of exclusion pairs by extracting the T8 frequencies
## Pivot wider the df and compute the final frequencies. Reduce row number to 186 (total 186 pairs)
pairs_interaction_Tend_freq <- pairs_freq_wide %>%
    group_by(Set, PairID, Community, Isolate1, Isolate2) %>%
    select(group_cols(), Isolate1InitialODFreq, Isolate1MeasuredFreqTend) %>%
    pivot_wider(names_from = Isolate1InitialODFreq, values_from = Isolate1MeasuredFreqTend) %>%
    mutate(PairIsolate1MeasuredFreqTend = paste(round(`5`, 2), round(`50`, 2), round(`95`, 2), sep = "_")) %>%
    select(-`5`, -`50`, -`95`)

if (FALSE) {
    ## Update the column `Isolate1Win`, which indicates the competitive exclusion that cannot be specified by fitness functions
    pairs_interaction_Tend_freq[pairs_interaction_Tend_freq$PairIsolate1MeasuredFreqTend == "1_1_1","Isolate1Win"] <- TRUE
    pairs_interaction_Tend_freq[pairs_interaction_Tend_freq$PairIsolate1MeasuredFreqTend == "0_0_0","Isolate1Win"] <- FALSE
    pairs_interaction_Tend_freq[pairs_interaction_Tend_freq$PairIsolate1MeasuredFreqTend == "1_NA_1","Isolate1Win"] <- TRUE
    pairs_interaction_Tend_freq[pairs_interaction_Tend_freq$PairIsolate1MeasuredFreqTend == "0_NA_0","Isolate1Win"] <- FALSE

    pairs_interaction_Tend_freq <- pairs_interaction_Tend_freq %>%
        select(Community, Isolate1, Isolate2, PairIsolate1MeasuredFreqTend, Isolate1Win) %>%
        arrange(Community, Isolate1, Isolate2, PairIsolate1MeasuredFreqTend, Isolate1Win)

}


# Determine the interaction types by fitness functions ----
# Table for determining interaction types
interaction_type_three <- tibble(
    FromRare = rep(c(1, -1, 0), each = 9),
    FromMedium = rep(rep(c(1, -1, 0), each = 3), 3),
    FromAbundant = rep(c(1, -1, 0), 9),
    InteractionType = NA,
    InteractionTypeFiner = NA
)
interaction_type_two <- tibble(
    FromRare = rep(c(1, -1, 0), each = 3),
    FromMedium = rep(NA, 9),
    FromAbundant = rep(c(1, -1, 0), 3),
    InteractionType = NA,
    InteractionTypeFiner = NA
)

interaction_type <- bind_rows(interaction_type_three, interaction_type_two)

## Assign interaction types to combinations of frequency changes signs
interaction_type$InteractionType[c(1,14, 10, 13, 28, 31, 32)] <- "exclusion"
interaction_type$InteractionType[c(2, 3, 5, 8, 9, 23, 26, 11, 12, 15, 17, 20, 29:30, 33:35)] <- "coexistence"
interaction_type$InteractionType[c(27, 36)] <- "neutrality"

## Assign finer interaction types to combinations of frequency changes signs
interaction_type$InteractionTypeFiner[c(1, 14, 28, 32)] <- "competitive exclusion"
interaction_type$InteractionTypeFiner[c(10, 13, 31)] <- "mutual exclusion"
interaction_type$InteractionTypeFiner[c(2, 3, 5, 8, 9, 23, 26, 29, 30, 33)] <- "stable coexistence"
interaction_type$InteractionTypeFiner[c(11, 12, 15, 17, 20, 34:35)] <- "frequency-dependent coexistence"
interaction_type$InteractionTypeFiner[c(27,36)] <- "neutrality"

interaction_type <- interaction_type %>%  mutate(FreqFunc = paste(FromRare, FromMedium, FromAbundant, sep = "_"))




# Combine the frequency data with the isolate data
pairs_interaction_fitness <- pairs_freq_wide %>%
    group_by(Set, PairID, Community, Isolate1, Isolate2) %>%
    select(group_cols(), Isolate1InitialODFreq, DifferenceTiniTend) %>%
    pivot_wider(names_from = Isolate1InitialODFreq, values_from = DifferenceTiniTend) %>%
    mutate(FreqFunc = paste(`5`, `50`, `95`, sep = "_")) %>%
    select(-`5`, -`50`, -`95`) %>%
    # Join two dfs: frequency changes and Tend frequency
    left_join(pairs_interaction_Tend_freq) %>%
    # Join the fitness and interaction tables
    left_join(interaction_type, by = "FreqFunc") %>%
    # Direction of links in exclusion pairs; from winner to loser
    mutate(From = case_when(
        FreqFunc %in% c("1_1_1", "1_NA_1") ~ Isolate1, # Isolate1 wins
        FreqFunc %in% c("-1_-1_-1", "-1_NA_-1") ~ Isolate2, # Isolate2 wins
        TRUE ~ Isolate1
    )) %>%
    mutate(To = case_when(
        FreqFunc %in% c("1_1_1", "1_NA_1") ~ Isolate2,
        FreqFunc %in% c("-1_-1_-1", "-1_NA_-1") ~ Isolate1,
        TRUE ~ Isolate2
    )) %>%
    ungroup()

pairs_interaction <- pairs_interaction_fitness %>%
    select(Set, Community, Isolate1, Isolate2, InteractionType, InteractionTypeFiner, From, To)



# Print the interaction tables ----
ff <- function(x) {
    x[x==1] <- "+"
    x[x==-1] <- "-"
    x[x==0] <- "0"
    return(x)
}

pairs_interaction_table <- pairs_interaction_fitness %>%
    filter(!is.na(FreqFunc)) %>%
    group_by(Set, FreqFunc, InteractionType, InteractionTypeFiner) %>%
    summarize(Count = n()) %>%
    ungroup() %>%
    separate(col = FreqFunc, into = c("FromRare", "FromMedium", "FromAbundant"), sep = "_") %>%
    mutate(across(starts_with("From"), ff)) %>%
    drop_na(InteractionType)


write_csv(pairs_interaction_fitness, "~/Dropbox/lab/emergent-coexistence/data/temp/pairs_interaction_fitness.csv")
write_csv(pairs_interaction, "~/Dropbox/lab/emergent-coexistence/data/output/pairs_interaction.csv")
write_csv(pairs_interaction_table, "~/Dropbox/lab/emergent-coexistence/data/output/pairs_interaction_table.csv")
