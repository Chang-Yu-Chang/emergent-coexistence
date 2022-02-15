#' Determine the pairwise interaction by the result of competition in three inital frequencies.
#' This script has the following content
#' 1. Read and plot the data of frequency changes
#' 2. Compute the sign of frequency changes
#' 3. Determine the interaction types by the fitness function
library(tidyverse)
library(data.table)
communities <- read_csv(here::here("data/output/communities.csv"))

# Frequency changes of pairs ----
pairs_freq <- read_csv(here::here("data/output/pairs_freq.csv")) %>%
    #mutate(Community = ordered(Community, levels = communities$Community)) %>%
    arrange(Time, Community, Isolate1, Isolate2, Isolate1InitialODFreq, Isolate2InitialODFreq)

## R Function for plotting frequencies changes
plot_frequency_change <- function (pairs_freq_comm) {
    comm <- unique(pairs_freq$Community)
    pairs_freq_comm %>%
        mutate(
            Isolate1 = ordered(Isolate1, 1:12),
            Isolate2 = ordered(Isolate2, 1:12),
            Isolate1InitialODFreq = factor(Isolate1InitialODFreq),
            Isolate2InitialODFreq = factor(Isolate2InitialODFreq),
            RawDataType = factor(RawDataType, levels = c("OD", "ODtoCFU", "CFU", "Sanger"))
        ) %>%
        ggplot(aes(x = Time, y = Isolate1MeasuredFreq, color = Isolate1InitialODFreq, group = Isolate1InitialODFreq)) +
        geom_point(aes(pch = RawDataType)) + geom_line() +
        geom_segment(aes(x = Time, xend = Time,
                         y = Isolate1MeasuredFreq - ErrorIsolate1MeasuredFreq,
                         yend = Isolate1MeasuredFreq + ErrorIsolate1MeasuredFreq)) +
        facet_grid(Isolate1 ~ Isolate2) +
        scale_y_continuous(limits = c(0, 1)) +
        theme_bw() +
        theme(legend.position = "top") +
        ggtitle(comm)
}

## Plot frequencies changes
p_pairs_freq_change_list <- rep(list(NA), length(communities$Community)) %>% setNames(communities$Community)
for (i in 1:length(communities$Community)) p_pairs_freq_change_list[[i]] <- pairs_freq %>% filter(Community == communities$Community[i]) %>% plot_frequency_change()

## Spread the df
pairs_freq_T0 <- pairs_freq %>%
    filter(Time == "T0") %>%
    mutate(Isolate1MeasuredFreqT0 = Isolate1MeasuredFreq,
           ErrorIsolate1MeasuredFreqT0 = ErrorIsolate1MeasuredFreq,
           RawDataTypeT0 = RawDataType) %>%
    select(Community, Isolate1, Isolate2, Isolate1InitialODFreq, Isolate2InitialODFreq,
           Isolate1MeasuredFreqT0, ErrorIsolate1MeasuredFreqT0, RawDataTypeT0)

pairs_freq_T8 <- pairs_freq %>%
    # cross-community and random networks are done in 3 transfers
    filter(Time == "T8" | Time == "T3") %>%
    mutate(Isolate1MeasuredFreqT8 = Isolate1MeasuredFreq,
           ErrorIsolate1MeasuredFreqT8 = ErrorIsolate1MeasuredFreq,
           Isolate1MeasuredFreq,
           RawDataTypeT8 = RawDataType) %>%
    select(Community, Isolate1, Isolate2, Isolate1InitialODFreq, Isolate2InitialODFreq,
           Isolate1MeasuredFreqT8, ErrorIsolate1MeasuredFreqT8, RawDataTypeT8)

pairs_freq_spread <- pairs_freq_T0 %>% left_join(pairs_freq_T8) %>%
    # Drop the incomplete pairs
    filter(!is.na(Isolate1MeasuredFreqT0), !is.na(Isolate1MeasuredFreqT8))


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

b = 10000
p <- 0.05 # Significant value

pairs_freq_spread$DifferenceT8T0 <- NA
pairs_freq_spread$DifferenceT8T0pvalue <- NA

for (i in 1:nrow(pairs_freq_spread)) {
    if (!is.na(pairs_freq_spread$Isolate1MeasuredFreqT8)) {
        T0_sd <- ifelse(is.na(pairs_freq_spread$ErrorIsolate1MeasuredFreqT0[i]), 0, pairs_freq_spread$ErrorIsolate1MeasuredFreqT0[i])
        T8_sd <- ifelse(is.na(pairs_freq_spread$ErrorIsolate1MeasuredFreqT8[i]), 0, pairs_freq_spread$ErrorIsolate1MeasuredFreqT8[i])

        # Draw from normal distribution
        df_temp <- data.frame(
            T0 = rnorm(b, mean = pairs_freq_spread$Isolate1MeasuredFreqT0[i], sd = T0_sd),
            T8 = rnorm(b, mean = pairs_freq_spread$Isolate1MeasuredFreqT8[i], sd = T8_sd))

        # p value
        temp <- sum((df_temp$T8 - df_temp$T0) > 0) / b
        pairs_freq_spread$DifferenceT8T0pvalue[i] <- min(temp, 1-temp)

        # Growing or decreasing
        # 1 means that the difference is growing from T0 to T8; 0 for non-significance; -1 for decreasing
        pairs_freq_spread$DifferenceT8T0[i] <- ifelse(temp > 0.5, 1, -1)
        #cat(i); cat(" ")
    }
}

## 0 for non-significance
pairs_freq_spread$DifferenceT8T0[pairs_freq_spread$DifferenceT8T0pvalue > p] <- 0



## Spread the df and paste the frequency changes. Reduce the row number to 186 (total 186 pairs)
pairs_interaction_fitness <- pairs_freq_spread %>%
    select(Community, Isolate1, Isolate2, Isolate1InitialODFreq, DifferenceT8T0) %>%
    group_by(Community, Isolate1, Isolate2) %>%
    pivot_wider(names_from = Isolate1InitialODFreq, values_from = DifferenceT8T0) %>%
    # Filter out those that are missing in the random networks (20211108)
    filter(!(is.na(`5`) & is.na(`95`) & str_detect(Community, "Ass"))) %>%
    mutate(FreqFunc = paste(`5`, `50`, `95`, sep = "_")) %>%
    select(-`5`, -`50`, -`95`)

# Extract the T8 frequencies ----
## Spread the df and compute the final frequencies. Reduce row number to 186 (total 186 pairs)
pairs_interaction_T8_freq <-
    pairs_freq_spread %>%
    select(Community, Isolate1, Isolate2, Isolate1InitialODFreq, Isolate1MeasuredFreqT8) %>%
    group_by(Community, Isolate1, Isolate2) %>%
    pivot_wider(names_from = Isolate1InitialODFreq, values_from = Isolate1MeasuredFreqT8) %>%
    mutate(PairIsolate1MeasuredFreqT8 = paste(round(`5`, 2), round(`50`, 2), round(`95`, 2), sep = "_"))


### Update the column `Isolate1Win`, which indicates the competitive exclusion that cannot be specified by fitness functions
pairs_interaction_T8_freq[pairs_interaction_T8_freq$PairIsolate1MeasuredFreqT8 == "1_1_1","Isolate1Win"] <- TRUE
pairs_interaction_T8_freq[pairs_interaction_T8_freq$PairIsolate1MeasuredFreqT8 == "0_0_0","Isolate1Win"] <- FALSE
pairs_interaction_T8_freq[pairs_interaction_T8_freq$PairIsolate1MeasuredFreqT8 == "1_NA_1","Isolate1Win"] <- TRUE
pairs_interaction_T8_freq[pairs_interaction_T8_freq$PairIsolate1MeasuredFreqT8 == "0_NA_0","Isolate1Win"] <- FALSE

pairs_interaction_T8_freq <- pairs_interaction_T8_freq %>%
    select(Community, Isolate1, Isolate2, PairIsolate1MeasuredFreqT8, Isolate1Win) %>%
    arrange(Community, Isolate1, Isolate2, PairIsolate1MeasuredFreqT8, Isolate1Win)


# Match two dfs: frequence changes and T8 frequenies ----
pairs_interaction_fitness <- pairs_interaction_fitness %>%
    left_join(pairs_interaction_T8_freq, by = c("Community", "Isolate1", "Isolate2"))

# Determine the interaction types by fitness functions ----
## Table for determining interaction types
interaction_type <- tibble(
    FromRare = rep(c(1, -1, 0), each = 9),
    FromMedium = rep(rep(c(1, -1, 0), each = 3), 3),
    FromAbundant = rep(c(1, -1, 0), 9),
    InteractionType = NA,
    InteractionTypeFiner = NA)

## Assign interaction types to combinations of frequency changes signs
interaction_type$InteractionType[c(1,14, 10, 13)] <- "exclusion"
interaction_type$InteractionType[c(2, 3, 5, 8, 9, 23, 26, 11, 12, 15, 17, 20)] <- "coexistence"
interaction_type$InteractionType[27] <- "neutrality"

## Assign finer interaction types to combinations of frequency changes signs
interaction_type$InteractionTypeFiner[c(1,14)] <- "competitive exclusion"
interaction_type$InteractionTypeFiner[c(10, 13)] <- "mutual exclusion"
interaction_type$InteractionTypeFiner[c(2, 3, 5, 8, 9, 23, 26)] <- "stable coexistence"
interaction_type$InteractionTypeFiner[c(11, 12, 15, 17, 20)] <- "frequency-dependent coexistence"
interaction_type$InteractionTypeFiner[27] <- "neutrality"

## Table for those with only two frequencies
interaction_type_twof <- tibble(
    FromRare = c(1,1,-1,-1),
    FromMedium = rep(NA, 4),
    FromAbundant = c(1,-1,1,-1),
    InteractionType = c("exclusion", "coexistence", "exclusion", "exclusion"),
    InteractionTypeFiner = c("competitive exclusion", "stable coexistence", "mutual exclusion", "competitive exclusion")
)

interaction_type <- bind_rows(interaction_type, interaction_type_twof) %>%
    mutate(FreqFunc = paste(FromRare, FromMedium, FromAbundant, sep = "_"))

## Join the fitness and interaction tables
pairs_interaction_fitness <- pairs_interaction_fitness %>%
    left_join(interaction_type, by = "FreqFunc") %>%
    as.data.table()

## Direction of links; from winner to loser
pairs_interaction_fitness[FreqFunc %in% c("1_1_1"), c("InteractionType", "InteractionTypeFiner", "From", "To") := list("exclusion", "competitive exclusion", Isolate1, Isolate2)]
pairs_interaction_fitness[FreqFunc %in% c("1_NA_1"), c("InteractionType", "InteractionTypeFiner", "From", "To") := list("exclusion", "competitive exclusion", Isolate1, Isolate2)]
pairs_interaction_fitness[FreqFunc %in% c("-1_-1_-1"), c("InteractionType", "InteractionTypeFiner", "From", "To") := list("exclusion", "competitive exclusion", Isolate2, Isolate1)]
pairs_interaction_fitness[FreqFunc %in% c("-1_NA_-1"), c("InteractionType", "InteractionTypeFiner", "From", "To") := list("exclusion", "competitive exclusion", Isolate2, Isolate1)]
pairs_interaction_fitness[!(FreqFunc %in% c("1_1_1", "1_NA_1", "-1_-1_-1", "-1_NA_-1")), c("From", "To") := list(Isolate1, Isolate2)]


### C11R2 isolate 3 and 6 freq 5:95 has contamination so its FreqFunc is not complete
pairs_interaction_fitness[Community == "C11R2" & Isolate1 == 3 & Isolate2 == 6, c("InteractionType", "InteractionTypeFiner", "From", "To") := list("exclusion", "competitive exclusion", Isolate2, Isolate1)]
### C2R8 isolate 3 and 4 freq 50:50 has contamination so its FreqFunc is not complete
pairs_interaction_fitness[Community == "C2R8" & Isolate1 == 3 & Isolate2 == 4, c("InteractionType", "InteractionTypeFiner", "From", "To") := list("coexistence", "stable coexistence", Isolate1, Isolate2)]


# Update interaction types by the final frequencies ----
pairs_interaction_fitness[Isolate1Win == TRUE, c("InteractionType", "InteractionTypeFiner", "From", "To") := list("exclusion", "competitive exclusion", Isolate1, Isolate2)]
pairs_interaction_fitness[Isolate1Win == FALSE, c("InteractionType", "InteractionTypeFiner", "From", "To") := list("exclusion", "competitive exclusion", Isolate2, Isolate1)]
# Select for the columns for pairs_interaction
pairs_interaction <- pairs_interaction_fitness %>%
    select(Community, Isolate1, Isolate2, InteractionType, InteractionTypeFiner, From, To)


# Interaction tables ----
# Fitness function: signs of frequencies changes in pairs
temp_df <- pairs_interaction_fitness %>%
    filter(!is.na(FreqFunc)) %>%
    group_by(FreqFunc, InteractionType, InteractionTypeFiner) %>%
    summarize(Count = n())

# Define the interaction table
signs <- 1:-1

## Full 27 combinations of fitness functions
pairs_interaction_table <- tibble(
    FromRare = rep(signs, each = 9),
    FromMedium = rep(rep(signs, each = 3), 3),
    FromAbundant = rep(signs, 9),
)

pairs_interaction_table$FreqFunc <-
    paste(round(pairs_interaction_table$FromRare, 2),
          round(pairs_interaction_table$FromMedium, 2),
          round(pairs_interaction_table$FromAbundant, 2), sep = "_")

## Merge
pairs_interaction_table <-
    pairs_interaction_table %>%
    left_join(temp_df) %>%
    select(-FreqFunc) %>%
    as_tibble()

ff <- function(x) {
    x[x==1] <- "+"
    x[x==-1] <- "-"
    x[x==0] <- "0"
    return(x)
}

pairs_interaction_table$FromRare <- ff(pairs_interaction_table$FromRare)
pairs_interaction_table$FromMedium <- ff(pairs_interaction_table$FromMedium)
pairs_interaction_table$FromAbundant <- ff(pairs_interaction_table$FromAbundant)
pairs_interaction_table <- pairs_interaction_table %>% drop_na(InteractionType)


write_csv(pairs_interaction, file = here::here("data/temp/pairs_interaction.csv"))
write_csv(interaction_type, file = here::here("data/temp/interaction_type.csv"))
write_csv(pairs_interaction_fitness, file = here::here("data/temp/pairs_interaction_fitness.csv"))
write_csv(pairs_interaction_table, file = here::here("data/output/pairs_interaction_table.csv"))
