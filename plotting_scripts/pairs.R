# CASEU and colony counting

library(tidyverse)
library(data.table)
library(cowplot)
source(here::here("plotting_scripts/network_functions.R"))

isolates <- read_csv(here::here("data/output/isolates.csv"), col_types = cols())
pairs <- read_csv(here::here("data/output/pairs.csv"), col_types = cols()) %>% mutate(InteractionType = ifelse(InteractionType == "neutrality", "coexistence", InteractionType))
#pairs_freq <- read_csv(here::here("data/output/pairs_freq.csv"), col_types = cols())
#pairs_example_outcomes_finer <- read_csv(here::here("data/output/pairs_example_outcomes_finer.csv"), col_types = cols())
communities <- read_csv(here::here("data/output/communities.csv"), col_types = cols()) %>%
    filter(Assembly == "self_assembly") %>%
    arrange(CommunitySize) %>%
    mutate(CommunityLabel = 1:13) %>%
    select(Community, CommunityLabel, CommunitySize, CommunityPairSize)
community_factor <- communities$Community
communities_size <- communities$CommunitySize


# CFU results
pairs_cfu <- read_csv(here::here("data/temp/pairs_CFU_freq_uncertainty.csv"), col_types = cols()) %>%
    #mutate(Isolate1 = as.character(Isolate1), Isolate2 = as.character(Isolate2)) %>%
    # Configurate the variable names
    mutate(Isolate1InitialODFreq = Isolate1Freq,
           Isolate2InitialODFreq = Isolate2Freq,
           Isolate1MeasuredFreq = Isolate1CFUFreq,
           ErrorIsolate1MeasuredFreq = ErrorIsolate1CFUFreq) %>%
    select(Community, Isolate1, Isolate2, Isolate1InitialODFreq, Isolate2InitialODFreq, Time, Isolate1MeasuredFreq, ErrorIsolate1MeasuredFreq, RawDataType, Contamination)

# CASEU results from pilot
switch_pairwise_column <- function (df, bypair = T) {
    if (any(is.factor(df$Isolate1))) df$Isolate1 <- as.numeric(df$Isolate1); df$Isolate2 <- as.numeric(df$Isolate2)
    if ("Isolate1FreqPredicted" %in% colnames(df)) {
        if (bypair == T) {
            temp_index <- df$Isolate1 > df$Isolate2
            df[temp_index, c("Isolate1", "Isolate2", "Isolate1Freq", "Isolate2Freq", "Isolate1FreqPredicted", "Isolate2FreqPredicted")] <-
                df[temp_index, c("Isolate2", "Isolate1", "Isolate2Freq", "Isolate1Freq", "Isolate2FreqPredicted", "Isolate1FreqPredicted")]

            df %>% arrange(Isolate1, Isolate2, Isolate1Freq) %>% return()
        } else if (bypair == F) {
            temp_index <- df$Isolate1Freq == 5
            df[temp_index, c("Isolate1", "Isolate2", "Isolate1Freq", "Isolate2Freq", "Isolate1FreqPredicted", "Isolate2FreqPredicted")] <-
                df[temp_index, c("Isolate2", "Isolate1", "Isolate2Freq", "Isolate1Freq", "Isolate2FreqPredicted", "Isolate1FreqPredicted")]

            df %>% arrange(Isolate1Freq, Isolate1, Isolate2) %>% return()
        }
    } else {

        if (bypair == T) {
            temp_index <- df$Isolate1 > df$Isolate2
            df[temp_index, c("Isolate1", "Isolate2", "Isolate1Freq", "Isolate2Freq")] <-
                df[temp_index, c("Isolate2", "Isolate1", "Isolate2Freq", "Isolate1Freq")]

            df %>% arrange(Isolate1, Isolate2, Isolate1Freq) %>% return()
        } else if (bypair == F) {
            temp_index <- df$Isolate1Freq == 5
            df[temp_index, c("Isolate1", "Isolate2", "Isolate1Freq", "Isolate2Freq")] <-
                df[temp_index, c("Isolate2", "Isolate1", "Isolate2Freq", "Isolate1Freq")]

            df %>% arrange(Isolate1Freq, Isolate1, Isolate2) %>% return()
        }
    }
}
CASEU_pilot2 <- read_csv(here::here("data/temp/CASEU_pilot2.csv")) %>%
    filter(Treatment == "C11R1") %>%
    mutate(Isolate1 = as.numeric(Isolate1), Isolate2 = as.numeric(Isolate2)) %>%
    switch_pairwise_column(bypair = T) %>%
    mutate(
        Community = Treatment,
        Isolate1Freq = Isolate1Freq * 100,
        Isolate2Freq = Isolate2Freq * 100,
        Time = "T8",
        Isolate1CFUFreq = Isolate1FreqPredicted,
        ErrorIsolate1CFUFreq = NA,
        RawDataType = "Sanger"
    ) %>%
    select(Community, Isolate1, Isolate2, Isolate1Freq, Isolate2Freq, Time, Isolate1CFUFreq, ErrorIsolate1CFUFreq, RawDataType)

CASEU_pilot3 <- read_csv(here::here("data/temp/CASEU_pilot3.csv")) %>%
    mutate(
        Isolate1Freq = Isolate1Freq * 100,
        Isolate2Freq = Isolate2Freq * 100,
        Time = "T8",
        Isolate1CFUFreq = Isolate1FreqPredicted,
        ErrorIsolate1CFUFreq = NA,
        RawDataType = "Sanger"
    ) %>%
    select(Community, Isolate1, Isolate2, Isolate1Freq, Isolate2Freq, Time, Isolate1CFUFreq, ErrorIsolate1CFUFreq, RawDataType)

CASEU_pilot4 <- read_csv(here::here("data/temp/CASEU_pilot4.csv")) %>%
    mutate(
        Isolate1Freq = Isolate1Freq * 100,
        Isolate2Freq = Isolate2Freq * 100,
        Time = "T8",
        Isolate1CFUFreq = Isolate1FreqPredicted,
        ErrorIsolate1CFUFreq = NA,
        RawDataType = "Sanger"
    ) %>%
    select(Community, Isolate1, Isolate2, Isolate1Freq, Isolate2Freq, Time, Isolate1CFUFreq, ErrorIsolate1CFUFreq, RawDataType)


pairs_caseu_pilot <- bind_rows(CASEU_pilot2, CASEU_pilot3, CASEU_pilot4) %>%
    rename(Isolate1InitialODFreq = Isolate1Freq, Isolate2InitialODFreq = Isolate2Freq, Isolate1MeasuredFreq = Isolate1CFUFreq)


# CASEU result from the 6 plates batch
pairs_caseu_sixplates <- read_csv("~/Dropbox/lab/invasion-network/data/temp/CASEU_six_plates_result.csv", col_types = cols()) %>%
    select(Community, Isolate1, Isolate2, Isolate1Freq, Isolate2Freq, Time, Isolate1FreqPredicted) %>%
    rename(Isolate1InitialODFreq = Isolate1Freq, Isolate2InitialODFreq = Isolate2Freq, Isolate1MeasuredFreq = Isolate1FreqPredicted) %>%
    mutate(RawDataType = "Sanger") %>%
    # Remove the contaminant
    filter(!(Community == "C11R2" & Isolate1 == 13)) %>%
    filter(!(Community == "C11R2" & Isolate2 == 13))
#sum(is.na(pairs_caseu_sixplates$Isolate1MeasuredFreq)) # 19 pairs with short trace

## Switch isolate order
temp_index <- pairs_caseu_sixplates$Isolate1 > pairs_caseu_sixplates$Isolate2
temp <- pairs_caseu_sixplates[temp_index,] %>% rename(Isolate1 = Isolate2, Isolate2 = Isolate1, Isolate1InitialODFreq = Isolate2InitialODFreq, Isolate2InitialODFreq = Isolate1InitialODFreq)
pairs_caseu_sixplates <- bind_rows(pairs_caseu_sixplates[!temp_index,], temp) %>%
    arrange(Community, Isolate1, Isolate2, Isolate1InitialODFreq)




# Use the casue pilot data to fill the short-trace pairs in the caseu six plates
## Only three fulfilled in this way, leaving 16 pair-frequency out of 372 in total that do not have the caseu data
temp <- pairs_caseu_pilot %>%
    filter(Isolate1InitialODFreq %in% c(5, 95)) %>%
    select(Community, Isolate1, Isolate2, Isolate1InitialODFreq, temp = Isolate1MeasuredFreq)

pairs_caseu_sixplates_short_trace <- pairs_caseu_sixplates %>%
    filter(is.na(Isolate1MeasuredFreq)) %>%
    left_join(temp) %>%
    mutate(Isolate1MeasuredFreq = ifelse(!is.na(temp), temp, NA)) %>%
    select(-temp)

pairs_caseu <- pairs_caseu_sixplates %>%
    filter(!is.na(Isolate1MeasuredFreq)) %>%
    bind_rows(pairs_caseu_sixplates_short_trace) %>%
    arrange(Community, Isolate1, Isolate2, Isolate1InitialODFreq)




# Compare caseu and cfu
pairs_freq <- bind_rows(pairs_cfu, pairs_caseu)

pairs_ID <- pairs_freq %>%
    distinct(Community, Isolate1, Isolate2) %>%
    mutate(Pair = 1:n())

pairs_freq %>%
    filter(Isolate1InitialODFreq != 50) %>%
    filter(Community == "C1R7") %>%
    left_join(pairs_ID) %>%
    mutate(Isolate1InitialODFreq = factor(Isolate1InitialODFreq, c(95,50,5))) %>%
    ggplot(aes(x = Time, y = Isolate1MeasuredFreq, color = Isolate1InitialODFreq, group = Isolate1InitialODFreq)) +
    geom_point() +
    geom_line() +
    scale_y_continuous(breaks = c(0, .5, 1), limits = c(0,1)) +
    #scale_x_discrete(labels = c(0,8)) +
    scale_color_manual(values = frequency_color, label = c("95%", "50%", "5%")) +
    facet_wrap(.~Pair) +
    theme_classic() +
    theme(
        panel.spacing = unit(2, "mm"), strip.text.x = element_blank(),
        panel.border = element_rect(color = 1, fill = NA, size = 1),
        panel.grid.minor.y = element_blank(),
        #axis.title = element_text(size = 10), axis.text = element_text(color = 1, size = 8),
        axis.title = element_blank(), axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_rect(fill = "white"),
        plot.background = element_blank()
    )




# CASEU pairwise result

## Spread the df
pairs_freq_T0 <- pairs_freq %>%
    filter(Time == "T0") %>%
    filter(Isolate1InitialODFreq != 50) %>%
    mutate(Isolate1MeasuredFreqT0 = Isolate1MeasuredFreq,
           ErrorIsolate1MeasuredFreqT0 = ErrorIsolate1MeasuredFreq,
           RawDataTypeT0 = RawDataType) %>%
    select(Community, Isolate1, Isolate2, Isolate1InitialODFreq, Isolate2InitialODFreq,
           Isolate1MeasuredFreqT0, ErrorIsolate1MeasuredFreqT0, RawDataTypeT0)

pairs_freq_T8 <- pairs_freq %>%
    filter(Time == "T7") %>%
    mutate(Isolate1MeasuredFreqT8 = Isolate1MeasuredFreq,
           ErrorIsolate1MeasuredFreqT8 = ErrorIsolate1MeasuredFreq,
           Isolate1MeasuredFreq,
           RawDataTypeT8 = RawDataType) %>%
    select(Community, Isolate1, Isolate2, Isolate1InitialODFreq, Isolate2InitialODFreq,
           Isolate1MeasuredFreqT8, ErrorIsolate1MeasuredFreqT8, RawDataTypeT8)

pairs_freq_spread <- pairs_freq_T0 %>% left_join(pairs_freq_T8) %>%
    # Drop the incomplete pairs
    filter(!is.na(Isolate1MeasuredFreqT0), !is.na(Isolate1MeasuredFreqT8)) %>%
    mutate(DifferenceT8T0 = Isolate1MeasuredFreqT8 - Isolate1MeasuredFreqT0) %>%
    mutate(DifferenceT8T0 = ifelse(DifferenceT8T0 > 0, 1, -1))
# 1 means that the difference is growing from T0 to T8; -1 for decreasing

## Spread the df and paste the frequency changes. Reduce the row number to 186 (total 186 pairs)
pairs_interaction_fitness <- pairs_freq_spread %>%
    select(Community, Isolate1, Isolate2, Isolate1InitialODFreq, DifferenceT8T0) %>%
    group_by(Community, Isolate1, Isolate2) %>%
    pivot_wider(names_from = Isolate1InitialODFreq, values_from = DifferenceT8T0) %>%
    mutate(FreqFunc = paste(`5`, `95`, sep = "_")) %>%
    select(-`5`, -`95`) %>%
    # remove NA
    filter(!str_detect(FreqFunc, "NA"))

# Extract the T8 frequencies ----
## Spread the df and compute the final frequencies. Reduce row number to 186 (total 186 pairs)
pairs_interaction_T8_freq <-
    pairs_freq_spread %>%
    select(Community, Isolate1, Isolate2, Isolate1InitialODFreq, Isolate1MeasuredFreqT8) %>%
    group_by(Community, Isolate1, Isolate2) %>%
    pivot_wider(names_from = Isolate1InitialODFreq, values_from = Isolate1MeasuredFreqT8) %>%
    mutate(PairIsolate1MeasuredFreqT8 = paste(round(`5`, 2), round(`95`, 2), sep = "_"))


### Update the column `Isolate1Win`, which indicates the competitive exclusion that cannot be specified by fitness functions
pairs_interaction_T8_freq[pairs_interaction_T8_freq$PairIsolate1MeasuredFreqT8 == "1_1","Isolate1Win"] <- TRUE
pairs_interaction_T8_freq[pairs_interaction_T8_freq$PairIsolate1MeasuredFreqT8 == "0_0","Isolate1Win"] <- FALSE
pairs_interaction_T8_freq <- pairs_interaction_T8_freq %>%
    select(Community, Isolate1, Isolate2, PairIsolate1MeasuredFreqT8, Isolate1Win) %>%
    arrange(Community, Isolate1, Isolate2, PairIsolate1MeasuredFreqT8, Isolate1Win)

# Match two dfs: frequency changes and T8 frequenies ----
pairs_interaction_fitness <- pairs_interaction_fitness %>%
    left_join(pairs_interaction_T8_freq, by = c("Community", "Isolate1", "Isolate2"))

# Determine the interaction types by fitness functions ----
## Table for determining interaction types
interaction_type <- tibble(
    FromRare = rep(c(1, -1), each = 2),
    FromAbundant = rep(c(1, -1), 2),
    InteractionType = NA,
    InteractionTypeFiner = NA)

## Assign interaction types to combinations of frequency changes signs
interaction_type$InteractionType[c(1,3,4)] <- "exclusion"
interaction_type$InteractionType[c(2)] <- "coexistence"
interaction_type$InteractionTypeFiner[c(1,4)] <- "competitive exclusion"
interaction_type$InteractionTypeFiner[c(2)] <- "coexistence"
interaction_type$InteractionTypeFiner[c(3)] <- "mutual exclusion"
interaction_type <- interaction_type %>% mutate(FreqFunc = paste(FromRare, FromAbundant, sep = "_"))

## Join the fitness and interaction tables
pairs_interaction_fitness <- pairs_interaction_fitness %>%
    left_join(interaction_type, by = "FreqFunc") %>%
    as.data.table()
## Direction of links; from winner to loser
pairs_interaction_fitness[FreqFunc %in% c("1_1"), c("InteractionType", "InteractionTypeFiner", "From", "To") := list("exclusion", "competitive exclusion", Isolate1, Isolate2)]
pairs_interaction_fitness[FreqFunc %in% c("-1_-1"), c("InteractionType", "InteractionTypeFiner", "From", "To") := list("exclusion", "competitive exclusion", Isolate2, Isolate1)]


# Update interaction types by the final frequencies ----
pairs_interaction_fitness[Isolate1Win == TRUE, c("InteractionType", "InteractionTypeFiner", "From", "To") := list("exclusion", "competitive exclusion", Isolate1, Isolate2)]
pairs_interaction_fitness[Isolate1Win == FALSE, c("InteractionType", "InteractionTypeFiner", "From", "To") := list("exclusion", "competitive exclusion", Isolate2, Isolate1)]
# Select for the columns for pairs_interaction
pairs_interaction <- pairs_interaction_fitness %>%
    select(Community, Isolate1, Isolate2, InteractionType, InteractionTypeFiner, From, To)


pairs_interaction_caseu <- pairs_interaction %>%
    as_tibble() %>%
    rename(InteractionTypeCASEU = InteractionType, InteractionTypeFinerCASEU = InteractionTypeFiner) %>%
    select(Community, Isolate1, Isolate2, InteractionTypeCASEU)

pairs %>%
    filter(Assembly == "self_assembly") %>%
    select(Community, Isolate1, Isolate2, InteractionType) %>%
    left_join(pairs_interaction_caseu) %>%
    mutate(ResultMatch = InteractionType == InteractionTypeCASEU) %>%
    group_by(ResultMatch) %>%
    count(ResultMatch)


pairs_interaction_caseu %>% count(InteractionTypeCASEU)
pairs %>% filter(Assembly == "self_assembly") %>% count(InteractionType)





# Shit CASEU and CFU result does not correspond...
# Good news: the current figure is using the approach we discussed: the majority of the results are CFU, with correction, and the morphologically similar strains are filled in by CASEU
# Bad news: caseu and cfu result do not quit match... even quantitatively. There are coexistence and exclusion, but caseu and cfu gives difference results in half of the pairs
# Solution? We do not need to be worried about caseu outcome. Only use the CFU data
# Bad news two? I did not do caseu for 50:50








if (FALSE) {
    pairs_freq %>%
        mutate(Isolate1InitialODFreq = factor(Isolate1InitialODFreq, c(95,50,5))) %>%
        ggplot(aes(x = Time, y = Isolate1MeasuredFreq, color = Isolate1InitialODFreq, group = Isolate1InitialODFreq)) +
        geom_point(size = 2) +
        geom_line(size = 1) +
        scale_y_continuous(breaks = c(0, .5, 1), limits = c(0,1)) +
        scale_x_discrete(labels = c(0,8)) +
        scale_color_manual(values = frequency_color, label = c("95%", "50%", "5%")) +
        facet_wrap(.~Pair) +
        theme_bw() +
        theme(panel.spacing = unit(2, "mm"), strip.text.x = element_blank(),
              panel.border = element_rect(color = 1, fill = NA, size = 1),
              panel.grid.minor.y = element_blank(),
              #axis.title = element_text(size = 10), axis.text = element_text(color = 1, size = 8),
              axis.title = element_blank(), axis.text = element_blank(),
              axis.ticks = element_blank(),
              panel.background = element_rect(fill = "white"),
              plot.background = element_blank()) +
        guides(color = "none") +
        labs(x = "Time (days)", y = "Frequency")

}





if (FALSE) {

    ## Fill in ambiguous pairs with CASEU pilot2 result. Note that if a pair has both CFU and CASEU result, CASEU will overwrite the CFU result
    for (i in 1:nrow(CASEU_pilot2)) {
        index_row <- which(pairs_freq$Time == "T8" &
                               pairs_freq$Community == CASEU_pilot2$Community[i] &
                               pairs_freq$Isolate1 == CASEU_pilot2$Isolate1[i] &
                               pairs_freq$Isolate2 == CASEU_pilot2$Isolate2[i] &
                               pairs_freq$Isolate1InitialODFreq == CASEU_pilot2$Isolate1Freq[i] &
                               pairs_freq$Isolate2InitialODFreq == CASEU_pilot2$Isolate2Freq[i])

        pairs_freq[index_row, c("Isolate1MeasuredFreq", "ErrorIsolate1MeasuredFreq","RawDataType")] <-
            CASEU_pilot2[i, c("Isolate1CFUFreq", "ErrorIsolate1CFUFreq", "RawDataType")]
    }

    ## Fill in ambiguous pairs with CASEU pilot3 result
    for (i in 1:nrow(CASEU_pilot3)) {
        index_row <- which(pairs_freq$Time == "T8" &
                               pairs_freq$Community == CASEU_pilot3$Community[i] &
                               pairs_freq$Isolate1 == CASEU_pilot3$Isolate1[i] &
                               pairs_freq$Isolate2 == CASEU_pilot3$Isolate2[i] &
                               pairs_freq$Isolate1InitialODFreq == CASEU_pilot3$Isolate1Freq[i] &
                               pairs_freq$Isolate2InitialODFreq == CASEU_pilot3$Isolate2Freq[i])

        pairs_freq[index_row, c("Isolate1MeasuredFreq", "ErrorIsolate1MeasuredFreq","RawDataType")] <-
            CASEU_pilot3[i, c("Isolate1CFUFreq", "ErrorIsolate1CFUFreq", "RawDataType")]
    }


    ## Fill in ambiguous pairs with CASEU pilot4 result
    for (i in 1:nrow(CASEU_pilot4)) {
        index_row <- which(pairs_freq$Time == "T8" &
                               pairs_freq$Community == CASEU_pilot4$Community[i] &
                               pairs_freq$Isolate1 == CASEU_pilot4$Isolate1[i] &
                               pairs_freq$Isolate2 == CASEU_pilot4$Isolate2[i] &
                               pairs_freq$Isolate1InitialODFreq == CASEU_pilot4$Isolate1Freq[i] &
                               pairs_freq$Isolate2InitialODFreq == CASEU_pilot4$Isolate2Freq[i])

        pairs_freq[index_row, c("Isolate1MeasuredFreq", "ErrorIsolate1MeasuredFreq","RawDataType")] <-
            CASEU_pilot4[i, c("Isolate1CFUFreq", "ErrorIsolate1CFUFreq", "RawDataType")]
    }
    pairs_freq <- pairs_freq %>% as_tibble()
}


# pairs_freq <- bind_rows(pairs_freq, pairs_freq_random)
# write_csv(pairs_freq, file = here::here("data/output/pairs_freq.csv"))













