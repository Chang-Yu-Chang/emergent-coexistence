# CASEU and colony counting

library(tidyverse)
library(data.table)
library(cowplot)
library(gridExtra)
source(here::here("plotting_scripts/network_functions.R"))

isolates <- read_csv(here::here("data/output/isolates.csv"), col_types = cols())
sequences_abundance <- read_csv(here::here("data/output/sequences_abundance.csv"), col_types = cols())
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
    select(Community, Isolate1, Isolate2, InteractionType, InteractionTypeFiner, From, To) %>%
    mutate(From = ifelse(InteractionType == "coexistence", Isolate1, From)) %>%
    mutate(To = ifelse(InteractionType == "coexistence", Isolate2, To))








#==============================================================================================================================================================================================
# Species abundance
## Panel B. isolate abundance
temp <- sequences_abundance %>%
    filter(AlignmentType == "local") %>%
    filter(AllowMismatch == Inf) %>%
    filter(BasePairMismatch <= 4) %>%
    mutate(Community = ordered(Community,  communities$Community)) %>%
    group_by(Community) %>%
    mutate(RankRelativeAbundance = rank(-RelativeAbundance))

isolates_abundance <- isolates %>%
    filter(Assembly == "self_assembly") %>%
    left_join(temp) %>%
    select(Community, Isolate, Fermenter, RelativeAbundance, RankRelativeAbundance)


#==============================================================================================================================================================================================
# Pairwise outcome by caseu only ----
## NOTE: CASEU has some missing data (n=186-172=14 pairs) because sanger returns short trace data for pairs
pairs_interaction_caseu <- pairs_interaction %>%
    as_tibble() %>%
    #rename(InteractionTypeCASEU = InteractionType, InteractionTypeFinerCASEU = InteractionTypeFiner) %>%
    #select(Community, Isolate1, Isolate2, InteractionTypeCASEU, From, To)
    select(Community, Isolate1, Isolate2, InteractionType, From, To)

pairs_interaction_caseu %>%
    group_by(Community) %>%
    count()

## Hierarchy
isolates_tournament_caseu <- communities %>%
    select(comm = Community, everything()) %>%
    rowwise() %>%
    mutate(pairs_comm = pairs_interaction_caseu %>% filter(Community == comm) %>% list()) %>%
    mutate(tournaments_comm = pairs_comm %>% tournament_rank() %>% list()) %>%
    select(Community = comm, tournaments_comm) %>%
    unnest(cols = tournaments_comm) %>%
    # Join abundance data
    left_join(isolates_abundance)

## Rank versus Abundance
lm1 <- isolates_tournament_caseu %>% glm(RelativeAbundance ~ Rank, data = .)
r_square1 <- with(summary(lm1), 1 - deviance/null.deviance) %>% round(2)
p_value1 <- coef(summary(lm1))[2,4] %>% round(2)

p_caseu1 <- isolates_tournament_caseu %>%
    ggplot() +
    geom_smooth(aes(x = Rank, y = RelativeAbundance), method = "lm", formula = "y~x") +
    geom_point(aes(x = Rank, y = RelativeAbundance), shape = 21, size = 3, stroke = 1, alpha = 0.7) +
    annotate("text", x = 6, y = 0.8, label = paste0("R-square=", r_square1, "\np=", p_value1)) +
    scale_x_continuous(breaks = 1:12) +
    scale_y_continuous(limits = c(0,1)) +
    theme_classic() +
    theme(legend.position = "top", legend.title = element_blank(), panel.border = element_rect(color = 1, fill = NA)) +
    labs(x = "Rank", y = "Relative abundance") +
    ggtitle("CASEU")


## Rank versus ranked abundance
lm2 <- isolates_tournament_caseu %>% glm(RankRelativeAbundance ~ Rank, data = .)
r_square2 <- with(summary(lm2), 1 - deviance/null.deviance) %>% round(3)
p_value2 <- coef(summary(lm2))[2,4] %>% round(3)

p_caseu2 <- isolates_tournament_caseu %>%
    ggplot() +
    geom_smooth(aes(x = Rank, y = RankRelativeAbundance), method = "lm", formula = "y~x") +
    geom_point(aes(x = Rank, y = RankRelativeAbundance), shape = 21, size = 3, stroke = 1, alpha = 0.7, position = position_jitter(height = 0.08, width = 0.08)) +
    annotate("text", x = 2, y = 6, label = paste0("R-square=", r_square2, "\np=", p_value2)) +
    scale_x_continuous(breaks = 1:12) +
    scale_y_continuous(breaks = 1:12) +
    theme_classic() +
    theme(legend.position = "top", legend.title = element_blank(), panel.border = element_rect(color = 1, fill = NA)) +
    labs(x = "Rank", y = "Relative abundance") +
    ggtitle("CASEU")




# Pairwise outcome by cfu ----
pairs_interaction_cfu <- pairs %>%
    filter(Assembly == "self_assembly") %>%
    select(Community, Isolate1, Isolate2, InteractionType, From, To)
pairs_interaction_cfu %>%
    group_by(Community) %>%
    count()


##
isolates_tournament_cfu <- communities %>%
    select(comm = Community, everything()) %>%
    rowwise() %>%
    mutate(pairs_comm = pairs %>% filter(Community == comm) %>% list()) %>%
    mutate(tournaments_comm = pairs_comm %>% tournament_rank() %>% list()) %>%
    select(Community = comm, tournaments_comm) %>%
    unnest(cols = tournaments_comm)


## Hierarchy
isolates_tournament_cfu <- communities %>%
    select(comm = Community, everything()) %>%
    rowwise() %>%
    mutate(pairs_comm = pairs %>% filter(Community == comm) %>% list()) %>%
    mutate(tournaments_comm = pairs_comm %>% tournament_rank() %>% list()) %>%
    select(Community = comm, tournaments_comm) %>%
    unnest(cols = tournaments_comm) %>%
    # Join abundance data
    left_join(isolates_abundance)

## Rank versus Abundance
lm1 <- isolates_tournament_cfu %>% glm(RelativeAbundance ~ Rank, data = .)
r_square1 <- with(summary(lm1), 1 - deviance/null.deviance) %>% round(2)
p_value1 <- coef(summary(lm1))[2,4] %>% round(3)

p_cfu1 <- isolates_tournament_cfu %>%
    ggplot() +
    geom_smooth(aes(x = Rank, y = RelativeAbundance), method = "lm", formula = "y~x") +
    geom_point(aes(x = Rank, y = RelativeAbundance), shape = 21, size = 3, stroke = 1, alpha = 0.7) +
    annotate("text", x = 6, y = 0.8, label = paste0("R-square=", r_square1, "\np=", p_value1)) +
    scale_x_continuous(breaks = 1:12) +
    scale_y_continuous(limits = c(0,1)) +
    theme_classic() +
    theme(legend.position = "top", legend.title = element_blank(), panel.border = element_rect(color = 1, fill = NA)) +
    labs(x = "Rank", y = "Relative abundance") +
    ggtitle("CFU")


## Rank versus ranked abundance
lm2 <- isolates_tournament_cfu %>% glm(RankRelativeAbundance ~ Rank, data = .)
r_square2 <- with(summary(lm2), 1 - deviance/null.deviance) %>% round(2)
p_value2 <- coef(summary(lm2))[2,4] %>% round(3)

p_cfu2 <- isolates_tournament_cfu %>%
    ggplot() +
    geom_smooth(aes(x = Rank, y = RankRelativeAbundance), method = "lm", formula = "y~x") +
    geom_point(aes(x = Rank, y = RankRelativeAbundance), shape = 21, size = 3, stroke = 1, alpha = 0.7, position = position_jitter(height = 0.08, width = 0.08)) +
    annotate("text", x = 2, y = 6.5, label = paste0("R-square=", r_square2, "\np=", p_value2)) +
    scale_x_continuous(breaks = 1:12) +
    scale_y_continuous(breaks = 1:12) +
    theme_classic() +
    theme(legend.position = "top", legend.title = element_blank(), panel.border = element_rect(color = 1, fill = NA)) +
    labs(x = "Rank", y = "Relative abundance") +
    ggtitle("CFU")


p <- plot_grid(p_caseu1, p_caseu2, p_cfu1, p_cfu2, nrow = 2, labels = LETTERS[1:4], axis = "tb", align = "hv")
ggsave(here::here("plots/FigS2.png"), p, width = 10, height = 8)





#==============================================================================================================================================================================================
# Pair CASEU and CFU result matching
pairs %>%
    filter(Assembly == "self_assembly") %>%
    select(Community, Isolate1, Isolate2, InteractionType) %>%
    left_join(pairs_interaction_caseu %>% rename(InteractionTypeCASEU = InteractionType)) %>%
    mutate(ResultMatch = InteractionType == InteractionTypeCASEU) %>%
    group_by(ResultMatch) %>%
    count(ResultMatch)

# Clean the frequency data
temp1 <- pairs_freq %>%
    select(Community, Isolate1, Isolate2, Isolate1InitialODFreq, Time, Isolate1MeasuredFreq, RawDataType) %>%
    mutate(Time = str_replace(Time, "T", "") %>% as.numeric()) %>%
    mutate(Isolate1InitialODFreq = factor(Isolate1InitialODFreq))
temp2 <- temp1 %>%
    filter(Time == 0) %>%
    mutate(RawDataType = "CFU")
temp3 <- temp1 %>%
    filter(Time == 0) %>%
    mutate(RawDataType = "CASEU")
temp4 <- temp1 %>%
    filter(RawDataType == "CFU", Time == 8)
temp5 <- temp1 %>%
    filter(RawDataType == "Sanger", Time == 7) %>%
    mutate(RawDataType = "CASEU", Time = 8)
pairs_freq_cleaned <- bind_rows(temp2, temp3, temp4, temp5)


# Clean the pairwise outcome data
temp <- pairs_interaction_caseu %>%
    select(InteractionType, Community, Isolate1, Isolate2) %>%
    mutate(RawDataType = "CASEU")

pairs_ID <- pairs %>%
    filter(Assembly == "self_assembly") %>%
    select(InteractionType, Community, Isolate1, Isolate2) %>%
    arrange(InteractionType) %>%
    mutate(PairID = factor(1:n())) %>%
    select(-InteractionType)


pairs_example_freq <- pairs %>%
    filter(Assembly == "self_assembly") %>%
    select(InteractionType, Community, Isolate1, Isolate2) %>%
    mutate(RawDataType = "CFU") %>%
    bind_rows(temp) %>%
    mutate(InteractionType = factor(InteractionType, c("exclusion", "coexistence"))) %>%
    arrange(RawDataType, InteractionType) %>%
    right_join(pairs_freq_cleaned) %>%
    left_join(pairs_ID) %>%
    select(PairID, RawDataType, InteractionType, Isolate1InitialODFreq, Time, Isolate1MeasuredFreq, Community, Isolate1, Isolate2)

## Waffle plot for caseu result
plot_example_freq <- function(pairs_freq) {
    # Extract params
    comm <- unique(pairs_freq$Community)
    isolate1 <- unique(pairs_freq$Isolate1)
    isolate2 <- unique(pairs_freq$Isolate2)
    interaction_type_finer <- pairs %>%
        filter(Community == comm, Isolate1 == isolate1, Isolate2 == isolate2) %>%
        pull(InteractionTypeFiner)

    #
    pairs_freq %>%
        mutate(Isolate1InitialODFreq = factor(Isolate1InitialODFreq, c(95,50,5))) %>%
        ggplot() +
        geom_rect(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, aes(fill = InteractionType), alpha = .1) +
        geom_point(size = 2, aes(x = Time, y = Isolate1MeasuredFreq, color = Isolate1InitialODFreq, group = Isolate1InitialODFreq)) +
        geom_line(size = 1, aes(x = Time, y = Isolate1MeasuredFreq, color = Isolate1InitialODFreq, group = Isolate1InitialODFreq)) +
        scale_y_continuous(breaks = c(0, .5, 1), limits = c(0,1)) +
        scale_color_manual(values = frequency_color, label = c("95%", "50%", "5%")) +
        scale_fill_manual(values = assign_interaction_color()) +
        theme_bw() +
        theme(panel.spacing = unit(2, "mm"),
              panel.border = element_rect(color = 1, fill = NA, size = 1),
              panel.grid = element_blank(),
              panel.grid.minor.y = element_blank(),
              axis.title = element_blank(), axis.text = element_blank(),
              axis.ticks = element_blank(),
              plot.background = element_blank(),
              plot.margin = margin(0,0,0,0, "mm")) +
        guides(color = "none", fill = "none") +
        labs(x = "Time", y = "Frequency")
}
temp_list <- pairs_example_freq %>%
    # remove cfu result
    filter(RawDataType == "CASEU") %>%
    arrange(InteractionType, PairID) %>%
    group_split(InteractionType, PairID) %>%
    lapply(plot_example_freq)


## Grid layout
m <- matrix(c(1:186, rep(NA, 4)), nrow = 10)
p_waffle <- arrangeGrob(grobs = temp_list, layout_matrix = m)
p_waffle <- plot_grid(p_waffle, NULL, rel_widths = c(3, 1), scale = c(.9, 1)) + paint_white_background()


make_legend_fill <- function() {
    temp1 <- pairs_interaction_caseu %>%
        group_by(InteractionType) %>%
        count() %>%
        ungroup() %>%
        mutate(Fraction = n/sum(n)) %>%
        mutate(InteractionType = factor(InteractionType, c("exclusion", "coexistence"))) %>%
        arrange(InteractionType)

    temp <- pairs %>%
        ggplot() +
        geom_tile(aes(x = Isolate1, y = Isolate2, fill = InteractionTypeFiner), height = .8, width = .8, alpha = .9) +
        scale_fill_manual(values = assign_interaction_color(),
                          breaks = c("exclusion", "coexistence"),
                          labels = paste0(temp1$InteractionType, " (", round(temp1$Fraction, 3) * 100,"%)")) +
        theme(legend.title = element_blank(),
              legend.position = "right",
              legend.spacing.y = unit("2", "mm"),
              legend.text = element_text(size = 12)
        ) +
        guides(fill = guide_legend(byrow = T)) +
        paint_white_background()
    return(get_legend(temp))

}
p_legend_fill <- make_legend_fill()

p <- p_waffle
ss = .3
p <- ggdraw(p_waffle) +
    draw_plot(p_legend_fill, x = .86, y = .7, width = ss/2, height = ss/2, hjust = 0.5, vjust = .5) +
    #draw_plot(p_legend_color, x = .77, y = .2, width = ss/2, height = ss/2, hjust = 0.5, vjust = .5) +
    theme(panel.background = element_blank(), plot.background = element_rect(color = NA, fill = "white"),
          plot.margin = unit(c(0,0,0,0), "mm"))

ggsave(here::here("plots/FigS3.png"), p, width = 13, height = 5)


# Plot side-by-side comparison
plot_sidebyside_freq <- function(tb) {
    # Update outcome data. If the T8 result is missing, set it to NA
    temp4 <- tb %>%
        filter(Time == 8) %>%
        select(RawDataType, InteractionType, Isolate1MeasuredFreq) %>%
        mutate(InteractionType = factor(InteractionType, c("exclusion", "coexistence", "NA"))) %>%
        mutate(RawDataType = factor(RawDataType, c("CFU", "CASEU")))
    temp4$InteractionType[is.na(temp4$Isolate1MeasuredFreq)] <- "NA"
    temp4$InteractionType[is.na(temp4$InteractionType)] <- "NA"
    temp4 <- distinct(temp4, RawDataType, InteractionType)


    tb %>%
        mutate(RawDataType = factor(RawDataType, c("CFU", "CASEU"))) %>%
        ggplot() +
        geom_rect(data = temp4, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, aes(fill = InteractionType, linetype = RawDataType), color = 1, alpha = .4, size = 1) +
        geom_point(size = 2, aes(x = Time, y = Isolate1MeasuredFreq, color = Isolate1InitialODFreq, group = Isolate1InitialODFreq)) +
        geom_line(size = 1, aes(x = Time, y = Isolate1MeasuredFreq, color = Isolate1InitialODFreq, group = Isolate1InitialODFreq)) +
        scale_y_continuous(breaks = c(0, .5, 1), limits = c(0,1)) +
        scale_color_manual(values = frequency_color, label = c("95%", "50%", "5%")) +
        scale_fill_manual(values = c(assign_interaction_color(), "NA" = "grey90")) +
        scale_linetype_manual(values = c("CFU" = 1, "CASEU" = 2), name = "") +
        facet_grid(.~RawDataType) +
        theme_bw() +
        theme(panel.spacing = unit(0, "mm"),
              strip.background = element_blank(),
              strip.text = element_blank(),
              #panel.border = element_rect(color = 1, fill = NA, size = 1),
              panel.border = element_blank(),
              panel.grid = element_blank(),
              panel.grid.minor.y = element_blank(),
              axis.title = element_blank(), axis.text = element_blank(),
              axis.ticks = element_blank(),
              plot.background = element_blank(),
              plot.margin = margin(0,0,0,0, "mm")) +
        guides(color = "none", fill = "none", linetype = "none") +
        labs(x = "Time", y = "Frequency")
}
temp_list <- pairs_example_freq %>%
    #filter(PairID %in% 1:4) %>%
    as_tibble() %>%
    arrange(PairID) %>%
    group_split(PairID) %>%
    lapply(plot_sidebyside_freq)


tb <- pairs_example_freq %>%
    filter(PairID %in% 53)

## Grid layout
m <- matrix(c(1:186, rep(NA, 4)), nrow = 10)
#p_waffle <- arrangeGrob(grobs = temp_list, layout_matrix = m)
p_waffle <- plot_grid(p_waffle, NULL, rel_widths = c(5, 1), scale = c(.9, 1)) + paint_white_background()


make_legend_fill_simple <- function() {
    temp1 <- pairs_interaction_caseu %>%
        group_by(InteractionType) %>%
        count() %>%
        ungroup() %>%
        mutate(Fraction = n/sum(n)) %>%
        mutate(InteractionType = factor(InteractionType, c("exclusion", "coexistence"))) %>%
        arrange(InteractionType)

    temp <- pairs %>%
        ggplot() +
        geom_tile(aes(x = Isolate1, y = Isolate2, fill = InteractionTypeFiner), height = .8, width = .8, alpha = .9) +
        scale_fill_manual(values = assign_interaction_color(),
                          breaks = c("exclusion", "coexistence"),
                          labels = paste0(temp1$InteractionType)) +
        theme(legend.title = element_blank(),
              legend.position = "right",
              legend.spacing.y = unit("2", "mm"),
              legend.text = element_text(size = 12)
        ) +
        guides(fill = guide_legend(byrow = T)) +
        paint_white_background()
    return(get_legend(temp))

}
p_legend_fill <- make_legend_fill_simple()
make_legend_linetype <- function() {
    tb <- pairs_example_freq %>%
        filter(PairID %in% 53)

    # Update outcome data. If the T8 result is missing, set it to NA
    temp4 <- tb %>%
        filter(Time == 8) %>%
        select(RawDataType, InteractionType, Isolate1MeasuredFreq) %>%
        mutate(InteractionType = factor(InteractionType, c("exclusion", "coexistence", "NA"))) %>%
        mutate(RawDataType = factor(RawDataType, c("CFU", "CASEU")))
    temp4$InteractionType[is.na(temp4$Isolate1MeasuredFreq)] <- "NA"
    temp4$InteractionType[is.na(temp4$InteractionType)] <- "NA"
    temp4 <- distinct(temp4, RawDataType, InteractionType)


    temp <- tb %>%
        mutate(RawDataType = factor(RawDataType, c("CFU", "CASEU"))) %>%
        ggplot() +
        geom_rect(data = temp4, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, aes(fill = InteractionType, linetype = RawDataType), color = 1, alpha = .4, size = 1) +
        geom_point(size = 2, aes(x = Time, y = Isolate1MeasuredFreq, color = Isolate1InitialODFreq, group = Isolate1InitialODFreq)) +
        geom_line(size = 1, aes(x = Time, y = Isolate1MeasuredFreq, color = Isolate1InitialODFreq, group = Isolate1InitialODFreq)) +
        scale_y_continuous(breaks = c(0, .5, 1), limits = c(0,1)) +
        scale_color_manual(values = frequency_color, label = c("95%", "50%", "5%")) +
        scale_fill_manual(values = c(assign_interaction_color(), "NA" = "grey90")) +
        scale_linetype_manual(values = c("CFU" = 1, "CASEU" = 2), name = "") +
        facet_grid(.~RawDataType) +
        theme_bw() +
        theme(panel.spacing = unit(0, "mm"),
              strip.background = element_blank(),
              strip.text = element_blank(),
              #panel.border = element_rect(color = 1, fill = NA, size = 1),
              panel.border = element_blank(),
              panel.grid = element_blank(),
              panel.grid.minor.y = element_blank(),
              axis.title = element_blank(), axis.text = element_blank(),
              axis.ticks = element_blank(),
              plot.background = element_blank(),
              plot.margin = margin(0,0,0,0, "mm")) +
        guides(color = "none", fill = "none") +
        labs(x = "Time", y = "Frequency")
    return(get_legend(temp))

}
p_legend_linetype <- make_legend_linetype()


ss = .3
p <- ggdraw(p_waffle) +
    draw_plot(p_legend_fill, x = .86, y = .7, width = ss/2, height = ss/2, hjust = 0.5, vjust = .5) +
    draw_plot(p_legend_linetype, x = .86, y = .2, width = ss/2, height = ss/2, hjust = 0.5, vjust = .5) +
    #draw_plot(p_legend_color, x = .77, y = .2, width = ss/2, height = ss/2, hjust = 0.5, vjust = .5) +
    theme(panel.background = element_blank(), plot.background = element_rect(color = NA, fill = "white"),
          plot.margin = unit(c(0,0,0,0), "mm"))

ggsave(here::here("plots/FigS4.png"), p, width = 15, height = 5)





pairs %>%
    filter(Assembly == "self_assembly") %>%
    select(InteractionType, Community, Isolate1, Isolate2) %>%
    mutate(RawDataType = "CFU") %>%
    bind_rows(temp) %>%
    left_join(pairs_ID) %>%
    select(PairID, InteractionType, RawDataType) %>%
    pivot_wider(names_from = RawDataType, values_from = InteractionType) %>%
    mutate(ResultMatch = CFU == CASEU) %>%
    group_by(ResultMatch) %>%
    count()



##
pairs_interaction_caseu %>% count(InteractionType)
pairs %>% filter(Assembly == "self_assembly") %>% count(InteractionType)

















