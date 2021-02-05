#' OD-CFU conversion for pairwise competition
#' This script calculate the mean of conversion


# Calculate the mean of measurement ----
## Step1: Calculate epsilon for monoculture from T8 ----
# Isolates
isolates_epsilon <- fread(here::here("data/temp/isolates_epsilon.csv"))

# Pairs
pairs_competition <- fread(here::here("data/temp/pairs_competition.csv")) %>%
  mutate(Isolate1 = ordered(Isolate1, 1:12), Isolate2 = ordered(Isolate2, 1:12))%>%
  mutate(ColonyCount2 = ColonyCount - ColonyCount1)

## Step2: Calculate cell frequency at T0 from OD ----
match_isolates_pairs <- function (pairs_df, isolates_df, variable_name = "Epsilon") {
  if (!("Isolate" %in% names(isolates_df))) {
    stop("Isolate df does not have a column of Isolate")
  } else if (!("Isolate1" %in% names(pairs_df)) | !("Isolate1" %in% names(pairs_df))) {
    stop("Pair df does not have columns of Isolate1 and Isolate2")
  }

  # Congifurate isolate df
  temp <- isolates_df %>%
    mutate(Isolate1 = ordered(Isolate, levels = 1:12), Isolate2 = ordered(Isolate, levels = 1:12)) %>%
    data.table::as.data.table()
  temp[,paste0("Isolate", 1:2)] <- temp[,c("Isolate", "Isolate")]
  temp[,paste0(variable_name, 1:2)] <- temp[,rep(variable_name, 2), with = FALSE] # temp is a data.table object
  temp <- as_tibble(temp) %>% mutate(Isolate1 = ordered(Isolate1, levels = 1:12), Isolate2 = ordered(Isolate2, levels = 1:12))

  # Match isolate 1 and isolate 2 by join isolates df to pairs df
  pairs_df %>%
    mutate(Isolate1 = ordered(Isolate1, levels = 1:12), Isolate2 = ordered(Isolate2, levels = 1:12)) %>%
    left_join(temp[,c("Community", "Isolate1", paste0(variable_name, 1))], by = c("Community", "Isolate1")) %>%
    left_join(temp[,c("Community", "Isolate2", paste0(variable_name, 2))], by = c("Community", "Isolate2")) %>%
    invisible() %>%
    return()
}

# Match epsilon from isolates table to pairs table. Use invnet function for it
pairs_eplison <- match_isolates_pairs(pairs_competition, isolates_epsilon, variable_name = "Epsilon")

# Cell frequency. Assume that the the inocula at T0 was standardized to OD = 0.1.
pairs_CFU_freq_T0 <- pairs_eplison %>%
  mutate(Isolate1CFUFreqT0 = (Isolate1Freq * .1 * Epsilon1) / ((Isolate1Freq * .1 * Epsilon1) + (Isolate2Freq * .1 * Epsilon2)))

## Step3: CFU relative fraction at T8 ----
# Calculate the cell frequency at T8. The data is from CFU plating.
pairs_CFU_freq_T8 <- pairs_eplison %>%
  select(Community, Isolate1, Isolate2, Isolate1Freq, Isolate2Freq, ColonyCount, ColonyCount1, ColonyCount2) %>%
  mutate(Isolate1CFUFreqT8 = ColonyCount1 / (ColonyCount1+ColonyCount2))

# Merge T0 and T8 CFU frequencies ----
pairs_CFU_freq <- pairs_CFU_freq_T0 %>%
  left_join(pairs_CFU_freq_T8, by = c("Community", "Isolate1", "Isolate2", "Isolate1Freq", "Isolate2Freq")) %>%
  select(Community, Isolate1, Isolate2, Isolate1Freq, Isolate2Freq, Isolate1CFUFreqT0, Isolate1CFUFreqT8) %>%
  gather("Time", "Isolate1CFUFreq", 6:7)


fwrite(pairs_CFU_freq, file = here::here("data/temp/pairs_CFU_freq.csv"))
