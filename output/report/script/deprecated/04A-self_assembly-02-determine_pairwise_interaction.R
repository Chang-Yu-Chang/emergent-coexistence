#' Determine the pairwise interactions based on frequency changes
#' 1. Complete the raw pairs data by joining the full pair list
#' 2. Pivot wider by isolates. Plot the frequency changes across transfers
#' 3. Pivot wider by transfer. Compute the freqence sign changes
#' 4. Pivot wider by initial frequencies. Compute the pairwise interactions

# Read the pairwise competition data; this may take a few seconds
df_pair <- fread(root$find_file("output/report/data/temp/df_pair.txt")) %>% as_tibble()

# Prepare the pairs list ----
#' For transfers not T0, the rows of excluded isolate in a exclusion pair are missing, so
#' we have to add the missing rows back

# Only take T0
df_temp <- df_pair %>%
  filter(Transfer == 0, Type == "consumer") %>%
  select(-Transfer, Abundance)

# Filter for isolate1 initial frequencies
df_temp2 <- df_temp %>% group_by(Assembly, Community, Well, Type) %>% slice(1) %>%
  mutate(Isolate1InitialFreq = Abundance) %>%
  select(-ID, -Abundance)
df_temp <- df_temp %>% left_join(df_temp2, by = c("Assembly", "Community", "Well", "Type")) %>% select(-Abundance)

# Make the pair list for all transfers
pairs_list <- rep(list(df_temp), length(unique(df_pair$Transfer))) %>% setNames(0:(length(.)-1)) %>%
  rbindlist(idcol = "Transfer") %>%
  mutate(Transfer = as.numeric(Transfer)) %>%
  mutate(Isolate1InitialFreq = factor(Isolate1InitialFreq)) %>%
  as_tibble()


# Pivot the df wider by isolate, so that each row is a pair-transfer ----
# Compute isolate 1 frequencies across transfers
df_pair_freq_transfer <- df_pair %>%
  filter(Type == "consumer") %>%
  # For T10, the rows of excluded isolate in a exclusion pair are missing, so
  right_join(pairs_list, by = c("Assembly", "Community", "Well", "Transfer", "Type", "ID")) %>%
  # Replace abundance that are NA into 0. THese isolates are excluded in the competition
  replace_na(list(Abundance = 0)) %>%
  arrange(Community, Well, Transfer) %>%
  # Assign the isolate order within each pair
  mutate(temp = rep(c("Isolate1", "Isolate2"), nrow(.)/2)) %>%
  group_by(Assembly, Community, Transfer, Type) %>%
  # Pivot wider by isolate, so isolate 1 and 2 are in separate columns
  pivot_wider(names_from = temp, values_from = c(ID, Abundance)) %>%
  # Compute the relative abundances
  mutate(temp_column = Abundance_Isolate1 + Abundance_Isolate2,
    Abundance_Isolate1 = Abundance_Isolate1 / temp_column,
    Abundance_Isolate2 = Abundance_Isolate2 / temp_column) %>%
  arrange(Community, Well, Transfer) %>%
  ungroup() %>%
  # Only keep the relative abundance of isolate 1
  mutate(Isolate1 = ID_Isolate1, Isolate2 = ID_Isolate2) %>%
  select(Assembly, Community, Well, Isolate1, Isolate2, Isolate1InitialFreq, Transfer, Abundance_Isolate1)


# Pivot wider by transfer. Compute the frequencies changes in each pair from T0 to T10 ----
df_pair_freq <- df_pair_freq_transfer %>%
  # To determine pairwise interactions, we only need T0 and T10
  filter(Transfer %in% c(0, 10)) %>%
  mutate(Transfer = paste0("T", Transfer)) %>%
  # Pivot wider by transfer
  pivot_wider(names_from = Transfer, values_from = Abundance_Isolate1) %>%
  arrange(Community, Well, T0) %>%
  # Fitness change
  mutate(SignFreq = (T10 - T0) > 0) %>%
  select(-Well) %>%
  arrange(Community, Isolate1, Isolate2)



# Pivot wider by initial frequnencies. Determine the pairwise interaction ----
interaction_type <- tibble(
  `0.05` = c(rep(c(T,F), each = 2), NA),
  `0.95` = c(rep(c(T,F), 2), NA),
  InteractionType = c("exclusion", "coexistence", "bistability", "exclusion", "undefined")
)

df_pair_interaction <- df_pair_freq %>%
  # Only use two initial frequencies to determine the interactions
  filter(T0 %in% c(0.05, 0.95)) %>%
  select(Assembly, Community, Isolate1, Isolate2, Isolate1InitialFreq, SignFreq) %>%
  # Pivot wider by isolate 1's initial frquency
  pivot_wider(names_from = Isolate1InitialFreq, values_from = SignFreq) %>%
  # Determine pairwise interactions by the combination of sign changes
  left_join(interaction_type, by = c("0.05", "0.95")) %>%
  mutate(From = Isolate1, To = Isolate2)

# Switch arrow sign when isolate2 outcompetes isolate1
temp_index <- which(df_pair_interaction$`0.05` == F & df_pair_interaction$`0.95` == F)
df_pair_interaction[temp_index, c("From", "To")] <- df_pair_interaction[temp_index, c("To", "From")]
df_pair_interaction <- df_pair_interaction %>% select(-`0.95`, -`0.05`)




















# Old code ----
if (FALSE) {
  file_list <- list.files(root$find_file("output/report/data/"))

  df_pair <-
    grep(file_list, pattern = "self_assembly-pair", value = T) %>%
    as.list() %>%
    lapply(function(x) {
      community <- x %>% gsub("self_assembly-pair-", "", .) %>% gsub("-T\\d{2}.txt", "", .)

      fread(paste0(root$find_file("output/report/data/"), x)) %>%
        mutate(Community = community) %>%
        return()
    }) %>% rbindlist() %>%
    mutate(ID = factor(ID)) %>%
    as_tibble()

  # Specify the pairs ----
  temp_df <- df_pair %>%
    filter(Transfer == 0, Type == "consumer") %>%
    mutate(Isolate = rep(paste0("Isolate", 1:2), nrow(.)/2))

  temp_df2 <- temp_df[rep(1:nrow(temp_df), max(df_pair$Transfer) + 1),] %>%
    mutate(Transfer = rep(0:max(df_pair$Transfer), each = nrow(temp_df))) %>%
    select(-Abundance)

  ## Pivot wider df_pair
  df_pair_wider <- df_pair %>%
    filter(Type == "consumer") %>%
    right_join(temp_df2) %>%
    pivot_wider(names_from = "Isolate", values_from = c(ID, Abundance)) %>%
    # Replace NA in extinct species
    replace_na(list(Abundance_Isolate1 = 0, Abundance_Isolate2 = 0)) %>%
    # Relative abundance
    mutate(RelativeAbundance_Isolate1 = Abundance_Isolate1 / (Abundance_Isolate1 + Abundance_Isolate2))

  # Plot relative abundances ----
  communities_name <- unique(df_pair$Community) %>% as.character()
  p_pair <- rep(list(NA), length(communities_name))

  ## Relative abundance changes
  for (i in 1:length(p_pair)) {
    p_pair[[i]] <-
      df_pair_wider %>%
      filter(Community == communities_name[i]) %>%
      {if (communities_name[i] != "W14a") filter(., Transfer <= 10) else {.}} %>%
      #  filter(Well %in% paste0("W", 0:200)) %>%
      ggplot(aes(x = Transfer, y = RelativeAbundance_Isolate1, color = Well)) +
      geom_point(size = .5) + geom_line() +
      #    scale_x_continuous(breaks = seq(0,10,2)) +
      scale_y_continuous(breaks = seq(0,1,.5)) +
      guides(color = F) +
      facet_grid(ID_Isolate1 ~ ID_Isolate2) +
      theme_bw() +
      ggtitle(paste0("Self-assembled community ", communities_name[i]))

  }


  # Compute the frequency changes in the frequency-pairs ----
  df_pair_freq_change <- df_pair_wider %>%
    filter(Type == "consumer") %>%
    # Only compare the frequency changes in T0 and T10
    filter(Transfer %in% c(0, 10)) %>%
    select(Community, Well, ID_Isolate1, ID_Isolate2, Transfer, RelativeAbundance_Isolate1)

  ## Pair specific ID
  df_pair_freq_change$Pair <- paste0(df_pair_freq_change$ID_Isolate1, "_", df_pair_freq_change$ID_Isolate2)

  ## Pivot wider the relative abundance of isolate 1 in T0 and T10
  df_pair_freq_change <- df_pair_freq_change  %>%
    select(Community, Well, Pair, Transfer, RelativeAbundance_Isolate1) %>%
    group_by(Community) %>%
    # Pivot wider
    pivot_wider(names_from = Transfer, values_from = RelativeAbundance_Isolate1, names_prefix = "RelativeAbundance_Isolate1_T") %>%
    select(Community, Pair, RelativeAbundance_Isolate1_T0, RelativeAbundance_Isolate1_T10) %>%
    # Sign of frequency changes for each initial frequencies
    mutate(FrequencyChange = (RelativeAbundance_Isolate1_T10 - RelativeAbundance_Isolate1_T0) > 0) %>%
    select(Community, Pair, RelativeAbundance_Isolate1_T0, FrequencyChange) %>%
    # Pivot wider again to make each row a specific pair
    pivot_wider(names_from = RelativeAbundance_Isolate1_T0, values_from = FrequencyChange,
      names_prefix = "InitialFrequency_Isolate1_") %>%
    select(Community, Pair, InitialFrequency_Isolate1_0.05, InitialFrequency_Isolate1_0.5, InitialFrequency_Isolate1_0.95)


  # Determine the pairwise interactions by fitness functions ----

  ## Since there is no frequency changes in the model, there are only 2^3 = 8 combinations of fitness function
  fitness_function_example <- data.frame(
    InitialFrequency_Isolate1_0.05 = rep(c(T,F), each = 4),
    InitialFrequency_Isolate1_0.5 = rep(rep(c(T, F), each = 2), 2),
    InitialFrequency_Isolate1_0.95 = rep(c(T, F), 4),
    InteractionType = c("exclusion", rep("coexistence", 3), "exclusion", "coexistence", "exclusion", "exclusion"),
    InteractionTypeFiner = 1:8,
    From = paste0("Isolate", c(1, 1, 1, 1, 1, 1, 1, 2)),
    To = paste0("Isolate", c(2, 2,2,2,2,2,2, 1)),
    stringsAsFactors = F
  )

  ## Pairwise interactions
  df_pair_interaction <- df_pair_freq_change %>%
    left_join(fitness_function_example) %>%
    select(Community, Pair, InteractionType, From, To) %>%
    {.}



  ## Fraction of pairwise competition
  df_pair_interaction_count <- df_pair_interaction %>%
    # Count each type of pairwise competition
    group_by(Community, InteractionType) %>%
    summarize(Count = n()) %>%
    # Undefined interactions because of the pairs of non-growers
    replace_na(list(InteractionType = "undefined")) %>%
    ungroup() %>%
    # Pivot wider by interaction types
    pivot_wider(names_from = "InteractionType", values_from = Count) %>%
    # Replace 0 count
    replace_na(list(undefined = 0)) %>%
    # Compute the fraction of coexistence
    mutate(FractionCoexistence = coexistence / (coexistence + exclusion + undefined))



  # Plot the interacton types
  ## Fraction of pairwise coexistence
  p_coexistence <- df_pair_interaction_count %>%
    ggplot() +
    geom_bar(aes(x = Community, y = FractionCoexistence), stat = "identity", color = "black") +
    scale_y_continuous(limits = c(0,1)) +
    theme_bw()

  ## Interaction frequency
  p_interaction_frequency <- df_pair_interaction %>%
    # Count each type of pairwise competition
    group_by(Community, InteractionType) %>%
    summarize(Count = n()) %>%
    # Undefined interactions because of the pairs of non-growers
    replace_na(list(InteractionType = "undefined")) %>%
    # Frequency of each interaction Type
    group_by(Community) %>%
    mutate(Frequency = Count / sum(Count)) %>%
    ggplot() +
    geom_bar(aes(x = Community, y = Frequency, fill = InteractionType), stat = "identity", color = "black") +
    #  scale_y_continuous(limits = c(0, 1)) +
    theme_bw() +
    theme(legend.position = "top")

}






