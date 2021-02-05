#' Compute the 1) competitive scores of isolates and
#' 2) degree of hierarchy of communities

# Read data
pairs_freq <- fread(root$find_file("data/temp/pairs_freq.csv")) %>%
  mutate(Community = ordered(Community, levels = communities_name))


# Competitive scores of isolates----
# Extract frequencies data at T8
temp_df <- pairs_freq %>%
  filter(Time == "T8") %>%
  mutate(Isolate2MeasuredFreq = 1 - Isolate1MeasuredFreq) %>%
  mutate(Freq1 = Isolate1MeasuredFreq, Freq2 = Isolate2MeasuredFreq) %>%
  mutate(InitFreq1 = Isolate1InitialODFreq) %>%
  group_by(Community, InitFreq1)

# Isolates
temp_df1 <- temp_df %>%
  select(Community, InitFreq1, Isolate1, Isolate2) %>%
  pivot_longer(cols = starts_with("Isolate"), names_to = "temp", values_to = "Isolate") %>%
  select(-temp)

# Frequencies
temp_df2 <- temp_df %>%
  select(Community, InitFreq1, Freq1, Freq2) %>%
  pivot_longer(cols = starts_with("Freq"), names_to = "temp", values_to = "Frequency") %>%
  select(-temp)

# Competitive score by matching the isolates and frequencis df
isolates_score <- temp_df1 %>%
  ungroup() %>%
  mutate(Frequency = temp_df2$Frequency) %>%
  group_by(Community, Isolate) %>%
  summarize(Count = n(),
    CompScoreMean = mean(Frequency, na.rm = T),
    CompScoreSd = sd(Frequency, na.rm = T))


# Degree of hierarchy of communities
communities_hierarchy <- pairs_freq %>%
  filter(Time == "T8") %>%
  # When s_i > s_j
  mutate(Freq = ifelse(Isolate1MeasuredFreq > 0.5, Isolate1MeasuredFreq, 1 - Isolate1MeasuredFreq)) %>%
  select(Community, Isolate1, Isolate2, Time, Isolate1InitialODFreq, Isolate1MeasuredFreq, Freq) %>%
  select(Community, Isolate1, Isolate2, Isolate1InitialODFreq, Freq) %>%
  # Sum up the competitive frequencies
  group_by(Community) %>%
  summarize(HierarchyScore = mean(Freq, na.rm = T)) %>%
  as_tibble()











