#' Filter the coexistence pairs

community_name <- c("C1R2", "C1R4", "C1R6", "C1R7", "C2R6", "C2R8", "C4R1", "C7R1", "C8R4", "C10R2", "C11R1", "C11R2", "C11R5")


# Pairs that coexist in three initial frequencies ----
# Read data
isolate_ID_match <- fread(root$find_file("data/temp/isolates_ID_match.csv")) %>%
  mutate(Community = ordered(Community, levels = community_name),
         Isolate = ordered(Isolate, 1:12))

isolates_RDP <- fread(root$find_file("data/temp/isolates_RDP.csv"))
pairs_CFU_freq_uncertainty <- fread(root$find_file("data/temp/pairs_CFU_freq_uncertainty.csv")) %>%
  mutate(Community = ordered(Community, levels = community_name),
    Isolate1 = ordered(Isolate1, 1:12), Isolate2 = ordered(Isolate2, 1:12))


# Pairs that coexist at T8 in all three initial frequencies
pairs_coext_all <- pairs_CFU_freq_uncertainty %>%
  arrange(Community, Isolate1, Isolate2, Isolate1Freq, Isolate2Freq, Time) %>%
  filter(Time == "T8") %>%
  filter(!(Isolate1CFUFreq == 1 | Isolate1CFUFreq == 0)) %>%
  group_by(Community, Isolate1, Isolate2) %>%
  summarize(Count = n()) %>%
  filter(Count == 3) %>% # Coexistence in all three initial frequencies
  select(-Count) %>% ungroup() %>%
  #  gather("temp", "Isolate", 2:3) %>%
  #  select(-temp)
  {.}

# Isolates that are in the coexitence pairs
isolates_coext_all <- pairs_coext_all %>%
  gather("temp", "Isolate", 2:3) %>%
  select(-temp) %>%
  left_join(isolate_ID_match, by = c("Community", "Isolate")) %>%
  mutate(Community = ordered(Community, levels = community_name),
    Isolate = ordered(Isolate)) %>%
  arrange(Community, Isolate) %>%
  distinct()

# Taxonomic identity of isolates
isolates_coext_all <- isolates_coext_all %>%
  left_join(isolates_RDP) %>%
  select(Community, Isolate, ExpID, Fermenter, Family, Genus)

# Taxonomic identy of pairs
pairs_coext_all <- pairs_coext_all %>%
  invnet::match_isolates_pairs(isolates_coext_all, variable_name = "Family") %>%
  invnet::match_isolates_pairs(isolates_coext_all, variable_name = "Genus") %>%
  invnet::match_isolates_pairs(isolates_coext_all, variable_name = "Fermenter")




# Pairs that coexist at 50-50 iniitial frequencies ----
isolates_ID_match <- fread(root$find_file("data/temp/isolates_ID_match.csv"))

pairs_coext_50_50 <-
  pairs_CFU_freq_uncertainty %>%
  filter(Time == "T8") %>%
  filter(Isolate1Freq == 50, Isolate1CFUFreq != 0, Isolate1CFUFreq != 1)


# Isolates that are in the coexitence pairs
isolates_coext_50_50 <- pairs_coext_50_50 %>%
  gather("temp", "Isolate", 2:3) %>%
  select(-temp) %>%
  select(Community, Isolate) %>%
  left_join(isolate_ID_match, by = c("Community", "Isolate")) %>%
  mutate(Community = ordered(Community, levels = community_name),
    Isolate = ordered(Isolate)) %>%
  arrange(Community, Isolate) %>%
  distinct() %>%
  left_join(isolates_RDP) %>%
  select(Community, Isolate, ExpID, Fermenter, Family, Genus)

pairs_coext_50_50 <- pairs_coext_50_50 %>%
  invnet::match_isolates_pairs(isolates_coext_50_50, variable_name = "Family") %>%
  invnet::match_isolates_pairs(isolates_coext_50_50, variable_name = "Genus") %>%
  invnet::match_isolates_pairs(isolates_coext_50_50, variable_name = "Fermenter")
{.}







