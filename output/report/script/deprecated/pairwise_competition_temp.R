# Align the colony count data

result <- fread(root$find_file("data/raw/pairwise_competition/result_pairwise_competition_keyin.csv"))
pairs_compt <- fread(root$find_file("data/raw/pairwise_competition/pairs1_compt_manual_check.csv"))
pairs_coext <- fread(root$find_file("data/raw/pairwise_competition/pairs1_coext_manual_check.csv"))

names(result)
names(pairs_compt)

pairs <- rbind(select(pairs_coext, Community, Isolate1, Isolate2, Isolate1Freq, Isolate2Freq, Experiment, ColonyCount, ColonyCount1, ColonyCount2),
  select(pairs_compt, Community, Isolate1, Isolate2, Isolate1Freq, Isolate2Freq, Experiment, ColonyCount, ColonyCount1, ColonyCount2))

result %>%
  select(-starts_with("ColonyCount")) %>%
  left_join(pairs) %>%
  fwrite(file = root$find_file("data/raw/pairwise_competition/result_pairwise_competition.csv"))
