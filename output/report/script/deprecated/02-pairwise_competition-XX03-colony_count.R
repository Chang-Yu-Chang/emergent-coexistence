#' Read the single isolate colony counts at T8

# Read colony counts at T8
colonies <- fread(root$find_file("data/raw/pairwise_competition/colony_count_dilution_factor.csv"))
colonies$ColonyCount[colonies$ColonyCount <= 0] <- 0

# Remove monoculture or pairs that have two dilution factors (two rows).
colonies <- distinct(colonies, Batch, Community, Isolate1, Isolate2, Isolate1Freq, .keep_all = T)

# Remove monoculture of community C11R1 at batch B2.
colonies <- filter(colonies, !(Community == "C11R1" & Batch == "B2" & Isolate1Freq == "monoculture"))
