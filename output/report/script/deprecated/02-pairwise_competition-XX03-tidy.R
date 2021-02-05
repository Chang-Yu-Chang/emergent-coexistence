#' pairwise-competition_03-tidy.R
#' Description: this script tidy up the pairwise and isolate data for later analyses and visualizations
#' Input objects: pairs1_coext, pairs1_compt, pairs3, isolates2
#' Output objects: pairs3_ferm, pairs3_ferm_ratio, isolates3, pairs1_coext_ferm1, pairs1_coext_ferm2

fermenter_family <- c("Enterobacteriaceae", "Aeromonadaceae")
non_fermenter_family <- c("Pseudomonadaceae", "Moraxellaceae", "Xanthomonadaceae", "Alcaligenaceae", "Comamonadaceae")

# Pairwise competitive outcomes and isolate fermenter identities----

# Need to compare the pairwise competitive outcomes with fermenter identity. Match the outcomes of pairwise competition to the fermenter and non-fermenter identities of isolates
pairs3_ferm <- 
  left_join(pairs3, mutate(isolates2, From=Isolate, Fermenter1=Fermenter) %>% select(Community, From, Fermenter1)) %>%
  left_join(mutate(isolates2, To=Isolate, Fermenter2=Fermenter) %>% select(Community, To, Fermenter2)) %>%
  filter(!is.na(Fermenter1), !is.na(Fermenter2))

# Pair types classified by fermenter or non-fermenter
pairs3_ferm$PairType <-
  apply(pairs3_ferm[,c("Fermenter1", "Fermenter2")], 1, function(x) sort(c(substr(x[1], 1, 1), substr(x[2], 1, 1))) %>% as.vector()) %>%
  apply(2, function(x) paste0(x[1], x[2]))

# Need another df for plotting the ratio with barplot. Calculate coexistence/exclusion ratio in fermenter and non-fermenter pairs
pairs3_ferm_ratio <- pairs3_ferm %>%
  filter(InteractionType %in% c("exclusion", "coexistence")) %>%
  group_by(PairType, InteractionType) %>%
  summarize(Count=n()) %>%
  group_by(PairType) %>%
  spread(InteractionType, Count) %>%
  mutate(Ratio=coexistence/exclusion)


# Pairwise competitive outcomes and update the isolate win/lose/draw ----

# Remove undefined pairs
pairs3_defined <- filter(pairs3, InteractionType%in%c("exclusion", "coexistence"))

# Calculate the isolate win/loss/draw table by using PlayerRatings::elo() 
df_list <- rep(list(NA), 13)

for(i in 1:length(community_name)) {
  df_temp <- filter(pairs3_defined, Community==community_name[i]) %>%
    mutate(Period=1, InteractionType=ifelse(InteractionType=="exclusion", 1, 1/2)) %>%
    select(Period, From, To, InteractionType) %>%
    PlayerRatings::elo() %>% `[[`(1) %>% # Calulate the table
    arrange(Player) %>% mutate(Isolate=Player) %>%
    select(Isolate, Rating, Games, Win, Draw) 
  # Normalize Competitive rank
  df_temp$Rating <- (df_temp$Rating - min(df_temp$Rating)) /(max(df_temp$Rating)-min(df_temp$Rating))
  
  df_list[[i]] <- cbind(Community=community_name[i], df_temp)
}

isolates3 <- do.call(rbind, df_list) 
isolates3 <- left_join(isolates2, isolates3, by=c("Community", "Isolate"))


# Write the updated isolate table 
if (write_file) write.csv(isolates3, file = "cleaned_data/isolates3.csv", row.names = F)
# Colony counts in pairwise competitive outcomes and isolate fermenter identities ----

# Match fermetation identity to the pairs
pairs1_coext_ferm <- 
  left_join(pairs1_coext, mutate(isolates2, From=Isolate, Fermenter1=Fermenter) %>% select(Community, From, Fermenter1)) %>%
  left_join(mutate(isolates2, To=Isolate, Fermenter2=Fermenter) %>% select(Community, To, Fermenter2)) %>%
  arrange(Community, Isolate1, Isolate2, Isolate1Freq, Isolate2Freq)

pairs1_coext_ferm <- filter(pairs1_coext_ferm, !is.na(Fermenter1), !is.na(Fermenter2))
pairs1_coext_ferm$PairType <- NULL
pairs1_coext_ferm$PairType[(pairs1_coext_ferm$Fermenter1=="Fermenter" & pairs1_coext_ferm$Fermenter2=="Fermenter")] <- "FF"
pairs1_coext_ferm$PairType[(pairs1_coext_ferm$Fermenter1=="Non-fermenter" & pairs1_coext_ferm$Fermenter2=="Non-fermenter")] <- "NN"
pairs1_coext_ferm$PairType[(pairs1_coext_ferm$Fermenter1=="Fermenter" & pairs1_coext_ferm$Fermenter2=="Non-fermenter") |
                             (pairs1_coext_ferm$Fermenter1=="Non-fermenter" & pairs1_coext_ferm$Fermenter2=="Fermenter")] <- "FN"


# Match only the stable coexistence. Each count is a pair-frequency
## Filter colony count df for only FN pairs
pairs1_coext_ferm1 <-  filter(pairs1_coext_ferm, PairType=="FN") %>%
  select(Community, Isolate1, Isolate2, Isolate1Freq, Isolate2Freq, Fermenter1, Fermenter2, ColonyCount1, ColonyCount2, ColonyCount, PairType)

## Filter competitive outcome df for only the stable coexistence
pairs3_stable_coexistence <- mutate(pairs2, Isolate1=From, Isolate2=To, Note=T) %>%
  filter(InteractionType%in%c("stable coexistence", "bistability (coexistence and exclusion)")) %>% 
  select(Community, Isolate1, Isolate2, Note)

## Match competitive outcomes to colony count
pairs1_coext_ferm1 <- pairs1_coext_ferm1 %>% left_join(pairs3_stable_coexistence, by = c("Community", "Isolate1", "Isolate2")) %>% filter(!is.na(Note))

## Need isolate1 in the pair always be fermenter. Switch the columns for isolate1 and 2
pairs1_coext_ferm1[pairs1_coext_ferm1$Fermenter1=="Non-fermenter", c("Isolate1", "Isolate2", "Isolate1Freq", "Isolate2Freq", "Fermenter1", "Fermenter2", "ColonyCount1", "ColonyCount2")] <- pairs1_coext_ferm1[pairs1_coext_ferm1$Fermenter1=="Non-fermenter", c("Isolate2", "Isolate1", "Isolate2Freq", "Isolate1Freq", "Fermenter2", "Fermenter1", "ColonyCount2", "ColonyCount1")]

## Filter for non-zero colony count and calculate the log ratios of fermenters over non-fermenters
pairs1_coext_ferm1 <- pairs1_coext_ferm1 %>%
  filter(ColonyCount1!=0, ColonyCount2!=0) %>% # Remove the 0 colony in bistability pairs
  select(-PairType, -Note) %>%
  mutate(logRatio=log(ColonyCount2/ColonyCount1, 10))

## Save df
if (write_file) write.csv(pairs1_coext_ferm1, file="cleaned_data/pairs_NF_ratio.csv", row.names = F)


# Match all coexistence (stable and bistabiliy)
## Filter colony count df for only FN pairs
pairs1_coext_ferm2 <-  filter(pairs1_coext_ferm, PairType=="FN") %>%
  select(Community, Isolate1, Isolate2, Isolate1Freq, Isolate2Freq, Fermenter1, Fermenter2, ColonyCount1, ColonyCount2, ColonyCount, PairType)

## Need isolate1 in the pair always be fermenter. Switch the columns for isolate1 and 2
pairs1_coext_ferm2[pairs1_coext_ferm2$Fermenter1=="Non-fermenter", c("Isolate1", "Isolate2", "Isolate1Freq", "Isolate2Freq", "Fermenter1", "Fermenter2", "ColonyCount1", "ColonyCount2")] <- pairs1_coext_ferm2[pairs1_coext_ferm2$Fermenter1=="Non-fermenter", c("Isolate2", "Isolate1", "Isolate2Freq", "Isolate1Freq", "Fermenter2", "Fermenter1", "ColonyCount2", "ColonyCount1")]

## Filter for non-zero colony count and calculate the log ratios of fermenters over non-fermenters
pairs1_coext_ferm2 <- pairs1_coext_ferm2 %>%
  filter(ColonyCount1!=0, ColonyCount2!=0) %>% # Remove the 0 colony in bistability pairs
  select(-PairType) %>%
  mutate(logRatio=log(ColonyCount2/ColonyCount1, 10))




# Pairwise competitive outcomes by family ----

# Need to compare the pairwise competitive outcomes with fermenter identity. Match the outcomes of pairwise competition to the fermenter and non-fermenter identities of isolates
pairs3_family <- 
  left_join(pairs3, mutate(isolates2, From=Isolate, Family1=Family) %>% select(Community, From, Family1)) %>%
  left_join(mutate(isolates2, To=Isolate, Family2=Family) %>% select(Community, To, Family2)) %>%
  filter(!is.na(Family1), !is.na(Family2))

# Pair types classified by family
pairs3_family$PairTypeFamily <-
  apply(pairs3_family[,c("Family1", "Family2")], 1, function(x) sort(c(substr(x[1], 1, 3), substr(x[2], 1, 3))) %>% as.vector()) %>%
  apply(2, function(x) paste0(x[1], x[2]))

# Need another df for plotting the ratio with barplot. Calculate coexistence/exclusion ratio in fermenter and non-fermenter pairs
pairs3_family_ratio <- pairs3_family %>%
  filter(InteractionType %in% c("exclusion", "coexistence")) %>%
  group_by(PairTypeFamily, InteractionType) %>%
  summarize(Count=n()) %>%
  group_by(PairTypeFamily) %>%
  spread(InteractionType, Count) %>%
  mutate(Ratio=coexistence/exclusion)

# Filter for only two fermenter families
#filter(pairs3_family_ratio, grepl("Ent|Aer", PairTypeFamily))

# Pairwise competitive outcomes by genus ----

# Need to compare the pairwise competitive outcomes with fermenter identity. Match the outcomes of pairwise competition to the fermenter and non-fermenter identities of isolates
pairs3_genus <- 
  left_join(pairs3, mutate(isolates2, From=Isolate, Genus1=Genus) %>% select(Community, From, Genus1)) %>%
  left_join(mutate(isolates2, To=Isolate, Genus2=Genus) %>% select(Community, To, Genus2)) %>%
  filter(!is.na(Genus1), !is.na(Genus2))

# Pair types classified by Genus
pairs3_genus$PairTypeGenus <-
  apply(pairs3_genus[,c("Genus1", "Genus2")], 1, function(x) sort(c(substr(x[1], 1, 3), substr(x[2], 1, 3))) %>% as.vector()) %>%
  apply(2, function(x) paste0(x[1], x[2]))


# Need another df for plotting the ratio with barplot. Calculate coexistence/exclusion ratio in fermenter and non-fermenter pairs
pairs3_genus_ratio <- pairs3_genus %>%
  filter(InteractionType %in% c("exclusion", "coexistence")) %>%
  group_by(PairTypeGenus, InteractionType) %>%
  summarize(Count=n()) %>%
  group_by(PairTypeGenus) %>%
  spread(InteractionType, Count) %>%
  mutate(Ratio=coexistence/exclusion)
