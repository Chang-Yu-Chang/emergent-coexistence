#' Import the data of pairwise competition outcomes



# Name community
community_name <- c("C1R2", "C1R4", "C1R6", "C1R7", "C2R6", "C2R8", "C4R1", "C7R1", "C8R4", "C10R2", "C11R1", "C11R2", "C11R5")
community_size <- c(4, 5, 5, 7, 4, 4, 3, 4, 3, 3, 9, 13, 5)

# Read isolate data
isolates2 <- fread("cleaned_data/isolates2.csv") %>%
  select(-Rating, -Games, -Win, -Draw, -Loss) %>% # The competitive results of isolates are updated in the Isolates section.
  mutate(Community=factor(Community, levels=community_name))

# Read pairs data
pairs2_cell <- fread("cleaned_data/pairs2_cell.csv")
pairs2 <- fread("cleaned_data/pairs2_cell_interaction.csv") %>%
  mutate(Community = factor(Community, community_name)) %>%
  select(Community, From, To, InteractionType) %>%
  arrange(Community)

# For pairs3, Consider bistability (coexistence and exclusion) as coexistence
pairs3 <- pairs2
pairs3$InteractionType[pairs3$InteractionType %in% c("bistability (coexistence and exclusion)", "stable coexistence")] <- "coexistence"


# Read the colony counting results
pairs1_coext <- fread("cleaned_data/pairs1_coext_manual_check.csv") %>% filter(CoexistenceOD==1 | is.na(CoexistenceOD)) %>% select(-Hard, -DilutionFactor, -CoexistenceOD)
pairs1_compt <- fread("cleaned_data/pairs1_compt_manual_check.csv")
