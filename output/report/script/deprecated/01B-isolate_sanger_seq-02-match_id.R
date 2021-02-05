#' Match 16S Sanger sequences with Djordje's isolate ID

# Reformat the data
isolates_16S <- data.table::fread("data/temp/isolates_16S.csv")


if (FALSE) {
  #' Load Djordje's isolate ID sheet from a csv file community_isolates_201806/seq_isolate_ids_for_R.csv:
  #'  `id` is then turned into `ExpID`
  #'  `seq` is the sanger sequence ID
  #'  `genus` is predicted from RDP


# Load isolate IDs (lab id, sequence and predicted genus) that is saved in a csv file.
# Djordje's isolate IDs
my.ids <- data.table::fread(here::here("data/raw/Sanger/", "seq_isolate_ids_for_R.csv")) %>% filter(grepl("\\.", id))

# Extract community and replicate of the isolates from ids
my.ids <- my.ids %>%
  dplyr::mutate(Community = paste0("C", strsplit(my.ids$id, "\\.") %>% sapply("[", 1) %>% as.numeric(),
                                   "R", strsplit(my.ids$id, "\\.") %>% sapply("[", 2) %>% as.numeric())) %>%
  setNames(c("ExpID", "ID", "Genus", "Community")) %>%
  filter(!is.na(ID), !is.na(Community)) %>%
  dplyr::select(-Genus)

# Subset Djordje's isolates; remove Jean's isolates in these Sanger sequencing batches
community_names <- unique(my.ids$Community)
isolates_RDP <-
  right_join(isolates_16S, my.ids, by = "ID") %>%
  dplyr::select(ExpID, ID, Community, Sequence) %>%
  mutate(Sequence = as.character(Sequence),
    Community = ordered(Community, levels = community_names))

#' There is an isolate called "azo", so it returns `NA` when coerced into `numeric` object.
#' Ignore this one that will be manually corrected
isolates_RDP$Community[isolates_RDP$ExpID == "azo.a1"] <- "C11R5"
}








