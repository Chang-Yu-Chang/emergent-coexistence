# Clean up the simulation data
library(tidyverse)
library(tidygraph)
library(igraph)
library(ggraph)
library(cowplot)
source(here::here("plotting_scripts/network_functions.R"))

output_dir <- "~/Dropbox/lab/invasion-network/simulation/data/raw9/"
input_independent <- read_csv(paste0(output_dir,"input_independent.csv"), col_types = cols())
input_pairs <- read_csv(paste0(output_dir, "input_pairs.csv"), col_types = cols())
input_row <- input_independent[1,]
sa <- input_independent$sa[1]
ma <- input_independent$ma[1]
sal <- tibble(Family = paste0("F", c(rep(0, sa), rep(1, sa))), Species = paste0("S", 0:(sa * 2 - 1)))
mal <- tibble(Class = paste0("T", c(rep(0, ma), rep(1, ma), rep(2, ma))), Resource = paste0("R", 0:(ma * 3 - 1)))
read_wide_file <- function(x) {
    tt <- read_csv(x, col_types = cols()) %>%
        # Remove abundance=0
        mutate_all(~replace(., .==0, NA)) %>%
        pivot_longer(cols = starts_with("W"), names_to = "Well", values_to = "Abundance", values_drop_na = T)
    if ("...1" %in% colnames(tt)) tt <- tt %>% rename(Family = ...1, Species = ...2)
    return(tt)
}
paint_white_background <- function(x) theme(plot.background = element_rect(color = NA, fill = "white"))

b = 1000


"
ran 1000 permutations and leave lab. Check the result
"


# Community abundance ----
df_comm_init <- paste0(output_dir, "selfAssembly-1_init.csv") %>%
    read_csv(col_types = cols()) %>%
    mutate_all(~replace(., .==0, NA)) %>%
    pivot_longer(cols = starts_with("W"), names_to = "Well", values_to = "Abundance", values_drop_na = T) %>%
    mutate(Transfer = 0, Time = 0) %>%
    mutate(Family = factor(Family, c("F0", "F1"))) %>%
    mutate(Species = factor(Species, paste0("S", 0:999))) %>%
    mutate(Community = factor(Well, paste0("W", 0:999)), .keep = "unused")
df_comm_timepoint <- list.files(output_dir) %>%
    str_subset("_T\\d+t\\d+") %>%
    str_subset("selfAssembly") %>%
    lapply(function(x) {
        temp <- str_replace(x, "selfAssembly-1_", "") %>%
            str_replace(".csv", "") %>%
            str_split("t") %>%
            unlist()
        paste0(output_dir, x) %>%
            read_csv(col_types = cols()) %>%
            mutate_all(~replace(., .==0, NA)) %>%
            pivot_longer(cols = starts_with("W"), names_to = "Well", values_to = "Abundance", values_drop_na = T) %>%
            rename(Family = ...1, Species = ...2) %>%
            mutate(Transfer = str_replace(temp[1], "T", ""), Time = temp[2]) %>%
            mutate(across(c("Transfer", "Time"), as.numeric))
    }) %>%
    bind_rows() %>%
    mutate(Family = factor(Family, c("F0", "F1"))) %>%
    mutate(Species = factor(Species, paste0("S", 0:999))) %>%
    mutate(Community = factor(Well, paste0("W", 0:999)), .keep = "unused")

df_communities_abundance <- bind_rows(df_comm_init, df_comm_timepoint) %>%
    rename(ID = Species)
df_communities <- df_communities_abundance %>%
    filter(Transfer == max(Transfer), Time == max(Time)) %>%
    arrange(Community) %>%
    group_by(Community) %>%
    summarize(Richness = n())

write_csv(df_communities, here::here("data/output/df_communities.csv"))
write_csv(df_communities_abundance, here::here("data/output/df_communities_abundance.csv"))

# Pairwise outcome of community pairs ----
determine_interaction <- function(pairs_init, pairs_end) {
    temp <- bind_rows(pairs_init, pairs_end) %>%
        select(-Family, -Abundance) %>%
        pivot_wider(names_from = Time, values_from = RelativeAbundance, names_prefix = "RelativeAbundance_") %>%
        # Fill abundance = NA with 0
        replace_na(list(RelativeAbundance_Tend = 0)) %>%
        # Frequency changes
        mutate(FrequencyChange = ifelse(RelativeAbundance_Tend - RelativeAbundance_Tinit > 0, "increase", "decrease")) %>%
        select(-RelativeAbundance_Tinit) %>%
        group_by(Community, Well) %>%
        mutate(Isolate = c(1,2)) %>%
        pivot_wider(names_from = Isolate, values_from = c(Species, FrequencyChange, RelativeAbundance_Tend), names_sep = "") %>%
        ungroup() %>%
        select(-FrequencyChange2, -RelativeAbundance_Tend2, -Well) %>%
        # Frequency changes in each pair
        mutate(Pair = rep(paste0("P", 1:(n()/2)), each = 2), Replicate = rep(1:2, n()/2)) %>%
        pivot_wider(names_from = Replicate, values_from = c(FrequencyChange1, RelativeAbundance_Tend1), names_prefix = "Replicate") %>%
        # Interactions
        mutate(Outcome = with(., case_when(
            (FrequencyChange1_Replicate1 == "increase" & FrequencyChange1_Replicate2 == "increase") ~ "win",
            (FrequencyChange1_Replicate1 == "decrease" & FrequencyChange1_Replicate2 == "decrease") ~ "lose",
            (FrequencyChange1_Replicate1 == "increase" & FrequencyChange1_Replicate2 == "decrease" & RelativeAbundance_Tend1_Replicate1 > 0.5) ~ "draw and Species 1 dominant",
            (FrequencyChange1_Replicate1 == "increase" & FrequencyChange1_Replicate2 == "decrease" & RelativeAbundance_Tend1_Replicate1 <= 0.5) ~ "draw and Species 2 dominant",
            (FrequencyChange1_Replicate1 == "decrease" & FrequencyChange1_Replicate2 == "increase") ~ "mutual",
            (is.na(FrequencyChange1_Replicate1) | is.na(FrequencyChange1_Replicate2)) ~ "no-growth",
        )))

    # Coexistence pairs
    df_coexistence <- temp %>% filter(str_detect(Outcome, "draw")) %>%
        mutate(InteractionType = "coexistence") %>%
        mutate(across(starts_with("Species"), as.character)) %>%
        mutate(temp = ifelse(Outcome == "draw and Species 2 dominant", Species2, NA),
               Species2 = ifelse(Outcome == "draw and Species 2 dominant", Species1, Species2),
               Species1 = ifelse(Outcome == "draw and Species 2 dominant", temp, Species1)) %>%
        select(Community, Species1, Species2, Pair, InteractionType)

    # Exclusion pairs
    df_exclusion <- temp %>% filter(Outcome == "win" | Outcome == "lose") %>%
        mutate(InteractionType = "exclusion") %>%
        mutate(across(starts_with("Species"), as.character)) %>%
        mutate(temp = ifelse(Outcome == "lose", Species2, NA),
               Species2 = ifelse(Outcome == "lose", Species1, Species2),
               Species1 = ifelse(Outcome == "lose", temp, Species1)) %>%
        select(Community, Species1, Species2, Pair, InteractionType)

    # No-growth pairs
    df_nogrowth <- temp %>% filter(Outcome == "no-growth") %>%
        mutate(InteractionType = "no-growth") %>%
        select(Community, Species1, Species2, Pair, InteractionType)
    # Mutual exclusion
    df_mutual <- temp %>% filter(Outcome == "mutual") %>%
        mutate(InteractionType = "mutual exclusion") %>%
        select(Community, Species1, Species2, Pair, InteractionType)

    bind_rows(df_coexistence, df_exclusion, df_nogrowth, df_mutual) %>%
        return()
}
df_cp_init <- input_pairs %>%
    filter(str_detect(init_N0, "communityPairs")) %>%
    pull(init_N0) %>%
    paste0(output_dir, .) %>%
    lapply(function(x) {
        read_wide_file(x) %>%
            mutate(Well = factor(Well, paste0("W", 0:999)), .keep = "unused") %>%
            mutate(Species = factor(Species, sal$Species)) %>%
            # Remove rare species (relative abundance <0.01)
            group_by(Well) %>%
            mutate(RelativeAbundance = Abundance/sum(Abundance)) %>%
            #filter(RelativeAbundance > 0.01) %>%
            arrange(Well, Species) %>%
            # Community, or network
            mutate(Community = str_replace(x, paste0(output_dir, "communityPairs_"), "")  %>% str_replace("-1_init.csv", ""))
    }) %>%
    bind_rows() %>%
    mutate(Time = "Tinit")

df_cp_end <- input_pairs %>%
    filter(str_detect(init_N0, "communityPairs")) %>%
    pull(init_N0) %>% str_replace("_init.csv", "_end.csv") %>%
    #str_subset("W[1-9]-") %>%  # Subset for part of result done
    paste0(output_dir, .) %>%
    lapply(function(x) {
        read_wide_file(x) %>%
            mutate(Well = factor(Well, paste0("W", 0:1000)), .keep = "unused") %>%
            mutate(Species = factor(Species, sal$Species)) %>%
            # Remove rare species (relative abundance <0.01)
            group_by(Well) %>%
            mutate(RelativeAbundance = Abundance/sum(Abundance)) %>%
            # filter(RelativeAbundance > 0.01) %>%
            arrange(Well, Species) %>%
            # Community, or network
            mutate(Community = str_replace(x, paste0(output_dir, "communityPairs_"), "")  %>% str_replace("-1_end.csv", ""))
    }) %>%
    bind_rows() %>%
    mutate(Time = "Tend")

df_pairs <- determine_interaction(df_cp_init, df_cp_end) %>%
    left_join(rename_with(sal, ~paste0(., 1))) %>%
    left_join(rename_with(sal, ~paste0(., 2))) %>%
    mutate(PairConspecific = with(., case_when(
        (Family1 == Family2) ~ "conspecific",
        (Family1 != Family2) ~ "heterospecific"
    ))) %>%
    mutate(Community = factor(Community, paste0("W", 0:999))) %>%
    arrange(Community) %>%
    # Change variable names so conforming to function requirement
    rename(ID1 = Species1, ID2 = Species2)

df_isolates_ID <- bind_rows(select(df_pairs, Community, ID = ID1), select(df_pairs, Community, ID = ID2)) %>%
    distinct(Community, ID) %>%
    arrange(Community, ID) %>%
    group_by(Community) %>%
    mutate(Isolate = 1:n()) %>%
    ungroup()

## Pairwise outcome
df_pairs <- df_pairs %>%
    group_by(Community) %>%
    left_join(rename_with(df_isolates_ID, ~ paste0(., "1"), !contains("Community"))) %>%
    left_join(rename_with(df_isolates_ID, ~ paste0(., "2"), !contains("Community"))) %>%
    mutate(From = Isolate1, To = Isolate2) %>%
    ungroup()

## Pairwise frequency
df_pairs_freq <- bind_rows(df_cp_init, df_cp_end) %>%
    rename(ID = Species) %>%
    select(-Family, -Abundance) %>%
    pivot_wider(names_from = Time, values_from = RelativeAbundance, names_prefix = "RelativeAbundance_") %>%
    # Fill abundance = NA with 0
    replace_na(list(RelativeAbundance_Tend = 0)) %>%
    select(-RelativeAbundance_Tinit) %>%
    group_by(Community, Well) %>%
    mutate(Isolate = c(1,2)) %>%
    pivot_wider(names_from = Isolate, values_from = c(ID, RelativeAbundance_Tend), names_sep = "") %>%
    left_join(rename_with(df_isolates_ID, ~ paste0(., "1"), !contains("Community"))) %>%
    left_join(rename_with(df_isolates_ID, ~ paste0(., "2"), !contains("Community"))) %>%
    rename(Isolate1MeasuredFreq = RelativeAbundance_Tend1) %>%
    ungroup() %>%
    select(-Well, -RelativeAbundance_Tend2)

write_csv(df_pairs, here::here("data/output/df_pairs.csv"))
write_csv(df_pairs_freq, here::here("data/output/df_pairs_freq.csv"))

# Isolate ----
df_isolates_tournament <- df_communities_abundance %>%
    filter(Time == max(Time), Transfer == max(Transfer)) %>%
    arrange(Community) %>%
    group_by(Community) %>%
    summarize(Richness = n()) %>%
    select(comm = Community, everything()) %>%
    rowwise() %>%
    mutate(pairs_comm = df_pairs %>% filter(Community == comm) %>% list()) %>%
    mutate(tournaments_comm = pairs_comm %>% tournament_rank() %>% list()) %>%
    select(Community = comm, tournaments_comm) %>%
    unnest(cols = tournaments_comm)

df_isolates <- df_isolates_ID %>%
    left_join(df_isolates_tournament) %>%
    ungroup()




# Hierarchy measures ----
## Hierarchy following pairs
compute_hierarchy1 <- function(pairs_mock) {
    pairs_temp <- pairs_mock %>%
        select(Isolate1, Isolate2, InteractionType, From, To)
    isolates_tournament <- tournament_rank(pairs_temp) %>% select(Isolate, Score)

    pairs_temp %>%
        left_join(rename_with(isolates_tournament, ~ paste0(., "1")), by = "Isolate1") %>%
        left_join(rename_with(isolates_tournament, ~ paste0(., "2")), by = "Isolate2") %>%
        filter(InteractionType == "exclusion") %>%
        mutate(WinnerScore = ifelse(From == Isolate1, Score1, Score2),
               LoserScore = ifelse(From == Isolate1, Score2, Score1)) %>%
        mutate(FollowRank = (WinnerScore > LoserScore) %>% factor(c(T,F))) %>%
        count(FollowRank, .drop = F, name = "Count") %>%
        mutate(FractionFollowRank = Count / sum(Count)) %>%
        filter(FollowRank == T) %>%
        pull(FractionFollowRank) %>%
        return()
}
randomize_pairs1 <- function(x) {
    # Shuffle pairs
    rng <- order(runif(nrow(x),0,1))
    x$InteractionType <- x$InteractionType[rng]
    # Shuffle dominance
    rng <- sample(1:nrow(x), size = nrow(x)/2, replace = F)
    temp <- x$From[rng]
    x$From[rng] <- x$To[rng]
    x$To[rng] <- temp

    return(x)
}

### Permutation
communities_randomized_list1 <- rep(list(NA), nrow(df_communities))
tt <- proc.time()
for (i in 1:length(communities_randomized_list1)) {
    temp_tt <- proc.time()
    pairs_temp <- df_pairs %>% filter(Community == df_communities$Community[i])
    isolates_temp <- tibble(Community = unique(pairs_temp$Community), Isolate = sort(unique(c(pairs_temp$Isolate1, pairs_temp$Isolate2))))

    communities_randomized_list1[[i]] <- tibble(Community = df_communities$Community[i], Replicate = 1:b) %>%
        mutate(pairs_comm = list(pairs_temp)) %>%
        mutate(pairs_randomized = map(pairs_comm, randomize_pairs1)) %>%
        rowwise() %>%
        mutate(h1 = compute_hierarchy1(pairs_randomized))
    cat("\n", df_communities$Community[i] %>% as.character())
    cat("\n", (proc.time() - temp_tt)[3], "seconds")
}

communities_hierarchy_randomized1 <- bind_rows(communities_randomized_list1) %>%
    select(Community, Replicate, h1) %>%
    mutate(Community = factor(Community))
### obv
communities_hierarchy1 <- df_pairs %>%
    nest(pairs_comm = -Community) %>%
    rowwise() %>%
    mutate(h1 = compute_hierarchy1(pairs_comm)) %>%
    select(Community, h1)


## Higgins 2017
compute_comp_score2 <- function(pairs_comm, target_isolate) {
    pairs_comm %>%
        filter(Isolate1 == target_isolate | Isolate2 == target_isolate) %>%
        mutate(Freq = ifelse(Isolate1 == target_isolate, Isolate1MeasuredFreq, 1-Isolate1MeasuredFreq)) %>%
        pull(Freq) %>%
        mean()
}
compute_hierarchy2 <- function(isolates_mock, pairs_mock) {
    ranking <- isolates_mock %>%
        rowwise() %>%
        mutate(CompetitiveScore = compute_comp_score2(pairs_mock, Isolate)) %>%
        arrange(desc(CompetitiveScore)) %>%
        # Ranking by competitive score
        pull(Isolate)

    pairs_mock %>%
        mutate(Freq1 = Isolate1MeasuredFreq, Freq2 = 1-Freq1) %>%
        select(Isolate1, Isolate2, Freq1, Freq2) %>%
        mutate(Pair = 1:n()) %>%
        pivot_longer(cols = c(-Pair), names_to = ".value", names_pattern = "(.+)[12]") %>%
        mutate(Isolate = ordered(Isolate, ranking)) %>%
        group_by(Pair) %>%
        arrange(Pair, Isolate) %>%
        slice(1) %>%
        pull(Freq) %>%
        mean() %>%
        return()

}
randomize_pairs2 <- function(x) {
    # Shuffle pairs
    rng <- order(runif(nrow(x),0,1))
    x$Isolate1MeasuredFreq <- x$Isolate1MeasuredFreq[rng]
    # Shuffle dominance
    rng <- sample(1:nrow(x), size = nrow(x)/2, replace = F)
    x$Isolate1MeasuredFreq[rng] <- 1 - x$Isolate1MeasuredFreq[rng]

    return(x)
}

### Permutation
communities_randomized_list2 <- rep(list(NA), nrow(df_communities))
tt <- proc.time()
for (i in 1:length(communities_randomized_list2)) {
    temp_tt <- proc.time()
    pairs_temp <- df_pairs_freq %>%
        filter(Community == df_communities$Community[i]) %>%
        select(Community, Isolate1, Isolate2, Isolate1MeasuredFreq)
    isolates_temp <- tibble(Community = unique(pairs_temp$Community), Isolate = sort(unique(c(pairs_temp$Isolate1, pairs_temp$Isolate2))))

    communities_randomized_list2[[i]] <- tibble(Community = df_communities$Community[i], Replicate = 1:b) %>%
        mutate(pairs_comm = list(pairs_temp), isolates_comm = list(isolates_temp)) %>%
        mutate(pairs_randomized = map(pairs_comm, randomize_pairs2)) %>%
        rowwise() %>%
        mutate(h2 = compute_hierarchy2(isolates_comm, pairs_randomized))
    cat("\n", df_communities$Community[i] %>% as.character())
    cat("\n", (proc.time() - temp_tt)[3], "seconds")
}

communities_hierarchy_randomized2 <- bind_rows(communities_randomized_list2) %>%
    select(Community, Replicate, h2) %>%
    mutate(Community = factor(Community))

### obv
communities_hierarchy2 <- df_pairs_freq %>%
    rename(comm = Community) %>%
    nest(pairs_comm = -comm) %>%
    rowwise() %>%
    mutate(isolates_comm = list(filter(df_isolates, Community == comm))) %>%
    mutate(h2 = compute_hierarchy2(isolates_comm, pairs_comm)) %>%
    select(Community = comm, h2)


# Join data
df_communities_hierarchy <- left_join(communities_hierarchy1, communities_hierarchy2) %>%
    pivot_longer(-Community, names_to = "Metric", values_to = "HierarchyScore")
df_communities_hierarchy_randomized <- left_join(communities_hierarchy_randomized1, communities_hierarchy_randomized2) %>%
    pivot_longer(-c(Community, Replicate), names_to = "Metric", values_to = "HierarchyScore")
write_csv(df_communities_hierarchy, here::here("data/output/df_communities_hierarchy.csv"))
write_csv(df_communities_hierarchy_randomized, here::here("data/output/df_communities_hierarchy_randomized.csv"))







# Make networks ----
## obv
temp <- df_communities %>%
    rename(comm = Community) %>%
    rowwise() %>%
    mutate(isolate_comm = filter(df_isolates, Community == comm) %>% list()) %>%
    mutate(pairs_comm = filter(df_pairs, Community == comm) %>% list()) %>%
    mutate(Network = make_network(isolate_comm, pairs_comm) %>% list())
net_simulated_list <- temp$Network

## permutation
net_simulated_randomized_list <- rep(list(rep(list(NA), b)), length(net_simulated_list))
names(net_simulated_randomized_list) <- df_communities$Community

tt <- proc.time()
for (i in 1:length(net_simulated_list)) {
    temp_tt <- proc.time()
    for (b_loop_index in 1:b) {
        set.seed(b_loop_index)
        net_simulated_randomized_list[[i]][[b_loop_index]] <- randomize_network(net_simulated_list[[i]])
        if (b_loop_index %% 100 == 0) cat("\n boostrap =", b_loop_index)
    }
    # Print
    cat("\n\n", df_communities$Community[i])
    cat("\n", (proc.time() - temp_tt)[3], "seconds")
    if (i == length(net_simulated_list)) cat("\n\n total time:", (proc.time() - tt)[3], "seconds\n\n")
}

save(net_simulated_list, file = "~/Dropbox/lab/invasion-network/data/output/network_simulated.Rdata")
save(net_simulated_randomized_list, file = "~/Dropbox/lab/invasion-network/data/output/network_simulated_randomized.Rdata")



# Motif ----
# load("~/Dropbox/lab/invasion-network/data/output/network_simulated.Rdata")
# load("~/Dropbox/lab/invasion-network/data/output/network_simulated_randomized.Rdata")

## Obv
df_motif <- df_communities %>%
    mutate(Network = net_simulated_list) %>%
    rowwise() %>%
    mutate(temp = list(tibble(Motif = 1:7, Count = count_motif(Network))), .keep = "unused") %>%
    unnest(cols = c(temp)) %>%
    group_by(Community) %>%
    mutate(Fraction = Count / sum(Count))

## Permutation
df_motif_randomized <- df_communities %>%
    mutate(NetworkList = net_simulated_randomized_list) %>%
    rowwise() %>%
    mutate(temp = list(lapply(NetworkList, function(x) {
        tibble(Motif = 1:7, Count = count_motif(x))
    }) %>% bind_rows(.id = "Replicate")), .keep = "unused") %>%
    unnest(cols = c(temp)) %>%
    group_by(Replicate, Community) %>%
    mutate(Fraction = Count / sum(Count))

## Find percentiles
df_motif_percentile <- df_motif_randomized %>%
    group_by(Community, Motif) %>%
    arrange(desc(Count)) %>%
    slice(ceiling(b * 0.05), floor(b * 0.95)) %>%
    mutate(Percentile = c("p5", "p95")) %>%
    select(Community, Motif, Percentile, Count) %>%
    pivot_wider(names_from = Percentile, values_from = Count)

df_motif <- df_motif %>%
    left_join(df_motif_percentile) %>%
    mutate(Sign = case_when(Count > p95 ~ "top",
                            Count < p5 ~ "bottom",
                            Count < p95 & Count > p5 ~ "n.s."))


write_csv(df_motif, here::here("data/output/df_motif.csv"))
write_csv(df_motif_randomized, here::here("data/output/df_motif_randomized.csv"))





# Diagonal analysis ----
adj_from_net <- function(net) {
  # Get adjacent matrix
  net_m <- get.adjacency(net, attr = "InteractionType", sparse = F)

  # Exclusion: win or lose
  temp_index <- which(net_m=="exclusion", arr.ind = T) %>% as.data.frame()
  for(i in 1:nrow(temp_index)) net_m[temp_index$col[i], temp_index$row[i]] <- "lose"
  net_m[net_m=="exclusion"] <- "win"

  # Bistability
  temp_index2 <- which(net_m=="bistability", arr.ind = T) %>% as.data.frame()
  for(i in 1:nrow(temp_index2)) net_m[temp_index2$col[i], temp_index2$row[i]] <- "bistability"
  net_m[net_m=="bistability"] <- "bistability"

  # Neutrality
  temp_index2 <- which(net_m=="neutrality", arr.ind = T) %>% as.data.frame()
  for(i in 1:nrow(temp_index2)) net_m[temp_index2$col[i], temp_index2$row[i]] <- "neutrality"
  net_m[net_m=="neutrality"] <- "neutrality"

  # Diagonal
  diag(net_m) <- "self"

  # Undefined
  temp_index <- which(net_m=="undefined", arr.ind = T) %>% as.data.frame()
  for(i in 1:nrow(temp_index)) net_m[temp_index$col[i], temp_index$row[i]] <- "undefined"

  # NA
  net_m[net_m == "" | is.na(net_m)] <- NA

  return(net_m)
}
## Find the fraction of coexistence as a function of distance to diagonal
diag_distance <- function (net, observation = F) {
  # Convert network to matrix
  temp_matrix <- adj_from_net(net)

  # Re-order the matrix axis by the isolates' competitive rank
  temp_rank <- tournament_rank(net)$Isolate
  temp_matrix <- temp_matrix[temp_rank, temp_rank]

  # Order the matrix axis by competitive score
  t1 <- which(temp_matrix == "coexistence", arr.ind = T) %>%
    as_tibble() %>%
    filter(col > row) %>%
    mutate(DistanceToDiagonal = abs(row - col)) %>%
    group_by(DistanceToDiagonal) %>%
    summarize(CountCoexistence = n())

  # Total count for each distance
  t2 <- which(temp_matrix == "win" | temp_matrix == "lose" | temp_matrix == "coexistence", arr.ind = T) %>%
      as_tibble() %>%
      filter(col > row) %>%
      mutate(DistanceToDiagonal = abs(row - col)) %>%
      group_by(DistanceToDiagonal) %>%
      summarize(TotalCount = n())

  #
  if (observation) {
      full_join(t1, t2, by = "DistanceToDiagonal") %>%
          replace_na(list(CountCoexistence = 0)) %>%
          return()
  } else {
      return(t1)
  }
}

## Observation
df_diag <- df_communities %>%
    mutate(Network = net_simulated_list) %>%
    rowwise() %>%
    mutate(temp = list(diag_distance(Network)), .keep = "unused") %>%
    unnest(cols = c(temp))

## Permutation
# Count the distance to diagonal in randomized networks
net_simulated_diag_randomized_list <- rep(list(NA), length(net_simulated_list))
names(net_simulated_diag_randomized_list) <- names(net_simulated_list)

df_diag_randomized <- df_communities %>%
    mutate(NetworkList = net_simulated_randomized_list) %>%
    rowwise() %>%
    mutate(temp = list(lapply(NetworkList, diag_distance) %>% bind_rows(.id = "Replicate")), .keep = "unused") %>%
    unnest(cols = c(temp))

# Save the result
write_csv(df_diag, file = here::here("data/output/df_diag.csv"))
write_csv(df_diag_randomized, file = here::here("data/output/df_diag_randomized.csv"))


if (FALSE) {
df_diag_randomized %>%
    group_by(DistanceToDiagonal, Replicate) %>%
    summarize(Count = sum(CountCoexistence)) %>%
    ggplot(aes(x = DistanceToDiagonal, y = Count, group = DistanceToDiagonal)) +
    geom_boxplot() +
    geom_jitter(height = 0, width = 0.1) +
    geom_point(data = df_diag %>% group_by(DistanceToDiagonal) %>% summarize(Count = sum(CountCoexistence)),
               aes(x = DistanceToDiagonal, y = Count, group = DistanceToDiagonal), color = "red") +
    theme_classic()

}














