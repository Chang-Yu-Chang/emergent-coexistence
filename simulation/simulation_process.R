# Clean up the simulation data
library(tidyverse)
library(tidygraph)
library(igraph)
library(ggraph)
library(cowplot)
source(here::here("plotting_scripts/network_functions.R"))

output_dir <- "~/Dropbox/lab/invasion-network/simulation/data/raw12/"
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

b = 1000


# Community abundance ----
## Self-assembled communities
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
    mutate(Assembly = "self_assembly") %>%
    rename(ID = Species) %>%
    select(Assembly, everything()) %>%
    # Remove rare species (relative abundance <0.01)
    group_by(Community, Transfer, Time) %>%
    mutate(RelativeAbundance = Abundance/sum(Abundance)) %>%
    filter(RelativeAbundance > 0.01) %>%
    ungroup()


df_communities_abundance %>%
    filter(Community %in% paste0("W", 0:9)) %>%
    filter(Transfer == max(Transfer), Time == max(Time)) %>%
    mutate(Community = factor(Community, c("W1", "W0", paste0("W", 2:20)))) %>%
    mutate(Family = ifelse(Family == "F0", "fermenter", ifelse(Family == "F1", "respirator", Family))) %>%
    ggplot() +
    geom_col(aes(x = Community, y = Abundance, fill = Family, group = ID),
             color = 1, size = .5, position = "fill") +
    scale_fill_manual(values = fermenter_color) +
    scale_x_discrete(label = 1:20) +
    scale_y_continuous(breaks = c(0, .5, 1), expand = c(0,0)) +
    #facet_grid(Time~.) +
    theme_classic() +
    theme(panel.border = element_rect(color = 1, fill = NA),
          axis.title = element_text(size = 15, color = 1),
          axis.text = element_text(size = 10, color = 1),
          legend.text = element_text(size = 15, color = 1),
          legend.title = element_text(size = 15, color = 1),
          legend.position = "top") +
    guides(alpha = "none", color = "none") +
    labs(fill = "")

write_csv(df_communities_abundance, here::here("data/output/df_communities_abundance.csv"))



# ## Pool networks
# df_pool_init <- paste0(output_dir, "poolNetwork-1_init.csv") %>%
#     read_csv(col_types = cols()) %>%
#     mutate(Transfer = 0, Time = 0) %>%
#     mutate(Family = factor(Family, c("F0", "F1"))) %>%
#     mutate(Species = factor(Species, paste0("S", 0:999))) %>%
#     mutate(Community = factor(Community, paste0("W", 0:999)), .keep = "unused")


# Pairwise outcome ----
determine_interaction <- function(pairs_init, pairs_end) {
    temp <- bind_rows(pairs_init, pairs_end) %>%
        ungroup() %>%
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

## Community pairs
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
    mutate(Time = "Tinit") %>%
    ungroup()
df_cp_end <-
    # input_pairs %>%
    # filter(str_detect(init_N0, "communityPairs")) %>%
    # pull(init_N0) %>% str_replace("_init.csv", "_end.csv") %>%
    list.files(output_dir, pattern = "communityPairs") %>% str_subset("end.csv") %>%
    #str_subset("W[1-9]") %>%  # Subset for part of result done
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
    mutate(Time = "Tend") %>%
    ungroup()


df_cp_freq <- bind_rows(df_cp_init, df_cp_end) %>%
    group_by(Community, Well) %>%
    select(Community, Well, Species, RelativeAbundance, Time) %>%
    pivot_wider(names_from = Time, values_from = RelativeAbundance, names_prefix = "RelativeAbundance_") %>%
    # Fill abundance = NA with 0
    replace_na(list(RelativeAbundance_Tend = 0)) %>%
    group_by(Community, Well) %>%
    mutate(Isolate = c(1,2)) %>%
    pivot_wider(names_from = Isolate, values_from = c(Species, RelativeAbundance_Tinit, RelativeAbundance_Tend), names_sep = "") %>%
    ungroup() %>%
    mutate(Assembly = "self_assembly") %>%
    rename(ID1 = Species1, ID2 = Species2) %>%
    select(Assembly, everything(), -Well)

write_csv(df_cp_freq, here::here("data/output/df_cp_freq.csv"))

if (FALSE) {
## Pool pairs
df_pp_init <- input_pairs %>%
    filter(str_detect(init_N0, "poolPairs")) %>%
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
            mutate(Community = str_replace(x, paste0(output_dir, "poolPairs_"), "")  %>% str_replace("-1_init.csv", ""))
    }) %>%
    bind_rows() %>%
    mutate(Time = "Tinit") %>%
    ungroup()
df_pp_end <- input_pairs %>%
    filter(str_detect(init_N0, "poolPairs")) %>%
    pull(init_N0) %>% str_replace("_init.csv", "_end.csv") %>%
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
            mutate(Community = str_replace(x, paste0(output_dir, "poolPairs_"), "")  %>% str_replace("-1_end.csv", ""))
    }) %>%
    bind_rows() %>%
    mutate(Time = "Tend") %>%
    ungroup()

}


## List of communities
temp1 <- df_cp_init %>%
    distinct(Community, Species) %>%
    arrange(Community) %>%
    group_by(Community) %>%
    summarize(Richness = n()) %>%
    mutate(Assembly = "self_assembly")
# temp2 <- df_pp_init %>%
#     distinct(Community, Species) %>%
#     group_by(Community) %>%
#     summarize(Richness = n()) %>%
#     mutate(Assembly = "random_assembly")
df_communities <- #bind_rows(temp1, temp2) %>%
    bind_rows(temp1) %>%
    select(Assembly, Community, Richness) %>%
    mutate(Assembly = factor(Assembly, c("self_assembly", "random_assembly"))) %>%
    mutate(Community = factor(Community, paste0("W", 0:999))) %>%
    arrange(Assembly, Community) %>%
    filter(Assembly == "self_assembly")

write_csv(df_communities, here::here("data/output/df_communities.csv"))

## Isolate ID
temp1 <- df_cp_init %>%
    distinct(Community, Species) %>%
    left_join(sal) %>%
    rename(ID = Species) %>%
    mutate(Assembly = "self_assembly") %>%
    select(Assembly, Community, Family, ID)
# temp2 <- df_pp_init %>%
#     distinct(Community, Species) %>%
#     left_join(sal) %>%
#     rename(ID = Species) %>%
#     mutate(Assembly = "random_assembly") %>%
#     select(Assembly, Community, Family, ID)
df_isolates_ID <- #bind_rows(temp1, temp2) %>%
    bind_rows(temp1) %>%
    mutate(Assembly = factor(Assembly, c("self_assembly", "random_assembly"))) %>%
    #filter(ID %in% unique(c(df_pp_init$Species, df_cp_init$Species))) %>%
    filter(ID %in% unique(c(df_cp_init$Species))) %>%
    group_by(Assembly, Community) %>%
    arrange(Assembly, Community) %>%
    mutate(Isolate = 1:n())




## Determine pairwise result
df_pairs <- bind_rows(
    determine_interaction(df_cp_init, df_cp_end) %>% mutate(Assembly = "self_assembly"),
    #determine_interaction(df_pp_init, df_pp_end) %>% mutate(Assembly = "random_assembly")
) %>%
    left_join(rename_with(sal, ~paste0(., 1))) %>%
    left_join(rename_with(sal, ~paste0(., 2))) %>%
    mutate(PairConspecific = with(., case_when(
        (Family1 == Family2) ~ "conspecific",
        (Family1 != Family2) ~ "heterospecific"
    ))) %>%
    mutate(Community = factor(Community, paste0("W", 0:999))) %>%
    mutate(Assembly = factor(Assembly, c("self_assembly", "random_assembly"))) %>%
    arrange(Assembly, Community) %>%
    # Change variable names so conforming to function requirement
    rename(ID1 = Species1, ID2 = Species2) %>%
    filter(ID1 %in% df_isolates_ID$ID & ID2 %in% df_isolates_ID$ID) %>%
    left_join(rename_with(df_isolates_ID, ~ paste0(., "1"), !contains("Community") & !contains("Assembly"))) %>%
    left_join(rename_with(df_isolates_ID, ~ paste0(., "2"), !contains("Community") & !contains("Assembly"))) %>%
    mutate(From = Isolate1, To = Isolate2) %>%
    ungroup() %>%
    select(Assembly, everything())


## Pairwise frequency
df_cp_freq <- bind_rows(df_cp_init, df_cp_end) %>%
    rename(ID = Species) %>%
    select(-Family, -Abundance) %>%
    pivot_wider(names_from = Time, values_from = RelativeAbundance, names_prefix = "RelativeAbundance_") %>%
    # Fill abundance = NA with 0
    replace_na(list(RelativeAbundance_Tend = 0)) %>%
    select(-RelativeAbundance_Tinit) %>%
    group_by(Community, Well) %>%
    mutate(Isolate = c(1,2)) %>%
    pivot_wider(names_from = Isolate, values_from = c(ID, RelativeAbundance_Tend), names_sep = "") %>%
    filter(ID1 %in% df_isolates_ID$ID, ID2 %in% df_isolates_ID$ID) %>%
    left_join(rename_with(df_isolates_ID, ~ paste0(., "1"), !contains("Community"))) %>%
    left_join(rename_with(df_isolates_ID, ~ paste0(., "2"), !contains("Community"))) %>%
    rename(Isolate1MeasuredFreq = RelativeAbundance_Tend1) %>%
    ungroup() %>%
    select(-Well, -RelativeAbundance_Tend2) %>%
    mutate(Assembly = "self_assembly")

# df_pp_freq <- bind_rows(df_pp_init, df_pp_end) %>%
#     rename(ID = Species) %>%
#     select(-Family, -Abundance) %>%
#     pivot_wider(names_from = Time, values_from = RelativeAbundance, names_prefix = "RelativeAbundance_") %>%
#     # Fill abundance = NA with 0
#     replace_na(list(RelativeAbundance_Tend = 0)) %>%
#     select(-RelativeAbundance_Tinit) %>%
#     group_by(Community, Well) %>%
#     mutate(Isolate = c(1,2)) %>%
#     pivot_wider(names_from = Isolate, values_from = c(ID, RelativeAbundance_Tend), names_sep = "") %>%
#     filter(ID1 %in% df_isolates_ID$ID, ID2 %in% df_isolates_ID$ID) %>%
#     left_join(rename_with(df_isolates_ID, ~ paste0(., "1"), !contains("Community"))) %>%
#     left_join(rename_with(df_isolates_ID, ~ paste0(., "2"), !contains("Community"))) %>%
#     rename(Isolate1MeasuredFreq = RelativeAbundance_Tend1) %>%
#     ungroup() %>%
#     select(-Well, -RelativeAbundance_Tend2) %>%
#     mutate(Assembly = "random_assembly")

df_pairs_freq <-
    bind_rows(df_cp_freq) %>%
    #bind_rows(df_cp_freq, df_pp_freq) %>%
    mutate(Assembly = factor(Assembly, c("self_assembly", "random_assembly"))) %>%
    select(Assembly, everything())


write_csv(df_pairs, here::here("data/output/df_pairs.csv"))
write_csv(df_pairs_freq, here::here("data/output/df_pairs_freq.csv"))

# Isolate ----
df_isolates_tournament <- df_communities %>%
    filter(Assembly == "self_assembly") %>%
    select(assem = Assembly, comm = Community, everything()) %>%
    rowwise() %>%
    mutate(pairs_comm = df_pairs %>% filter(Assembly == assem, Community == comm) %>% list()) %>%
    mutate(tournaments_comm = pairs_comm %>% tournament_rank() %>% list()) %>%
    select(Assembly = assem, Community = comm, tournaments_comm) %>%
    unnest(cols = tournaments_comm)

df_isolates <- df_isolates_ID %>%
    left_join(df_isolates_tournament) %>% # some species in community are in the transient to extinction but picked. They have low abundance
    ungroup()


write_csv(df_isolates, here::here("data/output/df_isolates.csv"))


# Hierarchy measures ----
## Hierarchy following pairs
if (FALSE) {
### Permutation
communities_randomized_list1 <- rep(list(NA), nrow(df_communities))
tt <- proc.time()
for (i in 1:length(communities_randomized_list1)) {
    temp_tt <- proc.time()
    assem <- df_communities$Assembly[i]
    comm <- df_communities$Community[i]

    pairs_temp <- df_pairs %>% filter(Assembly == assem, Community == comm)
    isolates_temp <- tibble(Community = unique(pairs_temp$Community), Isolate = sort(unique(c(pairs_temp$Isolate1, pairs_temp$Isolate2))))

    communities_randomized_list1[[i]] <-
        tibble(Assembly = assem, Community = comm, Replicate = 1:b) %>%
        mutate(pairs_comm = list(pairs_temp)) %>%
        mutate(pairs_randomized = map(pairs_comm, randomize_pairs1)) %>%
        rowwise() %>%
        mutate(h1 = compute_hierarchy1(pairs_randomized))
    cat("\n", assem, comm %>% as.character())
    cat("\n", (proc.time() - temp_tt)[3], "seconds")
}
communities_hierarchy_randomized1 <- bind_rows(communities_randomized_list1) %>%
    select(Community, Replicate, h1) %>%
    mutate(Community = factor(Community))

}
### obv
communities_hierarchy1 <- df_pairs %>%
    nest(pairs_comm = -c(Assembly, Community)) %>%
    rowwise() %>%
    mutate(h1 = compute_hierarchy1(pairs_comm)) %>%
    select(Assembly, Community, h1)


## Higgins 2017
if (FALSE) {

### Permutation
communities_randomized_list2 <- rep(list(NA), nrow(df_communities))
tt <- proc.time()
for (i in 1:length(communities_randomized_list2)) {
    temp_tt <- proc.time()
    assem <- df_communities$Assembly[i]
    comm <- df_communities$Community[i]
    pairs_temp <- df_pairs_freq %>%
        filter(Assembly == assem, Community == comm) %>%
        select(Community, Isolate1, Isolate2, Isolate1MeasuredFreq)
    isolates_temp <- tibble(Community = unique(pairs_temp$Community), Isolate = sort(unique(c(pairs_temp$Isolate1, pairs_temp$Isolate2))))

    communities_randomized_list2[[i]] <-
        tibble(Assembly == assem, Community = comm, Replicate = 1:b) %>%
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
}

### obv
communities_hierarchy2 <- df_pairs_freq %>%
    rename(assem = Assembly, comm = Community) %>%
    nest(pairs_comm = -c(assem, comm)) %>%
    rowwise() %>%
    mutate(isolates_comm = list(filter(df_isolates, Assembly == assem, Community == comm))) %>%
    mutate(h2 = compute_hierarchy2(isolates_comm, pairs_comm)) %>%
    select(Assembly = assem, Community = comm, h2)


# Join data
df_communities_hierarchy <- left_join(communities_hierarchy1, communities_hierarchy2) %>%
    pivot_longer(-c(Assembly, Community), names_to = "Metric", values_to = "HierarchyScore")
#df_communities_hierarchy_randomized <- left_join(communities_hierarchy_randomized1, communities_hierarchy_randomized2) %>%
#    pivot_longer(-c(Community, Replicate), names_to = "Metric", values_to = "HierarchyScore")
write_csv(df_communities_hierarchy, here::here("data/output/df_communities_hierarchy.csv"))
#write_csv(df_communities_hierarchy_randomized, here::here("data/output/df_communities_hierarchy_randomized.csv"))









# Make networks ----
## obv
temp <- df_communities %>%
    filter(Assembly == "self_assembly") %>%
    rename(assem = Assembly, comm = Community) %>%
    rowwise() %>%
    mutate(isolate_comm = filter(df_isolates, Assembly == assem, Community == comm) %>% list()) %>%
    mutate(pairs_comm = filter(df_pairs, Assembly == assem, Community == comm) %>% list()) %>%
    mutate(Network = make_network(isolate_comm, pairs_comm) %>% list())
net_simulated_list <- temp$Network %>% set_names(paste0(temp$assem, "_", temp$comm))

## permutation
net_simulated_randomized_list <- rep(list(rep(list(NA), b)), length(net_simulated_list))
names(net_simulated_randomized_list) <- df_communities$Community[1:20]

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
load("~/Dropbox/lab/invasion-network/data/output/network_simulated.Rdata")
load("~/Dropbox/lab/invasion-network/data/output/network_simulated_randomized.Rdata")

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
## Observation
df_diag <- df_communities %>%
    mutate(Network = net_simulated_list) %>%
    rowwise() %>%
    mutate(Diagonal = count_diag_coexistence(Network) %>% mutate(Community) %>% list()) %>%
    select(-Richness, -Network, -Community) %>%
    unnest(Diagonal) %>%
    bind_rows() %>%
    select(Assembly, Community, everything())

# Permutation
df_diag_randomized_list <- rep(list(NA), length(net_simulated_list))
names(df_diag_randomized_list) <- names(net_simulated_list)
tt <- proc.time()
for (i in 1:length(net_simulated_list)) {
    temp_tt <- proc.time()
    df_diag_randomized_list[[i]] <- lapply(net_simulated_randomized_list[[i]], count_diag_coexistence) %>%
        bind_rows(.id = "Replicate")
    # Print
    cat("\n\n", df_communities$Community[i])
    cat("\n", (proc.time() - temp_tt)[3], "seconds")
    if (i == length(net_simulated_list)) cat("\n\n total time:", (proc.time() - tt)[3], "seconds")
}
df_diag_randomized <- bind_rows(df_diag_randomized_list, .id = "Community") %>%
    separate(Community, into = c("Assembly", "Community"), sep = "_W") %>%
    mutate(Community = paste0("W", Community))
#
write_csv(df_diag, file = here::here("data/output/df_diag.csv"))
write_csv(df_diag_randomized, file = here::here("data/output/df_diag_randomized.csv"))




# Number of cliques or components ----
df_component <- df_communities %>%
    mutate(Network = net_simulated_list) %>%
    rowwise() %>%
    mutate(Component = count_component(Network)) %>%
    select(Assembly, Community, Component)

df_component_randomized <- df_communities %>%
    mutate(NetworkList = net_simulated_randomized_list) %>%
    rowwise() %>%
    mutate(temp = list(tibble(Replicate = 1:b, Component = sapply(NetworkList, count_component)))) %>%
    unnest(cols = c(temp)) %>%
    group_by(Replicate, Community) %>%
    select(Assembly, Community, Replicate, Component)


write_csv(df_component, file = here::here("data/output/df_component.csv"))
write_csv(df_component_randomized, file = here::here("data/output/df_component_randomized.csv"))











