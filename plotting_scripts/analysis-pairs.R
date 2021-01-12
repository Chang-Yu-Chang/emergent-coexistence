# Pairs from top-down communities

suppressWarnings(suppressMessages(library(tidyverse)))
suppressWarnings(suppressMessages(library(data.table)))
#options(dplyr.summarise.inform = FALSE)
suppressWarnings(suppressMessages(library(cowplot)))
suppressWarnings(suppressMessages(library(ggpubr)))
suppressWarnings(suppressMessages(library(tidygraph)))
suppressWarnings(suppressMessages(library(ggraph)))
suppressWarnings(suppressMessages(library(igraph)))
source("network_functions.R")

read_pair_data <- function (output_dir, pattern) {
    #output_dir <- "../data/raw/simulation/"
    #pattern <- "pair_from_random_species"
    df_pair_data <- list.files(output_dir, pattern = "_composition.txt", full.names = T) %>%
        grep(pattern, ., value = T) %>%
        lapply(fread) %>%
        rbindlist %>%
        filter(Type == "consumer") %>%
        group_by(exp_id, Well, Transfer) %>%
        mutate(RelativeAbundance = Abundance / sum(Abundance))
    return(df_pair_data)
}
read_pair_list <- function (output_dir, pattern) {
    #output_dir <- "../data/raw/simulation/"
    #pattern <- "pair_from_random_species"
    df_pair_initial <- list.files(output_dir, pattern = "\\d.txt", full.names = T) %>%
        grep(pattern, ., value = T) %>%
        lapply(function(x) {
            temp <- strsplit(x, "/") %>% unlist
            exp_id <- sub(".txt", "", temp[length(temp)])
            fread(x) %>% mutate(exp_id = exp_id)
        }) %>%
        rbindlist %>%
        filter(Type == "consumer") %>%
        select(exp_id, Well, Pair, InitialFrequency, Type, ID)

    temp <- rep(list(df_pair_initial), 11)
    for (i in 1:11) temp[[i]] <- mutate(temp[[i]], Transfer = i-1)
    df_pair_list <- rbindlist(temp)
    return(df_pair_list)
}
reshape_pair_data <- function(pair_data, pair_list) {
    # pair_data <- df_pair_data
    # pair_list <- df_pair_list
    df_pair_data2 <- pair_list %>%
        left_join(pair_data) %>%
        replace_na(list(Abundance = 0)) %>%
        mutate(Isolate = rep(c("Isolate1", "Isolate2"), nrow(.)/2)) %>%
        mutate(temp = rep(c("RelativeAbundance1", "RelativeAbundance2"), nrow(.)/2)) %>%
        select(-Abundance) %>%
        group_by(exp_id, Well, Pair, InitialFrequency, Transfer) %>%
        #pivot_wider(names_from = c("temp", values_from = "RelativeAbundance") %>%
        pivot_wider(names_from = c("Isolate", "temp"), values_from = c("ID", "RelativeAbundance")) %>%
        dplyr::rename(Isolate1 = ID_Isolate1_RelativeAbundance1, Isolate2 = ID_Isolate2_RelativeAbundance2,
                  RelativeAbundance1 = RelativeAbundance_Isolate1_RelativeAbundance1, RelativeAbundance2 = RelativeAbundance_Isolate2_RelativeAbundance2)
    return(df_pair_data2)
}
determine_pair_outcome <- function(pair_data_reshaped) {
    df_frequency_change <- pair_data_reshaped %>%
        filter(Transfer %in% c(0, 10)) %>%
        select(-RelativeAbundance2, -Type) %>%
        group_by(exp_id, Well, Pair, InitialFrequency) %>%
        pivot_wider(names_from = "Transfer", names_prefix = "T", values_from = "RelativeAbundance1") %>%
        mutate(FrequencyChange = ifelse((T10 - T0)>0, "T", "F")) %>%  # TRUE = Isolate1 increases, FALSE = Isolate1 decreases
        ungroup() %>%
        select(exp_id, Pair, Isolate1, Isolate2, FrequencyChange) %>%
        group_by(exp_id, Pair, Isolate1, Isolate2) %>%
        summarise(FrequencyChangePattern = paste0(FrequencyChange, collapse = "-")) %>%
        {.}
    map_frequency_pattern <- tibble(FrequencyChangePattern = c("T-T-T", "F-F-F", "T-F-F", "T-T-F", "T-F-T", "F-T-F", "F-F-T", "F-T-T"),
                                    InteractionType = c("exclusion", "exclusion", "coexistence", "coexistence", "coexistence", "coexistence", "mutual exclusion", "mutual exclusion"))

    df_frequency_change$From = NA
    df_frequency_change$To = NA
    df_frequency_change$InteractionType = NA

    for (i in 1:nrow(df_frequency_change)) {
        if (df_frequency_change$FrequencyChangePattern[i] == "T-T-T") {
            from <- df_frequency_change$Isolate1[i]
            to <- df_frequency_change$Isolate2[i]
            interaction_type <- "exclusion"
        } else if (df_frequency_change$FrequencyChangePattern[i] == "F-F-F") {
            from <- df_frequency_change$Isolate2[i]
            to <- df_frequency_change$Isolate1[i]
            interaction_type <- "exclusion"
        } else if (df_frequency_change$FrequencyChangePattern[i] %in% c("T-F-F", "T-T-F", "T-F-T", "F-T-F")) {
            from <- df_frequency_change$Isolate1[i]
            to <- df_frequency_change$Isolate2[i]
            interaction_type <- "coexistence"
        } else {
            from <- NA
            to <- NA
            interaction_type <- NA
        }

        df_frequency_change$From[i] <- from
        df_frequency_change$To[i] <- to
        df_frequency_change$InteractionType[i] <- interaction_type
    }

    pair_outcome <- df_frequency_change %>%
        select(exp_id, Pair, Isolate1, Isolate2, InteractionType, From, To) %>%
        ungroup()
    return(pair_outcome)
}
determine_community_motif <- function(pair_from_set) {
    isolate_list <- sort(unique(c(pair_from_set$Isolate1, pair_from_set$Isolate2)))
    nodes <- tibble(Isolate = 1:length(isolate_list), ID = isolate_list)
    edges <- pair_from_set %>%
        mutate(from = match(From, nodes$ID), to = match(To, nodes$ID)) %>%
        select(InteractionType, from, to)
    edges_coext <- edges[edges$InteractionType == "coexistence",]
    edges_coext[,c("from", "to")] <- edges_coext[,c("to", "from")] # Add the mutual edges for coexistence links
    edges <- rbind(edges, edges_coext)

    graph <- tbl_graph(nodes = nodes, edges = filter(edges, !is.na(InteractionType)), directed = T)
    motif_counts <- tibble(Motif = 1:7, Count = count_motif(graph))

    return(motif_counts)
}

plot_pair_line <- function(pair_data_reshape, exp_id) {
    # df_pair_data_reshaped <- pair_data_reshaped
    # exp_id = "simple_medium-pair_from_random_species_1-3"
    pair_data_reshape %>%
        filter(exp_id == exp_id) %>%
        mutate(Isolate1 = factor(Isolate1)) %>%
        mutate(Isolate2 = factor(Isolate2)) %>%
        ggplot(aes(x = Transfer, y = RelativeAbundance1, color = InitialFrequency)) +
        geom_point() +
        geom_line() +
        facet_grid(Isolate1~Isolate2) +
        theme_cowplot()
}
plot_pair_outcome <- function(pair_outcome, normalize_sum = F) {
    pair_outcome <- pair_outcome %>%
        group_by(seed, Community, InteractionType) %>%
        summarize(Count = n())
    if (normalize_sum) pair_outcome <- pair_outcome %>% group_by(seed, Community) %>% mutate(Count = Count / sum(Count))

    pair_outcome %>%
        group_by(seed, InteractionType) %>%
        summarize(MeanCount = mean(Count)) %>%
        ggplot(aes(x = InteractionType, y = MeanCount), color = 1) +
        geom_boxplot() +
        geom_jitter() +
        theme_cowplot()
}
plot_motif_count <- function(motif_count, normalize_sum = F) {
    if (normalize_sum) {motif_count <- motif_count %>% group_by(seed, Community) %>% mutate(Count = Count / sum(Count))}
    motif_count %>%
        group_by(seed, Motif) %>%
        summarize(MeanCount = mean(Count)) %>%
        ggplot(aes(x = Motif, y = MeanCount, group = Motif)) +
        geom_boxplot() +
        geom_jitter() +
        scale_x_continuous(breaks = 1:7) +
        guides(color = F) +
        theme_cowplot()
}


#

# list_tt <- c("simple_medium1", "simple_medium2", "simple_medium3", "simple_medium5", "simple_medium6", "simple_medium7",
#              "rich_medium1", "rich_medium2", "rich_medium3", "rich_medium5", "rich_medium6", "rich_medium7")
list_tt <- c("rich_medium5", "rich_medium6", "rich_medium7")

list_ss <- c("pair_from_random_species", "pair_from_top_down_community")
list_df_pair <- rep(list(NA), length(list_tt) * length(list_ss))
list_df_motif <- rep(list(NA), length(list_tt) * length(list_ss))

aggregate_result <- function (tt, ss) {

    # Random species
    df_pair_data <- read_pair_data("~/Dropbox/invasion-network/data/raw/simulation/", paste0(tt, "-", ss))
    df_pair_list <- read_pair_list("~/Dropbox/invasion-network/data/raw/simulation/", paste0(tt, "-", ss))
    df_pair_data_reshaped <- reshape_pair_data(df_pair_data, df_pair_list)

    df_pair_outcome <- determine_pair_outcome(df_pair_data_reshaped)
    df_motif_count <- df_pair_outcome %>%
        split.data.frame(f = .$exp_id) %>%
        lapply(determine_community_motif) %>%
        rbindlist(idcol = "exp_id")
    return(list(df_pair_outcome=df_pair_outcome, df_motif_count=df_motif_count))
}

for (i in 1:length(list_tt)) {
    for (j in 1:length(list_ss)) {
        temp <- aggregate_result(list_tt[i], list_ss[j])
        list_df_pair[[length(list_ss)*(i-1)+j]] <- temp$df_pair_outcome
        list_df_motif[[length(list_ss)*(i-1)+j]] <- temp$df_motif_count
        cat("\n", list_tt[i], "\t", list_ss[j])
    }
}

#
list_treatments <- tibble(
    folder_id = c(rep("simple_medium", 8), rep("rich_medium", 8)),
    #exp_id = c(paste0("simple_medium", 1:8), paste0("rich_medium", 1:8)),
    Medium = c(paste0("simple_medium", 1:8), paste0("rich_medium", 1:8)),
    rich_medium = c(rep(F, 8), rep(T, 8)),
    n_inoc = c(rep(10^3, 8), rep(10^6, 8)),
    n_propagation = c(rep(20, 8), rep(1,8)),
    dilution = c(rep(0.001, 16)),
    l = rep(c(0, 0, 0, 0, 0.5, 0.5, 0.5, 0.5), 2),
    q = rep(c(0, 0, 0.8, 0.8, 0, 0, 0.8, 0.8), 2),
    sn = c(rep(c(2100, 700), 8)),
    sf = c(rep(c(1, 3), 8)),
    Sgen = 0,
    rn = c(rep(c(90, 30), 8)),
    rf = c(rep(c(1, 3), 8)),
    n_wells = 10,
    power_alpha = 0.01,
    sampling = "Binary_Gamma",
    n_transfer = 10,
    n_transfer_selection = 10,
    save_function = F,
    composition_lograte = 1
)


df_pair <- rbindlist(list_df_pair)
df_motif <- rbindlist(list_df_motif)

df_pair <- df_pair %>%
    separate(col = exp_id, into = c("Medium", "temp", "Seed"), sep = "-") %>%
    mutate(Treatment = ifelse(grepl("random", temp), "random_species", "top_down_community")) %>%
    mutate(Community = sub("pair_from_random_species_", "", temp) %>% sub("pair_from_top_down_community_", "", .))

df_motif <- df_motif %>%
    separate(col = exp_id, into = c("Medium", "temp", "Seed"), sep = "-") %>%
    mutate(Treatment = ifelse(grepl("random", temp), "random_species", "top_down_community")) %>%
    mutate(Community = sub("pair_from_random_species_", "", temp) %>% sub("pair_from_top_down_community_", "", .))

# input_synthetic <- fread(paste0("../data/raw/simulation/mapping_files/input_synthetic_simple_medium.csv")) %>%
#     filter(grepl(tt, exp_id)) %>%
#     filter(grepl(ss, exp_id)) %>%
#     mutate(Community = sub(paste0(tt, "-", ss, "_"), "", exp_id) %>% sub("-\\d+$", "", .))

df_pair_aggregated <- df_pair %>%
    group_by(Medium, Treatment, Seed, Community, InteractionType) %>%
    summarize(Count = n()) %>%
    group_by(Medium, Treatment, Seed, Community) %>% mutate(Count = Count / sum(Count)) %>%
    group_by(Medium, Treatment, Seed, InteractionType) %>%
    summarize(MeanCount = mean(Count)) %>%
    left_join(list_treatments)

df_motif_aggregated <- df_motif %>%
    group_by(Medium, Treatment, Seed, Community) %>% mutate(Count = Count / sum(Count)) %>%
    group_by(Medium, Treatment, Seed, Motif) %>%
    summarize(MeanCount = mean(Count)) %>%
    left_join(list_treatments)


p1 <- df_pair_aggregated %>%
    ggplot(aes(x = InteractionType, y = MeanCount), color = 1) +
    geom_boxplot() +
    geom_jitter(shape = 21) +
    facet_grid(Medium ~ Treatment) +
    theme_cowplot()


p2 <-  df_motif_aggregated %>%
    ggplot(aes(x = Motif, y = MeanCount, group = Motif)) +
    geom_boxplot() +
    geom_jitter(shape = 21) +
    scale_x_continuous(breaks = 1:7) +
    guides(color = F) +
    facet_grid(Medium ~ Treatment) +
    theme_cowplot()

# p1
# p2
# ggsave("../plots/pairwise.png", plot = p1, width = 5, height = 10)
# ggsave("../plots/motif.png", plot = p2, width = 7, height = 10)

# df <- fread("~/Downloads/pnas.1712211114.sd02.csv")
# df %>%
#     filter(species.id %in% c("l.conif", "p.rufa.acer", "s.comm")) %>%
#     view

# Side by side
p3 <- df_pair_aggregated %>%
    #filter(Medium %in% paste0("simple_medium", 5:6)) %>%
    #filter(!is.na(InteractionType)) %>%
    ggboxplot(x = "Treatment", y = "MeanCount",
              color = "Treatment", palette = "jco",
              add = "jitter",
              facet.by = c("sf", "InteractionType"), short.panel.labs = FALSE) +
    theme(axis.text.x = element_blank())  +
    stat_compare_means(label = "p.format")

p4 <- df_motif_aggregated %>%
    #filter(Medium %in% paste0("simple_medium", 5:6)) %>%
    ggboxplot(x = "Treatment", y = "MeanCount",
              color = "Treatment", palette = "jco",
              add = "jitter",
              facet.by = c("sf", "Motif"), short.panel.labs = FALSE) +
    theme(axis.text.x = element_blank())  +
    stat_compare_means(label = "p.format")

ggsave("../plots/pairwise_rich.png", plot = p3, width = 7, height = 7)
ggsave("../plots/motif_rich.png", plot = p4, width = 7, height = 7)



# p1 <- df_pair_outcome %>%
#     left_join(input_synthetic) %>%
#     #filter(Community %in% 1:10, seed %in% 1:20) %>%
#     plot_pair_outcome(normalize_sum = T)
#
# p2 <- df_motif_count %>%
#     left_join(input_synthetic) %>%
#     #filter(Community %in% 1:3, seed %in% 1:6) %>%
#     plot_motif_count(normalize_sum = T)

#p <- plot_grid(p1, p2, align = "h", axis = "lr", nrow = 1)
#p + ggtitle(paste0(tt, "-", ss))



if(FALSE) {


    # Top-down communities
    df_pair_data <- read_pair_data("../data/raw/simulation/", paste0(tt, "-pair_from_top_down_community"))
    df_pair_list <- read_pair_list("../data/raw/simulation/", paste0(tt, "-pair_from_top_down_community"))
    df_pair_data_reshaped <- reshape_pair_data(df_pair_data, df_pair_list)

    df_pair_outcome <- determine_pair_outcome(df_pair_data_reshaped)
    df_motif_count <- df_pair_outcome %>%
        split.data.frame(f = .$exp_id) %>%
        lapply(determine_community_motif) %>%
        rbindlist(idcol = "exp_id")

    p1 <- df_pair_outcome %>%
        left_join(input_synthetic) %>%
        plot_pair_outcome(normalize_sum = T)

    p2 <- df_motif_count %>%
        left_join(input_synthetic) %>%
        plot_motif_count(normalize_sum = T)

    p <- plot_grid(p1, p2, align = "h", axis = "lr", nrow = 1)
    ggsave("../plots/Fig1_top_down.png", width = 10, height = 5)

}

if (FALSE) {
    #
    read_community_list <- function (community_composition) {
        community_composition %>%
            filter(Transfer == max(Transfer), Type == "consumer") %>%
            mutate(Community = paste0("Community", sub("W", "", Well))) %>%
            select(Community, Type, ID, Abundance)
    }
    read_pair_from_commmunity_list <- function(pair_from_community_list) {
        pair_from_community_list %>%
            filter(Type == "consumer") %>%
            select(Well, Community, Pair, InitialFrequency, ID) %>%
            group_by(Well) %>%
            arrange(Well) %>%
            mutate(Isolate = 1:2) %>%
            pivot_wider(names_from = Isolate, names_prefix = "Isolate", values_from = ID) %>%
            ungroup()
    }
    read_pair_from_community_competition <- function(pair_from_community_data, pair_from_community_list) {
        pair_from_community_data %>%
            filter(Type == "consumer") %>%
            left_join(select(pair_from_community_list, Well, Community, Pair, InitialFrequency), by = "Well") %>%
            #left_join(select(pair_from_community_list, Community, Pair, InitialFrequency), by = "Well") %>%
            group_by(Community, Pair, InitialFrequency, Transfer) %>%
            mutate(ID = factor(ID)) %>%
            mutate(TotalAbundance = sum(Abundance), RelativeAbundance = Abundance/TotalAbundance) %>%
            select(Community, Pair, InitialFrequency, Transfer, ID, RelativeAbundance) %>%
            ungroup()
    }
    determine_pair_from_community_outcome <- function(pairs_consumer, pair_list) {
        frequency_changes <- pairs_consumer %>%
            left_join(pair_list, by = c("Community", "Pair", "InitialFrequency")) %>%
            filter(ID == Isolate1) %>% select(-ID) %>%
            group_by(Community, Pair) %>%
            pivot_wider(names_from = c(Transfer), names_prefix = "T", values_from = RelativeAbundance) %>%
            mutate(FrequencyChange = ifelse((T5 - T0)>0, "T", "F")) %>%  # TRUE = Isolate1 increases, FALSE = Isolate1 decreases
            select(Community, Pair, Isolate1, Isolate2, FrequencyChange) %>%
            group_by(Community, Pair, Isolate1, Isolate2) %>%
            summarise(FrequencyChangePattern = paste0(FrequencyChange, collapse = "-"))

        map_frequency_pattern <- tibble(FrequencyChangePattern = c("T-T-T", "F-F-F", "T-F-F", "T-T-F", "T-F-T", "F-T-F", "F-F-T", "F-T-T"),
                                        InteractionType = c("exclusion", "exclusion", "coexistence", "coexistence", "coexistence", "coexistence", "mutual exclusion", "mutual exclusion"))

        frequency_changes$From = NA
        frequency_changes$To = NA
        frequency_changes$InteractionType = NA

        for (i in 1:nrow(frequency_changes)) {
            if (frequency_changes$FrequencyChangePattern[i] == "T-T-T") {
                from <- frequency_changes$Isolate1[i]
                to <- frequency_changes$Isolate2[i]
                interaction_type <- "exclusion"
            } else if (frequency_changes$FrequencyChangePattern[i] == "F-F-F") {
                from <- frequency_changes$Isolate2[i]
                to <- frequency_changes$Isolate1[i]
                interaction_type <- "exclusion"
            } else if (frequency_changes$FrequencyChangePattern[i] %in% c("T-F-F", "T-T-F", "T-F-T", "F-T-F")) {
                from <- frequency_changes$Isolate1[i]
                to <- frequency_changes$Isolate2[i]
                interaction_type <- "coexistence"
            } else {
                from <- NA
                to <- NA
                interaction_type <- NA
            }

            frequency_changes$From[i] <- from
            frequency_changes$To[i] <- to
            frequency_changes$InteractionType[i] <- interaction_type
        }

        frequency_changes %>%
            select(Community, Pair, Isolate1, Isolate2, InteractionType, From, To) %>%
            ungroup() %>%
            return()
    }
    determine_community_motif <- function(pairs_from_community) {
        # Make network
        isolate_list <- sort(unique(c(pairs_from_community$Isolate1, pairs_from_community$Isolate2)))
        nodes <- tibble(Isolate = 1:length(isolate_list), ID = isolate_list)
        edges <- pairs_from_community %>%
            mutate(from = match(From, nodes$ID), to = match(To, nodes$ID)) %>%
            select(InteractionType, from, to)
        edges_coext <- edges[edges$InteractionType == "coexistence",]
        edges_coext[,c("from", "to")] <- edges_coext[,c("to", "from")] # Add the mutual edges for coexistence links
        edges <- rbind(edges, edges_coext)

        #
        #if (any(is.na(pairs_from_community$InteractionType))) {
        #    motif_counts <- tibble(Motif = 1:7, Count = rep(NA, 7))
        #} else {
        graph <- tbl_graph(nodes = nodes, edges = filter(edges, !is.na(InteractionType)), directed = T)
        motif_counts <- tibble(Motif = 1:7, Count = count_motif(graph))
        #}

        return(motif_counts)
    }


    df_community_list <- fread(paste0("../data/raw/simulation/simple_medium-top_down_community-", i, "_composition.txt")) %>%
        read_community_list()

    # Pairs from trios
    df_pair_from_community_list <- fread(paste0("../data/raw/simulation/simple_medium-pair_from_top_down_community_", comm, "-", seed1, "_composition", ".txt")) %>%
        read_pair_from_commmunity_list()

    df_pair_from_community_competition <- fread(paste0("../data/raw/simulation/pair-from_top_down_community-", comm_seed, "-community", comms[j], "_composition.txt")) %>%
        read_pair_from_community_competition(df_pair_from_community_list)

    df_pair_from_community_outcome <- determine_pair_from_community_outcome(df_pair_from_community_competition, df_pair_from_community_list)

    df_community_motif <- determine_community_motif(df_pair_from_community_outcome) %>%
        mutate(Community = paste0("Community", comms[j])) %>%
        filter(Count != 0) %>%
        select(Community, Motif, Count)


    fread("../data/raw/simulation/simple_medium-pair_from_random_species_1-1.txt") %>%
        read_pair_from_community_competition(pair_from_community_list)









}

