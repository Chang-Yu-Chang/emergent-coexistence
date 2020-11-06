# Simulation. Random trios in speices pool

library(cowplot)
library(tidygraph)
library(ggraph)
library(tidyverse)
library(data.table)
source("network_functions.R")
source("analysis-pair-culturable_random.R")
source("analysis-trio-culturable_random.R")
source("analysis-pair-from_top_down_community.R")

# Read data
input_independent <- fread("../data/raw/simulation/mapping_files/input_independent.csv")
input_independent_trios <- input_independent %>% filter(grepl("trio-culturable_isolates", exp_id)) %>%
    filter(seed == 2)
interaction_color <- assign_interaction_color()

# Panel A: cartoon and temporal dynamics of pairs and trios
df_trio_competition <- fread(paste0("../data/raw/simulation/trio-culturable_isolates-2_composition.txt"))
df_trio_competition %>%
    filter(Well == "W0", Type == "consumer") %>% 
    mutate(ID = factor(ID)) %>% 
    group_by(Well, Transfer) %>% 
    mutate(RelativeAbundance = Abundance / sum(Abundance)) %>% 
    ggplot() +
    #geom_line() + geom_point() +
    geom_bar(aes(x = Transfer, y = RelativeAbundance, fill = ID, group = ID), stat = "identity") +
    scale_x_continuous(breaks = 1:10, expand = c(0,0)) +
    scale_y_continuous(expand = c(0,0)) +
    theme_bw() 

temp_list <- rep(list(NA), nrow(input_independent_pairs))
for (i in 1:nrow(input_independent_pairs)) {
    cat("\nexp_id = ", input_independent_pairs$exp_id[i])
    cat(",\tseed = ", input_independent_pairs$seed[i])
    current_seed <- input_independent_trios$seed[i]
    df_pair_list <- fread(paste0("../data/raw/simulation/pair-culturable-", current_seed, ".txt")) %>% 
        read_pair_list()
    
    df_pair_competition <- 
        fread(paste0("../data/raw/simulation/pair-culturable_isolates-", current_seed, "_composition.txt")) %>% 
        read_pair_competition(df_pair_list)
    
    df_pair_outcome <- determine_pair_outcome(df_pair_competition, df_pair_list)
    
    temp_list[[i]] <- df_pair_outcome
}
df_pair_outcomes <- bind_rows(temp_list, .id = "Seed")

p1 <- ggdraw() + draw_image("../plots/cartoons/Fig1A.png")

# Panel B: trio motif
temp_list <- rep(list(NA), nrow(input_independent_trios))
for (i in 1:nrow(input_independent_trios)) {
    cat("\nexp_id = ", input_independent_trios$exp_id[i])
    cat(",\tseed = ", input_independent_trios$seed[i])
    current_seed <- input_independent_trios$seed[i]
    
    # Trio
    df_trio_list <- fread(paste0("../data/raw/simulation/trio-culturable-", current_seed, ".txt")) %>% 
        read_trio_list()
    df_trio_competition <- fread(paste0("../data/raw/simulation/trio-culturable_isolates-", current_seed, "_composition.txt")) %>% 
        read_trio_competition(df_trio_list)
    df_trio_outcome <- determine_trio_outcome(df_trio_competition)
    
    # Pairs from trios
    df_pair_from_trio_list <- fread(paste0("../data/raw/simulation/pair-culturable_from_trio-", current_seed, ".txt")) %>% 
        read_pair_from_trio_list()
    df_pair_from_trio_competition <- fread(paste0("../data/raw/simulation/pair-culturable_from_trio-", current_seed, "_composition.txt")) %>% 
        read_pair_from_trio_competition(df_pair_from_trio_list)
    df_pair_from_trio_outcome <- determine_pair_from_trio_outcome(df_pair_from_trio_competition, df_pair_from_trio_list)
    df_trio_motif <- df_pair_from_trio_outcome %>% 
        split.data.frame(f=.$Trio) %>% 
        lapply(determine_trio_motif) %>% 
        bind_rows(.id = "Trio") %>% 
        filter(Count != 0) %>% 
        select(Trio, Motif)
    
    temp_list[[i]] <- df_trio_motif
    
}
df_trio_motif_aggregate <- bind_rows(temp_list, .id = "Seed") %>% 
    left_join(df_trio_outcome) %>% 
    mutate(Coexistence = ifelse(Richness == 3, "trio coexists", "trio does not coexist"))
trio_counts <- df_trio_motif_aggregate %>%
    group_by(Seed) %>% summarize(Count = n())

p2 <- df_trio_motif_aggregate %>% 
    filter(!is.na(Richness)) %>%
    ggplot() +
    geom_bar(aes(x = Motif, fill = Coexistence), stat = "count", color = 1) +
    geom_text(data = trio_counts, aes(label = paste0("n=", Count)), x = Inf, y = Inf, hjust = 1, vjust = 2) +
    scale_x_continuous(limits = c(0,8), breaks = 1:7, expand = c(0,0)) +
    scale_fill_manual(values = c("trio coexists" = "#557BAA", "trio does not coexist" = "#DB7469")) +
    facet_wrap(Seed~., scales = "free", ncol = 1) +
    theme_cowplot() +
    theme(legend.position = c(.4, .8), strip.background = element_blank(), strip.text = element_blank()) +
    guides(fill = guide_legend(title = ""))

p <- plot_grid(p1, p2, nrow = 1, rel_widths = c(4, 4), axis = "tblr", align = "h")
ggsave("../plots/Fig1.png", plot = p, width = 8, height = 4)


if (FALSE) {
    ps1 <- df_trio_motif_aggregate %>% 
        filter(!is.na(Richness)) %>%
        ggplot() +
        geom_bar(aes(x = Motif), stat = "count", fill = NA, color = 1) +
        geom_text(data = trio_counts, aes(label = paste0("n=", Count)), x = Inf, y = Inf, hjust = 1, vjust = 2) +
        scale_x_continuous(limits = c(0,8), breaks = 1:7, expand = c(0,0)) +
        scale_fill_manual(values = c("trio coexists" = "#557BAA", "trio does not coexist" = "#DB7469")) +
        facet_wrap(Seed~., scales = "free", ncol = 1) +
        theme_cowplot() +
        theme(legend.position = c(.4, .8), strip.background = element_blank(), strip.text = element_blank()) +
        guides(fill = guide_legend(title = ""))
    ggsave("../plots/FigS1.png", plot = ps1, width = 4, height = 4)
}


#----





# Panel A: random pairs
temp_list <- rep(list(NA), nrow(input_independent_pairs))
for (i in 1:nrow(input_independent_pairs)) {
    cat("\nexp_id = ", input_independent_pairs$exp_id[i])
    cat(",\tseed = ", input_independent_pairs$seed[i])
    df_pair_list <- fread(paste0("../data/raw/simulation/pair-culturable-", i, ".txt")) %>% 
        read_pair_list()
    
    df_pair_competition <- 
        fread(paste0("../data/raw/simulation/pair-culturable_isolates-", i, "_composition.txt")) %>% 
        read_pair_competition(df_pair_list)
    
    df_pair_outcome <- determine_pair_outcome(df_pair_competition, df_pair_list)
    
    temp_list[[i]] <- df_pair_outcome
}
df_pair_outcomes <- bind_rows(temp_list, .id = "Seed")

p1 <- df_pair_outcomes %>% 
    mutate(InteractionType = ifelse(is.na(InteractionType), "no-growth", InteractionType)) %>% 
    group_by(Seed, InteractionType) %>% 
    summarise(Count = n()) %>% 
    #filter(Seed %in% c(1,3)) %>% 
    ggplot() +
    geom_bar(aes(x = Seed, y = Count, fill = InteractionType), stat = "identity", color = 1) +
    scale_x_discrete(expand = c(0,0)) +
    scale_y_continuous(expand = c(0,0)) +
    scale_fill_manual(values = interaction_color) +
    facet_wrap(Seed~., scales = "free_x", ncol = 1) +
    theme_cowplot() +
    theme(legend.position = "right", legend.title = element_blank(), strip.background = element_blank(), strip.text = element_blank()) +
    panel_border(color = 1) +
    ggtitle("Random culturable pairs")
p1
ggsave("../plots/Fig2A.png", plot = p1, width = 5, height = 5)

if(FALSE) {
    p1 <- df_pair_outcomes %>% 
        group_by(Seed, InteractionType) %>% 
        summarise(Count = n()) %>% 
        ggplot() +
        geom_bar(aes(x = Seed, y = Count, fill = InteractionType), stat = "identity") +
        scale_x_discrete(expand = c(0,0)) +
        scale_y_continuous(expand = c(0,0)) +
        theme_cowplot() +
        theme(legend.position = "top", legend.title = element_blank()) +
        panel_border(color = 1) +
        ggtitle("Random pairs of culturable isolates")
}


# Panel B: random trios
temp_list <- rep(list(NA), nrow(input_independent_trios))
for (i in 1:nrow(input_independent_trios)) {
    cat("\nexp_id = ", input_independent_trios$exp_id[i])
    cat(",\tseed = ", input_independent_trios$seed[i])
    
    # Trio
    df_trio_list <- fread(paste0("../data/raw/simulation/trio-culturable-", i, ".txt")) %>% 
        read_trio_list()
    df_trio_competition <- fread(paste0("../data/raw/simulation/trio-culturable_isolates-", i, "_composition.txt")) %>% 
        read_trio_competition(df_trio_list)
    df_trio_outcome <- determine_trio_outcome(df_trio_competition)
    
    # Pairs from trios
    df_pair_from_trio_list <- fread(paste0("../data/raw/simulation/pair-culturable_from_trio-", i, ".txt")) %>% 
        read_pair_from_trio_list()
    df_pair_from_trio_competition <- fread(paste0("../data/raw/simulation/pair-culturable_from_trio-", i, "_composition.txt")) %>% 
        read_pair_from_trio_competition(df_pair_from_trio_list)
    df_pair_from_trio_outcome <- determine_pair_from_trio_outcome(df_pair_from_trio_competition, df_pair_from_trio_list)
    df_trio_motif <- df_pair_from_trio_outcome %>% 
        split.data.frame(f=.$Trio) %>% 
        lapply(determine_trio_motif) %>% 
        bind_rows(.id = "Trio") %>% 
        filter(Count != 0) %>% 
        select(Trio, Motif)
    
    temp_list[[i]] <- df_trio_motif
    
}
df_trio_motif_aggregate <- bind_rows(temp_list, .id = "Seed") %>% 
    left_join(df_trio_outcome) %>% 
    mutate(Coexistence = ifelse(Richness == 3, "trio coexists", "trio does not coexist"))
trio_counts <- df_trio_motif_aggregate %>%
    group_by(Seed) %>% summarize(Count = n())

p2 <- df_trio_motif_aggregate %>% 
    filter(!is.na(Richness)) %>%
    ggplot() +
    geom_bar(aes(x = Motif, fill = Coexistence), stat = "count", color = 1) +
    geom_text(data = trio_counts, aes(label = paste0("n=", Count)), x = -Inf, y = Inf, hjust = -1, vjust = 2) +
    scale_x_continuous(limits = c(0,8), breaks = 1:7, expand = c(0,0)) +
    scale_fill_manual(values = c("trio coexists" = "#557BAA", "trio does not coexist" = "#DB7469")) +
    facet_wrap(Seed~., scales = "free", ncol = 1) +
    theme_cowplot() +
    theme(legend.position = c(.05, .9), strip.background = element_blank(), strip.text = element_blank()) +
    guides(fill = guide_legend(title = "")) +
    ggtitle("Random culturable trios")
p2
#p2 <- plot_grid(p_trio_motifs, p_trio_motifs_coexist, nrow = 1, axis = "tb", align = "vh")

ggsave("../plots/Fig2B.png", plot = p2, width = 5, height = 5)

if (FALSE) {
    
    p_trio_motifs <- df_trio_motif_aggregate %>%
        ggplot() +
        geom_bar(aes(x = Motif), stat = "count", color = 1, fill = NA) +
        annotate("text", x = -Inf, y = Inf, label = paste0("n=", nrow(df_trio_motif_aggregate)), hjust = -0.5, vjust = 2) +
        scale_x_continuous(limits = c(0,8), breaks = 1:7, expand = c(0,0)) +
        theme_cowplot() +
        theme(legend.position = "top") +
        
    ps1 <- df_trio_motif_aggregate %>% 
        left_join(df_trio_outcome) %>% 
        ggplot() +
        geom_jitter(aes(x = Motif, y = Richness), shape = 21, width = 0.2, height = 0.2, size = 3) +
        scale_x_continuous(limits = c(1,7), breaks = 1:7) +
        scale_y_continuous(limits = c(0.5,3.5), breaks = 1:3) +
        facet_grid(.~Seed, scales = "free_y") + 
        theme_bw() +
        theme(panel.grid.minor = element_blank())
    ps1
    
#    ggsave("../plots/Fig2B.png", plot = p2, width = 10, height = 4)
}

# Panel C: top-down assembled communities
community_motif_list <- rep(list(NA), nrow(input_independent_community))
community_pair_list <- rep(list(NA), nrow(input_independent_community))
names(community_motif_list) <- input_independent_community$seed
names(community_pair_list) <- input_independent_community$seed
for (i in 1:nrow(input_independent_community)) {
    cat("\nexp_id = ", input_independent_community$exp_id[i])
    cat(",\tseed = ", input_independent_community$seed[i])
    comm_seed <- input_independent_community$seed[i]
    input_independent_pair_from_comm <- input_independent %>% filter(grepl("pair-from_top_down", exp_id), seed == comm_seed)
    n_comms <- nrow(input_independent_pair_from_comm)
#    if (i == 1) n_comms = 9 else n_comms <- 10
    comms <- input_independent_pair_from_comm %>% 
        filter(seed == i) %>% 
        pull(exp_id) %>% 
        gsub(paste0("pair-from_top_down_community-", comm_seed, "-community"), "", .)
    
    temp_list <- rep(list(NA), n_comms)
    temp_list2 <- rep(list(NA), n_comms)
    
    for (j in 1:length(comms)) {
#        if (i == 1 & j >=10) next
        cat("\n community = Community", comms[j])
        # Community
        df_community_list <- fread(paste0("../data/raw/simulation/community-top_down-", i, "_composition.txt")) %>% 
            read_community_list()
        
        # Pairs from trios
        df_pair_from_community_list <- fread(paste0("../data/raw/simulation/pair-from_top_down_community-", comm_seed, "-community", comms[j], ".txt")) %>% 
            read_pair_from_commmunity_list()
        
        df_pair_from_community_competition <- fread(paste0("../data/raw/simulation/pair-from_top_down_community-", comm_seed, "-community", comms[j], "_composition.txt")) %>% 
            read_pair_from_community_competition(df_pair_from_community_list)
        
        df_pair_from_community_outcome <- determine_pair_from_community_outcome(df_pair_from_community_competition, df_pair_from_community_list)
        
        df_community_motif <- determine_community_motif(df_pair_from_community_outcome) %>% 
            mutate(Community = paste0("Community", comms[j])) %>% 
            filter(Count != 0) %>% 
            select(Community, Motif, Count)
        
        temp_list[[j]] <- df_community_motif
        temp_list2[[j]] <- df_pair_from_community_outcome
    }
    community_motif_list[[i]] <- rbindlist(temp_list)
    community_pair_list[[i]] <- rbindlist(temp_list2)
}

community_pair <- rbindlist(community_pair_list, idcol = "Seed")
community_motif <- rbindlist(community_motif_list, idcol = "Seed")


p3 <- community_pair %>% 
    mutate(InteractionType = ifelse(is.na(InteractionType), "no-growth", InteractionType)) %>% 
    group_by(Seed, InteractionType) %>% 
    summarise(Count = n()) %>% 
    ggplot() +
    geom_bar(aes(x = Seed, y = Count, fill = InteractionType), stat = "identity", color = 1) +
    scale_x_discrete(expand = c(0,0)) +
    scale_y_continuous(expand = c(0,0)) +
    scale_fill_manual(values = interaction_color) +
    facet_wrap(Seed~., scales = "free", ncol = 1) +
    theme_cowplot() + 
    theme(legend.position = "right", legend.title = element_blank(), strip.background = element_blank(), strip.text = element_blank()) +
    panel_border(color = 1) +
    ggtitle("Pairs in top-down communities")

p4 <- community_motif %>% 
    mutate(Community = gsub("Community", "", Community)) %>% 
    ggplot() +
    geom_bar(aes(x = Motif, y = Count, fill = Community), stat = "identity", color = 1) +
    scale_x_continuous(limits = c(0,8), breaks = 1:7, expand = c(0,0)) +
    facet_wrap(Seed~., scales = "free", ncol = 1) +
    theme_cowplot() +
    theme(legend.position = "right", strip.background = element_blank(), strip.text = element_blank()) +
    ggtitle("")

p_random <- plot_grid(p1, p2, nrow = 1, rel_widths = c(4, 8), axis = "tblr", align = "h")
p_community <- plot_grid(p3, p4, nrow = 1, rel_widths = c(4, 8), axis = "tblr", align = "h")
p <- plot_grid(p_random, p_community, nrow = 1, axis = "tb")
ggsave("../plots/Fig2C.png", plot = p_random, width = 8, height = 8)
ggsave("../plots/Fig2D.png", plot = p_community, width = 8, height = 8)
ggsave("../plots/Fig2.png", plot = p, width = 16, height = 8)











seed_subset <- c(3)
n_seeds <- length(seed_subset)

p1 <- df_pair_outcomes %>% 
    mutate(InteractionType = ifelse(is.na(InteractionType), "no-growth", InteractionType)) %>% 
    filter(Seed %in% seed_subset) %>% 
    group_by(Seed, InteractionType) %>% 
    summarise(Count = n()) %>% 
    mutate(Medium = ifelse(Seed == 1, "rich medium", "single supplied resource")) %>% 
    #filter(Seed %in% c(1,3)) %>% 
    ggplot() +
    geom_bar(aes(x = Seed, y = Count, fill = InteractionType), stat = "identity", color = 1) +
    scale_x_discrete(expand = c(0,0)) +
    scale_y_continuous(expand = c(0,0)) +
    scale_fill_manual(values = interaction_color) +
    facet_wrap(Medium~., scales = "free_x", ncol = 1) +
    theme_cowplot() +
    theme(legend.position = "right", legend.title = element_blank()) +
    panel_border(color = 1) +
    ggtitle("Random culturable pairs")


p2 <- df_trio_motif_aggregate %>% 
    filter(!is.na(Richness)) %>%
    filter(Seed %in% seed_subset) %>% 
    mutate(Medium = ifelse(Seed == 1, "rich medium", "single supplied resource")) %>% 
    ggplot() +
    geom_bar(aes(x = Motif, fill = Coexistence), stat = "count", color = 1) +
    geom_text(data = mutate(filter(trio_counts, Seed %in% seed_subset), Medium = ifelse(Seed == 1, "rich medium", "single supplied resource")), aes(label = paste0("n=", Count)), x = -Inf, y = Inf, hjust = -1, vjust = 2) +
    scale_x_continuous(limits = c(0,8), breaks = 1:7, expand = c(0,0)) +
    scale_fill_manual(values = c("trio coexists" = "#557BAA", "trio does not coexist" = "#DB7469")) +
    facet_wrap(Medium~., scales = "free", ncol = 1) +
    theme_cowplot() +
    theme(legend.position = c(.5, .9)) +
    guides(fill = guide_legend(title = "")) +
    ggtitle("Random culturable trios")

p3 <- community_pair %>% 
    mutate(InteractionType = ifelse(is.na(InteractionType), "no-growth", InteractionType)) %>% 
    filter(Seed %in% seed_subset) %>% 
    group_by(Seed, InteractionType) %>% 
    summarise(Count = n()) %>% 
    mutate(Medium = ifelse(Seed == 1, "rich medium", "single supplied resource")) %>% 
    ggplot() +
    geom_bar(aes(x = Seed, y = Count, fill = InteractionType), stat = "identity", color = 1) +
    scale_x_discrete(expand = c(0,0)) +
    scale_y_continuous(expand = c(0,0)) +
    scale_fill_manual(values = interaction_color) +
    facet_wrap(Medium~., scales = "free", ncol = 1) +
    theme_cowplot() + 
    theme(legend.position = "right", legend.title = element_blank()) +
    panel_border(color = 1) +
    ggtitle("Pairs in top-down communities")

p4 <- community_motif %>% 
    mutate(Community = gsub("Community", "", Community)) %>% 
    filter(Seed %in% seed_subset) %>% 
    mutate(Medium = ifelse(Seed == 1, "rich medium", "single supplied resource")) %>% 
    ggplot() +
    geom_bar(aes(x = Motif, y = Count, fill = Community), stat = "identity", color = 1) +
    scale_x_continuous(limits = c(0,8), breaks = 1:7, expand = c(0,0)) +
    facet_wrap(Medium~., scales = "free", ncol = 1) +
    theme_cowplot() +
    theme(legend.position = "right") +
    ggtitle("")

p_random <- plot_grid(p1, p2, nrow = 1, rel_widths = c(4, 8), axis = "tblr", align = "h")
p_community <- plot_grid(p3, p4, nrow = 1, rel_widths = c(4, 8), axis = "tblr", align = "h")
p <- plot_grid(p_random, p_community, nrow = 1, axis = "tb")
ggsave("../plots/Fig2C_subset.png", plot = p_random, width = 8, height = 4*n_seeds)
ggsave("../plots/Fig2D_subset.png", plot = p_community, width = 8, height = 4*n_seeds)
ggsave("../plots/Fig2_subset.png", plot = p, width = 16, height = 4*n_seeds)

if (FALSE) {
    plot_community_temporal <- function(community_composition) {
        community_composition %>%
            filter(Well %in% paste0("W", 0:5)) %>%
            filter(Type == "consumer") %>% 
            mutate(ID = factor(ID)) %>% 
            group_by(Well, Transfer) %>% 
            mutate(RelativeAbundance = Abundance /sum(Abundance)) %>% 
            ggplot() +
            geom_bar(aes(x = Transfer, y = RelativeAbundance, fill = ID), stat = "identity", color = 1) +
            facet_wrap(Well~.) +
            scale_x_continuous(expand = c(0,0)) +
            scale_y_continuous(expand = c(0,0)) +
            theme_bw() +
            theme(legend.position = "none")
    }
    fread(paste0("../data/raw/simulation/community-top_down-2_composition.txt")) %>% 
        plot_community_temporal()
    
    
    fread(paste0("../data/raw/simulation/pair-from_top_down_community-1-community1_composition.txt")) %>% 
        #    filter(ID == 65) %>% 
        #    filter(Type == "consumer") 
        #filter(Well %in% paste0("W", c(91:92)))
        plot_community_temporal()
    
    fread(paste0("../data/raw/simulation/pair-from_top_down_community-", comm_seed, "-community", comms[j], "_composition.txt")) %>% 
        read_pair_from_community_competition(df_pair_from_community_list)
    
    
    filter(Type == "consumer") %>%
        left_join(select(pair_from_community_list, Well, Community, Pair, InitialFrequency), by = "Well") %>% 
        group_by(Community, Pair, InitialFrequency, Transfer) %>% 
        mutate(ID = factor(ID)) %>% 
        mutate(TotalAbundance = sum(Abundance), RelativeAbundance = Abundance/TotalAbundance) %>% 
        select(Community, Pair, InitialFrequency, Transfer, ID, RelativeAbundance) %>% 
        ungroup()
    
    plot_pair_temporal <- function(pair_composition) {
        pair_composition %>%
            #filter(Well %in% paste0("W", 0:100)) %>%
            #filter(Type == "consumer") %>% 
            #mutate(ID = factor(ID)) %>% 
            group_by(Community, Pair, InitialFrequency, Transfer) %>% 
            #mutate(RelativeAbundance = Abundance /sum(Abundance)) %>% 
            ggplot() +
            geom_bar(aes(x = Transfer, y = RelativeAbundance, fill = ID), stat = "identity", color = 1) +
            facet_grid(InitialFrequency~Pair) +
            scale_x_continuous(expand = c(0,0)) +
            scale_y_continuous(expand = c(0,0)) +
            theme_bw() +
            theme(legend.position = "none")
    }
    
    
    df_pair_from_community_competition %>% 
        filter(Community == 10) %>% 
        filter(Pair %in% paste0("Pair", 1:10)) %>% 
        plot_pair_temporal()
    
    
    fread("../data/raw/simulation/community-top_down-2_composition.txt") %>% 
        plot_community_temporal()
    
    
    # Panel XX: motif distribution, compared to randomized network
    b = 100
    
    cat("\n Randomizing the empirical graphs")
    temp_list <- rep(list(rep(list(NA), b)), length(graph_list))
    names(temp_list) <- names(graph_list)
    for (j in 1:length(graph_list)){
        cat("\ngraph:", names(graph_list)[j], "\n")
        for (i in 1:b) {
            temp_list[[j]][[i]] <- count_motif(randomize_network(graph_list[[j]]))
            if (i%%10 == 0) cat(i, " ")
        }
    }
    motif_counts <- temp_list %>%
        lapply(function(x) {
            lapply(x, function(y) {tibble(Motif = factor(1:7), Count = y)}) %>%
                rbindlist(idcol = "Seed")
        }) %>% 
        bind_rows(.id = "Community")
    
    motif_counts_p95 <- motif_counts %>% 
        group_by(Community, Motif) %>% 
        filter(Count >= quantile(Count, 0.95)) %>% 
        distinct(Community, Motif, Count) %>% 
        arrange(Community, Motif, Count) %>% 
        slice_min(Count) %>% 
        mutate(Percentile = "p95")
    motif_counts_p05 <- motif_counts %>% 
        group_by(Community, Motif) %>% 
        filter(Count <= quantile(Count, 0.05)) %>% 
        distinct(Community, Motif, Count) %>% 
        arrange(Community, Motif, Count) %>% 
        slice_max(Count) %>% 
        mutate(Percentile = "p05")
    motif_counts_percentile <- bind_rows(motif_counts_p05, motif_counts_p95) %>% 
        mutate(Community = ordered(Community, levels = community_names_ordered_by_size))
    
    
    colors <- c("observed" = "red", "random [5th and 95th percentiles]" = "black")
    
    p3 <- summary_network_motifs %>% 
        mutate(Community = ordered(Community, levels = community_names_ordered_by_size)) %>% 
        ggplot(aes(x = Motif, y = Count)) +
        geom_point(data = motif_counts_percentile, aes(x = Motif, y = Count, color = "random [5th and 95th percentiles]")) +
        geom_segment(data = pivot_wider(motif_counts_percentile, names_from = Percentile, values_from = Count), 
            aes(x = Motif, xend = Motif, y = p05, yend = p95, color = "random [5th and 95th percentiles]")) +
        geom_point(aes(color = "observed")) +
        scale_color_manual(values = colors) +
        facet_wrap(Community~., scales = "free_y", nrow = 2) +
        theme_cowplot() + 
        theme(legend.position = "bottom") +
        panel_border(color = "black") +
        labs(color = "")
    p3
    
    ggsave("../plots/Fig2C.png", plot = p3, width = 14, height = 4)
    
    
    # Panel XX: pooled networks
    motif_counts_aggregated <- motif_counts %>% 
        group_by(Seed, Motif) %>% 
        summarize(Count = sum(Count))
    
    motif_counts_aggregated_p95 <- motif_counts_aggregated %>% 
        group_by(Motif) %>% 
        filter(Count >= quantile(Count, 0.95)) %>% 
        distinct(Motif, Count) %>% 
        arrange(Motif, Count) %>% 
        slice_min(Count) %>% 
        mutate(Percentile = "p95")
    motif_counts_aggregated_p05 <- motif_counts_aggregated %>% 
        group_by(Motif) %>% 
        filter(Count <= quantile(Count, 0.05)) %>% 
        distinct(Motif, Count) %>% 
        arrange(Motif, Count) %>% 
        slice_min(Count) %>% 
        mutate(Percentile = "p05")
    motif_counts_aggregated_percentile <- bind_rows(motif_counts_aggregated_p05, motif_counts_aggregated_p95) 
    
    summary_network_motifs_aggregated <- summary_network_motifs %>% 
        group_by(Motif) %>% 
        summarize(Count = sum(Count))
    
    
    plot_example_motifs <- function(node_size=5) {
        temp_list <- rep(list(NA), 7)
        temp_id <- c(11, 7, 8, 12, 13, 14, 15)
        
        for (i in 1:7) {
            g <- as_tbl_graph(igraph::graph.isocreate(size = 3, temp_id[i]))
            layout <- create_layout(g, layout = 'circle')
            g <- activate(g, edges) %>% mutate(InteractionType = ifelse(edge_is_mutual(), "coexistence", "exclusion"))
            temp_list[[i]] <- g %>% activate(nodes) %>% mutate(x = layout$x, y = layout$y, graph = paste0(i))
        }
        
        merged_graph <- bind_graphs(temp_list)
        
        plot_competitive_network(merged_graph, layout = "example_motif", node_size = node_size) + 
            facet_nodes(~graph, nrow = 1)
    }
    p_example <- plot_example_motifs(node_size = 3)
    
    
    colors <- c("observed" = "red", "random [5th and 95th percentiles]" = "black")
    p4 <- summary_network_motifs_aggregated %>% 
        ggplot(aes(x = Motif, y = Count)) +
        geom_point(data = motif_counts_aggregated_percentile, aes(x = Motif, y = Count, color = "random [5th and 95th percentiles]"), size = 3) +
        geom_segment(data = pivot_wider(motif_counts_aggregated_percentile, names_from = Percentile, values_from = Count), 
            aes(x = Motif, xend = Motif, y = p05, yend = p95, color = "random [5th and 95th percentiles]")) +
        geom_point(aes(color = "observed"), size = 3) +
        scale_color_manual(values = colors) +
        facet_grid(.~Motif, scales = "free_x") +
        theme_cowplot() + 
        theme(legend.position = "bottom", strip.text = element_blank(), strip.background = element_blank(),
            axis.text.x = element_blank()) +
        labs(color = "")
    
    p <- plot_grid(p_example, p4, ncol = 1, axis = "rl", align = "vh", rel_heights = c(2, 7))
    p
    ggsave("../plots/Fig2D.png", plot = p, width = 8, height = 5)
    
    
    
    # Combining the plots
    p <- plot_grid(p2, p3, ncol = 1, align = "v", rel_heights = c(1, 2))
    
    ggsave("../plots/Fig2.png", plot = p, width = 10, height = 5)
    
    
    
    
    
    
    graph_list[[11]] %>%
        plot_competitive_network()
    
    # BArplot
    motif_counts %>%
        mutate(Community = ordered(Community, level = community_names)) %>% 
        group_by(Community, Seed) %>% 
        mutate(TotalMotifCount = sum(Count)) %>% 
        group_by(Community, Seed, Motif) %>% 
        summarize(RelativeMotifCount = Count/TotalMotifCount) %>% 
        ggplot(aes(x = Seed, y = RelativeMotifCount, fill = Motif)) +
        geom_bar(position = "stack", stat = "identity") +
        scale_x_continuous(expand = c(0,0)) +
        scale_y_continuous(expand = c(0,0)) +
        facet_grid(Community~.) + 
        theme_bw()
    
    
}
 #










