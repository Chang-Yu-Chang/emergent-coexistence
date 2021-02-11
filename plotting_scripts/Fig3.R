# Simulation. Top-down and bottom-up networks

library(cowplot)
library(tidygraph)
library(ggraph)
library(tidyverse)
library(data.table)
source("network_functions.R")
source("analysis-pair-culturable_random.R")
source("analysis-trio-culturable_random.R")
source("analysis-pair-from_top_down_community.R")

input_independent <- fread("../data/raw/simulation/mapping_files/input_independent.csv")
input_independent_community <- input_independent %>% filter(grepl("community-top_down", exp_id)) %>%
    filter(seed == 2)
interaction_color <- assign_interaction_color()

# Panel A: top-down assembly cartoon
p1 <- ggdraw() + draw_image("../plots/cartoons/Fig2A.png")

# Panel B: motif count of top-down assembled communities, relative to trios in pool
community_motif_list <- rep(list(NA), nrow(input_independent_community))
community_pair_list <- rep(list(NA), nrow(input_independent_community))
names(community_motif_list) <- input_independent_community$seed
names(community_pair_list) <- input_independent_community$seed
for (i in 1:nrow(input_independent_community)) {
    cat("\nexp_id =", input_independent_community$exp_id[i])
    cat("\tseed =", input_independent_community$seed[i])
    comm_seed <- input_independent_community$seed[i]
    input_independent_pair_from_comm <- input_independent %>% filter(grepl("pair-from_top_down", exp_id), seed == comm_seed)
    comms <- input_independent_pair_from_comm %>%
        filter(seed == input_independent_community$seed[i]) %>%
        pull(exp_id) %>%
        gsub(paste0("pair-from_top_down_community-", comm_seed, "-community"), "", .)
    #comms <- c(1:4, 6:10)
    n_comms <- length(comms)

    temp_list <- rep(list(NA), n_comms)
    temp_list2 <- rep(list(NA), n_comms)

    for (j in 1:length(comms)) {
        cat("\ncommunity", comms[j])
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

# p2 <- community_motif %>%
#     mutate(Community = gsub("Community", "", Community)) %>%
#     group_by(Seed, Motif) %>%
#     summarize(Count = sum(Count)) %>%
#     ggplot() +
#     #geom_bar(aes(x = Motif, y = Count, fill = Community), stat = "identity", color = 1) +
#     geom_bar(aes(x = Motif, y = Count, group = Motif), stat = "identity", fill = NA, color = 1) +
#     scale_x_continuous(limits = c(0,8), breaks = 1:7, expand = c(0,0)) +
#     facet_wrap(Seed~., scales = "free", ncol = 1) +
#     theme_cowplot() +
#     theme(legend.position = "right", strip.background = element_blank(), strip.text = element_blank())
# # p <- plot_grid(p1, p2, nrow = 1, rel_widths = c(4, 4), axis = "tblr", align = "h")
# ggsave("../plots/Fig2B.png", plot = p2, width = 4, height = 4)

community_motif_count_relative <- community_motif %>% group_by(Motif) %>% summarize(Count = sum(Count)) %>% ungroup() %>% mutate(RelativeCount = Count / sum(Count))
motif_counts <- df_trio_motif_aggregate %>% group_by(Motif) %>% summarize(Count = n())
my_colors <- c("species pool" = "black", "top-down" = "red")
p2 <- df_trio_motif_aggregate %>%
    group_by(Motif) %>%
    summarize(Count = n()) %>%
    ungroup() %>%
    mutate(RelativeCount = Count / sum(Count)) %>%
    ggplot() +
    geom_point(aes(x = Motif, y = RelativeCount, color = "species pool"), size = 3, shape = 21) +
    geom_point(data = community_motif_count_relative, aes(x = Motif, y = RelativeCount, color = "top-down"), size = 3, shape = 21) +
    scale_x_continuous(limits = c(0,8), breaks = 1:7, expand = c(0,0)) +
    scale_color_manual(values = my_colors) +
    theme_bw() +
    theme(legend.position = "top", legend.title = element_blank()) +
    labs(x = "Motif", y = "Relative count")
ggsave("../plots/Fig2B.png", plot = p2, width = 4, height = 4)


# Panel C: overrepresentation versus trio coexistence
df_community_motif_counts <- community_motif %>%
    group_by(Seed, Motif) %>%
    summarise(Count = sum(Count)) %>%
    group_by(Seed) %>%
    mutate(RelativeMotifCount = Count / sum(Count)) %>%
    mutate(Experiment = "Top-down community")

motif_count_trio_coexistence <- tibble(Seed = 2,
    Motif = 2:7,
    RelativeMotifCountDiff = df_community_motif_counts$RelativeMotifCount - df_trio_motif_counts$RelativeMotifCount,
    FractionTrioCoexistence = df_trio_motif_counts$FractionTrioCoexistence)

p3 <- motif_count_trio_coexistence %>%
    ggplot() +
    geom_point(aes(x = RelativeMotifCountDiff, y = FractionTrioCoexistence), shape = 21, size = 3) +
    geom_text(aes(x = RelativeMotifCountDiff, y = FractionTrioCoexistence, label = Motif), vjust = 1, hjust = -1) +
    geom_vline(xintercept = 0, linetype = 2) +
    scale_y_continuous(limits = c(0,1)) +
    theme_cowplot() +
    labs(x = "RelativeMotifCount[community] - RelativeMotifCount[trio]", y = "Fraction of trio coexistence")

ggsave("../plots/Fig2C.png", plot = p3, width = 6, height = 4)


if (FALSE) {
    # Random pair outcome
    temp <- df_pair_from_trio_outcome %>%
        group_by(InteractionType) %>%
        summarize(Count = n()) %>%
        mutate(Experiment = "Species\npool")

    # Pairs in top-down network
    ps1 <- community_pair %>%
        mutate(InteractionType = ifelse(is.na(InteractionType), "no-growth", InteractionType)) %>%
        group_by(Seed, InteractionType) %>%
        summarise(Count = n()) %>% ungroup() %>%
        select(-Seed) %>%
        mutate(Experiment = "Top-down\ncommunities") %>%
        bind_rows(temp) %>%
        group_by(Experiment) %>%
        mutate(RelativeCount = Count/sum(Count)) %>%
        ggplot() +
        geom_bar(aes(x = Experiment, y = RelativeCount, fill = InteractionType), stat = "identity", color = 1) +
        #scale_x_discrete(expand = c(0,0)) +
        scale_y_continuous(expand = c(0,0)) +
        scale_fill_manual(values = interaction_color) +
        #facet_wrap(Seed~., scales = "free", ncol = 1) +
        theme_cowplot() +
        theme(legend.position = "right", legend.title = element_blank(), strip.background = element_blank(), strip.text = element_blank()) +
        panel_border(color = 1) +
        labs(x = "", y = "Fraction of pairs")

    ggsave("../plots/FigS3.png", plot = ps1, width = 5, height = 4)


}













