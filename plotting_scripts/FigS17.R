library(tidyverse)
library(cowplot)
source(here::here("analysis/00-metadata.R"))
source(here::here("simulation/01-generate_input.R"))

# Model parameters
input_parameters <- read_csv(here::here("simulation/01-input_parameters.csv"), col_types = cols())
input_poolPairs <- read_csv(here::here("simulation/03a-input_poolPairs.csv"), col_types = cols())
input_withinCommunityPairs <- read_csv(here::here("simulation/03b-input_withinCommunityPairs.csv"), col_types = cols())

# Community pairs
withinCommunityPairs_N_freq <- read_csv(paste0(folder_simulation, "aggregated/13-withinCommunityPairs_N_freq.csv"), col_types = cols()) %>%
    mutate(Time = factor(Time, c("init", "end"))) %>%
    mutate(Pair = factor(Pair, paste0("P", 1:1000))) %>%
    mutate(InitialFrequency = factor(InitialFrequency, c(5, 50, 95)))

communities_richness <- read_csv(paste0(folder_simulation, "aggregated/12-communities_richness.csv"), col_types = cols()) %>%
    mutate(Community = factor(Community, paste0("W", 0:(nrow(input_withinCommunityPairs)-1)))) %>%
    mutate(PairSize = choose(Richness, 2)) %>%
    arrange(desc(Richness)) %>%
    slice(1:20) %>%
    mutate(CommunityLabel = 1:20)
withinCommunityPairs_N_outcome <- read_csv(paste0(folder_simulation, "aggregated/13-withinCommunityPairs_N_outcome.csv"), col_types = cols()) %>%
    mutate(Community = factor(Community, paste0("W", 0:(nrow(input_withinCommunityPairs)-1)))) %>%
    # Top 20 richness
    filter(Community %in% communities_richness$Community)


p1 <- withinCommunityPairs_N_outcome %>%
    group_by(Community, InteractionType, .drop = F) %>%
    summarize(Count = n()) %>%
    arrange(Community, InteractionType) %>%
    mutate(Fraction = Count / sum(Count)) %>%
    left_join(communities_richness) %>%
    ggplot() +
    geom_col(aes(x = CommunityLabel, y = Fraction, fill = InteractionType), color = 1, width = 0.8, linewidth = 0.5) +
    #geom_text(aes(x = Community, label = Richness), y = 0.9) +
    annotate("text", x = 1:20, y = 1.15, label = communities_richness$Richness, size = 4) +
    annotate("text", x = 21, y = 1.15, label = c("n. of species"), size = 4, hjust = 0) +
    annotate("segment", x = .5, xend = 25, y = 1.1, yend = 1.1, color = "black") +
    annotate("text", x = 21, y = 1.05, label = c("n. of tested pairs"), size = 4, hjust = 0) +
    annotate("text", x = 1:20, y = 1.05, label = communities_richness$PairSize, size = 4) +
    scale_fill_manual(values = interaction_color) +
    scale_x_continuous(breaks = 1:20, expand = c(0.01, 0)) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 1.3), breaks = seq(0, 1, .2)) +
    coord_cartesian(xlim = c(0.5, 20.5), ylim = c(0, 1), clip = "off") +
    theme_classic() +
    theme(

        legend.key.size = unit(.5, "cm"),
        legend.spacing.y = unit(.3, "cm"),
        panel.border = element_rect(color = 1, fill = NA),
        plot.margin = unit(c(1,.5,.5,.5), "cm")
    )  +
    guides(color = "none", fill = guide_legend(title = "")) +
    labs(x = "simulated community", y = "fraction") +
    ggtitle("")


# Empirical results ----
pairs <- read_csv(paste0(folder_data, "output/pairs.csv"), col_types = cols())
communities <- read_csv(paste0(folder_data, "temp/00c-communities.csv"), show_col_types = F)

p2 <- pairs %>%
    filter(!is.na(FitnessFunction)) %>%
    group_by(Community, InteractionType) %>%
    count(name = "Count") %>%
    group_by(Community) %>% mutate(Fraction = Count / sum(Count), TotalCount = sum(Count)) %>%
    left_join(communities, by = "Community") %>%
    replace_na(list(InteractionType = "unknown")) %>%
    ungroup() %>%
    ggplot() +
    geom_col(aes(x = CommunityLabel, fill = InteractionType, y = Fraction), color = 1, width = .8, linewidth = .5) +
    annotate("text", x = 1:13, y = 1.15, label = communities$CommunitySize, size = 4) +
    annotate("text", x = 14, y = 1.15, label = "n. of species", size = 4, hjust = 0) +
    annotate("segment", x = .5, xend = 18, y = 1.1, yend = 1.1, color = "black") +
    geom_text(aes(x = CommunityLabel, y = 1.05, label = TotalCount), size = 4) +
    annotate("text", x = 14, y = 1.05, label = "n. of tested pairs", size = 4, hjust = 0) +
    scale_fill_manual(values = assign_interaction_color(), breaks = c("coexistence", "exclusion", "unknown")) +
    scale_x_continuous(breaks = 1:13, expand = c(0.01, 0)) +
    scale_y_continuous(breaks = seq(0,1,0.2), limit = c(0, 1.3), expand = c(0,0)) +
    coord_cartesian(xlim = c(0.5, 13.5), ylim = c(0, 1), clip = "off") +
    theme_classic() +
    theme(
        legend.text = element_text(size = 10),
        legend.title = element_blank(),
        legend.key.size = unit(.5, "cm"),
        legend.spacing.y = unit(.3, "cm"),
        legend.position = "right",
        panel.border = element_rect(color = 1, fill = NA),
        axis.text = element_text(color = 1, size = 10),
        axis.title = element_text(color = 1, size = 10),
        plot.margin = unit(c(1,.5,.5,.5), "cm")
    ) +
    guides(fill = guide_legend(byrow = TRUE)) +
    labs(x = "empirical community", y = "fraction") +
    ggtitle("")

p <- plot_grid(
    plot_grid(NULL, NULL, ncol = 2, labels = LETTERS[c(1,2)]),
    plot_grid(p1, p2, nrow = 2, align = "v", axis = "rl", scale = .9,
              labels = LETTERS[c(3,4)], hjust = 0, label_x = 0.01),
    ncol = 1, rel_heights = c(1,2)
) + paint_white_background()
#p <- p2
ggsave(here::here("plots/FigS17-pairwise_model_empirical.png"), p, width = 8, height = 10)



withinCommunityPairs_fraction <- withinCommunityPairs_N_outcome %>%
    mutate(Community = factor(Community, communities_richness$Community)) %>%
    mutate(InteractionType = factor(InteractionType, c("exclusion", "coexistence"))) %>%
    group_by(Community, InteractionType, .drop = F) %>%
    summarize(Count = n()) %>%
    drop_na() %>%
    arrange(Community, InteractionType) %>%
    mutate(Fraction = Count / sum(Count)) %>%
    left_join(communities_richness) %>%
    select(Community, InteractionType, Fraction) %>%
    filter(InteractionType == "exclusion") %>%
    ungroup()













if (FALSE) {
    library(tidyverse)
    library(cowplot)
    source(here::here("analysis/00-metadata.R"))
    source(here::here("simulation/01-generate_input.R"))

    # 0. parameters ----
    input_parameters <- read_csv(here::here("simulation/01-input_parameters.csv"), col_types = cols())
    input_communities <- read_csv(here::here("simulation/02b-input_communities.csv"), col_types = cols())
    input_poolPairs <- read_csv(here::here("simulation/03a-input_poolPairs.csv"), col_types = cols())
    input_withinCommunityPairs <- read_csv(here::here("simulation/03b-input_withinCommunityPairs.csv"), col_types = cols())

    mcrm_family_colors <- RColorBrewer::brewer.pal(10, "Set3") %>% setNames(paste0("F", 0:9))

    family_names <- paste0("family ", 1:10) %>% setNames(paste0("F", 0:9))
    resource_names <- LETTERS[1:10] %>% setNames(paste0("R", 0:9))
    n_timesteps <- input_communities$n_timesteps[1]
    n_timepoints <- input_communities$n_timepoints[1]


    # 1. Pool pairs ----
    poolPairs_N_freq <- read_csv(paste0(folder_simulation, "aggregated/13-poolPairs_N_freq.csv"), col_types = cols()) %>%
        mutate(Time = factor(Time, c("init", "end"))) %>%
        mutate(Pair = factor(Pair, paste0("P", 1:1000))) %>%
        mutate(InitialFrequency = factor(InitialFrequency, c(5, 50, 95)))

    # Barplot
    monocultureSets_richness <- read_csv(paste0(folder_simulation, "aggregated/03-monocultureSets_richness.csv"), col_types = cols()) %>%
        mutate(Community = factor(Community, paste0("W", 0:19))) %>%
        mutate(PairSize = choose(Richness, 2))
    poolPairs_N_outcome <- read_csv(paste0(folder_simulation, "aggregated/13-poolPairs_N_outcome.csv"), col_types = cols()) %>%
        mutate(Community = factor(Community, paste0("W", 0:19)))

    p1 <- poolPairs_N_outcome %>%
        group_by(Community, InteractionType, .drop = F) %>%
        summarize(Count = n()) %>%
        arrange(Community, InteractionType) %>%
        mutate(Fraction = Count / sum(Count)) %>%
        left_join(monocultureSets_richness) %>%
        ggplot() +
        geom_col(aes(x = Community, y = Fraction, fill = InteractionType), color = 1) +
        annotate("text", x = 1:20, y = 1.15, label = monocultureSets_richness$Richness, size = 4) +
        annotate("text", x = 21, y = 1.15, label = c("n. of species"), size = 4, hjust = 0) +
        annotate("segment", x = .5, xend = 25, y = 1.1, yend = 1.1, color = "black") +
        annotate("text", x = 21, y = 1.05, label = c("n. of tested pairs"), size = 4, hjust = 0) +
        annotate("text", x = 1:20, y = 1.05, label = monocultureSets_richness$PairSize, size = 4) +
        scale_fill_manual(values = interaction_color) +
        scale_x_discrete(breaks = paste0("W", 0:19)) +
        scale_y_continuous(expand = c(0, 0), limits = c(0, 1.2)) +
        coord_cartesian(xlim = c(0.5, 20.5), ylim = c(0, 1), clip = "off") +
        theme_classic() +
        theme(plot.margin = unit(c(2,.5,.5,.5), "cm"))  +
        guides(color = "none", fill = guide_legend(title = "")) +
        labs(x = "Set")
    p1
    #ggsave(here::here("simulation/plots/23-poolPairs.png"), p, width = 8, height = 4)

    # 2. Community pairs ----
    withinCommunityPairs_N_freq <- read_csv(paste0(folder_simulation, "aggregated/13-withinCommunityPairs_N_freq.csv"), col_types = cols()) %>%
        mutate(Time = factor(Time, c("init", "end"))) %>%
        mutate(Pair = factor(Pair, paste0("P", 1:1000))) %>%
        mutate(InitialFrequency = factor(InitialFrequency, c(5, 50, 95)))

    # # Line plot
    # withinCommunityPairs_N_freq %>%
    #     filter(Community == "W6") %>%
    #     ggplot(aes(x = Time, y = Frequency1, color = InitialFrequency, group = InitialFrequency)) +
    #     geom_line(linewidth = 1) +
    #     geom_point(size = 1, shape = 21) +
    #     scale_color_manual(values = c("F0" = "#8A89C0", "F1" = "#FFCB77")) +
    #     scale_linetype_manual(values = c("F0" = 1, "F1" = 2)) +
    #     #scale_y_log10() +
    #     facet_wrap(~Pair) +
    #     theme_classic() +
    #     theme(panel.border = element_rect(color = 1, fill = NA)) +
    #     guides(alpha = "none", color = "none") +
    #     labs()


    # Barplot ----
    communities_richness <- read_csv(paste0(folder_simulation, "aggregated/12-communities_richness.csv"), col_types = cols()) %>%
        mutate(Community = factor(Community, paste0("W", 0:19))) %>%
        mutate(PairSize = choose(Richness, 2))
    withinCommunityPairs_N_outcome <- read_csv(paste0(folder_simulation, "aggregated/13-withinCommunityPairs_N_outcome.csv"), col_types = cols()) %>%
        mutate(Community = factor(Community, paste0("W", 0:19)))

    p2 <- withinCommunityPairs_N_outcome %>%
        group_by(Community, InteractionType, .drop = F) %>%
        summarize(Count = n()) %>%
        arrange(Community, InteractionType) %>%
        mutate(Fraction = Count / sum(Count)) %>%
        left_join(communities_richness) %>%
        ggplot() +
        geom_col(aes(x = Community, y = Fraction, fill = InteractionType), color = 1) +
        #geom_text(aes(x = Community, label = Richness), y = 0.9) +
        annotate("text", x = 1:20, y = 1.15, label = communities_richness$Richness, size = 4) +
        annotate("text", x = 21, y = 1.15, label = c("n. of species"), size = 4, hjust = 0) +
        annotate("segment", x = .5, xend = 25, y = 1.1, yend = 1.1, color = "black") +
        annotate("text", x = 21, y = 1.05, label = c("n. of tested pairs"), size = 4, hjust = 0) +
        annotate("text", x = 1:20, y = 1.05, label = communities_richness$PairSize, size = 4) +
        scale_fill_manual(values = interaction_color) +
        scale_x_discrete(breaks = paste0("W", 0:19)) +
        scale_y_continuous(expand = c(0, 0), limits = c(0, 1.2)) +
        coord_cartesian(xlim = c(0.5, 20.5), ylim = c(0, 1), clip = "off") +
        theme_classic() +
        theme(plot.margin = unit(c(2,.5,.5,.5), "cm"))  +
        guides(color = "none", fill = guide_legend(title = "")) +
        labs(x = "Community")

    # Assemble panels ----
    p <- plot_grid(p1, p2, nrow = 2, align = "v", axis = "rl", labels = c("Pool pairs", "Within-community pairs"), hjust = 0, label_x = 0.01) + paint_white_background()
    ggsave(here::here("plots/FigS17-pairwise_outcomes.png"), p, width = 8, height = 8)

    #ggsave(here::here("simulation/plots/23-withinCommunityPairs.png"), p, width = 8, height = 4)


    poolPair_fraction <- poolPairs_N_outcome %>%
        group_by(Community, InteractionType, .drop = F) %>%
        summarize(Count = n()) %>%
        arrange(Community, InteractionType) %>%
        mutate(Fraction = Count / sum(Count)) %>%
        left_join(monocultureSets_richness) %>%
        select(Community, InteractionType, Fraction) %>%
        filter(InteractionType == "exclusion") %>%
        ungroup()

    withinCommunityPairs_fraction <- withinCommunityPairs_N_outcome %>%
        mutate(InteractionType = factor(InteractionType, c("exclusion", "coexistence"))) %>%
        group_by(Community, InteractionType, .drop = F) %>%
        summarize(Count = n()) %>%
        drop_na() %>%
        arrange(Community, InteractionType) %>%
        mutate(Fraction = Count / sum(Count)) %>%
        left_join(communities_richness) %>%
        select(Community, InteractionType, Fraction) %>%
        filter(InteractionType == "exclusion") %>%
        ungroup()

    t.test(poolPair_fraction$Fraction, withinCommunityPairs_fraction$Fraction) %>%
        broom::tidy()


    stderror <- function(x) sd(x)/sqrt(length(x))
    mean(poolPair_fraction$Fraction)
    stderror(poolPair_fraction$Fraction)

    mean(withinCommunityPairs_fraction$Fraction)
    stderror(withinCommunityPairs_fraction$Fraction)




}
