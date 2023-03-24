library(tidyverse)
library(cowplot)
library(tidygraph)
library(ggraph)
source(here::here("analysis/00-metadata.R"))

pairs <- read_csv(paste0(folder_data, "output/pairs.csv"), show_col_types = F)
communities <- read_csv(paste0(folder_data, "temp/00c-communities.csv"), show_col_types = F)
pairs_freq <- read_csv(paste0(folder_data, "temp/93-pairs_freq.csv"), show_col_types = F)
load(paste0(folder_data, "temp/95-communities_network.Rdata"))

# Clean up pairs data
pairs <- pairs %>%
    # Remove no-colony pairs
    unite(col = "Pair", Community, Isolate1, Isolate2, sep = "_", remove = F) %>%
    filter(!(Pair %in% pairs_no_colony)) %>%
    # Remove low-accuracy model pairs
    filter(AccuracyMean > 0.9)

# Figure 2A: cartoon----
pA <- ggdraw() + draw_image(here::here("plots/cartoons/Fig2A.png")) + paint_white_background()

# Figure 2B: example network C2R6 ----
## Main network
node_size = 15
net <- communities_network %>%
    filter(Community == "C2R6") %>%
    pull(Network) %>% `[[`(1) %>%
    activate(nodes) %>%
    mutate(x = c(0, 0, 1, 1), y = c(1, 0, 1, 0))
p1 <- net %>%
    ggraph(layout = "nicely") +
    geom_edge_link(aes(color = InteractionType), width = 1,
                   arrow = arrow(length = unit(node_size/2-3, "mm"), type = "closed", angle = 30, ends = "last"),
                   start_cap = circle(node_size-7, "mm"),
                   end_cap = circle(node_size-7, "mm")) +
    scale_edge_color_manual(values = interaction_color, breaks = c("coexistence", "exclusion")) +
    scale_x_continuous(limits = c(-0.2, 1.2), breaks = c(0, .5, 1)) +
    scale_y_continuous(limits = c(-0.2, 1.2), breaks = c(0, .5, 1)) +
    theme_void() +
    theme(
        legend.position = "none",
        legend.key.size = unit(1.3, "line"),
        panel.background = element_blank(),
        plot.margin = unit(c(0,0,0,0), "mm"),
        plot.background = element_rect(fill = NA, color = NA),
    ) +
    draw_image(here::here("plots/cartoons/Fig2B_1.png"), x = -0.5, y = 0.75, vjust = 0.25, hjust = 0, clip = "on", scale = .3) +
    draw_image(here::here("plots/cartoons/Fig2B_2.png"), x = 0.5, y = 0.75, vjust = 0.25, hjust = 0, clip = "on", scale = .3) +
    draw_image(here::here("plots/cartoons/Fig2B_3.png"), x = -0.5, y = -0.25, vjust = 0.25, hjust = 0, clip = "on", scale = .3) +
    draw_image(here::here("plots/cartoons/Fig2B_4.png"), x = 0.5, y = -0.25, vjust = 0.25, hjust = 0, clip = "on", scale = .3)

## frequency plots
pairs_example_freq <- pairs %>%
    filter(Community == "C2R6") %>%
    select(Community, starts_with("Isolate"), starts_with("Interaction")) %>%
    mutate(Pair = 1:n()) %>%
    left_join(pairs_freq, by = c("Community", "Isolate1", "Isolate2"), multiple = "all") %>%
    mutate(Isolate1InitialODFreq = factor(Isolate1InitialODFreq)) %>%
    mutate(InteractionType = factor(InteractionType, c("exclusion", "coexistence")))

plot_example_freq <- function(pairs_freq) {
    pairs_freq %>%
        mutate(Isolate1InitialODFreq = factor(Isolate1InitialODFreq, c(95,50,5))) %>%
        ggplot(aes(x = Time, y = Isolate1CFUFreqMean, color = Isolate1InitialODFreq, group = Isolate1InitialODFreq)) +
        geom_line(linewidth = 1) +
        geom_point(size = 1.5) +
        geom_segment(aes(x = Time, xend = Time,
                         y = Isolate1CFUFreqMean + 2*Isolate1CFUFreqSd,
                         yend = Isolate1CFUFreqMean - 2*Isolate1CFUFreqSd,
                         color = Isolate1InitialODFreq), linewidth = 1) +
        scale_y_continuous(breaks = c(0, .5, 1), limits = c(-0.2, 1.2)) +
        scale_x_discrete(labels = c("start", "end")) +
        scale_color_manual(values = frequency_color, label = c("95%", "50%", "5%")) +
        facet_wrap(.~Pair) +
        theme_bw() +
        theme(
            panel.spacing = unit(2, "mm"), strip.text.x = element_blank(),
            panel.border = element_rect(color = 1, fill = NA, linewidth = 1),
            panel.grid.minor.y = element_blank(),
            panel.background = element_rect(fill = "white"),
            axis.title = element_blank(), axis.text = element_blank(),
            axis.ticks = element_blank(),
            plot.background = element_blank()
        ) +
        guides(color = "none") +
        labs(x = "time", y = "frequency")
}
p_pairs_example_freq_list <- pairs_example_freq %>% group_split(Pair) %>% lapply(plot_example_freq)
p_pairs_example_freq_list[[1]] <- p_pairs_example_freq_list[[1]] + theme(axis.text = element_text(size = 8), axis.title = element_text(size = 10))
p_legend_color <- get_legend({p_pairs_example_freq_list[[1]] +
        theme(legend.background = element_blank(), legend.text = element_text(size = 10), legend.title = element_blank()) +
        guides(color = guide_legend(title = "Initial frequency"))})
ss <- .25
pB <- ggdraw(p1) +
    draw_plot(p_pairs_example_freq_list[[1]], x = -.02, y = .5, width = ss*1.5, height = ss*1.5, hjust = .5, vjust = .5) +
    draw_plot(p_pairs_example_freq_list[[2]], x = .5, y = .85, width = ss, height = ss*1.1, hjust = .5, vjust = .5) +
    draw_plot(p_pairs_example_freq_list[[3]], x = .35, y = .58, width = ss, height = ss*1.1, hjust = .5, vjust = .5) +
    draw_plot(p_pairs_example_freq_list[[4]], x = .65, y = .58, width = ss, height = ss*1.1, hjust = .5, vjust = .5) +
    draw_plot(p_pairs_example_freq_list[[5]], x = .5, y = .15, width = ss, height = ss*1.1, hjust = .5, vjust = .5) +
    draw_plot(p_pairs_example_freq_list[[6]], x = .95, y = .55, width = ss, height = ss*1.1, hjust = .5, vjust = .5) +
    #draw_plot(p_legend_color, x = .05, y = .9, width = ss/2, height = ss/2, hjust = 0.5, vjust = .5) +
    theme(panel.background = element_blank(), plot.background = element_rect(color = NA, fill = "white"),
          plot.margin = unit(c(0,3,3,5), "mm"))

# Figure 2C: pairwise outcomes per community ----
pC <- pairs %>%
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
    theme(legend.text = element_text(size = 10),
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
    labs(x = "community", y = "fraction")

# Assemble panels ----
p_bottom <- plot_grid(pB, pC, nrow = 1, labels = c("B", "C"), scale = c(0.9, 0.9), rel_widths = c(1, 1.8), axis = "b")
p <- plot_grid(pA, p_bottom, nrow = 2, labels = c("A", ""), scale = c(.95, .95), rel_heights = c(.8, 1)) + paint_white_background()
ggsave(here::here("plots/Fig2.png"), p, width = 10, height = 6)


# Stats ----

pairs <- read_csv(paste0(folder_data, "output/pairs.csv"), show_col_types = F)
communities <- read_csv(paste0(folder_data, "temp/00c-communities.csv"), show_col_types = F)
pairs_freq <- read_csv(paste0(folder_data, "temp/93-pairs_freq.csv"), show_col_types = F)

pairs <- pairs %>%
    # Remove no-colony pairs
    unite(col = "Pair", Community, Isolate1, Isolate2, sep = "_", remove = F) %>%
    filter(!(Pair %in% pairs_no_colony)) %>%
    # Remove low-accuracy model pairs
    filter(AccuracyMean > 0.9)
pairs_freq <- pairs_freq %>%
    left_join(distinct(pairs, PairID, Community, Isolate1, Isolate2)) %>%
    select(PairID, everything()) %>%
    filter(!is.na(PairID))

compute_freq <- function (n_incidence, n_total) {
    paste0(n_incidence, "/", n_total, "=", round(n_incidence/n_total*100, 1), "%")
}

# One species went extinct anyways
pairs_freq %>%
    group_by(PairID) %>%
    filter(Time == "T8") %>%
    select(PairID, Isolate1InitialODFreq, Isolate1CFUFreqMean) %>%
    pivot_wider(names_from = Isolate1InitialODFreq, names_prefix = "F", values_from = Isolate1CFUFreqMean) %>%
    filter((F5 == 0 & F50 == 0 & F95 == 0) | (F5 == 1 & F50 == 1 & F95 == 1)) %>%
    nrow() %>%
    compute_freq(nrow(pairs))

# Losing species declines in frequency but not extinct
temp_id <- pairs_freq %>%
    group_by(PairID) %>%
    filter(Time == "T8") %>%
    select(PairID, Isolate1InitialODFreq, Isolate1CFUFreqMean) %>%
    pivot_wider(names_from = Isolate1InitialODFreq, names_prefix = "F", values_from = Isolate1CFUFreqMean) %>%
    filter(!((F5 == 0 & F50 == 0 & F95 == 0) | (F5 == 1 & F50 == 1 & F95 == 1))) %>%
    pull(PairID)

pairs %>%
    filter(InteractionType == "exclusion") %>%
    filter(PairID %in% temp_id) %>%
    nrow() %>%
    compute_freq(nrow(pairs))

# Coexisting at 1) stable equilibrim or 2) one negative frequency-dependent equilibirum
unique(pairs$InteractionTypeFiner)
pairs %>%
    filter(InteractionType == "coexistence") %>%
    filter(InteractionTypeFiner %in% c("stable coexistence", "frequency-dependent coexistence")) %>%
    nrow() %>%
    compute_freq(nrow(filter(pairs, InteractionType == "coexistence")))

# Coexisting at 1) 5% or 95%, or 2) neutrality
pairs %>%
    filter(InteractionType == "coexistence") %>%
    filter(InteractionTypeFiner %in% c("coexistence at 95%", "coexistence at 5%",
                                       "2-freq neutrality", "3-freq neutrality")) %>%
    nrow() %>%
    compute_freq(nrow(filter(pairs, InteractionType == "coexistence")))

# Undetermined pairs
pairs %>%
    filter(InteractionType == "unknown") %>%
    nrow() %>%
    compute_freq(nrow(pairs))


























