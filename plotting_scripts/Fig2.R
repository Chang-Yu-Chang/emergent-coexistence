library(tidyverse)
library(cowplot)
library(tidygraph)
library(ggraph)
source(here::here("analysis/00-metadata.R"))

pairs <- read_csv(paste0(folder_data, "output/pairs.csv"), show_col_types = F)
communities <- read_csv(paste0(folder_data, "temp/00c-communities.csv"), show_col_types = F)
pairs_freq <- read_csv(paste0(folder_data, "temp/93a-pairs_freq.csv"), show_col_types = F)
load(paste0(folder_data, "temp/95-communities_network.Rdata"))

pairs_freq <- pairs_freq %>% left_join(pairs) %>% remove_ineligible_pairs()
pairs <- remove_ineligible_pairs(pairs)

# Figure 2A: cartoon----
pA <- ggdraw() + draw_image(here::here("plots/cartoons/Fig2A.png")) + paint_white_background()
#pA <- ggdraw()

# Figure 2B: example network C2R6 ----
## Main network
node_size = 15
net <- communities_network %>%
    filter(Community == "C1R2") %>%
    pull(Network) %>% `[[`(1) %>%
    activate(nodes) %>%
    mutate(x = c(0, 0, 1, 1), y = c(1, 0, 1, 0))
p1 <- net %>%
    ggraph(layout = "nicely") +
    geom_edge_link(aes(color = outcome), width = 1,
                   arrow = arrow(length = unit(node_size/2-3, "mm"), type = "closed", angle = 30, ends = "last"),
                   start_cap = circle(node_size-7, "mm"),
                   end_cap = circle(node_size-7, "mm")) +
    scale_edge_color_manual(values = outcome_colors) +
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
pairs_freq_comm <- pairs_freq %>%
    filter(Community == "C1R2")

add_color_legend <- function (p) {
    p +
        theme(legend.background = element_blank(),
              legend.text = element_text(size = 8),
              legend.title = element_blank()) +
        guides(color = guide_legend(title = "Initial frequency"))
}
remove_axes <- function (p) {
    p + theme(
        axis.title = element_blank(),
        axis.text = element_blank()
    )
}
plot_freq <- function(pairs_freq) {
    line_size = 1
    pairs_freq %>%
        flip_winner_species_freq %>%
        mutate(Isolate1InitialODFreq = factor(Isolate1InitialODFreq, c(95,50,5))) %>%
        mutate(Time = as.character(Time) %>% str_replace("T", "") %>% as.numeric()) %>%
        mutate(Nudge = case_when(Isolate1InitialODFreq == 5 ~ .3, Isolate1InitialODFreq == 50 ~.6, Isolate1InitialODFreq == 95 ~ .9)) %>%
        ggplot() +
        geom_rect(aes(fill = outcome), color = 1, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, alpha = .05, linewidth = .8) +
        geom_hline(linewidth = line_size/2, yintercept = c(0,1), linetype = 1, color = "black") +
        #geom_hline(data = pairs_freq_mean_three, aes(yintercept = Isolate1CFUFreq_mean_three), linewidth = line_size/2, color = "black", linetype = 2) +
        geom_line(aes(x = Time+Nudge, y = Isolate1CFUFreqMean, color = Isolate1InitialODFreq, group = Isolate1InitialODFreq), linewidth = line_size) +
        #geom_point(aes(x = Time+Nudge, y = Isolate1CFUFreqMean, color = Isolate1InitialODFreq, group = Isolate1InitialODFreq), size = line_size*2, shape = 21) +
        geom_segment(aes(x = Time+Nudge, xend = Time+Nudge, y = Isolate1CFUFreqPercentile5, yend = Isolate1CFUFreqPercentile95, color = factor(Isolate1InitialODFreq)), linewidth = line_size) +
        scale_x_continuous(breaks = c(0,8), limits = c(-3,11)) +
        scale_y_continuous(breaks = c(0, .5, 1), labels = c("0", "0.5", "1")) +
        scale_color_manual(values = frequency_color, label = c("95%", "50%", "5%")) +
        scale_fill_manual(values = outcome_colors, labels = outcome_labels) +
        theme_classic() +
        theme(
            panel.spacing = unit(2, "mm"),
            panel.border = element_blank(),
            panel.background = element_rect(fill = "white"),
            panel.grid.minor.y = element_blank(),
            axis.text = element_text(size = 8),
            axis.title = element_text(size = 10),
            legend.position = "right",
            legend.title = element_text(size = 15),
            legend.text = element_text(size = 10),
            legend.key.size = unit(5, "mm"),
            legend.spacing.y = unit(5, "mm"),
            strip.text = element_blank(),
            plot.background = element_rect(color = NA, fill = NA)
        ) +
        guides(
            fill = "none",
            color = "none"
        ) +
        labs(x = "transfer", y = "frequency")
}

p_pairs_freq_list <- pairs_freq_comm %>%  mutate(PairID = factor(PairID, 1:200)) %>% arrange(PairID) %>% group_split(PairID) %>% lapply(plot_freq)
for (i in 2:6) p_pairs_freq_list[[i]] <- remove_axes(p_pairs_freq_list[[i]])
#p_pairs_freq_list[[1]] <- p_pairs_freq_list[[1]] + theme(axis.text = element_text(size = 8), axis.title = element_text(size = 10))
p_legend_color <- get_legend(add_color_legend(p_pairs_freq_list[[1]]))
ss <- .25
pB <- ggdraw(p1) +
    draw_plot(p_pairs_freq_list[[1]], x = -.02, y = .5, width = ss*1.5, height = ss*1.5, hjust = .5, vjust = .5) +
    draw_plot(p_pairs_freq_list[[2]], x = .5, y = .85, width = ss, height = ss*1.1, hjust = .5, vjust = .5) +
    draw_plot(p_pairs_freq_list[[3]], x = .35, y = .58, width = ss, height = ss*1.1, hjust = .5, vjust = .5) +
    draw_plot(p_pairs_freq_list[[4]], x = .65, y = .58, width = ss, height = ss*1.1, hjust = .5, vjust = .5) +
    draw_plot(p_pairs_freq_list[[5]], x = .5, y = .15, width = ss, height = ss*1.1, hjust = .5, vjust = .5) +
    draw_plot(p_pairs_freq_list[[6]], x = .95, y = .55, width = ss, height = ss*1.1, hjust = .5, vjust = .5) +
    #draw_plot(p_legend_color, x = .05, y = .9, width = ss/2, height = ss/2, hjust = 0.5, vjust = .5) +
    theme(panel.background = element_blank(), plot.background = element_rect(color = NA, fill = "white"),
          plot.margin = unit(c(0,3,3,5), "mm"))

# Figure 2C: pairwise outcomes per community ----
pC <- pairs %>%
    group_by(Community, outcome) %>%
    count(name = "Count") %>%
    # Total count
    group_by(Community) %>% mutate(Fraction = Count / sum(Count), TotalCount = sum(Count)) %>%
    left_join(communities, by = "Community") %>%
    ungroup() %>%
    ggplot() +
    geom_col(aes(x = CommunityLabel, fill = outcome, y = Fraction), color = 1, width = .8, linewidth = .5) +
    annotate("text", x = 1:13, y = 1.15, label = communities$CommunitySize, size = 4) +
    annotate("text", x = 14, y = 1.15, label = "n. of species", size = 4, hjust = 0) +
    annotate("segment", x = .5, xend = 18, y = 1.1, yend = 1.1, color = "black") +
    geom_text(aes(x = CommunityLabel, y = 1.05, label = TotalCount), size = 4) +
    annotate("text", x = 14, y = 1.05, label = "n. of tested pairs", size = 4, hjust = 0) +
    scale_fill_manual(values = outcome_colors, labels = outcome_labels) +
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
    labs(x = "community", y = "fraction")

# Assemble panels ----
p_bottom <- plot_grid(pB, pC, nrow = 1, labels = c("B", "C"), scale = c(0.9, 0.9), rel_widths = c(1, 1.8), axis = "b")
p <- plot_grid(pA, p_bottom, nrow = 2, labels = c("A", ""), scale = c(.95, .95), rel_heights = c(.8, 1)) + paint_white_background()
ggsave(here::here("plots/Fig2.png"), p, width = 10, height = 6)


# Stat
pairs %>%
    group_by(outcome) %>%
    count(name = "Count") %>%
    ungroup() %>%
    mutate(Fraction = Count / sum(Count), TotalCount = sum(Count))




