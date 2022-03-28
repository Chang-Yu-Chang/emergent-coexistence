# Figure 1: species coexistence is an emergent property

library(tidyverse)
library(tidymodels)
library(cowplot)
library(ggraph)
library(tidygraph)
library(officer)
library(flextable)
library(ggprism)
library(ggsci)
library(gridExtra)
source(here::here("plotting_scripts/network_functions.R"))

sequences_abundance <- read_csv(here::here("data/temp/sequences_abundance.csv"), col_types = cols())
communities <- read_csv(here::here("data/output/communities.csv"), col_types = cols())
isolates <- read_csv(here::here("data/output/isolates.csv"), col_types = cols())
pairs <- read_csv(here::here("data/output/pairs.csv"), col_types = cols()) %>% mutate(InteractionType = ifelse(InteractionType == "neutrality", "coexistence", InteractionType))
pairs_freq <- read_csv(here::here("data/output/pairs_freq.csv"), col_types = cols())
load("~/Dropbox/lab/invasion-network/data/output/network.Rdata")
load("~/Dropbox/lab/invasion-network/data/output/network_randomized.Rdata")
community_factor <- c(communities %>% filter(str_detect(Community, "C\\d")) %>% arrange(CommunitySize) %>% pull(Community),
                      communities %>% filter(str_detect(Community, "Ass")) %>% pull(Community))
communities_size <- communities %>% mutate(Community = factor(Community, community_factor)) %>% arrange(Community) %>% pull(CommunitySize)


# Figure 1 ----------------------------------------------------------------------------------------
p <- ggdraw() + draw_image(here::here("plots/cartoons/Fig1.png")) + theme(plot.background = element_rect(fill = "white", color = NA))
ggsave(here::here("plots/Fig1.png"), p, width = 27, height = 15)

# Figure 2 ----------------------------------------------------------------------------------------
# Figure 2A: networks of one community. C7R1
## Main network
net <- net_list$C7R1 %>%
    activate(nodes) %>%
    mutate(x = c(0, 0, 1, 1), y = c(1, 0, 1, 0))
node_size = 15
p1 <- net %>%
    ggraph(layout = "nicely") +
    #geom_node_point(fill = "grey", size = node_size, shape = 21, colour = "black", stroke = node_size/5) +
    #geom_node_text(aes(label = Isolate)) +
    geom_edge_link(aes(color = InteractionType), width = 2,
                   arrow = arrow(length = unit(node_size/2-1, "mm"), type = "closed", angle = 30, ends = "last"),
                   start_cap = circle(node_size/2+1, "mm"),
                   end_cap = circle(node_size/2+1, "mm")) +
    scale_edge_color_manual(values = interaction_color) +
    scale_x_continuous(limits = c(-0.4, 1.4), breaks = c(0, .5, 1)) +
    scale_y_continuous(limits = c(-0.4, 1.4), breaks = c(0, .5, 1)) +
    theme_void() +
    #theme_bw() +
    theme(legend.position = "none", panel.background = element_blank(),
          plot.margin=unit(c(3,3,3,3),"mm"), plot.background = element_rect(fill = NA, color = NA)) +
    labs() +
    draw_image(here::here("plots/cartoons/Fig1B_1.png"), x = -0.5, y = 0.75, vjust = 0.25, hjust = 0, clip = "on", scale = .3) +
    draw_image(here::here("plots/cartoons/Fig1B_2.png"), x = -0.5, y = -0.25, vjust = 0.25, hjust = 0, clip = "on", scale = .3) +
    draw_image(here::here("plots/cartoons/Fig1B_3.png"), x = 0.5, y = 0.75, vjust = 0.25, hjust = 0, clip = "on", scale = .3) +
    draw_image(here::here("plots/cartoons/Fig1B_4.png"), x = 0.5, y = -0.25, vjust = 0.25, hjust = 0, clip = "on", scale = .3)


p_temp <- net %>%
    ggraph(layout = "nicely") +
    geom_edge_link(aes(color = InteractionType), width = 2,
                   arrow = arrow(length = unit(2, "mm"), type = "closed", angle = 30, ends = "last")) +
    scale_edge_color_manual(values = interaction_color) +
    theme_void() +
    theme(legend.text = element_text(size = 12), legend.position = "top", legend.title = element_blank(), legend.direction = "vertical", legend.background = element_blank())

p_legend_network <- get_legend(p_temp)

## frequency plots
pairs_example_freq <- pairs %>%
    filter(Community == "C7R1") %>%
    select(Community, starts_with("Isolate"), starts_with("Interaction")) %>%
    mutate(Pair = 1:n()) %>%
    left_join(pairs_freq, by = c("Community", "Isolate1", "Isolate2")) %>%
    mutate(Isolate1InitialODFreq = factor(Isolate1InitialODFreq)) %>%
    mutate(InteractionType = factor(InteractionType, c("exclusion", "coexistence")))
plot_example_freq <- function(pairs_freq) {
    pairs_freq %>%
        mutate(Isolate1InitialODFreq = factor(Isolate1InitialODFreq, c(95,50,5))) %>%
        ggplot(aes(x = Time, y = Isolate1MeasuredFreq, color = Isolate1InitialODFreq, group = Isolate1InitialODFreq)) +
        geom_point(size = 2) +
        geom_line(size = 1) +
        scale_y_continuous(breaks = c(0, .5, 1), limits = c(0,1)) +
        scale_x_discrete(labels = c(0,8)) +
        scale_color_manual(values = frequency_color, label = c("95%", "50%", "5%")) +
        facet_wrap(.~Pair) +
        theme_bw() +
        theme(panel.spacing = unit(2, "mm"), strip.text.x = element_blank(),
              panel.border = element_rect(color = 1, fill = NA, size = 1),
              panel.grid.minor.y = element_blank(),
              #axis.title = element_text(size = 10), axis.text = element_text(color = 1, size = 8),
              axis.title = element_blank(), axis.text = element_blank(),
              axis.ticks = element_blank(),
              panel.background = element_rect(fill = "white"),
              plot.background = element_blank()) +
        guides(color = "none") +
        labs(x = "Time (days)", y = "Frequency")
}
p_pairs_example_freq_list <- pairs_example_freq %>% group_split(Pair) %>% lapply(plot_example_freq)
p_pairs_example_freq_list[[1]] <- p_pairs_example_freq_list[[1]] + theme(axis.text = element_text(size = 10), axis.title = element_text(size = 10))
p_legend_inset <- {p_pairs_example_freq_list[[1]] +
        theme(legend.background = element_blank(), legend.text = element_text(size = 12), legend.title = element_text(size = 15)) +
        guides(color = guide_legend(title = "Initial frequencies"))} %>%
    cowplot::get_legend()

ss <- .2
pA <- ggdraw(p1) +
    draw_plot(p_pairs_example_freq_list[[1]], x = .05, y = .5, width = ss*1.5, height = ss*1.5, hjust = .5, vjust = .5) +
    draw_plot(p_pairs_example_freq_list[[2]], x = .5, y = .85, width = ss, height = ss, hjust = .5, vjust = .5) +
    draw_plot(p_pairs_example_freq_list[[3]], x = .4, y = .6, width = ss, height = ss, hjust = .5, vjust = .5) +
    draw_plot(p_pairs_example_freq_list[[4]], x = .6, y = .6, width = ss, height = ss, hjust = .5, vjust = .5) +
    draw_plot(p_pairs_example_freq_list[[5]], x = .5, y = .15, width = ss, height = ss, hjust = .5, vjust = .5) +
    draw_plot(p_pairs_example_freq_list[[6]], x = .85, y = .5, width = ss, height = ss, hjust = .5, vjust = .5) +
    draw_plot(p_legend_inset, x = .15, y = .9, width = ss/2, height = ss/2, hjust = 0, vjust = .5) +
    draw_plot(p_legend_network, x = .15, y = .15, width = ss/2, height = ss/2, hjust = 0.5, vjust = .5) +
    theme(panel.background = element_blank(), plot.background = element_rect(color = NA, fill = "white"),
          plot.margin = unit(c(10,0,0,10), "mm"))

ggsave(here::here("plots/Fig2A-example_network.png"), pA, width = 5, height = 5)

# Figure 2B: All network graphs
net_list <- net_list %>% `[`(as.character(community_factor))
## subset the self-assembly networks
plot_competitive_network_grey <- function(x, node_size, edge_width){
    plot_competitive_network(x, node_size = node_size, edge_width = edge_width) +
        theme(plot.background = element_rect(fill = "grey90", color = NA),
              panel.background = element_rect(fill = "grey90", color = NA))

}

p_net_list <- communities %>%
    filter(Assembly == "self_assembly") %>%
    select(Community, CommunitySize) %>%
    arrange(CommunitySize) %>%
    mutate(Network = net_list[1:13]) %>%
    mutate(CommunitySize = max(CommunitySize) / CommunitySize / 4) %>%
    rowwise() %>%
    mutate(p_net = plot_competitive_network_grey(Network, 0, CommunitySize) %>% list()) %>%
    pull(p_net)

pB_title <- ggdraw() +
    draw_label("Pairwise networks of 13 replicate communities",fontface = 'bold',x = 0,hjust = 0) +
    theme(plot.margin = margin(0, 0, 0, 7))
p1 <- plot_grid(plotlist = p_net_list, nrow = 1, scale = 1.3) + paint_white_background()
if (FALSE) {
    p1 <- plot_grid(p1_title,
                    plot_grid(plotlist = p_net_list[1:5], nrow = 1, labels = 1:5),
                    plot_grid(plotlist = p_net_list[6:10], nrow = 1, labels = 6:10),
                    plot_grid(plotlist = p_net_list[11:13], nrow = 1, labels = 11:13),
                    ncol = 1, rel_heights = c(.2,1,1,2), scale = .9) +
        paint_white_background()

}

## pairwise outcomes per community
pairs_count <- pairs %>%
    filter(Assembly == "self_assembly") %>%
    group_by(Community) %>%
    summarize(Count = n()) %>%
    mutate(Community = factor(Community, community_factor)) %>%
    arrange(Community) %>%
    mutate(CommunityLabel = factor(1:13))
p2 <- pairs %>%
    filter(Assembly == "self_assembly") %>%
    group_by(Community, InteractionType) %>%
    summarize(Count = n()) %>%
    mutate(Fraction = Count / sum(Count)) %>%
    ungroup() %>%
    mutate(Community = factor(Community, community_factor)) %>%
    arrange(Community) %>%
    mutate(CommunityLabel = rep(1:13, each = 2) %>% factor()) %>%
    ggplot() +
    geom_col(aes(x = CommunityLabel, fill = InteractionType, y = Fraction), color = 1, width = .8, size = .5) +
    geom_text(data = pairs_count, aes(x = CommunityLabel, y = .1, label = paste0("n=", Count)), vjust = -.5) +
    #geom_rect(xmin = -Inf, xmax = Inf, ymin = 0, ymax = 1, color = grey(0.1), fill = NA, size = .5) +
    scale_fill_manual(values = assign_interaction_color()) +
    #scale_x_discrete(labels = 1:13) +
    scale_y_continuous(breaks = c(0,.5,1), limit = c(0, 1), expand = c(0,0)) +
    facet_grid(.~factor(CommunityLabel, 1:13), scales = "free_x") +
    theme_classic() +
    theme(strip.background = element_blank(),
          strip.text = element_blank(),
          legend.text = element_text(size = 12),
          axis.text = element_text(color = 1, size = 12),
          axis.title = element_text(color = 1, size = 12),
          panel.spacing = unit(0, "mm"),
          legend.title = element_blank(), legend.position = "right") +
    guides(fill = guide_legend(reverse = T)) +
    labs(x = "Community", y = "Fraction")

pB <- plot_grid(p1, p2, ncol = 1, scale = .9, rel_heights = c(1, 3), axis = "lr", align = "v") + paint_white_background()
ggsave(here::here("plots/Fig2B-all_networks.png"), pB, width = 12, height = 3)


# Figure 2C: pairwise competition
pairs_example_outcomes_finer <- read_csv(here::here("data/output/pairs_example_outcomes_finer.csv"))
temp1 <- pairs %>%
    filter(Assembly == "self_assembly") %>%
    mutate(InteractionType = factor(InteractionType, c("exclusion", "coexistence"))) %>%
    mutate(InteractionTypeFiner = factor(InteractionTypeFiner, c("competitive exclusion", "mutual exclusion", "stable coexistence", "frequency-dependent coexistence", "neutrality"))) %>%
    group_by(InteractionType, InteractionTypeFiner) %>%
    summarize(Count = n()) %>%
    ungroup() %>%
    mutate(Fraction = Count / sum(Count))

temp2 <- temp1 %>%
    group_by(InteractionType) %>%
    summarize(Count = sum(Count), Fraction = sum(Fraction))

p1 <- temp1 %>%
    ggplot() +
    geom_bar(aes(x = InteractionType, y = Fraction, fill = InteractionTypeFiner), stat="identity", color = 1, alpha = .8) +
    geom_text(x = Inf, y = Inf, label = paste0("n = ", sum(temp2$Count)), vjust = 1, hjust = 2) +
    geom_text(data = temp2, aes(x = InteractionType, y = Fraction, label = paste0(round(Fraction, 3) * 100,"%")), nudge_y = .05, size = 5) +
    scale_fill_manual(values = assign_interaction_color(level = "finer"),
                      breaks = c("competitive exclusion", "mutual exclusion", "stable coexistence", "frequency-dependent coexistence", "neutrality"),
                      labels = paste0(temp1$InteractionTypeFiner, " (", round(temp1$Fraction, 3) * 100,"%)")) +
    scale_y_continuous(limits = c(0, 1), expand = c(0,0),
                       breaks = c(0, .5, 1), labels = paste0(c(0,50,100), "%")) +
    theme_classic() +
    theme(axis.title.x = element_blank(), axis.title.y = element_text(size = 12),
          axis.text.y = element_text(size = 12, color = 1),
          axis.text.x = element_text(size = 12, color = "black", angle = 15, vjust = 1, hjust = 1),
          legend.text = element_text(size = 12),
          legend.background = element_blank(),
          legend.position = "right") +
    guides() +
    labs(x = "", y = "Fraction", fill = "")

pC <- p1
ggsave(here::here("plots/Fig2C-pairwise_competition.png"), pC, width = 6, height = 4)


p_temp <- temp2 %>%
    ggplot() +
    geom_bar(aes(x = InteractionType, y = Fraction, fill = InteractionType), stat="identity", color = 1) +
    geom_text(x = Inf, y = Inf, label = paste0("n = ", sum(temp2$Count)), vjust = 1, hjust = 2) +
    geom_text(data = temp2, aes(x = InteractionType, y = Fraction, label = paste0(round(Fraction, 3) * 100,"%")), nudge_y = .05, size = 5) +
    scale_fill_manual(values = assign_interaction_color(),
                      breaks = c("competitive exclusion", "mutual exclusion", "stable coexistence", "frequency-dependent coexistence", "neutrality"),
                      labels = paste0(temp1$InteractionTypeFiner, " (", round(temp1$Fraction, 3) * 100,"%)")) +
    scale_y_continuous(limits = c(0, 1), expand = c(0,0),
                       breaks = c(0, .5, 1), labels = paste0(c(0,50,100), "%")) +
    theme_classic() +
    theme(axis.title.x = element_blank(), axis.title.y = element_text(size = 12),
          axis.text.y = element_text(size = 12, color = 1),
          axis.text.x = element_text(size = 12, color = "black", angle = 15, vjust = 1, hjust = 1),
          legend.text = element_text(size = 12),
          legend.background = element_blank(),
          legend.position = "right") +
    guides() +
    labs(x = "", y = "Fraction", fill = "")


ggsave(here::here("plots/Fig2C-talk.png"), p_temp, width = 3, height = 4)


if (FALSE) {
    temp <- pairs %>% filter(Assembly == "self_assembly") %>%
        mutate(InteractionType = ifelse(InteractionType == "neutrality", "coexistence", InteractionType)) %>%
        mutate(InteractionType = factor(InteractionType, c("exclusion", "coexistence"))) %>%
        group_by(InteractionType) %>% summarize(Count = n()) %>% ungroup() %>% mutate(Fraction = Count / sum(Count))
    p_pairs_interaction <- temp %>%
        ggplot() +
        geom_col(aes(x = InteractionType, y = Count, fill = InteractionType), color = 1, width = 0.7) +
        geom_text(x = Inf, y = Inf, label = paste0("n = ", sum(temp$Count)), vjust = 1, hjust = 1.5, size = 5) +
        geom_text(aes(x = InteractionType, y = Count, label = paste0(round(Fraction, 3) * 100,"%")), nudge_y = 10, size = 5
        ) +
        scale_fill_manual(values = assign_interaction_color(level = "simple"), breaks = ordered(c("coexistence", "exclusion"))) +
        scale_y_continuous(limits = c(0, 150), expand = c(0,0)) +
        theme_classic() +
        theme(axis.title.x = element_blank(), axis.title.y = element_text(size = 15),
              axis.text.x = element_text(size = 15, color = "black", angle = 15, vjust = 1, hjust = 1),
              axis.text.y = element_text(size = 15, color = "black"),
              legend.position = "none") +
        labs(x = "", y = "Number of pairs", fill = "")

    pC <- p_pairs_interaction

}

p_top <- plot_grid(pA, pC, nrow = 1, labels = c("A", "C"), scale = .9, rel_widths = c(1,1.5))
p <- plot_grid(p_top, pB, ncol = 1, labels = c("", "B"), scale = c(1, 1), rel_heights = c(1.5,1), axis = "lr", align = "v") + paint_white_background()
ggsave(here::here("plots/Fig2.png"), p, width = 12, height = 8)



# Figure 4 ----
## Figure 4A One example network
node_size = 5
edge_width = 1.5

p_net <- net_list$C7R1 %>%
    activate(nodes) %>%
    mutate(y = -Rank) %>%
    group_by(Rank) %>%
    mutate(x = {seq(0, 1, length.out = n() + 2) %>% `[`(c(-1, -length(.)))}) %>% # + rnorm(n(), 0, .5)) %>%
    ungroup() %>%
    activate(edges) %>%
    ggraph(layout = "nicely") +
    geom_edge_link(aes(color = InteractionType), width = edge_width,
                  arrow = arrow(length = unit(edge_width, "mm"), type = "closed", angle = 30, ends = "last"),
                  start_cap = circle(node_size/2+1, "mm"),
                  end_cap = circle(node_size/2+1, "mm")) +
    scale_edge_color_manual(values = assign_interaction_color(level = "matrix"),
                            breaks = c("exclusion", "exclusion violating rank"),
                            labels = c("exclusion following rank", "exclusion violating rank")) +
    scale_x_continuous(limits = c(.2, .8), expand = c(0,0)) +
    scale_y_continuous(limits = c(-5, 0), breaks = -4:-1, labels = 4:1) +
    theme_void() +
    #theme_bw() +
    theme(
        panel.grid.major.y = element_line(color = "grey90"),
        legend.title = element_blank(),
        strip.text = element_blank(),
        plot.margin=unit(c(0,0,0,0),"mm"),
        legend.position = "none",
        axis.text.y = element_text(color = 1, size = 10),
        axis.title.y = element_text (color = 1, size = 10, angle = 90),
        legend.text = element_text(size = 13)
    ) +
    guides(color = "none") +
    labs(y = "Rank")


pA <- p_net +
    draw_image(here::here("plots/cartoons/Fig1B_1.png"), x = 0, y = -1, vjust = 0.5, hjust = 0, clip = "on", scale = .8) +
    draw_image(here::here("plots/cartoons/Fig1B_3.png"), x = -.17, y = -2, vjust = 0.5, hjust = 0, clip = "on", scale = .8) +
    draw_image(here::here("plots/cartoons/Fig1B_4.png"), x = .17, y = -2, vjust = 0.5, hjust = 0, clip = "on", scale = .8) +
    draw_image(here::here("plots/cartoons/Fig1B_2.png"), x = 0, y = -4, vjust = .6, hjust = 0, clip = "on", scale = .8) +
    paint_white_background()

ggsave(here::here("plots/Fig4A-example.png"), pA, width = 3, height = 3)



## Figure 4B. Experiment
set.seed(1)
node_size = 2
edge_width = .8
plot_network_hierarchy <- function(net, n_rank = 12, n_break = 12) {
    graph_ranked <- net %>%
        activate(nodes) %>%
        select(Isolate, Rank, PlotRank) %>%
        activate(edges) %>%
        mutate(fromRank = .N()$PlotRank[match(from, .N()$Isolate)],
               toRank = .N()$PlotRank[match(to, .N()$Isolate)])

    n_nodes <- igraph::vcount(graph_ranked)
    n_exclusion_violation <- graph_ranked %>%
        activate(edges) %>%
        filter(fromRank > toRank, InteractionType == "exclusion") %>%
        igraph::ecount()

    if (n_exclusion_violation != 0) {
        pairs_exclusion_violation <- graph_ranked %>%
            activate(edges) %>%
            filter(fromRank > toRank, InteractionType == "exclusion") %>%
            mutate(temp = fromRank, fromRank = toRank, toRank = temp) %>%
            select(-temp) %>%
            mutate(InteractionType = "exclusion violating rank") %>%
            as_tibble

        graph_ranked <- graph_ranked %>%
            activate(edges) %>%
            filter(fromRank <= toRank) %>%
            bind_edges(pairs_exclusion_violation)
    }

    temp <- graph_ranked %>%
        # Node position
        activate(nodes) %>%
        mutate(y = -Rank) %>%
        group_by(Rank) %>%
        mutate(x = {seq(0, 1, length.out = n() + 2) %>% `[`(c(-1, -length(.)))}) %>% # + rnorm(n(), 0, .5)) %>%
        ungroup() %>%
        # Filter out coexistence edges
        activate(edges) %>%
        filter(InteractionType != "coexistence") %>%
        mutate(Temp = sample(c(-1, 1), size = n(), replace = T))
    strength_angle <- as_tibble(temp)$Temp * 0.05

    temp %>%
        ggraph(layout = "nicely") +
        geom_hline(yintercept = c(-n_rank:-1), color = "grey90") +
        #geom_segment(aes(x = 0.5, xend = 0.5, y = -Inf, yend = -1), color = "grey90") +
        geom_node_point(size = node_size, shape = 21, fill = "grey", stroke = node_size/5, color = "black") +
        geom_edge_arc(strength = strength_angle, alpha = 0.5,
                      aes(color = InteractionType), width = edge_width,
                      arrow = arrow(length = unit(edge_width, "mm"), type = "closed", angle = 30, ends = "last"),
                      start_cap = circle(node_size/2, "mm"),
                      end_cap = circle(node_size/2, "mm")) +
        scale_edge_color_manual(values = assign_interaction_color(level = "matrix"),
                                breaks = c("exclusion", "exclusion violating rank"),
                                labels = c("exclusion pair that follows rank", "exclusion pair that violates rank")) +
        scale_x_continuous(limits = c(0.1, 0.9), expand = c(0,0)) +
        scale_y_continuous(limits = c(-n_break-1, 0), breaks = -n_break:-1, labels = n_break:1) +
        theme_void() +
        #theme_bw() +
        theme(
            #panel.grid.major.y = element_line(color = "grey90"),
            #panel.grid.minor.y = element_line(color = "grey80"),
            #panel.grid.major.x = element_blank(),
            #panel.grid.minor.x = element_blank(),
            legend.position = "none",
            #panel.border = element_rect(color = "grey80", fill = NA),
            legend.title = element_blank(),
            axis.title = element_blank(),
            strip.text = element_blank(),
            plot.margin=unit(c(0,0,0,0),"mm")
        ) +
        labs(y = "Rank")



}
communities_net <- communities %>%
    mutate(Community = factor(Community, community_factor)) %>%
    arrange(Community) %>%
    mutate(Network = net_list) %>%
    filter(Assembly == "self_assembly") %>%
    select(Community, CommunitySize, Network) %>%
    rowwise() %>%
    mutate(p_net = plot_network_hierarchy(Network, n_rank = CommunitySize) %>% list())

p_net_hierarchy_list <- communities_net %>% pull(p_net)
p_net_hierarchy_list[[13]] <- p_net_hierarchy_list[[13]] +
    scale_y_continuous(limits = c(-12-1, 0), breaks = -12:-1, labels = 12:1, position = "right") +
    theme(axis.title.y = element_text(color = 1, size = 10, angle = 270),
          axis.text.y = element_text(color = 1, size = 10))
pB_title <- ggdraw() + draw_label("Experiment", fontface = 'bold', x = 0, hjust = 0) + theme(plot.margin = margin(5, 0, 5, 7))
pB_axistitle <- ggdraw() + draw_label("Community", fontface = 'plain', x = .5, hjust = .5) + theme(plot.margin = margin(5, 0, 5, 7))
p_temp <- plot_grid(plotlist = p_net_hierarchy_list,
                 rel_widths = c(communities_net$CommunitySize / max(communities_net$CommunitySize)),
                 labels = 1:13, label_x = 0.5, hjust = c(rep(.5, 12), 1),
                 nrow = 1, axis = "tb", align = "h") + paint_white_background()
pB <- plot_grid(pB_title, pB_axistitle, p_temp, ncol = 1, rel_heights = c(.1, .1, 1)) + paint_white_background()
#pA <- ggdraw(p_temp) + draw_plot(p_legend, x = 0.1, y = 0.1, .3, .3)
ggsave(here::here("plots/Fig4B-network_hierarchy_experiment.png"), pB, width = 10, height = 3)


## Figure 4C. Simulation
df_communities <- read_csv(here::here("data/output/df_communities.csv"), col_types = cols())
load("~/Dropbox/lab/invasion-network/data/output/network_simulated.Rdata")

df_communities_ID <- df_communities %>%
    filter(Assembly == "self_assembly") %>%
    group_by(Community) %>%
    filter(Richness > 3) %>% ungroup() %>%
    slice(1:10) %>%
    arrange(Richness) %>%
    mutate(CommunityLabel = factor(1:10))

set.seed(1)
node_size = 2
edge_width = .8
df_communities_net <- df_communities %>%
    filter(Assembly == "self_assembly") %>%
    select(Community, Richness) %>%
    mutate(Network = net_simulated_list[paste0("self_assembly_W", 0:19)]) %>%
    right_join(select(df_communities_ID, Community, CommunityLabel)) %>%
    arrange(Richness) %>%
    rowwise() %>%
    mutate(p_net = plot_network_hierarchy(Network, n_rank = Richness, n_break = 8) %>% list())

p_net_simulated_hierarchy_list <- df_communities_net %>% pull(p_net)
p_net_simulated_hierarchy_list[[10]] <- p_net_simulated_hierarchy_list[[10]] +
    scale_y_continuous(limits = c(-8-1, 0), breaks = -8:-1, labels = 8:1, position = "right") +
    theme(axis.title.y = element_text(color = 1, size = 10, angle = 270),
          axis.text.y = element_text(color = 1, size = 10))

pC_title <- ggdraw() + draw_label("Simulation", fontface = 'bold', x = 0, hjust = 0) + theme(plot.margin = margin(5, 0, 5, 7))
pC_axistitle <- ggdraw() + draw_label("Community", fontface = 'plain', x = .5, hjust = .5)
p_temp <- plot_grid(plotlist = p_net_simulated_hierarchy_list,
                rel_widths = df_communities_net$Richness / max(df_communities_net$Richness),
                labels = 1:10, label_x = 0.5, hjust = c(rep(.5, 9), 1),
                nrow = 1, axis = "tb", align = "h") + paint_white_background()
pC <- plot_grid(pC_title, pC_axistitle, p_temp, ncol = 1, rel_heights = c(.1, .1, 1)) + paint_white_background()
ggsave(here::here("plots/Fig4C-network_hierarchy_simulation.png"), pC, width = 10, height = 3)



# Figure 4D: Hierarchy
communities_hierarchy <- read_csv(here::here("data/output/communities_hierarchy.csv")) %>% left_join(select(communities, Assembly, Community))
df_communities_hierarchy <- read_csv(here::here("data/output/df_communities_hierarchy.csv"))
pD <- communities_hierarchy_obv <- bind_rows(
    mutate(communities_hierarchy, Treatment = "experiment"),
    mutate(df_communities_hierarchy, Treatment = "simulation")
    ) %>%
    mutate(Assembly = factor(Assembly, c("self_assembly", "random_assembly"))) %>%
    filter(Metric == "h1") %>%
    ggplot(aes(x = Treatment, y = HierarchyScore), position = "dodge2") +
    geom_boxplot(width = .5, lwd = .8, outlier.color = NA) +
    geom_point(shape = 1, size = 2, stroke = .8) +
    scale_y_continuous(limits = c(0,1.01), breaks = c(0, .25, .5, .75, 1)) +
    #facet_grid(.~Treatment, scales = "free_x") +
    theme_classic() +
    #scale_color_d3() +
    theme(panel.grid.major.y = element_line(color = "grey", linetype = 2),
          panel.spacing = unit(0, "mm"),
          #legend.position = c(0.8, 0.3),
          axis.text = element_text(size = 10, color = 1),
          axis.text.x = element_text(size = 10, color = "black", angle = 15, vjust = 1, hjust = 1),
          axis.title = element_text(size = 10, color = 1),
          legend.title = element_blank(),
          plot.title = element_text(size = 10, color = 1),
          panel.border = element_rect(fill = NA, color = 1, size = 1.5)) +
    labs(x = "", y = "Score") +
    guides(color = "none") +
    ggtitle("Hierarchy")

ggsave(here::here("plots/Fig4D-hierarchy.png"), pD, width = 4, height = 4)

#
p_top <- plot_grid(pA, pB, nrow = 1, rel_widths = c(1, 3), scale = c(.7, .9), labels = c("A", "B"))
p_bottom <- plot_grid(pC, pD, nrow = 1, rel_widths = c(4, 1), scale = .9, labels = c("C", "D"), axis = "b", align = "h")
p <- plot_grid(p_top, p_bottom, nrow = 2, rel_heights = c(1,1)) + paint_white_background()
p_legend <- {p_net_hierarchy_list[[13]] + theme(legend.position = "right", legend.text = element_text(size = 10), legend.background = element_rect(fill = NA, color = NA))} %>% get_legend()
p <- ggdraw(p) + draw_plot(p_legend,.4,.6,.1,.1)
ggsave(here::here("plots/Fig4.png"), p, width = 10, height = 6)



# Test on computing the hierarchy ----
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
## 2 species, one loop, h = 0
pairs_mock <- tibble(Isolate1 = c(1),
       Isolate2 = c(2),
       InteractionType = rep("exclusion", 1),
       From = c(2),
       To = c(1))
compute_hierarchy1(pairs_mock)
tournament_rank(pairs_mock)


## 3 species, one loop, h = 0
pairs_mock <- tibble(Isolate1 = c(1,1,2),
       Isolate2 = c(2,3,3),
       InteractionType = rep("exclusion", 3),
       From = c(1,3,2),
       To = c(2,1,3))

compute_hierarchy1(pairs_mock)
tournament_rank(pairs_mock)

## 4 species, Two intransitive loop . h = 0.5
pairs_mock <- tibble(Isolate1 = c(1,1,1,2,2,3),
       Isolate2 = c(2,3,4,3,4,4),
       InteractionType = rep("exclusion", 6),
       From = c(2,1,1,2,4,3),
       To = c(1,3,4,3,2,4))

compute_hierarchy1(pairs_mock)
tournament_rank(pairs_mock)

## 4 species, Two intransitive loop . h = 0.5
pairs_mock <- tibble(Isolate1 = c(1,1,1,2,2,3),
       Isolate2 = c(2,3,4,3,4,4),
       InteractionType = rep("exclusion", 6),
       From = c(2,1,1,2,4,3),
       To = c(1,3,4,3,2,4))

compute_hierarchy1(pairs_mock)
tournament_rank(pairs_mock)


#
pairs_mock <- tibble(Isolate1 = c(1,1,1,2,2,3),
       Isolate2 = c(2,3,4,3,4,4),
       InteractionType = c("exclusion", "coexistence", "exclusion", "coexistence", "exclusion", "coexistence"),
       From = c(2,1,1,2,4,3),
       To = c(1,3,4,3,2,4))

#
n_species = 100
p_exclusion = 0.5
p_flip = 0.5


draw_h <- function(n_species, p_exclusion, p_flip = 0.5) {
    pairs_mock <- tibble(Isolate1 = 1:n_species, Isolate2 = 1:n_species) %>%
        expand(Isolate1, Isolate2) %>%
        filter(Isolate1 < Isolate2) %>%
        mutate(InteractionType = sample(c("exclusion", "coexistence"), size = n(), prob = c(p_exclusion, 1-p_exclusion), replace = T)) %>%
        mutate(From = Isolate1, To = Isolate2)
    temp_index <- which(pairs_mock$InteractionType == "exclusion")
    exchange_index <- temp_index[rbinom(length(temp_index), size = 1, p_flip) %>% as.logical]
    temp <- pairs_mock$From[exchange_index]
    pairs_mock$From[exchange_index] <- pairs_mock$To[exchange_index]
    pairs_mock$To[exchange_index] <- temp
    tournament_rank(pairs_mock)
    compute_hierarchy1(pairs_mock)

}


df <- tibble(n_species = 2:15) %>%
    slice(rep(1:n(), each = 3)) %>%
    mutate(p_exclusion = rep(c(.1, .5, 1), 14)) %>%
    slice(rep(1:n(), each = 30)) %>%
    rowwise() %>%
    mutate(h = draw_h(n_species, p_exclusion))


p <- df %>%
    drop_na() %>%
    mutate(p_exclusion = factor(p_exclusion),
           n_species = factor(n_species)) %>%
    ggplot() +
    geom_boxplot(aes(x = n_species, y = h, color = p_exclusion), outlier.color = NA) +
    geom_point(aes(x = n_species, y = h, color = p_exclusion), stroke = .5, shape = 21, size = 1,
               position = position_jitterdodge()) +
    #scale_x_continuous(breaks = 1:10) +
    #facet_grid(.~p_exclusion) +
    theme_classic() +
    theme(panel.border = element_rect(fill = NA, color = 1))

ggsave(here::here("plots/Fig4S-simulated_h.png"), p, width = 10, height = 4)

draw_h(n_species=100, p_exclusion=1, p_flip = 0.5)










if (FALSE) {

    # Figure 4A: one network in matrix form for example
    net <- net_list$C1R2
    n_nodes <- igraph::vcount(net)
    net_rank <- net %>%  activate(nodes) %>% arrange(Rank) %>% pull(Rank)
    p_net_matrix <- net %>%
        activate(edges) %>%
        mutate(fromRank = .N()$PlotRank[match(from, .N()$Isolate)], toRank = .N()$PlotRank[match(to, .N()$Isolate)]) %>%
        filter(fromRank <= toRank) %>%
        # Add self links
        bind_edges(tibble(from = 1:n_nodes, to = 1:n_nodes, fromRank = 1:n_nodes, toRank = 1:n_nodes, InteractionType = "self")) %>%
        mutate(fromRank = factor(fromRank, n_nodes:1), toRank = factor(toRank)) %>%
        as_tibble() %>%
        ggplot() +
        geom_tile(aes(x = toRank, y = fromRank, fill = InteractionType), width = 0.9, height = 0.9) +
        scale_x_discrete(position = "top", expand = c(0,0), labels = paste0("rank ", net_rank)) +
        scale_y_discrete(position = "right", expand = c(0,0), labels = paste0("rank ", rev(net_rank))) +
        scale_fill_manual(values = c(assign_interaction_color(), "self" = "black"), breaks = c("exclusion", "coexistence")) +
        theme_classic() +
        theme(axis.line = element_blank(), axis.ticks = element_blank(),
              legend.background = element_blank(), legend.position = c(.2, .2), legend.spacing.y = unit(3, "mm"),
              axis.text = element_text(color = "black"),
              plot.margin = unit(c(5,5,1,1), "mm")) +
        guides(fill = guide_legend(title = "", keywidth = 5, keyheight = 5, default.unit = "mm", byrow = T)) +
        labs(x = "", y = "") +
        # Add isolate cartoons
        draw_image(here::here("plots/cartoons/Fig1B_1.png"), x = 0.5, y = 4, vjust = -0.6, hjust = 0, clip = "on", scale = 1) +
        draw_image(here::here("plots/cartoons/Fig1B_2.png"), x = 1.5, y = 4, vjust = -0.6, hjust = 0, clip = "on", scale = 1) +
        draw_image(here::here("plots/cartoons/Fig1B_3.png"), x = 2.5, y = 4, vjust = -0.6, hjust = 0, clip = "on", scale = 1) +
        draw_image(here::here("plots/cartoons/Fig1B_4.png"), x = 3.5, y = 4, vjust = -0.6, hjust = 0, clip = "on", scale = 1) +
        draw_image(here::here("plots/cartoons/Fig1B_1.png"), x = 4.5, y = 3.5, vjust = 0, hjust = -0.8, clip = "on", scale = 1) +
        draw_image(here::here("plots/cartoons/Fig1B_2.png"), x = 4.5, y = 2.5, vjust = 0, hjust = -0.8, clip = "on", scale = 1) +
        draw_image(here::here("plots/cartoons/Fig1B_3.png"), x = 4.5, y = 1.5, vjust = 0, hjust = -0.8, clip = "on", scale = 1) +
        draw_image(here::here("plots/cartoons/Fig1B_4.png"), x = 4.5, y = 0.5, vjust = 0, hjust = -0.8, clip = "on", scale = 1)
    pA <- plot_grid(p_net_matrix, scale = .8) + paint_white_background()
    ggsave(here::here("plots/Fig4A-example_matrix.png"), pA, width = 4, height = 4)

    # Figure 4B. all matrices, empirical
    plot_adjacent_matrix_margin <- function(x) {
        plot_adjacent_matrix(x) + theme(plot.margin = grid::unit(c(5,0,3,0), "mm"))
    }
    p_net_matrix <- communities %>%
        mutate(Community = factor(Community, community_factor)) %>%
        filter(Assembly == "self_assembly") %>%
        select(Community, Richness = CommunitySize) %>%
        arrange(Community) %>%
        mutate(Network = net_list[1:13]) %>%
        mutate(Richness = Richness / max(Richness)) %>%
        rowwise() %>%
        mutate(p_net = plot_adjacent_matrix_margin(Network) %>% list())

    pB_legend <- {plot_adjacent_matrix(net_list$C11R2, show.legend = T) +
            scale_fill_manual(values = assign_interaction_color(level = "matrix"), breaks = c("self", "exclusion", "coexistence", "exclusion violating rank")) +
            theme(legend.direction = "vertical", legend.title = element_blank()) +
            guides(fill = guide_legend(keywidth = 5, keyheight = 5, default.unit = "mm"))} %>%
        get_legend()

    pB_title <- ggdraw() + draw_label("Experiment", fontface = 'bold',x = 0,hjust = 0) + theme(plot.margin = margin(5, 0, 5, 7))
    p_temp <- plot_grid(pB_title,
                        plot_grid(plotlist = p_net_matrix$p_net[1:10], scale = p_net_matrix$Richness[1:10]/max(p_net_matrix$Richness[1:10]), nrow = 1, labels = 1:10) + theme(plot.background = element_rect(fill = "grey90", color = NA)),
                        plot_grid(plotlist = p_net_matrix$p_net[11:13], scale = p_net_matrix$Richness[11:13]/max(p_net_matrix$Richness[11:13]), nrow = 1, labels = 11:13)+ theme(plot.background = element_rect(fill = "grey90", color = NA)),
                        ncol = 1, rel_heights = c(.2,1,2)) +
        paint_white_background()

    pB <- ggdraw(p_temp) +
        draw_plot(pB_legend, x = .1, y = .2, width = .1, height = .1, hjust = .5, vjust = .5)
    ggsave(here::here("plots/Fig4B-matrices.png"), pB, width = 8, height = 4)

    # Figure 4C. all matrices, simulation
    df_communities <- read_csv(here::here("data/output/df_communities.csv"), col_types = cols())
    df_communities_abundance <- read_csv(here::here("data/output/df_communities_abundance.csv"), col_types = cols())
    load("~/Dropbox/lab/invasion-network/data/output/network_simulated.Rdata")
    df_communities_ID <- df_communities %>%
        filter(Assembly == "self_assembly") %>%
        #filter(Transfer == max(Transfer), Time == max(Time)) %>%
        group_by(Community) %>%
        #mutate(RelativeAbundance = Abundance / sum(Abundance)) %>%
        #filter(RelativeAbundance > 0.01) %>%
        #summarize(Richness = n()) %>%
        filter(Richness > 3) %>% ungroup() %>%
        slice(1:10) %>%
        arrange(Richness) %>%
        mutate(CommunityLabel = factor(1:10))
    p_net_simulated_matrix <- df_communities_ID %>%
        mutate(Network = net_simulated_list[paste0("self_assembly_", Community)]) %>%
        #mutate(Richness = Richness / max(Richness)) %>%
        rowwise() %>%
        mutate(p_net = plot_adjacent_matrix_margin(Network) %>% list())
    pC_title <- ggdraw() + draw_label("Simulation", fontface = 'bold',x = 0,hjust = 0) + theme(plot.margin = margin(5, 0, 5, 7))
    pC <- plot_grid(pC_title,
                    plot_grid(plotlist = p_net_simulated_matrix$p_net[1:5], scale = p_net_simulated_matrix$Richness[1:5]/max(p_net_simulated_matrix$Richness[1:10]), nrow = 1, labels = 1:5) + theme(plot.background = element_rect(fill = "grey90", color = NA)),
                    plot_grid(plotlist = p_net_simulated_matrix$p_net[6:10], scale = p_net_simulated_matrix$Richness[6:10]/max(p_net_simulated_matrix$Richness[6:10]), nrow = 1, labels = 6:10)+ theme(plot.background = element_rect(fill = "grey90", color = NA)),
                    ncol = 1, rel_heights = c(.2,1,1)) +
        paint_white_background()
    ggsave(here::here("plots/Fig4C-matrices_simulation.png"), pC, width = 8, height = 4)


    # Figure 4D: Hierarchy
    ## Experiments
    communities_hierarchy <- read_csv(here::here("data/output/communities_hierarchy.csv")) %>% left_join(select(communities, Assembly, Community))
    ## Simulation
    df_communities_hierarchy <- read_csv(here::here("data/output/df_communities_hierarchy.csv"))

    # Binding the data
    ## Observation
    communities_hierarchy_obv <- bind_rows(
        mutate(communities_hierarchy, Treatment = "experiment"),
        mutate(df_communities_hierarchy, Treatment = "simulation")
    ) %>%
        mutate(Assembly = factor(Assembly, c("self_assembly", "random_assembly")))

    # Plot
    pD <- communities_hierarchy_obv %>%
        ggplot(aes(x = Treatment, y = HierarchyScore, color = Metric), position = "dodge2") +
        geom_boxplot(width = .5, lwd = .8) +
        geom_point(shape = 1, size = 2, stroke = .8, position = position_jitterdodge(jitter.width = .2)) +
        scale_y_continuous(limits = c(0,1.01), breaks = c(0, .25, .5, .75, 1)) +
        facet_grid(.~Treatment, scales = "free_x") +
        theme_classic() +
        scale_color_d3() +
        theme(panel.grid.major.y = element_line(color = "grey", linetype = 2),
              panel.spacing = unit(0, "mm"),
              legend.position = c(0.8, 0.3),
              axis.text = element_text(size = 10, color = 1),
              axis.text.x = element_text(size = 15, color = "black", angle = 15, vjust = 1, hjust = 1),
              axis.title = element_text(size = 15, color = 1),
              legend.title = element_blank(),
              panel.border = element_rect(fill = NA, color = 1, size = 1.5)) +
        labs(x = "", y = "Score") +
        ggtitle("Hierarchy")
    ggsave(here::here("plots/Fig4D-hierarchy.png"), pD, width = 4, height = 4)

    # Figure 4E: nontransitivity
    ## titles
    p0 <- ggdraw() + draw_label("Nontransitivity", x = .5, hjust = 0.5) + theme(plot.margin = margin(0, 0, 5, 0))
    ## motif diagrams
    load(here::here("data/output/motif_list.Rdata"))
    p_motif_list <- lapply(motif_list, function(x) plot_competitive_network(x, node_size = 3, edge_width = 1))
    p1 <- plot_grid(plotlist = p_motif_list[c(1)], nrow = 1, greedy = T)

    ## Motif count
    plot_motif_count <- function(networks_motif_summary, networks_motif_randomized_summary, motif_subset) {
        networks_motif_randomized_summary %>%
            filter(Motif %in% motif_subset) %>%
            group_by(Motif, Replicate) %>%
            mutate(p5 = quantile(Count, 0.05), p95 = quantile(Count, 0.95)) %>%
            ggplot() +
            geom_vline(xintercept = 0, color = 1) +
            geom_hline(yintercept = 0, color = 1, linetype = 2) +
            geom_histogram(aes(y = Fraction, x = after_stat(count / max(count))), alpha = .3, color = 1, binwidth = 0.005) +
            geom_point(data = filter(networks_motif_summary, Motif %in% motif_subset), aes(x = 0, y = Fraction, color = "observed network"), pch = 1, size = 2, stroke = 2, inherit.aes = F) +
            scale_color_manual(values = c("observed network" = "red")) +
            facet_grid(.~Motif) +
            scale_y_continuous(limits = c(-0.001, 0.3), breaks = scales::pretty_breaks(n=2)) +
            scale_x_continuous(breaks = scales::pretty_breaks(n=2)) +
            theme_cowplot() +
            theme(panel.background = element_rect(color = 1, size = 1.5, fill = NA),
                  panel.spacing = unit(0, "mm"), strip.text = element_blank(),
                  legend.position = "top",
                  axis.title = element_blank(), axis.text = element_text(size = 8),
                  plot.background = element_rect(fill = "white", color = NA)) +
            guides(fill = "none", color = "none") +
            labs(x = "Probability density", y = "Fraction of motif")
    }
    networks_motif <- read_csv(here::here("data/output/networks_motif.csv")) %>% filter(str_detect(Community, "C\\d"))
    networks_motif_randomized <- read_csv(here::here("data/output/networks_motif_randomized.csv")) %>% filter(str_detect(Community, "C\\d")) %>% mutate(Community = ordered(Community, communities$Community))
    networks_motif_summary <- networks_motif %>%
        group_by(Motif) %>%
        summarize(Count = sum(Count)) %>%
        mutate(TotalCount = sum(Count), Fraction = Count / TotalCount)
    networks_motif_randomized_summary <- networks_motif_randomized %>%
        group_by(Replicate, Motif) %>%
        summarize(Count = sum(Count)) %>%
        mutate(TotalCount = sum(Count), Fraction = Count / TotalCount)
    ##
    df_motif <- read_csv(here::here("data/output/df_motif.csv"))
    df_motif_randomized <- read_csv(here::here("data/output/df_motif_randomized.csv"))
    df_motif_summary <- df_motif %>%
        filter(Assembly == "self_assembly") %>%
        group_by(Motif) %>%
        summarize(Count = sum(Count)) %>%
        mutate(TotalCount = sum(Count), Fraction = Count / TotalCount)
    df_motif_randomized_summary <- df_motif_randomized %>%
        filter(Assembly == "self_assembly") %>%
        group_by(Replicate, Motif) %>%
        summarize(Count = sum(Count)) %>%
        mutate(TotalCount = sum(Count), Fraction = Count / TotalCount)
    motif_obv <- bind_rows(
        mutate(networks_motif_summary, Treatment = "experiment"),
        mutate(df_motif_summary, Treatment = "simulation")
    )
    motif_pmt <- bind_rows(
        mutate(networks_motif_randomized_summary, Treatment = "experiment"),
        mutate(df_motif_randomized_summary, Treatment = "simulation")
    )

    p2 <- motif_pmt %>%
        filter(Motif %in% 1) %>%
        group_by(Motif, Replicate) %>%
        mutate(p5 = quantile(Count, 0.05), p95 = quantile(Count, 0.95)) %>%
        ggplot() +
        geom_vline(xintercept = 0, color = 1) +
        geom_hline(yintercept = 0, color = 1, linetype = 2) +
        geom_histogram(aes(y = Fraction, x = after_stat(count / max(count)), color = "permutation"), alpha = .3, binwidth = 0.005) +
        geom_point(data = filter(motif_obv, Motif %in% 1), aes(x = 0, y = Fraction, color = "observation"), pch = 1, size = 2, stroke = 2, inherit.aes = F) +
        scale_color_manual(values = c("observation" = "red", "permutation" = "black")) +
        facet_grid(.~Treatment) +
        scale_y_continuous(limits = c(-0.001, 0.3), breaks = scales::pretty_breaks(n=2)) +
        scale_x_continuous(breaks = scales::pretty_breaks(n=2)) +
        theme_classic() +
        theme(panel.background = element_rect(color = 1, size = 1.5, fill = NA),
              panel.spacing = unit(0, "mm"),
              legend.position = c(0.8, 0.8),
              legend.title = element_blank(),
              legend.background = element_blank(),
              axis.text = element_text(size = 8),
              plot.title = element_text(size = 15, face = "plain"),
              plot.background = element_rect(fill = "white", color = NA)) +
        guides(fill = "none") +
        labs(x = "Probability density", y = "Fraction") +
        ggtitle("Nontransitivity")

    #pE <- plot_grid(p0, p2, ncol = 1, rel_heights = c(.2, 3), axis = "lr", align = "v") + paint_white_background()
    pE <- p2
    ggsave(here::here("plots/Fig4E-motifs.png"), pE, width = 4, height = 4)


    #
    p_top <- plot_grid(pA, pB, nrow = 1, labels = LETTERS[1:2], rel_widths = c(1,2), scale = c(.8, .9))
    p_bottom_2 <- plot_grid(pD, pE, nrow = 1, labels = LETTERS[4:5], rel_widths = c(1,1), scale = c(.9, .9, .9), axis = "tb", align = "h")
    p_bottom <- plot_grid(pC, p_bottom_2, nrow = 1, labels = c("C", ""), rel_widths = c(1,1), scale = c(.9, 1))
    p <- plot_grid(p_top, p_bottom, nrow = 2) + paint_white_background()
    ggsave(here::here("plots/Fig4.png"), p, width = 14, height = 8)

}

if (FALSE){




    # Figure 3D: motif of all networks
    ## titles
    p0_1 <- ggdraw() + draw_label("Nontransitive", x = .5, hjust = 0.5) + theme(plot.margin = margin(0, 0, 5, 0))
    p0_2 <- ggdraw() + draw_label("Hierarchical", x = .5, hjust = 0.5) + theme(plot.margin = margin(0, 0, 5, 0))
    ## motif diagrams
    load(here::here("data/output/motif_list.Rdata"))
    p_motif_list <- lapply(motif_list, function(x) plot_competitive_network(x, node_size = 2, edge_width = 1))
    p1_1 <- plot_grid(plotlist = p_motif_list[c(1)], nrow = 1, greedy = T)
    p1_2 <- plot_grid(plotlist = p_motif_list[c(2,3,5)], nrow = 1, greedy = T)

    ## Motif count
    networks_motif <- read_csv(here::here("data/output/networks_motif.csv")) %>% filter(str_detect(Community, "C\\d"))
    networks_motif_randomized <- read_csv(here::here("data/output/networks_motif_randomized.csv")) %>% filter(str_detect(Community, "C\\d")) %>% mutate(Community = ordered(Community, communities$Community))
    networks_motif_summary <- networks_motif %>%
        group_by(Motif) %>%
        summarize(Count = sum(Count)) %>%
        mutate(TotalCount = sum(Count), Fraction = Count / TotalCount)
    networks_motif_randomized_summary <- networks_motif_randomized %>%
        group_by(Replicate, Motif) %>%
        summarize(Count = sum(Count)) %>%
        mutate(TotalCount = sum(Count), Fraction = Count / TotalCount)
    plot_motif_count <- function(networks_motif_summary, networks_motif_randomized_summary, motif_subset) {
        networks_motif_randomized_summary %>%
            filter(Motif %in% motif_subset) %>%
            group_by(Motif, Replicate) %>%
            mutate(p5 = quantile(Count, 0.05), p95 = quantile(Count, 0.95)) %>%
            ggplot() +
            geom_vline(xintercept = 0, color = 1) +
            geom_hline(yintercept = 0, color = 1, linetype = 2) +
            geom_histogram(aes(y = Fraction, x = after_stat(count / max(count))), alpha = .3, color = 1, binwidth = 0.005) +
            geom_point(data = filter(networks_motif_summary, Motif %in% motif_subset), aes(x = 0, y = Fraction, color = "observed network"), pch = 1, size = 2, stroke = 2, inherit.aes = F) +
            scale_color_manual(values = c("observed network" = "red")) +
            facet_grid(.~Motif) +
            scale_y_continuous(limits = c(-0.001, 0.3), breaks = scales::pretty_breaks(n=2)) +
            scale_x_continuous(breaks = scales::pretty_breaks(n=2)) +
            theme_cowplot() +
            theme(panel.background = element_rect(color = 1, size = 1.5, fill = NA),
                  panel.spacing = unit(0, "mm"), strip.text = element_blank(),
                  legend.position = "top",
                  axis.title = element_blank(), axis.text = element_text(size = 8),
                  plot.background = element_rect(fill = "white", color = NA)) +
            guides(fill = "none", color = "none") +
            labs(x = "Probability density", y = "Fraction of motif")
    }
    p2_1 <- plot_motif_count(networks_motif_summary, networks_motif_randomized_summary, 1)
    p2_2 <- plot_motif_count(networks_motif_summary, networks_motif_randomized_summary, c(2,3,5))

    ##
    p_left <- plot_grid(p0_1, p1_1, p2_1, ncol = 1, rel_heights = c(.4, 1,2.5), axis = "tblr", align = "v")
    p_right <- plot_grid(p0_2, p1_2, p2_2, ncol = 1, rel_heights = c(.4, 1,2.5), axis = "tblr", align = "v")
    p_upper <- plot_grid(p_left, p_right, nrow = 1, rel_widths = c(1.5,3), axis = "tblr", align = "h") + paint_white_background()
    pD_xlab <- ggdraw() + draw_label("Probability density", x = 0.5, hjust = .5) + theme(plot.margin = margin(0, 0, 0, 0))
    pD <- plot_grid(p_upper, pD_xlab, ncol = 1, rel_heights = c(1, .1), scale = c(1, .5)) + paint_white_background()
    ggsave(here::here("plots/Fig3D-motifs.png"), pD, width = 5, height = 4)


    # Figure 3E: diagonal analysis. One community for example
    networks_diag <- read_csv(here::here("data/output/networks_diag.csv"), col_types = cols())
    networks_diag_randomized <- read_csv(here::here("data/output/networks_diag_randomized.csv"), col_types = cols())

    ## Permutation
    networks_diag_randomized_sum <- networks_diag_randomized %>%
        filter(str_detect(Community, "C\\d")) %>%
        group_by(Replicate, RankDifference) %>%
        summarize(Count = sum(Count), TotalCount = sum(TotalCount)) %>%
        mutate(Fraction = Count/TotalCount)
    networks_diag_randomized_sample_size <- networks_diag_randomized_sum %>%
        group_by(RankDifference) %>%
        summarize(Count = n())

    ## Observation
    networks_diag_sum <- networks_diag %>%
        filter(str_detect(Community, "C\\d")) %>%
        group_by(RankDifference) %>%
        summarize(ObservedCount = sum(Count), ObservedTotalCount = sum(TotalCount)) %>%
        mutate(ObservedFraction = ObservedCount/ObservedTotalCount)

    ## Statistics
    stat_diag <- networks_diag_randomized_sum %>%
        group_by(RankDifference) %>%
        # Find percentile
        bind_rows(tibble(Replicate = 0, RankDifference = 0:11, Count = 0)) %>%
        left_join(networks_diag_sum) %>%
        arrange(RankDifference, desc(Fraction)) %>%
        mutate(Percentile = (1:n())/n()) %>%
        filter(Count <= ObservedCount) %>%
        group_by(RankDifference) %>%
        slice(1) %>%
        select(RankDifference, Percentile) %>%
        # Asterisk
        mutate(Significance = case_when(Percentile < 0.001 | Percentile > 0.999 ~ "***",
                                        Percentile < 0.01 | Percentile > 0.99 ~ "**",
                                        Percentile < 0.05 | Percentile > 0.95 ~ "*",
                                        Percentile > 0.05 & Percentile < 0.95 ~ "n.s."),
               Sign = case_when(Percentile < 0.05 ~ "top",
                                Percentile > 0.95 ~ "bottom",
                                Percentile > 0.05 & Percentile < 0.95 ~ "n.s."))

    pE <- networks_diag_sum %>%
        ggplot() +
        # Sample size
        #geom_text(data = networks_diag_randomized_sample_size, aes(x = RankDifference, y = Inf, label = paste0("n=", Count)), vjust = 3.5) +
        # Asterisk
        geom_text(data = stat_diag, aes(x = RankDifference, y = Inf, label = Significance), vjust = 2) +
        geom_rect(data = stat_diag, aes(xmin = RankDifference-0.5, xmax = RankDifference+0.5, fill = Sign), ymin = -Inf, ymax = Inf, alpha = .2) +
        # Random networks
        geom_boxplot(data = networks_diag_randomized_sum, aes(x = RankDifference, y = Fraction, group = RankDifference, color = "permutation"),
                     outlier.size = 1) +
        geom_jitter(data = networks_diag_randomized_sum, aes(x = RankDifference, y = Fraction, group = RankDifference, color = "permutation"),
                    size = .2, alpha = 0.5, width = .3, shape = 21, height = .01) +
        # Observed networks
        geom_point(aes(x = RankDifference, y = ObservedFraction, group = RankDifference, color = "observation"), size = 2) +
        geom_line(aes(x = RankDifference, y = ObservedFraction, color = "observation")) +
        scale_x_continuous(breaks = 0:11, expand = c(0,0)) +
        scale_y_continuous(limits = c(-0.02, 1.1), breaks = c(0, .5, 1)) +
        scale_color_manual(values = c("observation" = "red", "permutation" = "black"))+
        scale_fill_manual(values = c("top" = "blue", "bottom" = "red", "n.s." = "grey")) +
        theme_classic() +
        theme(legend.position = "top", legend.title = element_blank(),
              legend.text = element_text(size = 12),
              axis.title.x = element_text(size = 15, color = 1),
              panel.border = element_rect(fill = NA, color = 1, size = 1.5)) +
        guides(fill = "none") +
        labs(x = "|i-j|", y = "Fraction of pairwise coexitence")

    ggsave(here::here("plots/Fig3E-diagonal_analysis.png"), pE, width = 4, height = 4)


    #
    p_top <- plot_grid(pA, pB, nrow = 1, labels = LETTERS[1:2], rel_widths = c(1,2), scale = c(.8, .9))
    p_bottom <- plot_grid(pC, pD, pE, nrow = 1, labels = LETTERS[3:5], rel_widths = c(1,2,2), scale = c(.9, .9, .9))
    p <- plot_grid(p_top, p_bottom, nrow = 2) + paint_white_background()
    ggsave(here::here("plots/Fig3.png"), p, width = 12, height = 8)

}






temp <- net_list[[8]] %>%
    # Node position
    activate(nodes) %>%
    mutate(y = -Rank) %>%
    group_by(Rank) %>%
    mutate(x = {seq(0, 1, length.out = n() + 2) %>% `[`(c(-1, -length(.)))}) %>% # + rnorm(n(), 0, .5)) %>%
    ungroup() %>%
    # Filter out edges
    activate(edges) %>%
    filter(InteractionType == "exclusion") %>%
    mutate(Temp = sample(c(-1, 1), size = n(), replace = T))
strength_angle <- as_tibble(temp)$Temp * 0.01

temp %>%
    ggraph(layout = "nicely") +
    geom_segment(aes(x = 0.5, xend = 0.5, y = -Inf, yend = -1), color = "grey80") +
    geom_node_point(size = node_size, shape = 21, fill = "grey", stroke = node_size/5, color = "black") +
    geom_edge_arc(strength = strength_angle, force_flip = T, aes(color = InteractionType), width = edge_width,
                  arrow = arrow(length = unit(edge_width, "mm"), type = "closed", angle = 30, ends = "last"),
                  start_cap = circle(node_size/2, "mm"),
                  end_cap = circle(node_size/2, "mm")) +
    scale_edge_color_manual(values = interaction_color) +
    scale_x_continuous(limits = c(-0.1, 1.1), expand = c(0,0)) +
    scale_y_continuous(limits = c(-13, 0), breaks = -12:-1, labels = 12:1) +
    theme_void() +
    #theme_bw() +
    theme(
        panel.grid.major.y = element_line(color = "grey80"),
        #panel.grid.minor.y = element_line(color = "grey80"),
        #panel.grid.major.x = element_blank(),
        #panel.grid.minor.x = element_blank(),
        legend.position = "none",
        #panel.border = element_rect(color = "grey80", fill = NA),
        legend.title = element_blank(),
        axis.title = element_blank(),
        strip.text = element_blank(),
        plot.margin=unit(c(0,0,0,0),"mm")
    ) +
    labs(y = "Rank")



# Figure S1. glucose isolate abundance  ----------------------------------------------------------------------------------------------------
## Panel A. cartoon of self-assembly experiment and isolate abundance
p1 <- ggdraw() + draw_image(here::here("plots/cartoons/FigS1.png")) + theme(plot.background = element_rect(fill = "white", color = NA))
## Panel B. isolate abundance
temp <- sequences_abundance %>%
    filter(AlignmentType == "local") %>%
    filter(AllowMismatch == Inf) %>%
    filter(BasePairMismatch <= 4) %>%
    mutate(Community = ordered(Community,  communities$Community))
color_sets <- tibble(Color = c("yellow", "deepskyblue3", "blue", "darkorchid2", "firebrick", "orange2", "grey"),
                     Family = c("Aeromonadaceae", "Enterobacteriaceae", "Moraxellaceae", "Pseudomonadaceae","Comamonadaceae","Alcaligenaceae", "Sphingobacteriaceae"))
p2 <- temp %>%
    mutate(Community = factor(Community, community_factor)) %>%
    arrange(Community) %>%
    ggplot() +
    geom_bar(aes(x = Community, y = RelativeAbundance, fill = Family), size = .3, color = "grey30", position = "stack", stat = "identity") +
    theme_bw() +
    scale_fill_manual(values = setNames(color_sets$Color, color_sets$Family)) +
    scale_x_discrete(labels = 1:13) +
    scale_y_continuous(breaks = c(0, .5, 1), expand = c(0,0), limits = c(0, 1)) +
    theme(#axis.text.x = element_text(angle = 90, vjust = 0.5),
          panel.grid = element_blank(),
          axis.title = element_text(size = 12),
          panel.border = element_rect(color = 1, size = 1)) +
    labs(y = "Relative abundance")
## Stats
temp %>%
    group_by(Community) %>%
    summarize(Total = sum(RelativeAbundance)) %>%
    summarize(Mean = mean(Total))

#p <- plot_grid(p1, p2, nrow = 1, rel_widths = c(1,2), labels = LETTERS[1:2], scale = c(1, .8)) + theme(plot.background = element_rect(fill = "white", color = NA))
p <- p2
ggsave(here::here("plots/FigS1-glucose_community_bar.png"), p, width = 5, height = 3)





# Figure S2. Raw pairwise frequency plots  ----------------------------------------------------------------------------------------------------
pairs_example_freq <- pairs %>%
    filter(str_detect(Community, "C\\d+")) %>%
    left_join(isolates %>% select(ID, Rank, PlotRank, Community) %>% rename_with(~ paste0(., 1), -Community)) %>%
    left_join(isolates %>% select(ID, Rank, PlotRank, Community) %>% rename_with(~ paste0(., 2), -Community)) %>%
    select(Community, starts_with("Isolate"), starts_with("Interaction"), starts_with("PlotRank"), starts_with("Rank")) %>%
    left_join(pairs_freq, by = c("Community", "Isolate1", "Isolate2")) %>%
    mutate(Isolate1InitialODFreq = factor(Isolate1InitialODFreq)) %>%
    mutate(InteractionType = factor(InteractionType, c("exclusion", "coexistence")))
plot_example_freq <- function(pairs_freq) {
    # Extract params
    comm <- unique(pairs_freq$Community)
    isolate1 <- unique(pairs_freq$Isolate1)
    isolate2 <- unique(pairs_freq$Isolate2)
    interaction_type <- pairs %>%
        filter(Community == comm, Isolate1 == isolate1, Isolate2 == isolate2) %>%
        pull(InteractionType)

    #
    pairs_freq %>%
        mutate(Isolate1InitialODFreq = factor(Isolate1InitialODFreq, c(95,50,5))) %>%
        ggplot(aes(x = Time, y = Isolate1MeasuredFreq, color = Isolate1InitialODFreq, group = Isolate1InitialODFreq)) +
        geom_point(size = 2) +
        geom_line(size = 1) +
        scale_y_continuous(breaks = c(0, .5, 1), limits = c(0,1)) +
        scale_color_manual(values = frequency_color, label = c("95%", "50%", "5%")) +
        theme_bw() +
        theme(panel.spacing = unit(2, "mm"),
              panel.border = element_rect(color = 1, fill = NA, size = 1),
              panel.grid = element_blank(),
              panel.grid.minor.y = element_blank(),
              axis.title = element_blank(), axis.text = element_blank(),
              axis.ticks = element_blank(),
              #panel.background = element_rect(fill = "white"),
              panel.background = element_rect(fill = alpha(ifelse(interaction_type == "coexistence", "#557BAA",
                                                                  ifelse(interaction_type == "exclusion", "#DB7469", NA)), 0.5)),
              plot.background = element_blank()) +
        guides(color = "none") +
        labs(x = "Time", y = "Frequency") +
        #        ggtitle(paste0(isolate1, "-", isolate2)) +
        NULL
}
plot_community_freq <- function(pairs_example_freq) {
    comm <- unique(pairs_example_freq$Community)
    # Make figures
    temp_freq <- pairs_example_freq %>%
        filter(Community == comm) %>%
        #mutate(Isolate1 = ordered(Isolate1, temp), Isolate2 = ordered(Isolate2, temp)) %>%
        rowwise() %>%
        mutate(PlotRankOrder = paste0(min(PlotRank1, PlotRank2), "_", max(PlotRank1, PlotRank2))) %>%
        mutate(PlotRankOrder = ordered(PlotRankOrder, paste0(rep(1:12, each = 12), "_", rep(1:12, 12)))) %>%
        arrange(PlotRankOrder)
    temp_list <- temp_freq %>%
        as_tibble() %>%
        arrange(PlotRankOrder) %>%
        group_split(PlotRankOrder) %>%
        lapply(plot_example_freq)

    # Grid
    n_isolate <- communities$CommunitySize[communities$Community == comm]
    m <- matrix(NA, n_isolate-1, n_isolate-1)
    m[lower.tri(m, diag = T)] <- 1:choose(n_isolate,2)
    m <- t(m)
    arrangeGrob(grobs = temp_list, layout_matrix = m)

}
p_list <- pairs_example_freq %>%
    mutate(Community = factor(Community, community_factor)) %>%
    group_split(Community) %>%
    lapply(plot_community_freq)

p_top <- plot_grid(plotlist = p_list[1:10], ncol = 3, scale = communities_size[1:10]/6, labels = 1:10) + paint_white_background()
p_body <- plot_grid(plotlist = p_list[11:12], ncol = 1, scale = communities_size[11:13]/10, labels = 11:12) + paint_white_background()
p_bottom <- plot_grid(plotlist = p_list[13], ncol = 1, scale = .9, labels = 13) + paint_white_background()
ggsave(here::here("plots/FigS2-1-frequency_plots.png"), p_top, width = 9, height = 12)
ggsave(here::here("plots/FigS2-2-frequency_plots.png"), p_body, width = 9, height = 12)
ggsave(here::here("plots/FigS2-3-frequency_plots.png"), p_bottom, width = 10, height = 10)
#pS16 <- plot_grid(p_top, p_bottom, nrow = 2, rel_heights = c(1.5,1)) + paint_white_background()
#ggsave(here::here("plots/Fig1S16-frequency_plots.png"), pS16, width = 20, height = 20)

## Number of pairs using CASEU
pairs_freq %>%
    filter(str_detect(Community, "C\\d+R\\d+")) %>%
    filter(Time == "T8") %>%
    #group_by(Community, Isolate1, Isolate2) %>%
    unite("Pair", Community, Isolate1, Isolate2, sep = "_") %>%
    select(Pair, Isolate1InitialODFreq, RawDataType) %>%
    pivot_wider(names_from = Isolate1InitialODFreq, values_from = RawDataType, names_prefix = "f") %>%
    mutate(Method = ifelse(f5 == "Sanger" | f50 == "Sanger" | f95 == "Sanger", "sanger",
                           ifelse(f5 == "CFU" & f50 == "CFU" & f95 == "CFU", "cfu", NA))) %>%
    group_by(Method) %>%
    summarize(Count = n())


# Figure S3. Finer pairwise competition outcomes  ----------------------------------------------------------------------------------------------------
pairs_example_outcomes_finer <- read_csv(here::here("data/output/pairs_example_outcomes_finer.csv"))
## The frequencies of coexistence vs. exclusion
p1 <- pairs %>%
    filter(Assembly == "self_assembly") %>%
    mutate(InteractionTypeFiner = factor(InteractionTypeFiner, c("competitive exclusion", "mutual exclusion", "stable coexistence", "frequency-dependent coexistence", "neutrality"))) %>%
    group_by(InteractionTypeFiner) %>%
    summarize(Count = n()) %>%
    ungroup() %>%
    mutate(Fraction = Count / sum(Count)) %>%
    ggplot() +
    geom_col(aes(x = InteractionTypeFiner, y = Count, fill = InteractionTypeFiner), color = 1) +
    geom_text(x = -Inf, y = Inf, label = paste0("n = ", sum(temp$Count)), vjust = 1, hjust = -0.1) +
    geom_text(aes(x = InteractionTypeFiner, y = Count, label = paste0(round(Fraction, 3) * 100,"%")), nudge_y = 6) +
    scale_fill_manual(values = assign_interaction_color(level = "finer")) +
    scale_x_discrete(breaks = c("competitive exclusion", "stable coexistence", "neutrality", "mutual exclusion", "frequency-dependent coexistence"),
                     labels = c("competitive\nexclusion", "stable\ncoexistence", "neutrality", "mutual\nexclusion", "frequency-dependent\ncoexistence")) +
    scale_y_continuous(limits = c(0, 150), expand = c(0,0)) +
    theme_classic() +
    theme(axis.title.x = element_blank(), axis.title.y = element_text(size = 10),
          axis.text.x = element_text(size = 10, color = "black"), axis.text.y = element_text(color = "black"),
          legend.position = "none") +
    labs(x = "", y = "Number of pairs", fill = "")
## Pairs example dynamics
p2 <- pairs_example_outcomes_finer %>%
    left_join(pairs_freq, by = c("Community", "Isolate1", "Isolate2")) %>%
    mutate(Isolate1InitialODFreq = factor(Isolate1InitialODFreq)) %>%
    mutate(InteractionTypeFiner = factor(InteractionTypeFiner, c("competitive exclusion", "mutual exclusion", "stable coexistence", "frequency-dependent coexistence", "neutrality"))) %>%
    ggplot(aes(x = Time, y = Isolate1MeasuredFreq, color = Isolate1InitialODFreq, group = Isolate1InitialODFreq)) +
    geom_point(size = 2) +
    geom_line(size = 1) +
    scale_y_continuous(breaks = c(0, .5, 1), limits = c(0,1)) +
    scale_color_manual(values = frequency_color, label = c("95%", "50%", "5%")) +
    facet_grid(.~InteractionTypeFiner) +
    theme_bw() +
    theme(panel.spacing = unit(2, "mm"), strip.text.x = element_blank(),
          panel.border = element_rect(color = 1, fill = NA, size = 1)) +
    guides(color = "none") +
    labs(x = "Time", y = "Frequency")

p <- plot_grid(p1, p2, ncol = 1, axis = "lf", align = "h", rel_heights = c(2,1))
ggsave(here::here("plots/FigS3-pairwise_outcomes_finer.png"), p, width = 7, height = 4)


# Figure S4. pairwise outcomes per community  ----------------------------------------------------------------------------------------------------
pairs_count <- pairs %>%
    filter(Assembly == "self_assembly") %>%
    group_by(Community) %>%
    summarize(Count = n())
p <- pairs %>%
    mutate(Community = factor(Community, communities$Community)) %>%
    filter(Assembly == "self_assembly") %>%
    ggplot() +
    geom_bar(aes(x = Community, fill = InteractionType), color = 1, position = position_fill(), size = .5) +
    geom_text(data = pairs_count, aes(x = Community, y = 1, label = Count), vjust = -.5) +
    geom_rect(xmin = -Inf, xmax = Inf, ymin = 0, ymax = 1, color = grey(0.1), fill = NA, size = .5) +
    scale_fill_manual(values = assign_interaction_color()) +
    scale_y_continuous(breaks = c(0,.5,1), limit = c(0, 1.15), expand = c(0,0)) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
          axis.line.y = element_blank(),
          legend.title = element_blank(), legend.position = "top") +
    labs(y = "Fraction")
ggsave(here::here("plots/FigS4-pairs_counts.png"), p, width = 4, height = 3)



# Figure S5. Hierarchy metrics per community ----------------------------------------------------------------------------------------------------
communities_hierarchy <- read_csv(here::here("data/output/communities_hierarchy.csv")) %>% mutate(Community = factor(Community, community_factor)) %>% arrange(Community)
communities_hierarchy_randomized <- read_csv(here::here("data/output/communities_hierarchy_randomized.csv")) %>% mutate(Community = factor(Community, community_factor)) %>% arrange(Community, Replicate)

communities_hierarchy <- communities_hierarchy %>%
    rename(HierarchyScoreObv = HierarchyScore) %>%
    right_join(communities_hierarchy_randomized) %>%
    group_by(Community, Metric) %>%
    arrange(desc(HierarchyScore)) %>%
    mutate(Percentile = (1:n())/n()) %>%
    filter(HierarchyScoreObv > HierarchyScore) %>%
    slice(1) %>%
    select(Community, Metric, HierarchyScoreObv, Percentile) %>%
    mutate(Significance = ifelse(Percentile < 0.05 | Percentile > 0.95, T, F))

p1 <- communities_hierarchy %>%
    filter(Metric == "h1") %>%
    ggplot() +
    geom_rect(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, aes(fill = Significance), alpha = .5) +
    geom_histogram(data = filter(communities_hierarchy_randomized, Metric == "h1"), aes(y = HierarchyScore), color = 1, fill = "white") +
    geom_hline(aes(yintercept = HierarchyScoreObv), color = "red") +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 3)) +
    scale_y_continuous(breaks = c(0.5, 0.75, 1)) +
    scale_fill_manual(values = c("TRUE" = "#DB7469", "FALSE" = "white")) +
    #facet_wrap(Community ~., ncol = 2, scales = "free_x", dir = "v") +
    facet_wrap(Community~., ncol = 4, scales = "free_x", dir = "v") +
    theme_classic() +
    theme(panel.border = element_rect(color = 1, fill = NA, size = 1.5)) +
    guides(fill = "none") +
    labs(x = "Probability density", y = "Hierarchy score") +
    ggtitle("Rank-following pairs")

p2 <- communities_hierarchy %>%
    filter(Metric == "h2") %>%
    ggplot() +
    geom_rect(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, aes(fill = Significance), alpha = .5) +
    geom_histogram(data = filter(communities_hierarchy_randomized, Metric == "h2"), aes(y = HierarchyScore), color = 1, fill = "white") +
    geom_hline(aes(yintercept = HierarchyScoreObv), color = "red") +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 3)) +
    scale_y_continuous(breaks = c(0.5, 0.75, 1)) +
    scale_fill_manual(values = c("TRUE" = "#DB7469", "FALSE" = "white")) +
    #facet_wrap(Community ~., ncol = 2, scales = "free_x", dir = "v") +
    facet_wrap(Community~., ncol = 4, scales = "free_x", dir = "v") +
    theme_classic() +
    theme(panel.border = element_rect(color = 1, fill = NA, size = 1.5)) +
    guides(fill = "none") +
    labs(x = "Probability density", y = "Hierarchy score") +
    ggtitle("Final relative frequency")

p <- plot_grid(p1, p2, ncol = 2, scale = .95, labels = c("A", "B")) + paint_white_background()
ggsave(here::here("plots/FigS5-hierarchy.png"), p, width = 10, height = 8)
if (FALSE) {

    plot_top <- function(x) {
        communities_hierarchy %>%
            filter(Metric == x) %>%
            ggplot(aes(x = Metric, y = HierarchyScore)) +
            geom_boxplot(width = .3) +
            geom_jitter(shape = 1, size = 2, width = .1) +
            scale_x_discrete (position = "top", labels = c("1" = "Higgins et al 2017", "2" = "Rank-following pairs", "3" = "Fraction of transitive motifs")) +
            scale_y_continuous(limits = c(0,1), breaks = c(0, 0.5, 1)) +
            facet_wrap(.~Metric, scales = "free", nrow = 1) +
            theme_classic() +
            theme(axis.title.x = element_blank(), strip.text = element_blank(),
                  axis.text.x = element_text(size = 13), axis.ticks.x = element_blank(),
                  panel.border = element_rect(fill = NA, color = 1, size = 1.5)) +
            labs(x = "Metrics", y = "Score")
    }
    plot_bottom <- function(x) {
        communities_hierarchy %>%
            filter(Metric == x) %>%
            ggplot() +
            geom_rect(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, aes(fill = Significance), alpha = .5) +
            geom_histogram(data = communities_hierarchy_randomized, aes(y = HierarchyScore1), color = 1, fill = "white") +
            geom_hline(aes(yintercept = HierarchyScore), color = "red") +
            geom_text(aes(x = Inf, y = -Inf, label = paste0("p=", Percentile)), vjust = -1, hjust = 1.5) +
            scale_x_continuous(breaks = scales::pretty_breaks(n = 3)) +
            scale_y_continuous(breaks = c(0.5, 0.75, 1)) +
            scale_fill_manual(values = c("TRUE" = "#DB7469", "FALSE" = "white")) +
            facet_wrap(Community ~., ncol = 2, scales = "free_x", dir = "v") +
            theme_classic() +
            theme(panel.grid.major = element_blank(), legend.position = "none", panel.border = element_rect(color = 1, fill = NA, size = 1.5)) +
            labs(x = "Count", y = "Hierarchy score")
    }
    p_top <- plot_grid(plot_top(1), plot_top(2), plot_top(3), nrow = 1, labels = LETTERS[1:3], align = "hv")
    p_bottom <- plot_grid(plot_bottom(1), plot_bottom(2), plot_bottom(3), nrow = 1)

    p <- plot_grid(p_top, p_bottom, ncol = 1, rel_heights = c(1, 4))

}



# Figure S6. total motif count and example motifs  ----------------------------------------------------------------------------------------------------
## motif diagram
load(here::here("data/output/motif_list.Rdata"))
p_motif_list <- lapply(motif_list, function(x) plot_competitive_network(x, node_size = 3))
p1 <- plot_grid(plotlist = p_motif_list, nrow = 1)

## motif counts
networks_motif <- read_csv(here::here("data/output/networks_motif.csv")) %>% filter(str_detect(Community, "C\\d"))
networks_motif_randomized <- read_csv(here::here("data/output/networks_motif_randomized.csv")) %>% filter(str_detect(Community, "C\\d")) %>% mutate(Community = ordered(Community, communities$Community))
networks_motif_randomized_total_percentile <- networks_motif_randomized %>%
    group_by(Replicate, Motif) %>%
    summarize(Count = sum(Count)) %>%
    group_by(Motif) %>%
    arrange(Motif, Count) %>%
    slice(c(1000 * 0.05, 1000 * 0.95)) %>%
    mutate(Percentile = c("p5", "p95")) %>%
    select(Motif, Count, Percentile)
networks_motif_total <- networks_motif %>%
    group_by(Motif) %>%
    summarize(Count = sum(Count)) %>%
    left_join(pivot_wider(networks_motif_randomized_total_percentile, names_from = Percentile, values_from = Count)) %>%
    mutate(Sign = case_when(Count > p95 ~ "top",
                            Count < p5 ~ "bottom",
                            Count < p95 & Count > p5 ~ "n.s."))

p2 <- networks_motif_randomized %>%
    group_by(Motif, Replicate) %>%
    summarize(Count = sum(Count)) %>%
    group_by(Motif) %>%
    # 5% and 95% percentiles in randomized networks
    #mutate(p5 = quantile(Count, 0.05), p95 = quantile(Count, 0.95), ColoredTails = ifelse(Count <= p5, "tail", ifelse(Count >= p95, "head", "body"))) %>%
    ggplot() +
    geom_rect(data = networks_motif_total, aes(fill = Sign), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, alpha = .2) +
    geom_vline(xintercept = 0, color = 1) +
    geom_histogram(aes(y = Count), binwidth = 2, color = 1, fill = "white") +
    geom_point(data = networks_motif_total, x = 0, aes(y = Count, color = "observed network"), pch = 1, size = 2, stroke = 2, inherit.aes = F) +
    scale_color_manual(values = c("observed network" = "red")) +
    facet_grid(.~Motif, scales = "free_x") +
    scale_fill_manual(values = c("top" = "blue", "bottom" = "red", "n.s." = "grey")) +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 2)) +
    theme_cowplot() +
    theme(panel.background = element_rect(color = 1, size = 1), panel.spacing = unit(0, "mm"),
          strip.background = element_rect(color = NA, fill = NA, size = 1)) +
    guides(color = "none", fill = "none") +
    labs(x = "Probability density", y = "Motif count")

p <- plot_grid(p1, p2, ncol = 1, rel_heights = c(1,3), axis = "lr", align = "v", labels = LETTERS[1:2]) + paint_white_background()
ggsave(here::here("plots/FigS6-total_motif_counts.png"), p, width = 10, height = 5)


# Figure S7. matrix, network, and motif per community ----------------------------------------------------------------------------------------------------
## Matrix and graph
networks_motif <- read_csv(here::here("data/output/networks_motif.csv")) %>% filter(str_detect(Community, "C\\d"))
networks_motif_randomized <- read_csv(here::here("data/output/networks_motif_randomized.csv")) %>% filter(str_detect(Community, "C\\d")) %>% mutate(Community = ordered(Community, communities$Community))
p_net_matrix_list <- lapply(net_list, function(x) plot_adjacent_matrix(x) + theme(plot.margin = grid::unit(c(5,0,3,0), "mm")))
p_net_list <- lapply(net_list, function(x) plot_competitive_network(x, node_size = 1, edge_width = 1) + theme(plot.background = element_rect(fill = NA)))
p_list <- rep(list(NA), length(p_net_list))
for (i in 1:length(net_list)) p_list[[i]] <- ggdraw(p_net_matrix_list[[i]]) + draw_plot(plot = p_net_list[[i]], x = -.1, y = -.1, width = 0.7, height = 0.7)
## Motif count
plot_motif_count <- function (x = 1) {
    motif_randomized_subset <- networks_motif_randomized_percentile %>%
        filter(Community %in% communities$Community[x]) %>%
        mutate(Community = factor(Community, communities$Community[x]))
    motif_community_subset <- networks_motif %>%
        filter(Community %in% communities$Community[x]) %>%
        mutate(Community = factor(Community, communities$Community[x]))

    ggplot() +
        # 5% and 95% percentiles in randomized networks
        geom_point(data = motif_randomized_subset, aes(x = Motif, y = Count, group = Motif, color = "randomized network")) +
        geom_segment(data = motif_randomized_subset %>% pivot_wider(id_cols = c(Community, Motif), names_from = Percentile, values_from = Count),
                     aes(x = Motif, xend = Motif, y = p5, yend = p95, color = "randomized network")) +
        # Observations
        geom_point(data = motif_community_subset, aes(x = Motif, y = Count, color = "observed network")) +
        scale_x_continuous(breaks = 1:7) +
        scale_color_manual(values = c("observed network" = "red", "randomized network" = "black"))+
        #facet_wrap(Community ~., scale = "free_y", nrow = 1)  +
        theme_classic() +
        theme(panel.background = element_rect(color = 1, size = 1), legend.position = "none")
}
p_motif_count_list <- rep(list(NA), length(p_net_list))
for (i in 1:length(p_net_list)) {
    if (i %in% c(1, 6)) p_motif_count_list[[i]] <- plot_motif_count(i) + theme(axis.title.x = element_blank())
    if (i %in% c(2:5, 7:10)) p_motif_count_list[[i]] <- plot_motif_count(i) + theme(axis.title = element_blank())
    if (i == 11) p_motif_count_list[[i]] <- plot_motif_count(i)
    if (i %in% 12:13) p_motif_count_list[[i]] <- plot_motif_count(i) + theme(axis.title.y = element_blank())
}
## Example matrix
plot_example_matrix <- function(graph) {
    graph_ranked <- graph %>%
        activate(nodes) %>%
        select(Isolate, PlotRank) %>%
        activate(edges) %>%
        mutate(
            fromRank = .N()$PlotRank[match(from, .N()$Isolate)],
            toRank = .N()$PlotRank[match(to, .N()$Isolate)])
    n_nodes <- igraph::vcount(graph_ranked)
    interaction_color <- assign_interaction_color(level = "matrix")

    graph_ranked %>%
        filter(fromRank <= toRank) %>%
        bind_edges(tibble(from = 1:n_nodes, to = 1:n_nodes, fromRank = 1:n_nodes, toRank = 1:n_nodes, InteractionType = "self")) %>%
        as_tibble() %>%
        ggplot() +
        geom_tile(aes(x = factor(toRank), y = ordered(fromRank, n_nodes:1), fill = InteractionType), width = 0.9, height = 0.9) +
        scale_x_discrete(position = "top", labels = c("top\nrank", rep("", n_nodes-2), "bottom\nrank")) +
        scale_y_discrete(position = "right", labels = c("bottom\nrank", rep("", n_nodes-2), "top\nrank")) +
        scale_fill_manual(breaks = c("exclusion", "coexistence", "self"), values = c(assign_interaction_color(), "self" = "black")) +
        theme_bw() +
        theme(axis.ticks = element_blank(), axis.title = element_blank(), legend.title = element_blank(),
              axis.text = element_text(size = 10, color = 1),
              panel.border = element_blank(), panel.grid = element_blank())
}
p_example_matrix <- plot_example_matrix(net_list[[1]]) + guides(fill = "none") + theme(plot.background = element_rect(fill = "grey90"), panel.background = element_rect(fill = "grey90"))
## Get legend for matrix
shared_legend_matrix <- get_legend(plot_example_matrix(net_list[[1]]) + theme(legend.text = element_text(size = 15), legend.justification = "right", plot.background = element_rect(fill = "grey90"), panel.background = element_rect(fill = "grey90")))
## Get legend for line
p_temp <- plot_motif_count(1) + theme(legend.position = "right", legend.title = element_blank(), legend.text = element_text(size = 15))
shared_legend_line <- cowplot::get_legend(p_temp)
p1 <- plot_grid(plotlist = c(list(p_example_matrix), list(shared_legend_matrix), list(NULL), list(shared_legend_line), rep(list(NULL), 1)), nrow = 1) + theme(plot.background = element_rect(fill = "white", color = NA))
p2 <- list(p_list[1:5], p_motif_count_list[1:5],
           p_list[6:10], p_motif_count_list[6:10],
           p_list[11:13], rep(list(NULL), 2),
           p_motif_count_list[11:13]) %>%
    unlist(recursive = F) %>%
    plot_grid(plotlist = ., labels = c(communities$Community[1:5], rep("", 5), communities$Community[6:10], rep("", 5), communities$Community[11:13], rep("", 7)),
              ncol = 5, axis = "tbrl", align = "v")
p <- plot_grid(p1, p2, ncol = 1, rel_heights = c(1,6)) + theme(plot.background = element_rect(fill = "white", color = NA))
ggsave(here::here("plots/FigS7-networks_matrix.png"), p, width = 10, height = 12)


# Figure S8. Diagonal analysis ----------------------------------------------------------------------------------------------------
networks_diag <- read_csv(here::here("data/output/networks_diag.csv"), col_types = cols())
networks_diag_randomized <- read_csv(here::here("data/output/networks_diag_randomized.csv"), col_types = cols())

## Each community
networks_diag_summary <- networks_diag %>%
    mutate(Community = factor(Community, community_factor)) %>%
    group_by(Community, RankDifference) %>%
    summarize(ObservedCount = sum(Count), ObservedTotalCount = sum(TotalCount)) %>%
    mutate(ObservedFraction = ObservedCount/ObservedTotalCount)
networks_diag_randomized_summary <- networks_diag_randomized %>%
    mutate(Community = factor(Community, community_factor)) %>%
    group_by(Community, Replicate, RankDifference) %>%
    summarize(Count = sum(Count), TotalCount = sum(TotalCount)) %>%
    left_join(select(networks_diag_sum, RankDifference)) %>%
    mutate(Fraction = Count/TotalCount)
## Statistics
temp <- tibble(Replicate = 0, RankDifference = 0:11, Count = 0) %>%
    slice(rep(1:12, 17)) %>%
    mutate(Community = rep(communities$Community, each = 12)) %>%
    mutate(Community = factor(Community, community_factor))
stat_diag <- networks_diag_randomized_summary %>%
    mutate(Community = factor(Community, community_factor)) %>%
    group_by(Community, RankDifference) %>%
    # find percentile
    #    bind_rows(temp) %>%
    left_join(networks_diag_summary) %>%
    arrange(RankDifference, desc(Count)) %>%
    mutate(Percentile = (1:n())/n()) %>%
    filter(Count < ObservedCount) %>%
    group_by(Community, RankDifference) %>%
    slice(1) %>%
    select(RankDifference, Percentile) %>%
    # Asterisk
    mutate(Significance = case_when(Percentile < 0.001 | Percentile > 0.999 ~ "***",
                                    Percentile < 0.01 | Percentile > 0.99 ~ "**",
                                    Percentile < 0.05 | Percentile > 0.95 ~ "*",
                                    Percentile > 0.05 & Percentile < 0.95 ~ "n.s."),
           Sign = case_when(Percentile < 0.05 ~ "top",
                            Percentile > 0.95 ~ "bottom",
                            Percentile > 0.05 & Percentile < 0.95 ~ "n.s."))

p <- networks_diag_summary %>%
    ggplot() +
    # Asterisk
    geom_text(data = stat_diag, aes(x = RankDifference, y = Inf, label = Significance), vjust = 2) +
    geom_rect(data = stat_diag, aes(xmin = RankDifference-0.5, xmax = RankDifference+0.5, fill = Sign), ymin = -Inf, ymax = Inf, alpha = .2) +
    # Random networks
    geom_boxplot(data = networks_diag_randomized_summary, aes(x = RankDifference, y = Fraction, group = RankDifference, color = "randomized network"),
                 outlier.size = 1) +
    geom_jitter(data = networks_diag_randomized_summary, aes(x = RankDifference, y = Fraction, group = RankDifference, color = "randomized network"),
                size = .1, alpha = 0.5, width = .3) +
    # Observed networks
    geom_point(aes(x = RankDifference, y = ObservedFraction, group = RankDifference, color = "observed network"), size = 2) +
    geom_line(aes(x = RankDifference, y = ObservedFraction, color = "observed network")) +
    scale_x_continuous(breaks = 0:11) +
    scale_y_continuous(limits = c(0, 1.1), breaks = c(0, .5, 1)) +
    scale_color_manual(values = c("observed network" = "red", "randomized network" = "black")) +
    scale_fill_manual(values = c("top" = "blue", "bottom" = "red", "n.s." = "grey")) +
    facet_wrap(Community~., ncol = 4) +
    theme_classic() +
    theme(legend.position = "top", legend.title = element_blank(),
          panel.border = element_rect(fill = NA, color = 1, size = 1.5),
          strip.background = element_blank()) +
    guides(fill = "none", color = "none") +
    labs(x = "|i-j|", y = "Fraction of pairwise coexitence")
ggsave(here::here("plots/FigS8-diagonal_analysis.png"), p, width = 15, height = 12)


# Figure S9 Cliques/Components in experiments and simulations ----------------------------------------------------------------------------------------------------
# Data
## Experiments
networks_component <- read_csv(here::here("data/output/networks_component.csv")) %>% left_join(select(communities, Assembly, Community))
networks_component_randomized <- read_csv(here::here("data/output/networks_component_randomized.csv")) %>% left_join(select(communities, Assembly, Community))
## Simulation
df_component <- read_csv(here::here("data/output/df_component.csv"))
df_component_randomized <- read_csv(here::here("data/output/df_component_randomized.csv"))

# Binding the data
## Observation
networks_component_obv <- bind_rows(
    mutate(networks_component, Treatment = "experiment"),
    mutate(df_component, Treatment = "simulation")
) %>%
    mutate(Assembly = factor(Assembly, c("self_assembly", "random_assembly"))) %>%
    group_by(Treatment, Assembly) %>%
    summarize(Component = sum(Component))
## Permutation
networks_component_pmt <- bind_rows(
    mutate(networks_component_randomized, Treatment = "experiment"),
    mutate(df_component_randomized, Treatment = "simulation")
) %>%
    mutate(Assembly = factor(Assembly, c("self_assembly", "random_assembly"))) %>%
    group_by(Treatment, Assembly, Replicate) %>%
    summarize(Component = sum(Component))

# Plot
p <- networks_component_pmt %>%
    ggplot() +
    geom_histogram(aes(x = Component, color = "permutation"), binwidth = 1, color = 1, fill = NA) +
    geom_vline(data = networks_component_obv, aes(color = "observation", xintercept = Component)) +
    scale_color_manual(values = c("observation" = "red", "permutation" = "black")) +
    facet_grid(Treatment~Assembly, scales = "free") +
    theme_classic() +
    theme(panel.border = element_rect(color = 1, fill = NA),
          legend.position = "right", legend.title = element_blank()) +
    #guides(color = "none") +
    labs(x = "Number of cliques", y = "Permutation")

ggsave(here::here("plots/FigS9-components.png"), p, width = 6, height = 4)


# Figures not used yet. ----
# Figure 01S0. Summary stat for random assembly and self-assembly in exp and sim ----------------------------------------------------------------------------------------------------
if (FALSE) {
    ## Treat across-community as random assembly
    p1 <- pairs %>%
        mutate(InteractionType = factor(InteractionType, c("exclusion", "coexistence"))) %>%
        mutate(Assembly = ifelse(Assembly == "across_community", "random_assembly", Assembly)) %>%
        mutate(Assembly = factor(Assembly, c("self_assembly", "random_assembly"))) %>%
        group_by(Assembly, InteractionType) %>% summarize(Count = n()) %>%
        mutate(Fraction = Count / sum(Count), TotalCount = sum(Count)) %>%
        ggplot() +
        geom_col(aes(x = Assembly, y = Count, fill = InteractionType), position = "fill", color = 1, width = 0.8) +
        geom_text(aes(label = paste0("n=", TotalCount), x = Assembly), y = Inf, vjust = 1.5) +
        scale_fill_manual(values = assign_interaction_color(level = "simple")) +
        scale_y_continuous(expand = c(0,0), breaks = c(0, .5, 1)) +
        theme_classic() +
        theme(axis.title.x = element_blank(), legend.position = "top",
              axis.text.x = element_text(size = 10, color = "black", angle = 15, vjust = 1, hjust = 1)) +
        guides(fill = "none") +
        labs(x = "", y  = "Fraction", fill = "")

    ## Treat across-community and random assembly separately
    p2 <- pairs %>%
        mutate(InteractionType = factor(InteractionType, c("exclusion", "coexistence"))) %>%
        #mutate(Assembly = ifelse(Assembly == "across_community", "random_assembly", Assembly)) %>%
        mutate(Assembly = factor(Assembly, c("self_assembly", "across_community", "random_assembly"))) %>%
        group_by(Assembly, InteractionType) %>% summarize(Count = n()) %>%
        mutate(Fraction = Count / sum(Count), TotalCount = sum(Count)) %>%
        ggplot() +
        geom_col(aes(x = Assembly, y = Count, fill = InteractionType), position = "fill", color = 1, width = 0.8) +
        geom_text(aes(label = paste0("n=", TotalCount), x = Assembly), y = Inf, vjust = 1.5) +
        scale_fill_manual(values = assign_interaction_color(level = "simple")) +
        scale_y_continuous(expand = c(0,0), breaks = c(0, .5, 1)) +
        theme_classic() +
        theme(axis.title.x = element_blank(), legend.position = "top",
              axis.text.x = element_text(size = 10, color = "black", angle = 15, vjust = 1, hjust = 1)) +
        guides(fill = "none") +
        labs(x = "", y  = "Fraction", fill = "")

}
# Pairwise coexistence
pairs <- read_csv(here::here("data/output/pairs.csv"), col_types = cols()) %>% mutate(InteractionType = ifelse(InteractionType == "neutrality", "coexistence", InteractionType))
df_pairs <- read_csv(here::here("data/output/df_pairs.csv"))

pairs_merged <- pairs %>%
    select(Assembly, InteractionType) %>%
    mutate(Assembly = ifelse(Assembly == "across_community", "random_assembly", Assembly)) %>%
    mutate(Assembly = factor(Assembly, c("self_assembly", "random_assembly"))) %>%
    mutate(Treatment = "experiment") %>%
    bind_rows(df_pairs %>% select(Assembly, InteractionType) %>% mutate(Treatment = "simulation")) %>%
    mutate(Assembly = factor(Assembly, c("self_assembly", "random_assembly")))
p1 <- pairs_merged %>%
    group_by(Treatment, Assembly, InteractionType) %>%
    summarize(Count = n()) %>%
    mutate(Fraction = Count / sum(Count), TotalCount = sum(Count)) %>%
    ggplot() +
    geom_col(aes(x = Assembly, y = Count, fill = InteractionType), position = "fill", color = 1, width = 0.8) +
    geom_text(aes(label = paste0("n=", TotalCount), x = Assembly), y = Inf, vjust = 2) +
    scale_fill_manual(values = assign_interaction_color(level = "simple")) +
    scale_y_continuous(expand = c(0,0), breaks = c(0, .5, 1)) +
    facet_wrap(Treatment~Assembly, nrow = 1, scales = "free_x") +
    #facet_grid(.~Treatment) +
    theme_classic() +
    theme(axis.title.x = element_blank(), legend.position = "right",
          axis.text.x = element_blank(),
          #axis.text.x = element_text(size = 10, color = "black", angle = 15, vjust = 1, hjust = 1),
          panel.border = element_rect(color = 1, fill = NA, size = 1)) +
    labs(x = "", y  = "Fraction", fill = "") +
    ggtitle("Pairwise outcomes")

## Stat
### Experiment
pairs_merged %>%
    mutate(InteractionType = factor(InteractionType, c("exclusion", "coexistence"))) %>%
    filter(Treatment == "experiment") %>%
    chisq_test(InteractionType ~ Assembly)
### Simulation
pairs_merged %>%
    mutate(InteractionType = factor(InteractionType, c("exclusion", "coexistence"))) %>%
    filter(Treatment == "simulation") %>%
    chisq_test(InteractionType ~ Assembly)

# Panel B. Hierarchy
# Data
## Experiments
communities_hierarchy <- read_csv(here::here("data/output/communities_hierarchy.csv")) %>% left_join(select(communities, Assembly, Community))
communities_hierarchy_randomized <- read_csv(here::here("data/output/communities_hierarchy_randomized.csv")) %>% left_join(select(communities, Assembly, Community))
## Simulation
df_communities_hierarchy <- read_csv(here::here("data/output/df_communities_hierarchy.csv"))
df_communities_hierarchy_randomized <- read_csv(here::here("data/output/df_communities_hierarchy_randomized.csv"))

# Binding the data
## Observation
communities_hierarchy_obv <- bind_rows(
    mutate(communities_hierarchy, Treatment = "experiment"),
    mutate(df_communities_hierarchy, Treatment = "simulation")
) %>%
    mutate(Assembly = factor(Assembly, c("self_assembly", "random_assembly")))
## Permutation
communities_hierarchy_pmt <- bind_rows(
    mutate(communities_hierarchy_randomized, Treatment = "experiment"),
    mutate(df_communities_hierarchy_randomized, Treatment = "simulation")
) %>%
    mutate(Assembly = factor(Assembly, c("self_assembly", "random_assembly")))

# Plot
p2 <- communities_hierarchy_obv %>%
    ggplot(aes(x = Metric, y = HierarchyScore)) +
    geom_boxplot(width = .5, lwd = .8) +
    geom_jitter(shape = 1, size = 2, width = .2, stroke = .8) +
    scale_y_continuous(limits = c(0,1.01), breaks = c(0, .25, .5, .75, 1)) +
    facet_wrap(Treatment~Assembly, scales = "free_x", nrow = 1) +
    theme_classic() +
    theme(panel.grid.major.y = element_line(color = "grey", linetype = 2),
          axis.text = element_text(size = 10, color = 1),
          axis.title = element_text(size = 15, color = 1),
          panel.border = element_rect(fill = NA, color = 1, size = 1.5)) +
    labs(x = "Hierarchy", y = "Score") +
    ggtitle("Hierarchy")


# Panel C. Clique
# Data
## Experiments
networks_component <- read_csv(here::here("data/output/networks_component.csv")) %>% left_join(select(communities, Assembly, Community))
networks_component_randomized <- read_csv(here::here("data/output/networks_component_randomized.csv")) %>% left_join(select(communities, Assembly, Community))
## Simulation
df_component <- read_csv(here::here("data/output/df_component.csv"))
df_component_randomized <- read_csv(here::here("data/output/df_component_randomized.csv"))

# Binding the data
## Observation
networks_component_obv <- bind_rows(
    mutate(networks_component, Treatment = "experiment"),
    mutate(df_component, Treatment = "simulation")
) %>%
    mutate(Assembly = factor(Assembly, c("self_assembly", "random_assembly"))) %>%
    group_by(Treatment, Assembly) %>%
    summarize(Component = sum(Component))
## Permutation
networks_component_pmt <- bind_rows(
    mutate(networks_component_randomized, Treatment = "experiment"),
    mutate(df_component_randomized, Treatment = "simulation")
) %>%
    mutate(Assembly = factor(Assembly, c("self_assembly", "random_assembly"))) %>%
    group_by(Treatment, Assembly, Replicate) %>%
    summarize(Component = sum(Component))

# Plot
p3 <- networks_component_pmt %>%
    ggplot() +
    geom_histogram(aes(x = Component, color = "permutation"), binwidth = 1, color = 1, fill = NA) +
    geom_vline(data = networks_component_obv, aes(color = "observation", xintercept = Component)) +
    scale_color_manual(values = c("observation" = "red", "permutation" = "black")) +
    facet_wrap(Treatment~Assembly, nrow = 1, scales = "free_x") +
    theme_classic() +
    theme(panel.border = element_rect(color = 1, fill = NA),
          legend.position = "right", legend.title = element_blank()) +
    labs(x = "Number of cliques", y = "Permutation") +
    ggtitle("Clique")

p <- plot_grid(p1, p2, p3, ncol = 1, labels = c("A", "B", "C"), axis = "lr", align = "v") + paint_white_background()
ggsave(here::here("plots/01S0-pairwsie_competition_assembly.png"), p, width = 10, height = 10)








# Figure 01S1. Analysis of self-assembly ----------------------------------------------------------------------------------------------------
# Figure A. Matrices
p_net_matrix_list <- lapply(net_list, function(x) plot_adjacent_matrix(x) + theme(plot.margin = grid::unit(c(5,0,3,0), "mm")))
## Subset the random-assembly networks
p_net_matrix_list <- p_net_matrix_list %>% `[`(communities %>% arrange(CommunitySize) %>% filter(str_detect(Community, "C\\d")) %>% pull(Community))
pA <- plot_grid(plotlist = p_net_matrix_list[1:13], scale = .9, nrow = 2, labels = names(p_net_matrix_list)) + theme(plot.background = element_rect(fill = "grey90", color = NA)) + paint_white_background()

# Figure B. pairwise outcome
pB <- pairs %>%
    filter(Assembly == "self_assembly") %>%
    group_by(InteractionType) %>%
    mutate(InteractionType = factor(InteractionType, c("exclusion", "coexistence"))) %>%
    summarize(Count = n()) %>%
    mutate(TotalCount = sum(Count), Fraction = Count/TotalCount) %>%
    ggplot() +
    geom_col(aes(x = InteractionType, y = Count, fill = InteractionType), color = 1, position = position_dodge(width = .8), width = .8) +
    geom_text(aes(x = InteractionType, y = Count, label = paste0(round(Fraction, 3) * 100,"%")), nudge_y = 5, size = 3) +
    geom_text(x = Inf, y = Inf, aes(label = paste0("n=", TotalCount)), vjust = 1, hjust = 1, size = 3) +
    #geom_text(aes(x = Assembly, y = Count, label = paste0(round(Fraction, 2)* 100, "%")), size = 5, position = position_dodge(width = .8)) +
    scale_fill_manual(values = assign_interaction_color(level = "simple")) +
    scale_y_continuous(expand = c(0,0), limits = c(0, 140), breaks = scales::pretty_breaks(n = 3)) +
    theme_classic() +
    theme(axis.title.x = element_blank(), axis.title.y = element_text(size = 10),
          axis.text.x = element_text(size = 10, color = "black", angle = 15, vjust = 1, hjust = 1),
          axis.text.y = element_text(size = 10, color = "black")) +
    labs(x = "", y = "Number of pairs", fill = "")

# Figure C: Hierarchy
communities_hierarchy <- read_csv(here::here("data/output/communities_hierarchy.csv"), col_types = cols())
pC <- communities_hierarchy %>%
    filter(str_detect(Community, "C\\d")) %>%
    ggplot(aes(x = Metric, y = HierarchyScore)) +
    geom_boxplot(width = .5, lwd = .8) +
    geom_jitter(shape = 1, size = 2, width = .2, stroke = .8) +
    scale_x_discrete (position = "bottom", labels = c("1" = "h2", "2" = "h1", "3" = "Fraction of transitive motifs")) +
    scale_y_continuous(limits = c(0,1.01), breaks = c(0, .25, .5, .75, 1)) +
    #facet_wrap(.~Metric, scales = "free", nrow = 1) +
    theme_classic() +
    theme(panel.grid.major.y = element_line(color = "grey", linetype = 2),
          axis.text = element_text(size = 10, color = 1),
          axis.title = element_text(size = 15, color = 1),
          panel.border = element_rect(fill = NA, color = 1, size = 1.5)) +
    labs(x = "Hierarchy", y = "Score")

# Figure D: motif of all networks
## motif diagram
load(here::here("data/output/motif_list.Rdata"))
p_motif_list <- lapply(motif_list, function(x) plot_competitive_network(x, node_size = 2, edge_width = 1))
p1 <- plot_grid(plotlist = p_motif_list, nrow = 1, scale = 1.3)

## Motif count
networks_motif <- read_csv(here::here("data/output/networks_motif.csv"))
networks_motif_randomized <- read_csv(here::here("data/output/networks_motif_randomized.csv"))

## Permutation
networks_motif_randomized_summary <- networks_motif_randomized %>%
    filter(str_detect(Community, "C\\d")) %>%
    select(Community, Replicate, Motif, Count) %>%
    group_by(Replicate, Motif) %>%
    summarize(Count = sum(Count)) %>%
    group_by(Replicate) %>%
    mutate(TotalCount = sum(Count)) %>%
    mutate(Fraction = Count/TotalCount)
networks_motif_randomized_summary_percentile <- networks_motif_randomized_summary %>%
    group_by(Motif) %>%
    arrange(Motif, Count) %>%
    slice(c(1000 * 0.05, 1000 * 0.95)) %>%
    mutate(Percentile = c("p5", "p95")) %>%
    select(Motif, Count, Percentile)
## Observation
networks_motif_summary <- networks_motif %>%
    filter(str_detect(Community, "C\\d")) %>%
    group_by(Motif) %>%
    summarize(Count = sum(Count)) %>%
    left_join(pivot_wider(networks_motif_randomized_summary_percentile, names_from = Percentile, values_from = Count)) %>%
    mutate(Sign = case_when(Count > p95 ~ "top",
                            Count < p5 ~ "bottom",
                            Count < p95 & Count > p5 ~ "n.s.")) %>%
    mutate(TotalCount = sum(Count), Fraction = Count/TotalCount)

p2 <- networks_motif_randomized_summary %>%
    group_by(Motif, Replicate) %>%
    mutate(p5 = quantile(Count, 0.05), p95 = quantile(Count, 0.95)) %>%
    ggplot() +
    geom_rect(data = networks_motif_summary, aes(fill = Sign), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, alpha = .2) +
    geom_vline(xintercept = 0, color = 1) +
    geom_hline(yintercept = 0, color = 1, linetype = 2) +
    geom_histogram(aes(y = Fraction, x = after_stat(count / max(count))), alpha = .3, color = 1, bins = 50) +
    geom_point(data = networks_motif_summary, aes(x = 0, y = Fraction, color = "observed network"), pch = 1, size = 2, stroke = 2, inherit.aes = F) +
    scale_fill_manual(values = c("top" = "blue", "bottom" = "red", "n.s." = "grey")) +
    scale_color_manual(values = c("observed network" = "red")) +
    facet_grid(.~Motif) +
    scale_x_continuous(breaks = c(0,0.5), labels = c("0", "0.5")) +
    theme_cowplot() +
    theme(panel.background = element_rect(color = 1, size = 1.5, fill = NA),
          panel.spacing = unit(0, "mm"), strip.background = element_rect(color = NA, fill = NA),
          legend.position = "top",
          axis.title = element_text(size = 10), axis.text = element_text(size = 8),
          plot.background = element_rect(fill = "white", color = NA)) +
    guides(fill = "none", color = "none") +
    labs(x = "Probability density", y = "Fraction of motif")

pD <- plot_grid(p1, p2, ncol = 1, rel_heights = c(1, 5), axis = "lr", align = "v") + paint_white_background()


# Figure E: diagonal analysis. One community for example
networks_diag <- read_csv(here::here("data/output/networks_diag.csv"), col_types = cols())
networks_diag_randomized <- read_csv(here::here("data/output/networks_diag_randomized.csv"), col_types = cols())
## Permutation
networks_diag_randomized_sum <- networks_diag_randomized %>%
    filter(str_detect(Community, "C\\d")) %>%
    group_by(Replicate, RankDifference) %>%
    summarize(Count = sum(Count), TotalCount = sum(TotalCount)) %>%
    mutate(Fraction = Count/TotalCount)

## Observation
networks_diag_sum <- networks_diag %>%
    filter(str_detect(Community, "C\\d")) %>%
    group_by(RankDifference) %>%
    summarize(ObservedCount = sum(Count), ObservedTotalCount = sum(TotalCount)) %>%
    mutate(ObservedFraction = ObservedCount/ObservedTotalCount)

## Statistics
stat_diag <- networks_diag_randomized_sum %>%
    group_by(RankDifference) %>%
    # Find percentile
    bind_rows(tibble(Replicate = 0, RankDifference = 0:11, Count = 0)) %>%
    left_join(networks_diag_sum) %>%
    arrange(RankDifference, desc(Fraction)) %>%
    mutate(Percentile = (1:n())/n()) %>%
    filter(Count <= ObservedCount) %>%
    group_by(RankDifference) %>%
    slice(1) %>%
    select(RankDifference, Percentile) %>%
    # Asterisk
    mutate(Significance = case_when(Percentile < 0.001 | Percentile > 0.999 ~ "***",
                                    Percentile < 0.01 | Percentile > 0.99 ~ "**",
                                    Percentile < 0.05 | Percentile > 0.95 ~ "*",
                                    Percentile > 0.05 & Percentile < 0.95 ~ "n.s."),
           Sign = case_when(Percentile < 0.05 ~ "top",
                            Percentile > 0.95 ~ "bottom",
                            Percentile > 0.05 & Percentile < 0.95 ~ "n.s."))

pE <- networks_diag_sum %>%
    ggplot() +
    # Asterisk
    geom_text(data = stat_diag, aes(x = RankDifference, y = Inf, label = Significance), vjust = 2) +
    geom_rect(data = stat_diag, aes(xmin = RankDifference-0.5, xmax = RankDifference+0.5, fill = Sign), ymin = -Inf, ymax = Inf, alpha = .2) +
    # Random networks
    geom_boxplot(data = networks_diag_randomized_sum, aes(x = RankDifference, y = Fraction, group = RankDifference, color = "permutation"),
                 outlier.size = 1) +
    geom_jitter(data = networks_diag_randomized_sum, aes(x = RankDifference, y = Fraction, group = RankDifference, color = "permutation"),
                size = .2, alpha = 0.5, width = .3, shape = 21, height = .01) +
    # Observed networks
    geom_point(aes(x = RankDifference, y = ObservedFraction, group = RankDifference, color = "observation"), size = 2) +
    geom_line(aes(x = RankDifference, y = ObservedFraction, color = "observation")) +
    scale_x_continuous(breaks = 0:11, expand = c(0,0)) +
    scale_y_continuous(limits = c(-0.02, 1.1), breaks = c(0, .5, 1)) +
    scale_color_manual(values = c("observation" = "red", "permutation" = "black"))+
    scale_fill_manual(values = c("top" = "blue", "bottom" = "red", "n.s." = "grey")) +
    theme_classic() +
    theme(legend.position = "top", legend.title = element_blank(),
          legend.text = element_text(size = 12),
          axis.title.x = element_text(size = 15, color = 1),
          panel.border = element_rect(fill = NA, color = 1, size = 1.5)) +
    guides(fill = "none") +
    labs(x = "|i-j|", y = "Fraction of pairwise coexitence")

#
p_top <- plot_grid(pA, pB, nrow = 1, rel_widths = c(2,1), scale = c(.9, .9), labels = c("A", "B"))
p_bottom <- plot_grid(pC, pD, pE, nrow = 1, labels = LETTERS[3:5], rel_widths = c(1,2,1.5), scale = c(.8, .9, .9), axis = "t", align = "h")
p <- plot_grid(p_top, p_bottom, nrow = 2, scale = c(1, 1), rel_heights = c(1.5,2)) + paint_white_background()
title <- ggdraw() + draw_label("Experiment: self-assembly networks", fontface='bold', x = 0, hjust = -0.1)
p <- plot_grid(title, p, ncol = 1, rel_heights = c(0.1, 1)) + paint_white_background()
ggsave(here::here("plots/01S1-self_assembly_networks.png"), p, width = 12, height = 8)





# Figure 01S2. Analysis of random assembly ----------------------------------------------------------------------------------------------------
# Figure A. Matrices
p_net_matrix_list <- lapply(net_list, function(x) plot_adjacent_matrix(x) + theme(plot.margin = grid::unit(c(5,0,3,0), "mm")))
## Subset the random-assembly networks
p_net_matrix_list <- p_net_matrix_list %>% `[`(communities %>% arrange(CommunitySize) %>% filter(str_detect(Community, "Ass")) %>% pull(Community))
pA <- plot_grid(plotlist = p_net_matrix_list[1:4], scale = .9, nrow = 1, labels = names(p_net_matrix_list)) + theme(plot.background = element_rect(fill = "grey90", color = NA)) + paint_white_background()

# Figure B. pairwise outcome
pB <- pairs %>%
    filter(Assembly != "self_assembly") %>%
    group_by(InteractionType) %>%
    mutate(InteractionType = factor(InteractionType, c("exclusion", "coexistence"))) %>%
    summarize(Count = n()) %>%
    mutate(TotalCount = sum(Count), Fraction = Count/TotalCount) %>%
    ggplot() +
    geom_col(aes(x = InteractionType, y = Count, fill = InteractionType), color = 1, position = position_dodge(width = .8), width = .8) +
    geom_text(aes(x = InteractionType, y = Count, label = paste0(round(Fraction, 3) * 100,"%")), nudge_y = 5, size = 3) +
    geom_text(x = Inf, y = Inf, aes(label = paste0("n=", TotalCount)), vjust = 1, hjust = 1, size = 3) +
    #geom_text(aes(x = Assembly, y = Count, label = paste0(round(Fraction, 2)* 100, "%")), size = 5, position = position_dodge(width = .8)) +
    scale_fill_manual(values = assign_interaction_color(level = "simple")) +
    scale_y_continuous(expand = c(0,0), limits = c(0, 75), breaks = scales::pretty_breaks(n = 3)) +
    theme_classic() +
    theme(axis.title.x = element_blank(), axis.title.y = element_text(size = 10),
          axis.text.x = element_text(size = 10, color = "black", angle = 15, vjust = 1, hjust = 1),
          axis.text.y = element_text(size = 10, color = "black")) +
    labs(x = "", y = "Number of pairs", fill = "")


# Figure C: Hierarchy
communities_hierarchy <- read_csv(here::here("data/output/communities_hierarchy.csv"), col_types = cols())
pC <- communities_hierarchy %>%
    filter(str_detect(Community, "Ass")) %>%
    ggplot(aes(x = Metric, y = HierarchyScore)) +
    geom_boxplot(width = .5, lwd = .8) +
    geom_jitter(shape = 1, size = 2, width = .2, stroke = .8) +
    scale_x_discrete (position = "bottom", labels = c("1" = "h2", "2" = "h1", "3" = "Fraction of transitive motifs")) +
    scale_y_continuous(limits = c(0,1.01), breaks = c(0, .25, .5, .75, 1)) +
    #facet_wrap(.~Metric, scales = "free", nrow = 1) +
    theme_classic() +
    theme(panel.grid.major.y = element_line(color = "grey", linetype = 2),
          axis.text = element_text(size = 10, color = 1),
          axis.title = element_text(size = 15, color = 1),
          panel.border = element_rect(fill = NA, color = 1, size = 1.5)) +
    labs(x = "Hierarchy", y = "Score")

# Figure D: motif of all networks
## motif diagram
load(here::here("data/output/motif_list.Rdata"))
p_motif_list <- lapply(motif_list, function(x) plot_competitive_network(x, node_size = 2, edge_width = 1))
p1 <- plot_grid(plotlist = p_motif_list, nrow = 1, scale = 1.3)

## Motif count
networks_motif <- read_csv(here::here("data/output/networks_motif.csv"))
networks_motif_randomized <- read_csv(here::here("data/output/networks_motif_randomized.csv"))

## Permutation
networks_motif_randomized_summary <- networks_motif_randomized %>%
    filter(str_detect(Community, "Ass")) %>%
    select(Community, Replicate, Motif, Count) %>%
    group_by(Replicate, Motif) %>%
    summarize(Count = sum(Count)) %>%
    group_by(Replicate) %>%
    mutate(TotalCount = sum(Count)) %>%
    mutate(Fraction = Count/TotalCount)
networks_motif_randomized_summary_percentile <- networks_motif_randomized_summary %>%
    group_by(Motif) %>%
    arrange(Motif, Count) %>%
    slice(c(1000 * 0.05, 1000 * 0.95)) %>%
    mutate(Percentile = c("p5", "p95")) %>%
    select(Motif, Count, Percentile)
## Observation
networks_motif_summary <- networks_motif %>%
    filter(str_detect(Community, "Ass")) %>%
    group_by(Motif) %>%
    summarize(Count = sum(Count)) %>%
    left_join(pivot_wider(networks_motif_randomized_summary_percentile, names_from = Percentile, values_from = Count)) %>%
    mutate(Sign = case_when(Count > p95 ~ "top",
                            Count < p5 ~ "bottom",
                            Count < p95 & Count > p5 ~ "n.s.")) %>%
    mutate(TotalCount = sum(Count), Fraction = Count/TotalCount)

p2 <- networks_motif_randomized_summary %>%
    group_by(Motif, Replicate) %>%
    mutate(p5 = quantile(Count, 0.05), p95 = quantile(Count, 0.95)) %>%
    ggplot() +
    geom_rect(data = networks_motif_summary, aes(fill = Sign), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, alpha = .2) +
    geom_vline(xintercept = 0, color = 1) +
    geom_hline(yintercept = 0, color = 1, linetype = 2) +
    geom_histogram(aes(y = Fraction, x = after_stat(count / max(count))), alpha = .3, color = 1, bins = 50) +
    geom_point(data = networks_motif_summary, aes(x = 0, y = Fraction, color = "observed network"), pch = 1, size = 2, stroke = 2, inherit.aes = F) +
    scale_fill_manual(values = c("top" = "blue", "bottom" = "red", "n.s." = "grey")) +
    scale_color_manual(values = c("observed network" = "red")) +
    facet_grid(.~Motif) +
    scale_x_continuous(breaks = c(0,0.5), labels = c("0", "0.5")) +
    theme_cowplot() +
    theme(panel.background = element_rect(color = 1, size = 1.5, fill = NA),
          panel.spacing = unit(0, "mm"), strip.background = element_rect(fill = NA, color = NA),
          legend.position = "top",
          axis.title = element_text(size = 10), axis.text = element_text(size = 8),
          plot.background = element_rect(fill = "white", color = NA)) +
    guides(fill = "none", color = "none") +
    labs(x = "Probability density", y = "Fraction of motif")

pD <- plot_grid(p1, p2, ncol = 1, rel_heights = c(1, 5), axis = "lr", align = "v") + paint_white_background()



# Figure E: diagonal analysis. One community for example
networks_diag <- read_csv(here::here("data/output/networks_diag.csv"), col_types = cols())
networks_diag_randomized <- read_csv(here::here("data/output/networks_diag_randomized.csv"), col_types = cols())
## Permutation
networks_diag_randomized_sum <- networks_diag_randomized %>%
    filter(str_detect(Community, "Ass")) %>%
    group_by(Replicate, RankDifference) %>%
    summarize(Count = sum(Count), TotalCount = sum(TotalCount)) %>%
    mutate(Fraction = Count/TotalCount)

## Observation
networks_diag_sum <- networks_diag %>%
    filter(str_detect(Community, "Ass")) %>%
    group_by(RankDifference) %>%
    summarize(ObservedCount = sum(Count), ObservedTotalCount = sum(TotalCount)) %>%
    mutate(ObservedFraction = ObservedCount/ObservedTotalCount)

## Statistics
stat_diag <- networks_diag_randomized_sum %>%
    group_by(RankDifference) %>%
    # Find percentile
    bind_rows(tibble(Replicate = 0, RankDifference = 0:11, Count = 0)) %>%
    left_join(networks_diag_sum) %>%
    arrange(RankDifference, desc(Fraction)) %>%
    mutate(Percentile = (1:n())/n()) %>%
    filter(Count <= ObservedCount) %>%
    group_by(RankDifference) %>%
    slice(1) %>%
    select(RankDifference, Percentile) %>%
    # Asterisk
    mutate(Significance = case_when(Percentile < 0.001 | Percentile > 0.999 ~ "***",
                                    Percentile < 0.01 | Percentile > 0.99 ~ "**",
                                    Percentile < 0.05 | Percentile > 0.95 ~ "*",
                                    Percentile > 0.05 & Percentile < 0.95 ~ "n.s."),
           Sign = case_when(Percentile < 0.05 ~ "top",
                            Percentile > 0.95 ~ "bottom",
                            Percentile > 0.05 & Percentile < 0.95 ~ "n.s."))

pE <- networks_diag_sum %>%
    ggplot() +
    # Asterisk
    geom_text(data = stat_diag, aes(x = RankDifference, y = Inf, label = Significance), vjust = 2) +
    geom_rect(data = stat_diag, aes(xmin = RankDifference-0.5, xmax = RankDifference+0.5, fill = Sign), ymin = -Inf, ymax = Inf, alpha = .2) +
    # Random networks
    geom_boxplot(data = networks_diag_randomized_sum, aes(x = RankDifference, y = Fraction, group = RankDifference, color = "permutation"),
                 outlier.size = 1) +
    geom_jitter(data = networks_diag_randomized_sum, aes(x = RankDifference, y = Fraction, group = RankDifference, color = "permutation"),
                size = .2, alpha = 0.5, width = .3, shape = 21, height = .01) +
    # Observed networks
    geom_point(aes(x = RankDifference, y = ObservedFraction, group = RankDifference, color = "observation"), size = 2) +
    geom_line(aes(x = RankDifference, y = ObservedFraction, color = "observation")) +
    scale_x_continuous(breaks = 0:11, expand = c(0,0)) +
    scale_y_continuous(limits = c(-0.02, 1.1), breaks = c(0, .5, 1)) +
    scale_color_manual(values = c("observation" = "red", "permutation" = "black"))+
    scale_fill_manual(values = c("top" = "blue", "bottom" = "red", "n.s." = "grey")) +
    theme_classic() +
    theme(legend.position = "top", legend.title = element_blank(),
          legend.text = element_text(size = 12),
          axis.title.x = element_text(size = 15, color = 1),
          panel.border = element_rect(fill = NA, color = 1, size = 1.5)) +
    guides(fill = "none") +
    labs(x = "|i-j|", y = "Fraction of pairwise coexitence")

#
p_top <- plot_grid(pA, pB, nrow = 1, rel_widths = c(2,1), scale = c(.9, .9), labels = c("A", "B"))
p_bottom <- plot_grid(pC, pD, pE, nrow = 1, labels = LETTERS[3:5], rel_widths = c(1,2,1.5), scale = c(.8, .9, .9), axis = "t", align = "h")
p <- plot_grid(p_top, p_bottom, nrow = 2, scale = c(1, 1), rel_heights = c(1.5,2))
title <- ggdraw() + draw_label("Experiment: random-assembly networks", fontface='bold', x = 0, hjust = -0.1)
p <- plot_grid(title, p, ncol = 1, rel_heights = c(0.1, 1)) + paint_white_background()
ggsave(here::here("plots/01S2-random_networks.png"), p, width = 12, height = 8)




# Table S1. Pairwise interaction tables ----
ft1 <- read_csv(here::here("data/output/pairs_interaction_table.csv")) %>%
    mutate(InteractionType = ifelse(InteractionType == "neutrality", "coexistence", InteractionType)) %>%
    arrange(InteractionType, InteractionTypeFiner, FromRare, FromMedium, FromAbundant) %>%
    setNames(c("From 5%", "From 50%", "From 95%", "Outcome", "Finer outcome", "Count")) %>%
    flextable() %>%
    width(j = 1:3, width = 1) %>%
    width(j = 5, width = 2.5)
save_as_image(ft1, here::here("plots/TableS1.png"))
## Count the totoal number
#read_csv(here::here("data/output/pairs_interaction_table.csv")) %>% pull(Count) %>% sum


# deprecated ----
if (FALSE) {
    # Figure 1S1: two hypotheses; two network structures
    load(here::here("data/output/motif_list.Rdata"))
    node_size = 5; text_size = 3
    p1 <- motif_list[[7]] %>%
        ggraph(layout = "circle") +
        geom_node_point(shape = 21, size = node_size, fill = "gray", stroke = node_size/5) +
        geom_edge_link(aes(color = InteractionType), width = node_size/5,
                       arrow = arrow(length = unit(node_size/2, "mm"), type = "closed", angle = 30, ends = "last"),
                       start_cap = circle(node_size/2+1, "mm"),
                       end_cap = circle(node_size/2+1, "mm")) +
        geom_text(x = 0.25, y = 0.433, label = "coexistence", angle = -30, vjust = -1.5, fontface = "bold", size = text_size, color = interaction_color["coexistence"]) +
        geom_text(x = -1.1, y = 0, label = "All-coexistence", angle = 90, size = text_size*1.2, color = 1) +
        scale_edge_color_manual(values = interaction_color) +
        scale_color_manual(values = c("black" = "black", "white" = NA)) +
        scale_shape_manual(values = c("real" = 21, "fake" = NA)) +
        scale_x_continuous(limits = c(-1.2, 1.2)) +
        scale_y_continuous(limits = c(-1.2, 1.2)) +
        theme_minimal() +
        theme(legend.position = "none", panel.grid = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
              plot.background = element_rect(fill = NA, color = NA), plot.margin = margin(0,0,0,0, unit = "pt")) +
        labs(x = "", y = "")

    p2 <- motif_list[[1]] %>%
        ggraph(layout = "circle") +
        geom_node_point(shape = 21, size = node_size, fill = "gray", stroke = node_size/5) +
        geom_edge_link(aes(color = InteractionType), width = node_size/5,
                       arrow = arrow(length = unit(node_size/2, "mm"), type = "closed", angle = 30, ends = "last"),
                       start_cap = circle(node_size/2+1, "mm"),
                       end_cap = circle(node_size/2+1, "mm")) +
        geom_text(x = 0.25, y = 0.433, label = "exclusion", angle = -30, vjust = -1.5, fontface = "bold", size = text_size, color = interaction_color["exclusion"]) +
        geom_text(x = -1.1, y = 0, label = "Nontransitivity", angle = 90, size = text_size*1.2, color = 1) +
        scale_edge_color_manual(values = interaction_color) +
        scale_color_manual(values = c("black" = "black", "white" = NA)) +
        scale_shape_manual(values = c("real" = 21, "fake" = NA)) +
        scale_x_continuous(limits = c(-1.2, 1.2)) +
        scale_y_continuous(limits = c(-1.2, 1.2)) +
        theme_minimal() +
        theme(legend.position = "none", panel.grid = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
              plot.background = element_rect(fill = NA, color = NA), plot.margin = margin(0,0,0,0, unit = "pt")) +
        labs(x = "", y = "")

    pS1 <- plot_grid(p1, p2, ncol = 1) + theme(panel.background = element_rect(color = NA, fill = "white"))
    ggsave(here::here("plots/Fig1S1-two_hypothese.png"), pS1, width = 1.5, height = 3)



    # Figure 1S2: nontransitive motif count
    networks_motif <- read_csv(here::here("data/output/networks_motif.csv"), col_types = cols()) %>% filter(str_detect(Community, "C\\d"))
    networks_motif_randomized <- read_csv(here::here("data/output/networks_motif_randomized.csv"), col_types = cols()) %>% filter(str_detect(Community, "C\\d")) %>% mutate(Community = ordered(Community, communities$Community))

    ## Diagram of nontransitive motif
    load(here::here("data/output/motif_list.Rdata"))
    p_motif1 <- plot_competitive_network(motif_list[[1]], node_size = 3)

    ##
    networks_motif_total <- networks_motif %>%
        group_by(Motif) %>%
        summarize(Count = sum(Count)) %>%
        filter(Motif == 1)

    p_motif_count <- networks_motif_randomized %>%
        filter(Motif == 1) %>%
        group_by(Motif, Replicate) %>%
        summarize(Count = sum(Count)) %>%
        group_by(Motif) %>%
        # 5% and 95% percentiles in randomized networks
        mutate(p5 = quantile(Count, 0.05), p95 = quantile(Count, 0.95), ColoredTails = ifelse(Count <= p5, "tail", ifelse(Count >= p95, "head", "body"))) %>%
        ggplot() +
        geom_vline(xintercept = 0, color = 1) +
        geom_hline(yintercept = 0, color = 1) +
        geom_histogram(aes(y = Count, x = after_stat(count / max(count)), fill = ColoredTails), alpha = .3, color = 1) +
        geom_point(data = filter(networks_motif, Motif == 1), aes(x = 0, y = Count, color = "observed network"), pch = 1, size = 2, stroke = 2, inherit.aes = F) +
        scale_fill_manual(values = c("head" = "#FF0000A0", "body" = "#A0A0A0A0", "tail" = "#FF0000A0"),
                          labels = c("head" = "top 5%", "body" = "middle", "tail" = "bottom 5%")) +
        scale_color_manual(values = c("observed network" = "red")) +
        coord_flip() +
        facet_grid(.~Motif) +
        scale_x_continuous(breaks = scales::pretty_breaks(n = 3)) +
        theme_cowplot() +
        theme(panel.background = element_rect(color = 1, size = 1.5, fill = NA),
              panel.spacing = unit(0, "mm"), strip.text = element_blank(),
              legend.position = "top",
              axis.title = element_text(size = 10), axis.text = element_text(size = 8),
              plot.background = element_rect(fill = "white", color = NA)) +
        #guides(fill = guide_legend(title = "randomized network"), color = guide_legend(title = "")) +
        guides(fill = "none", color = "none") +
        labs(x = "Probability density", y = "Count of nontransitive motif")

    p <- p_motif_count
    ggsave(here::here("plots/Fig1S2-nontransitive_motif_counts.png"), p, width = 2, height = 3)


    # Figure 1S3: hierarchy metrics
    communities_hierarchy <- read_csv(here::here("data/output/communities_hierarchy.csv")) %>% mutate(Community = factor(Community, communities$Community))
    p <- communities_hierarchy %>%
        pivot_longer(cols = matches("\\d"), names_to = c("Variable", "Metric"), names_pattern = "(.*)(\\d)") %>%
        pivot_wider(names_from = Variable) %>%
        mutate(Significance = as.logical(Significance)) %>%
        #filter(Metric == 1) %>%
        ggplot(aes(x = Metric, y = HierarchyScore)) +
        geom_boxplot(width = .5, lwd = .8) +
        geom_jitter(shape = 1, size = 2, width = .2, stroke = .8) +
        scale_x_discrete (position = "bottom", labels = c("1" = "Higgins2017", "2" = "Rank-following\npairs", "3" = "Fraction of transitive motifs")) +
        scale_y_continuous(limits = c(0.5,1), breaks = c(0.5, 0.75, 1)) +
        #facet_wrap(.~Metric, scales = "free", nrow = 1) +
        theme_classic() +
        theme(panel.grid.major.x = element_line(color = "grey", linetype = 2),
              axis.title.x = element_blank(), strip.text = element_blank(),
              axis.text = element_text(size = 8, color = 1),
              axis.ticks.x = element_blank(),
              panel.border = element_rect(fill = NA, color = 1, size = 1.5)) +
        labs(x = "Metrics", y = "Hierarchy score")

    ggsave(here::here("plots/Fig1S3-hierarchy_metrics.png"), p, width = 5, height = 5)


    # Figure 1S5: pairwise competition cartoon
    pS5 <- ggdraw() + draw_image(here::here("plots/cartoons/FigS2.png")) + theme(plot.background = element_rect(fill = "white", color = NA))
    ggsave(here::here("plots/Fig1S5-pairwise_experiment_cartoon.png"), pS5, width = 4, height = 3)

    # Figure 1S6: finer grain pairwise coexistence
    pairs_example_outcomes_finer <- read_csv(here::here("data/output/pairs_example_outcomes_finer.csv"))
    ## Plot pairs example dynamics
    p_pairs_example_outcomes_finer <- pairs_example_outcomes_finer %>%
        left_join(pairs_freq, by = c("Community", "Isolate1", "Isolate2")) %>%
        mutate(Isolate1InitialODFreq = factor(Isolate1InitialODFreq)) %>%
        mutate(InteractionTypeFiner = factor(InteractionTypeFiner, names(assign_interaction_color(level = "finer")))) %>%
        ggplot(aes(x = Time, y = Isolate1MeasuredFreq, color = Isolate1InitialODFreq, group = Isolate1InitialODFreq)) +
        geom_point(size = 2) +
        geom_line(size = 1) +
        scale_y_continuous(breaks = c(0, .5, 1), limits = c(0,1)) +
        scale_color_manual(values = frequency_color, label = c("95%", "50%", "5%")) +
        facet_grid(.~InteractionTypeFiner) +
        theme_bw() +
        theme(panel.spacing = unit(2, "mm"), strip.text.x = element_blank(),
              panel.border = element_rect(color = 1, fill = NA, size = 1)) +
        guides(color = "none") +
        labs(x = "Time", y = "Frequency")
    ## The frequencies of coexistence vs. exclusion
    temp <- pairs %>% filter(Assembly == "self_assembly") %>%
        mutate(InteractionTypeFiner = factor(InteractionTypeFiner, names(assign_interaction_color(level = "finer")))) %>%
        group_by(InteractionTypeFiner) %>% summarize(Count = n()) %>% ungroup() %>% mutate(Fraction = Count / sum(Count))
    p_pairs_interaction_finer <- temp %>%
        ggplot() +
        geom_col(aes(x = InteractionTypeFiner, y = Count, fill = InteractionTypeFiner), color = 1) +
        geom_text(x = -Inf, y = Inf, label = paste0("n = ", sum(temp$Count)), vjust = 1, hjust = -0.1) +
        geom_text(aes(x = InteractionTypeFiner, y = Count, label = paste0(round(Fraction, 3) * 100,"%")), nudge_y = 6) +
        scale_fill_manual(values = assign_interaction_color(level = "finer")) +
        scale_x_discrete(breaks = c("competitive exclusion", "stable coexistence", "neutrality", "mutual exclusion", "frequency-dependent coexistence"),
                         labels = c("competitive\nexclusion", "stable\ncoexistence", "neutrality", "mutual\nexclusion", "frequency-dependent\ncoexistence")) +
        scale_y_continuous(limits = c(0, 150), expand = c(0,0)) +
        theme_classic() +
        theme(axis.title.x = element_blank(), axis.title.y = element_text(size = 10),
              axis.text.x = element_text(size = 10, color = "black"), axis.text.y = element_text(color = "black"),
              legend.position = "none") +
        labs(x = "", y = "Number of pairs", fill = "")
    pS6 <- plot_grid(p_pairs_interaction_finer, p_pairs_example_outcomes_finer, ncol = 1, axis = "lf", align = "h", rel_heights = c(2,1))
    ggsave(here::here("plots/Fig1S6-pairwise_outcomes_finer.png"), pS6, width = 7, height = 4)



    # Figure 1S7: total motif count and example motifs
    ## motif diagram
    load(here::here("data/output/motif_list.Rdata"))
    p_motif_list <- lapply(motif_list, function(x) plot_competitive_network(x, node_size = 3))
    p2 <- plot_grid(plotlist = p_motif_list, nrow = 1)

    ## motif counts
    networks_motif <- read_csv(here::here("data/output/networks_motif.csv")) %>% filter(str_detect(Community, "C\\d"))
    networks_motif_randomized <- read_csv(here::here("data/output/networks_motif_randomized.csv")) %>% filter(str_detect(Community, "C\\d")) %>% mutate(Community = ordered(Community, communities$Community))
    networks_motif_randomized_total_percentile <- networks_motif_randomized %>%
        group_by(Replicate, Motif) %>%
        summarize(Count = sum(Count)) %>%
        group_by(Motif) %>%
        arrange(Motif, Count) %>%
        slice(c(1000 * 0.05, 1000 * 0.95)) %>%
        mutate(Percentile = c("p5", "p95")) %>%
        select(Motif, Count, Percentile)
    networks_motif_total <- networks_motif %>%
        group_by(Motif) %>%
        summarize(Count = sum(Count)) %>%
        left_join(pivot_wider(networks_motif_randomized_total_percentile, names_from = Percentile, values_from = Count)) %>%
        mutate(Sign = case_when(Count > p95 ~ "top",
                                Count < p5 ~ "bottom",
                                Count < p95 & Count > p5 ~ "n.s."))

    p3 <- networks_motif_randomized %>%
        group_by(Motif, Replicate) %>%
        summarize(Count = sum(Count)) %>%
        group_by(Motif) %>%
        # 5% and 95% percentiles in randomized networks
        #mutate(p5 = quantile(Count, 0.05), p95 = quantile(Count, 0.95), ColoredTails = ifelse(Count <= p5, "tail", ifelse(Count >= p95, "head", "body"))) %>%
        ggplot() +
        geom_rect(data = networks_motif_total, aes(fill = Sign), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, alpha = .2) +
        geom_vline(xintercept = 0, color = 1) +
        geom_histogram(aes(y = Count), binwidth = 2, color = 1, fill = "white") +
        geom_point(data = networks_motif_total, x = 0, aes(y = Count, color = "observed network"), pch = 1, size = 2, stroke = 2, inherit.aes = F) +
        scale_color_manual(values = c("observed network" = "red")) +
        facet_grid(.~Motif, scales = "free_x") +
        scale_fill_manual(values = c("top" = "blue", "bottom" = "red", "n.s." = "grey")) +
        scale_x_continuous(breaks = scales::pretty_breaks(n = 2)) +
        theme_cowplot() +
        theme(panel.background = element_rect(color = 1, size = 1), panel.spacing = unit(0, "mm"),
              strip.background = element_rect(color = NA, fill = NA, size = 1)) +
        guides(color = "none", fill = "none") +
        labs(x = "Probability density", y = "Motif count")

    #p_S4 <- plot_grid(p1, p2, p3, ncol = 1, rel_heights = c(4,1,2), axis = "lr", align = "v", labels = LETTERS[1:3]) + theme(plot.background = element_rect(fill = "white", color = NA))
    pS7 <- plot_grid(p2, p3, ncol = 1, rel_heights = c(1,3), axis = "lr", align = "v", labels = LETTERS[1:2]) + theme(plot.background = element_rect(fill = "white", color = NA))
    ggsave(here::here("plots/Fig1S7-total_motif_counts.png"), pS7, width = 10, height = 5)




    # Figure 1S8: network matrices and motifs
    ## Matrix and graph
    networks_motif <- read_csv(here::here("data/output/networks_motif.csv")) %>% filter(str_detect(Community, "C\\d"))
    networks_motif_randomized <- read_csv(here::here("data/output/networks_motif_randomized.csv")) %>% filter(str_detect(Community, "C\\d")) %>% mutate(Community = ordered(Community, communities$Community))
    networks_motif_randomized_percentile <- read_csv(here::here("data/output/networks_motif_randomized_percentile.csv")) %>% filter(str_detect(Community, "C\\d")) %>% mutate(Community = ordered(Community, communities$Community))
    p_net_matrix_list <- lapply(net_list, function(x) plot_adjacent_matrix(x) + theme(plot.margin = grid::unit(c(5,0,3,0), "mm")))
    p_net_list <- lapply(net_list, function(x) plot_competitive_network(x, node_size = 2) + theme(plot.background = element_rect(fill = NA)))
    p_list <- rep(list(NA), length(p_net_list))
    for (i in 1:length(net_list)) p_list[[i]] <- ggdraw(p_net_matrix_list[[i]]) + draw_plot(plot = p_net_list[[i]], x = -.1, y = -.1, width = 0.7, height = 0.7)
    ## Motif count
    plot_motif_count <- function (x = 1) {
        motif_randomized_subset <- networks_motif_randomized_percentile %>%
            filter(Community %in% communities$Community[x]) %>%
            mutate(Community = factor(Community, communities$Community[x]))
        motif_community_subset <- networks_motif %>%
            filter(Community %in% communities$Community[x]) %>%
            mutate(Community = factor(Community, communities$Community[x]))

        ggplot() +
            # 5% and 95% percentiles in randomized networks
            geom_point(data = motif_randomized_subset, aes(x = Motif, y = Count, group = Motif, color = "randomized network")) +
            geom_segment(data = motif_randomized_subset %>% pivot_wider(id_cols = c(Community, Motif), names_from = Percentile, values_from = Count),
                         aes(x = Motif, xend = Motif, y = p5, yend = p95, color = "randomized network")) +
            # Observations
            geom_point(data = motif_community_subset, aes(x = Motif, y = Count, color = "observed network")) +
            scale_x_continuous(breaks = 1:7) +
            scale_color_manual(values = c("observed network" = "red", "randomized network" = "black"))+
            #facet_wrap(Community ~., scale = "free_y", nrow = 1)  +
            theme_classic() +
            theme(panel.background = element_rect(color = 1, size = 1), legend.position = "none")
    }
    p_motif_count_list <- rep(list(NA), length(p_net_list))
    for (i in 1:length(p_net_list)) {
        if (i %in% c(1, 6)) p_motif_count_list[[i]] <- plot_motif_count(i) + theme(axis.title.x = element_blank())
        if (i %in% c(2:5, 7:10)) p_motif_count_list[[i]] <- plot_motif_count(i) + theme(axis.title = element_blank())
        if (i == 11) p_motif_count_list[[i]] <- plot_motif_count(i)
        if (i %in% 12:13) p_motif_count_list[[i]] <- plot_motif_count(i) + theme(axis.title.y = element_blank())
    }

    ## Example matrix
    plot_example_matrix <- function(graph) {
        graph_ranked <- graph %>%
            activate(nodes) %>%
            select(Isolate, PlotRank) %>%
            activate(edges) %>%
            mutate(
                fromRank = .N()$PlotRank[match(from, .N()$Isolate)],
                toRank = .N()$PlotRank[match(to, .N()$Isolate)])
        n_nodes <- igraph::vcount(graph_ranked)
        interaction_color <- assign_interaction_color(level = "matrix")

        graph_ranked %>%
            filter(fromRank <= toRank) %>%
            bind_edges(tibble(from = 1:n_nodes, to = 1:n_nodes, fromRank = 1:n_nodes, toRank = 1:n_nodes, InteractionType = "self")) %>%
            as_tibble() %>%
            ggplot() +
            geom_tile(aes(x = factor(toRank), y = ordered(fromRank, n_nodes:1), fill = InteractionType), width = 0.9, height = 0.9) +
            scale_x_discrete(position = "top", labels = c("top\nrank", rep("", n_nodes-2), "bottom\nrank")) +
            scale_y_discrete(position = "right", labels = c("bottom\nrank", rep("", n_nodes-2), "top\nrank")) +
            scale_fill_manual(breaks = c("exclusion", "coexistence", "self"), values = c(assign_interaction_color(), "self" = "black")) +
            theme_bw() +
            theme(axis.ticks = element_blank(), axis.title = element_blank(), legend.title = element_blank(),
                  axis.text = element_text(size = 10, color = 1),
                  panel.border = element_blank(), panel.grid = element_blank())
    }
    p_example_matrix <- plot_example_matrix(net_list[[1]]) + guides(fill = "none") + theme(plot.background = element_rect(fill = "grey90"), panel.background = element_rect(fill = "grey90"))
    ## Get legend for matrix
    shared_legend_matrix <- get_legend(plot_example_matrix(net_list[[1]]) + theme(legend.text = element_text(size = 15), legend.justification = "right", plot.background = element_rect(fill = "grey90"), panel.background = element_rect(fill = "grey90")))
    ## Get legend for line
    p_temp <- plot_motif_count(1) + theme(legend.position = "right", legend.title = element_blank(), legend.text = element_text(size = 15))
    shared_legend_line <- cowplot::get_legend(p_temp)
    p1 <- plot_grid(plotlist = c(list(p_example_matrix), list(shared_legend_matrix), list(NULL), list(shared_legend_line), rep(list(NULL), 1)), nrow = 1) + theme(plot.background = element_rect(fill = "white", color = NA))
    p2 <- list(p_list[1:5], p_motif_count_list[1:5],
               p_list[6:10], p_motif_count_list[6:10],
               p_list[11:13], rep(list(NULL), 2),
               p_motif_count_list[11:13]) %>%
        unlist(recursive = F) %>%
        plot_grid(plotlist = ., labels = c(communities$Community[1:5], rep("", 5), communities$Community[6:10], rep("", 5), communities$Community[11:13], rep("", 7)),
                  ncol = 5, axis = "tbrl", align = "v")
    pS8 <- plot_grid(p1, p2, ncol = 1, rel_heights = c(1,6)) + theme(plot.background = element_rect(fill = "white", color = NA))
    ggsave(here::here("plots/Fig1S8-networks_matrix.png"), pS8, width = 10, height = 12)


    # Figure 1S9: hierarchy score for communities, compared to 1000 randomized ones
    ## 1. Higgins2017
    ## 2. Violation of ranks
    ## 3. Fraction of transitive motifs
    communities_hierarchy <- read_csv(here::here("data/output/communities_hierarchy.csv")) %>% mutate(Community = factor(Community, communities$Community))
    communities_hierarchy_randomized <- read_csv(here::here("data/output/communities_hierarchy_randomized.csv")) %>% mutate(Community = factor(Community, communities$Community))
    networks_motif <- read_csv(here::here("data/output/networks_motif.csv")) %>% filter(Motif == 2) %>% filter(str_detect(Community, "C\\d"))
    networks_motif_randomized <- read_csv(here::here("data/output/networks_motif_randomized.csv")) %>% filter(Motif == 2) %>% filter(str_detect(Community, "C\\d")) %>% mutate(Community = ordered(Community, communities$Community))

    b = 1000
    networks_motif_pvalue <- networks_motif_randomized %>%
        group_by(Community) %>%
        # Add a bottom row to each community
        bind_rows(communities %>% select(Community) %>% mutate(Replicate = b+1, Motif = 2, Count = -1, Fraction = -.001) %>% filter(str_detect(Community, "C\\d"))) %>%
        arrange(Community, desc(Fraction)) %>%
        left_join(networks_motif %>% select(Community, Motif, FractionObserved = Fraction)) %>%
        mutate(Percentile = 0:b / b) %>%
        filter(FractionObserved > Fraction) %>%
        slice(1) %>%
        mutate(Significance = Percentile < 0.05) %>%
        select(Community, HierarchyScore3 = FractionObserved, Percentile3 = Percentile, Significance3 = Significance) %>%
        mutate(Community = factor(Community, communities$Community)) %>%
        arrange(Community)

    communities_hierarchy <- communities_hierarchy %>% left_join(networks_motif_pvalue) %>%
        pivot_longer(cols = matches("\\d"), names_to = c("Variable", "Metric"), names_pattern = "(.*)(\\d)") %>%
        pivot_wider(names_from = Variable) %>%
        mutate(Significance = as.logical(Significance))

    plot_top <- function(x) {
        communities_hierarchy %>%
            filter(Metric == x) %>%
            ggplot(aes(x = Metric, y = HierarchyScore)) +
            geom_boxplot(width = .3) +
            geom_jitter(shape = 1, size = 2, width = .1) +
            scale_x_discrete (position = "top", labels = c("1" = "Higgins et al 2017", "2" = "Rank-following pairs", "3" = "Fraction of transitive motifs")) +
            scale_y_continuous(limits = c(0,1), breaks = c(0, 0.5, 1)) +
            facet_wrap(.~Metric, scales = "free", nrow = 1) +
            theme_classic() +
            theme(axis.title.x = element_blank(), strip.text = element_blank(),
                  axis.text.x = element_text(size = 13), axis.ticks.x = element_blank(),
                  panel.border = element_rect(fill = NA, color = 1, size = 1.5)) +
            labs(x = "Metrics", y = "Score")
    }
    plot_bottom <- function(x) {
        communities_hierarchy %>%
            filter(Metric == x) %>%
            ggplot() +
            geom_rect(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, aes(fill = Significance), alpha = .5) +
            geom_histogram(data = communities_hierarchy_randomized, aes(y = HierarchyScore1), color = 1, fill = "white") +
            geom_hline(aes(yintercept = HierarchyScore), color = "red") +
            geom_text(aes(x = Inf, y = -Inf, label = paste0("p=", Percentile)), vjust = -1, hjust = 1.5) +
            scale_x_continuous(breaks = scales::pretty_breaks(n = 3)) +
            scale_y_continuous(breaks = c(0.5, 0.75, 1)) +
            scale_fill_manual(values = c("TRUE" = "#DB7469", "FALSE" = "white")) +
            facet_wrap(Community ~., ncol = 2, scales = "free_x", dir = "v") +
            theme_classic() +
            theme(panel.grid.major = element_blank(), legend.position = "none", panel.border = element_rect(color = 1, fill = NA, size = 1.5)) +
            labs(x = "Count", y = "Hierarchy score")
    }
    p_top <- plot_grid(plot_top(1), plot_top(2), plot_top(3), nrow = 1, labels = LETTERS[1:3], align = "hv")
    p_bottom <- plot_grid(plot_bottom(1), plot_bottom(2), plot_bottom(3), nrow = 1)

    if (FALSE) {
        p1 <- communities_hierarchy %>%
            filter(Metric == 1) %>%
            ggplot() +
            geom_rect(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, aes(fill = Significance), alpha = .5) +
            geom_histogram(data = communities_hierarchy_randomized, aes(y = HierarchyScore1), color = 1, fill = "white") +
            geom_hline(aes(yintercept = HierarchyScore), color = "red") +
            geom_text(aes(x = Inf, y = -Inf, label = paste0("p=", Percentile)), vjust = -1, hjust = 1.5) +
            scale_x_continuous(breaks = scales::pretty_breaks(n = 3)) +
            scale_y_continuous(breaks = c(0.5, 0.75, 1)) +
            scale_fill_manual(values = c("TRUE" = "#DB7469", "FALSE" = "white")) +
            facet_wrap(Community ~., ncol = 2, scales = "free_x", dir = "v") +
            theme_classic() +
            theme(panel.grid.major = element_blank(), legend.position = "none", panel.border = element_rect(color = 1, fill = NA, size = 1.5)) +
            labs(x = "Count", y = "Hierarchy score") +
            ggtitle("Higgins et al 2017")

        p2 <- communities_hierarchy %>%
            filter(Metric == 2) %>%
            ggplot() +
            geom_rect(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, aes(fill = Significance), alpha = .5) +
            geom_histogram(data = communities_hierarchy_randomized, aes(y = HierarchyScore2), color = 1, fill = "white") +
            geom_hline(aes(yintercept = HierarchyScore), color = "red") +
            geom_text(aes(x = Inf, y = -Inf, label = paste0("p=", Percentile)), vjust = -1, hjust = 1.5) +
            scale_x_continuous(breaks = scales::pretty_breaks(n = 3)) +
            scale_y_continuous(breaks = c(0, 0.5, 1)) +
            scale_fill_manual(values = c("TRUE" = "#DB7469", "FALSE" = "white")) +
            facet_wrap(Community ~., ncol = 2, scales = "free_x", dir = "v") +
            theme_classic() +
            theme(panel.grid.major = element_blank(), legend.position = "none", panel.border = element_rect(color = 1, fill = NA, size = 1.5)) +
            labs(x = "Count", y = "Hierarchy score") +
            ggtitle("Violation of ranks")

        p3 <- communities_hierarchy %>%
            filter(Metric == 3) %>%
            ggplot() +
            geom_rect(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, aes(fill = Significance), alpha = .5) +
            geom_histogram(data = networks_motif_randomized, aes(y = Fraction), color = 1, fill = "white") +
            geom_hline(aes(yintercept = HierarchyScore), color = "red") +
            geom_text(aes(x = Inf, y = Inf, label = paste0("p=", Percentile)), vjust = 1.5, hjust = 1.5) +
            scale_x_continuous(breaks = scales::pretty_breaks(n = 3)) +
            scale_y_continuous(breaks = c(0,0.5,1)) +
            scale_fill_manual(values = c("TRUE" = "#DB7469", "FALSE" = "white")) +
            facet_wrap(Community ~., ncol = 2, scales = "free_x", dir = "v") +
            theme_classic() +
            theme(panel.grid.major = element_blank(), legend.position = "none", panel.border = element_rect(color = 1, fill = NA, size = 1.5)) +
            labs(x = "Count", y = "Hierarchy score") +
            ggtitle("Fraction of transitive motifs")
    }
    pS9 <- plot_grid(p_top, p_bottom, ncol = 1, rel_heights = c(1, 4))

    ggsave(here::here("plots/Fig1S9-hierarchy.png"), pS9, width = 9, height = 13)



    # Figure 1S10: diagonal analysis
    networks_diag <- read_csv(here::here("data/output/networks_diag.csv"), col_types = cols())
    networks_diag_randomized <- read_csv(here::here("data/output/networks_diag_randomized.csv"), col_types = cols())
    # Overall
    networks_diag_sum <- networks_diag %>%
        group_by(DistanceToDiagonal) %>%
        summarize(ObservedCountCoexistenceSum = sum(CountCoexistence),
                  ObservedCountTotalSum = sum(TotalCount)) %>%
        mutate(ObservedFractionCoexistenceSum = ObservedCountCoexistenceSum/ObservedCountTotalSum)
    networks_diag_randomized_sum <- networks_diag_randomized %>%
        group_by(Replicate, DistanceToDiagonal) %>%
        summarize(CountCoexistenceSum = sum(CountCoexistence)) %>%
        left_join(select(networks_diag_sum, DistanceToDiagonal, CountTotalSum = ObservedCountTotalSum)) %>%
        mutate(FractionCoexistenceSum = CountCoexistenceSum/CountTotalSum)
    ## Statistics
    stat_diag <- networks_diag_randomized_sum %>%
        group_by(DistanceToDiagonal) %>%
        # find percentile
        bind_rows(tibble(Replicate = 0, DistanceToDiagonal = 1:11, CountCoexistenceSum = 0)) %>%
        left_join(networks_diag_sum) %>%
        arrange(DistanceToDiagonal, desc(CountCoexistenceSum)) %>%
        mutate(Percentile = (1:n())/n()) %>%
        filter(CountCoexistenceSum <= ObservedCountCoexistenceSum) %>%
        group_by(DistanceToDiagonal) %>%
        slice(1) %>%
        select(DistanceToDiagonal, Percentile) %>%
        # Asterisk
        mutate(Significance = case_when(Percentile < 0.001 | Percentile > 0.999 ~ "***",
                                        Percentile < 0.01 | Percentile > 0.99 ~ "**",
                                        Percentile < 0.05 | Percentile > 0.95 ~ "*",
                                        Percentile > 0.05 & Percentile < 0.95 ~ "n.s."),
               Sign = case_when(Percentile < 0.05 ~ "top",
                                Percentile > 0.95 ~ "bottom",
                                Percentile > 0.05 & Percentile < 0.95 ~ "n.s."))

    p1 <- networks_diag_sum %>%
        ggplot() +
        # Asterisk
        geom_text(data = stat_diag, aes(x = DistanceToDiagonal, y = Inf, label = Significance), vjust = 2) +
        geom_rect(data = stat_diag, aes(xmin = DistanceToDiagonal-0.5, xmax = DistanceToDiagonal+0.5, fill = Sign), ymin = -Inf, ymax = Inf, alpha = .2) +
        # Random networks
        geom_boxplot(data = networks_diag_randomized_sum, aes(x = DistanceToDiagonal, y = FractionCoexistenceSum, group = DistanceToDiagonal, color = "randomized network"),
                     outlier.size = 1) +
        geom_jitter(data = networks_diag_randomized_sum, aes(x = DistanceToDiagonal, y = FractionCoexistenceSum, group = DistanceToDiagonal, color = "randomized network"),
                    size = .1, alpha = 0.5, width = .3) +
        # Observed networks
        geom_point(aes(x = DistanceToDiagonal, y = ObservedFractionCoexistenceSum, group = DistanceToDiagonal, color = "observed network"), size = 2) +
        geom_line(aes(x = DistanceToDiagonal, y = ObservedFractionCoexistenceSum, color = "observed network")) +
        scale_x_continuous(breaks = 1:11, expand = c(0,0)) +
        scale_y_continuous(limits = c(0, 1.1), breaks = c(0, .5, 1)) +
        scale_color_manual(values = c("observed network" = "red", "randomized network" = "black")) +
        scale_fill_manual(values = c("top" = "blue", "bottom" = "red", "n.s." = "grey")) +
        theme_classic() +
        theme(legend.position = "right", legend.title = element_blank(),
              panel.border = element_rect(fill = NA, color = 1, size = 1.5),
              strip.background = element_blank()) +
        guides(fill = "none") +
        labs(x = "|i-j|", y = "Fraction of pairwise coexitence")

    # Each community
    networks_diag_sum <- networks_diag %>%
        mutate(Community = factor(Community, community_factor)) %>%
        group_by(Community, DistanceToDiagonal) %>%
        summarize(ObservedCountCoexistenceSum = sum(CountCoexistence),
                  ObservedCountTotalSum = sum(TotalCount)) %>%
        mutate(ObservedFractionCoexistenceSum = ObservedCountCoexistenceSum/ObservedCountTotalSum)
    networks_diag_randomized_sum <- networks_diag_randomized %>%
        mutate(Community = factor(Community, community_factor)) %>%
        group_by(Community, Replicate, DistanceToDiagonal) %>%
        summarize(CountCoexistenceSum = sum(CountCoexistence)) %>%
        left_join(select(networks_diag_sum, DistanceToDiagonal, CountTotalSum = ObservedCountTotalSum)) %>%
        mutate(FractionCoexistenceSum = CountCoexistenceSum/CountTotalSum)
    ## Statistics
    temp <- tibble(Replicate = 0, DistanceToDiagonal = 1:11, CountCoexistenceSum = 0) %>%
        slice(rep(1:11, 13)) %>%
        mutate(Community = rep(str_match(communities$Community, "C\\d+R\\d+"), each = 11)) %>%
        mutate(Community = factor(Community, community_factor))
    stat_diag <- networks_diag_randomized_sum %>%
        mutate(Community = factor(Community, community_factor)) %>%
        group_by(Community, DistanceToDiagonal) %>%
        # find percentile
        bind_rows(temp) %>%
        left_join(networks_diag_sum) %>%
        arrange(DistanceToDiagonal, desc(CountCoexistenceSum)) %>%
        mutate(Percentile = (1:n())/n()) %>%
        filter(CountCoexistenceSum <= ObservedCountCoexistenceSum) %>%
        group_by(Community, DistanceToDiagonal) %>%
        slice(1) %>%
        select(DistanceToDiagonal, Percentile) %>%
        # Asterisk
        mutate(Significance = case_when(Percentile < 0.001 | Percentile > 0.999 ~ "***",
                                        Percentile < 0.01 | Percentile > 0.99 ~ "**",
                                        Percentile < 0.05 | Percentile > 0.95 ~ "*",
                                        Percentile > 0.05 & Percentile < 0.95 ~ "n.s."),
               Sign = case_when(Percentile < 0.05 ~ "top",
                                Percentile > 0.95 ~ "bottom",
                                Percentile > 0.05 & Percentile < 0.95 ~ "n.s."))

    p2 <- networks_diag_sum %>%
        ggplot() +
        # Asterisk
        geom_text(data = stat_diag, aes(x = DistanceToDiagonal, y = Inf, label = Significance), vjust = 2) +
        geom_rect(data = stat_diag, aes(xmin = DistanceToDiagonal-0.5, xmax = DistanceToDiagonal+0.5, fill = Sign), ymin = -Inf, ymax = Inf, alpha = .2) +
        # Random networks
        geom_boxplot(data = networks_diag_randomized_sum, aes(x = DistanceToDiagonal, y = FractionCoexistenceSum, group = DistanceToDiagonal, color = "randomized network"),
                     outlier.size = 1) +
        geom_jitter(data = networks_diag_randomized_sum, aes(x = DistanceToDiagonal, y = FractionCoexistenceSum, group = DistanceToDiagonal, color = "randomized network"),
                    size = .1, alpha = 0.5, width = .3) +
        # Observed networks
        geom_point(aes(x = DistanceToDiagonal, y = ObservedFractionCoexistenceSum, group = DistanceToDiagonal, color = "observed network"), size = 2) +
        geom_line(aes(x = DistanceToDiagonal, y = ObservedFractionCoexistenceSum, color = "observed network")) +
        scale_x_continuous(breaks = 1:11) +
        scale_y_continuous(limits = c(0, 1.1), breaks = c(0, .5, 1)) +
        scale_color_manual(values = c("observed network" = "red", "randomized network" = "black")) +
        scale_fill_manual(values = c("top" = "blue", "bottom" = "red", "n.s." = "grey")) +
        facet_wrap(Community~., nrow = 3, scales = "free_x") +
        theme_classic() +
        theme(legend.position = "top", legend.title = element_blank(),
              panel.border = element_rect(fill = NA, color = 1, size = 1.5),
              strip.background = element_blank()) +
        guides(fill = "none", color = "none") +
        labs(x = "|i-j|", y = "Fraction of pairwise coexitence")

    p_top <- plot_grid(p1, NULL, nrow = 1)
    pS10 <- plot_grid(p_top, p2, ncol = 1, labels = LETTERS[1:2],
                      rel_heights = c(1,2), scale = c(1, 1)) + paint_white_background()
    ggsave(here::here("plots/Fig1S10-diagonal_analysis.png"), pS10, width = 12, height = 12)
    if (FALSE) {
        # Count
        pS10 <- networks_diag_sum %>%
            ggplot() +
            # Random networks
            geom_boxplot(data = networks_diag_randomized_sum, aes(x = DistanceToDiagonal, y = CountCoexistenceSum, group = DistanceToDiagonal, color = "randomized network"),
                         outlier.size = 1) +
            geom_jitter(data = networks_diag_randomized_sum, aes(x = DistanceToDiagonal, y = CountCoexistenceSum, group = DistanceToDiagonal, color = "randomized network"),
                        size = .1, alpha = 0.5, width = .3) +
            # Observed networks
            geom_point(aes(x = DistanceToDiagonal, y = ObservedCountCoexistenceSum, group = DistanceToDiagonal, color = "observed network"), size = 2) +
            geom_line(aes(x = DistanceToDiagonal, y = ObservedCountCoexistenceSum, color = "observed network")) +
            # Asterisk
            geom_text(data = stat_diag, aes(x = DistanceToDiagonal, y = Inf, label = Significance), vjust = 2) +
            scale_x_continuous(breaks = 1:11) +
            scale_y_continuous(limits = c(0, 40)) +
            scale_color_manual(values = c("observed network" = "red", "randomized network" = "black"))+
            theme_classic() +
            theme(legend.position = "top", legend.title = element_blank(),
                  panel.border = element_rect(fill = NA, color = 1, size = 1.5)) +
            #labs(x = "Distance to diagonal (|i-j|)", y = "Count of pairwise coexistence")
            labs(x = "Difference in rankings", y = "Count of pairwise coexistence")
    }



    # Figure 1S11: coexistence pair count for both random assembly and self-assembly
    temp <- pairs %>%
        mutate(Assembly = factor(Assembly, c("random_assembly", "across_community", "self_assembly"))) %>%
        filter(Assembly %in% c("random_assembly", "self_assembly")) %>%
        group_by(Assembly, InteractionType) %>% summarize(Count = n()) %>% mutate(Fraction = Count / sum(Count))
    n_size <- temp %>% summarize(Count = sum(Count))
    pS11 <- temp %>%
        ggplot() +
        geom_col(aes(x = Assembly, y = Count, fill = InteractionType), position = "fill", color = 1, width = 0.8) +
        scale_fill_manual(values = assign_interaction_color(level = "simple")) +
        scale_x_discrete(labels = c("random_assembly" = paste0("random assembly\nn=", pull(filter(n_size, Assembly == "random_assembly"), Count)),
                                    "self_assembly" = paste0("self assembly\nn=", pull(filter(n_size, Assembly == "self_assembly"), Count)))) +
        scale_y_continuous(expand = c(0,0), breaks = c(0, .5, 1)) +
        theme_classic() +
        theme(axis.title.x = element_blank(), legend.position = "top", axis.text.x = element_text(size = 10)) +
        labs(x = "", y  = "Fraction", fill = "")

    ggsave(here::here("plots/Fig1S11-pairwsie_competition_assembly.png"), pS11, width = 3, height = 3)

    ## Stat: does assembly explain for coexistence ratio?
    observed_stat <- pairs %>%
        filter(Assembly %in% c("random_assembly", "self_assembly")) %>%
        chisq_test(InteractionType ~ Assembly) %>% pull(statistic)
    null_stat <- pairs %>%
        filter(Assembly %in% c("random_assembly", "self_assembly")) %>%
        specify(InteractionType ~ Assembly, success = "coexistence") %>%
        hypothesize(null = "independence") %>%
        generate(reps = 1000, type = "permute") %>%
        calculate(stat = "Chisq", order = c("random_assembly", "self_assembly"))
    null_stat %>%
        get_p_value(obs_stat = observed_stat, direction = "two-sided")





    # Figure 1S12: fraction of pairwise coexistence across communities
    pairs_count <- pairs %>%
        filter(Assembly == "self_assembly") %>%
        group_by(Community) %>%
        summarize(Count = n())

    pS12 <- pairs %>%
        mutate(Community = factor(Community, communities$Community)) %>%
        filter(Assembly == "self_assembly") %>%
        ggplot() +
        geom_bar(aes(x = Community, fill = InteractionType), color = 1, position = position_fill(), size = .5) +
        geom_text(data = pairs_count, aes(x = Community, y = 1, label = Count), vjust = -.5) +
        geom_rect(xmin = -Inf, xmax = Inf, ymin = 0, ymax = 1, color = grey(0.1), fill = NA, size = .5) +
        scale_fill_manual(values = assign_interaction_color()) +
        scale_y_continuous(breaks = c(0,.5,1), limit = c(0, 1.15), expand = c(0,0)) +
        theme_classic() +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
              axis.line.y = element_blank(),
              legend.title = element_blank(), legend.position = "top") +
        labs(y = "Fraction")
    ggsave(here::here("plots/Fig1S12-pairs_counts.png"), pS12, width = 4, height = 3)


    # Figure 1S13: relative abundance within pairs
    pairs_freq_FN <- pairs_freq %>%
        filter(Isolate1InitialODFreq == 50) %>%
        filter(Time == "T8") %>%
        left_join(pairs) %>%
        filter(PairFermenter == "FN", InteractionType == "coexistence") %>%
        mutate(FractionFermenter = ifelse(Fermenter1, Isolate1MeasuredFreq, 1 - Isolate1MeasuredFreq),
               FractionRespirator = ifelse(Fermenter1, 1 - Isolate1MeasuredFreq, Isolate1MeasuredFreq)) %>%
        mutate(RF_Ratio = FractionRespirator / FractionFermenter) %>%
        select(Community, Isolate1, Isolate2, Isolate1InitialODFreq, Isolate1MeasuredFreq,
               PairFermenter, Fermenter1, Fermenter2,
               FractionFermenter, FractionRespirator, RF_Ratio)

    p1 <- pairs_freq_FN %>%
        ggplot() +
        geom_boxplot(aes(x = PairFermenter, y = RF_Ratio), color = 1) +
        geom_jitter(aes(x = PairFermenter, y = RF_Ratio), color = 1, shape = 1, size = 2, width = 0.2) +
        scale_y_log10(limits = c(0.01, 1), minor_breaks = rep(1:9, 4)*(10^rep(-2:1, each = 9)), guide = "prism_minor") +
        theme_classic() +
        theme(axis.title.x = element_blank(), axis.text.x = element_blank(),
              axis.ticks.y = element_line(size = .8), axis.ticks.length.y = unit(2, "mm")) +
        labs(y = "R/F")

    p2 <- pairs_freq_FN %>%
        arrange(FractionFermenter) %>%
        mutate(Pair = 1:n()) %>%
        select(Pair, FractionFermenter, FractionRespirator) %>%
        pivot_longer(cols = starts_with("Fraction"), names_prefix = "Fraction", names_to = "Fermenter", values_to = "Fraction") %>%
        mutate(Fermenter = tolower(Fermenter)) %>%
        ggplot() +
        geom_col(aes(x = Pair, y = Fraction, fill = Fermenter), color = 1) +
        scale_x_continuous(breaks = 1:22, expand = c(0,0)) +
        scale_y_continuous(breaks = c(0,.5, 1), expand = c(0,0)) +
        scale_fill_manual(values = fermenter_color) +
        theme_classic() +
        theme(legend.position = "right", legend.title = element_blank())
    pS13 <- plot_grid(p1, p2, nrow = 1, rel_widths = c(1,3), labels = LETTERS[1:2], axis = "tb", align = "h")
    ggsave(here::here("plots/Fig1S13-pair_relative_abundance.png"), pS13, width = 7, height = 3)



    # Figure 1S14: species coexist at least once. Count by the network components
    networks_component <- read_csv(here::here("data/output/networks_component.csv")) %>% mutate(Community = factor(Community, community_factor))
    networks_component_randomized <- read_csv(here::here("data/output/networks_component_randomized.csv")) %>% mutate(Community = factor(Community, community_factor))

    ## Number of clusters
    p1 <- networks_component_randomized %>%
        ggplot(aes(x = Community, y = NumberCluster, group = Community)) +
        geom_boxplot() +
        geom_point(position = position_jitter(height = 0)) +
        geom_point(data = networks_component, aes(x = Community, y = NumberCluster), color = "red") +
        theme_classic() +
        labs(y = "Number of cluster")

    ## Size of clusters
    observation <- networks_component %>%
        group_by(Community) %>%
        summarize(SdSizeCluster = sd(SizeCluster, na.rm = T)) %>%
        replace_na(list(SdSizeCluster = 0))
    permutation <- networks_component_randomized %>%
        group_by(Community, Replicate) %>%
        summarize(SdSizeCluster = sd(SizeCluster, na.rm = T)) %>%
        replace_na(list(SdSizeCluster = 0))
    p2 <- permutation %>%
        ggplot(aes(x = Community, y = SdSizeCluster, group = Community)) +
        geom_boxplot() +
        geom_point(position = position_jitter(height = 0)) +
        geom_point(data = observation, aes(x = Community, y = SdSizeCluster), color = "red") +
        theme_classic() +
        labs(y = "SD of cluster size")

    pS14 <- plot_grid(p1, p2, ncol = 1, labels = c("A", "B"))
    ggsave(here::here("plots/Fig1S14-network_component.png"), pS14, width = 10, height = 5)

    # Figure 1S15. node degree
    networks_degree <- read_csv(here::here("data/output/networks_degree.csv")) %>% mutate(Community = factor(Community, community_factor))
    networks_degree_randomized <- read_csv(here::here("data/output/networks_degree_randomized.csv")) %>% mutate(Community = factor(Community, community_factor))

    observation <- networks_degree %>%
        group_by(Community) %>%
        summarize(SdDegree = sd(Degree))

    permutation <- networks_degree_randomized %>%
        group_by(Community, Replicate) %>%
        summarize(SdDegree = sd(Degree))

    pS15 <- permutation %>%
        ggplot(aes(x = Community, y = SdDegree, group = Community)) +
        geom_boxplot() +
        geom_point(position = position_jitter(height = 0)) +
        geom_point(data = observation, aes(x = Community, y = SdDegree), color = "red") +
        theme_classic() +
        labs(y = "SD of node degree")

    if (FALSE) {
        pS15 <- networks_degree %>%
            ggplot() +
            geom_histogram(aes(x = Degree, fill = "observation"), binwidth = 1, alpha = .2) +
            geom_histogram(data = networks_degree_randomized, aes(x = Degree, y = after_stat(count)/1000, fill = "randomized network"), binwidth = 1, alpha = .2) +
            scale_x_continuous(breaks = 0:10) +
            scale_fill_manual(values = c("observation" = "red", "randomized network" = "black")) +
            theme_classic() +
            theme(legend.position = "top", legend.title = element_blank()) +
            labs(x = "Node degree", y = "Count")

    }

    ggsave(here::here("plots/Fig1S15-node_degree.png"), pS15, width = 10, height = 3)


    # Figure 1S16. Raw pair frequencies
    pairs_example_freq <- pairs %>%
        filter(str_detect(Community, "C\\d+")) %>%
        left_join(isolates %>% select(ID, Rank, PlotRank, Community) %>% rename_with(~ paste0(., 1), -Community)) %>%
        left_join(isolates %>% select(ID, Rank, PlotRank, Community) %>% rename_with(~ paste0(., 2), -Community)) %>%
        select(Community, starts_with("Isolate"), starts_with("Interaction"), starts_with("PlotRank"), starts_with("Rank")) %>%
        left_join(pairs_freq, by = c("Community", "Isolate1", "Isolate2")) %>%
        mutate(Isolate1InitialODFreq = factor(Isolate1InitialODFreq)) %>%
        mutate(InteractionType = factor(InteractionType, c("exclusion", "coexistence")))
    plot_example_freq <- function(pairs_freq) {
        # Extract params
        comm <- unique(pairs_freq$Community)
        isolate1 <- unique(pairs_freq$Isolate1)
        isolate2 <- unique(pairs_freq$Isolate2)
        interaction_type <- pairs %>%
            filter(Community == comm, Isolate1 == isolate1, Isolate2 == isolate2) %>%
            pull(InteractionType)

        #
        pairs_freq %>%
            mutate(Isolate1InitialODFreq = factor(Isolate1InitialODFreq, c(95,50,5))) %>%
            ggplot(aes(x = Time, y = Isolate1MeasuredFreq, color = Isolate1InitialODFreq, group = Isolate1InitialODFreq)) +
            geom_point(size = 2) +
            geom_line(size = 1) +
            scale_y_continuous(breaks = c(0, .5, 1), limits = c(0,1)) +
            scale_color_manual(values = frequency_color, label = c("95%", "50%", "5%")) +
            theme_bw() +
            theme(panel.spacing = unit(2, "mm"),
                  panel.border = element_rect(color = 1, fill = NA, size = 1),
                  panel.grid = element_blank(),
                  panel.grid.minor.y = element_blank(),
                  axis.title = element_blank(), axis.text = element_blank(),
                  axis.ticks = element_blank(),
                  #panel.background = element_rect(fill = "white"),
                  panel.background = element_rect(fill = alpha(ifelse(interaction_type == "coexistence", "#557BAA",
                                                                      ifelse(interaction_type == "exclusion", "#DB7469", NA)), 0.5)),
                  plot.background = element_blank()) +
            guides(color = "none") +
            labs(x = "Time", y = "Frequency") +
            #        ggtitle(paste0(isolate1, "-", isolate2)) +
            NULL
    }
    plot_community_freq <- function(pairs_example_freq) {
        comm <- unique(pairs_example_freq$Community)
        # Make figures
        temp_freq <- pairs_example_freq %>%
            filter(Community == comm) %>%
            #mutate(Isolate1 = ordered(Isolate1, temp), Isolate2 = ordered(Isolate2, temp)) %>%
            rowwise() %>%
            mutate(PlotRankOrder = paste0(min(PlotRank1, PlotRank2), "_", max(PlotRank1, PlotRank2))) %>%
            mutate(PlotRankOrder = ordered(PlotRankOrder, paste0(rep(1:12, each = 12), "_", rep(1:12, 12)))) %>%
            arrange(PlotRankOrder)
        temp_list <- temp_freq %>%
            as_tibble() %>%
            arrange(PlotRankOrder) %>%
            group_split(PlotRankOrder) %>%
            lapply(plot_example_freq)

        # Grid
        n_isolate <- communities$CommunitySize[communities$Community == comm]
        m <- matrix(NA, n_isolate-1, n_isolate-1)
        m[lower.tri(m, diag = T)] <- 1:choose(n_isolate,2)
        m <- t(m)
        arrangeGrob(grobs = temp_list, layout_matrix = m)

    }
    p_list <- pairs_example_freq %>%
        mutate(Community = factor(Community, community_factor)) %>%
        group_split(Community) %>%
        lapply(plot_community_freq)

    p_top <- plot_grid(plotlist = p_list[1:10], ncol = 3, scale = communities_size[1:10]/6, labels = 1:10) + paint_white_background()
    p_body <- plot_grid(plotlist = p_list[11:12], ncol = 1, scale = communities_size[11:13]/10, labels = 11:12) + paint_white_background()
    p_bottom <- plot_grid(plotlist = p_list[13], ncol = 1, scale = .9, labels = 13) + paint_white_background()
    ggsave(here::here("plots/Fig1S16-1-frequency_plots.png"), p_top, width = 9, height = 12)
    ggsave(here::here("plots/Fig1S16-2-frequency_plots.png"), p_body, width = 9, height = 12)
    ggsave(here::here("plots/Fig1S16-3-frequency_plots.png"), p_bottom, width = 10, height = 10)
    #pS16 <- plot_grid(p_top, p_bottom, nrow = 2, rel_heights = c(1.5,1)) + paint_white_background()
    #ggsave(here::here("plots/Fig1S16-frequency_plots.png"), pS16, width = 20, height = 20)

    ## Number of pairs using CASEU
    pairs_freq %>%
        filter(str_detect(Community, "C\\d+R\\d+")) %>%
        filter(Time == "T8") %>%
        #group_by(Community, Isolate1, Isolate2) %>%
        unite("Pair", Community, Isolate1, Isolate2, sep = "_") %>%
        select(Pair, Isolate1InitialODFreq, RawDataType) %>%
        pivot_wider(names_from = Isolate1InitialODFreq, values_from = RawDataType, names_prefix = "f") %>%
        mutate(Method = ifelse(f5 == "Sanger" | f50 == "Sanger" | f95 == "Sanger", "sanger",
                               ifelse(f5 == "CFU" & f50 == "CFU" & f95 == "CFU", "cfu", NA))) %>%
        group_by(Method) %>%
        summarize(Count = n())







    if (FALSE) {
        # Figure 1A: t1wo example networks.
        n_species = 4
        make_example_graph <- function(n_species, prob = c(1,1)) {
            # Probability = {"coexistence", "exclusion"}
            nodes <- tibble(Isolate = 1:n_species, x = 1:n_species, y = rep(0, n_species))
            edges <- tibble(From = 1:n_species, To = 1:n_species) %>%
                tidyr::expand(From, To) %>%
                filter(From < To) %>%
                mutate(InteractionType = sample(c("coexistence", "exclusion"), choose(n_species, 2), replace = T, prob = prob))
            edges_coexistence <- edges %>%
                filter(InteractionType == "coexistence") %>%
                mutate(temp = From, From = To, To = temp, .keep = "unused") %>% select(-temp)
            edges <- bind_rows(edges, edges_coexistence)
            example_graph <- tbl_graph(nodes, edges)
            return(example_graph)
        }
        plot_example_graph <- function(example_graph, node_size, edge_width) {
            example_graph %>%
                ggraph(layout = "nicely") +
                geom_edge_link(aes(color = InteractionType), edge_width = edge_width,
                               arrow = arrow(length = unit(edge_width, "mm"), type = "closed", angle = 30, ends = "last"),
                               start_cap = circle(node_size/2, "mm"), end_cap = circle(node_size/2, "mm")) +
                geom_node_point(fill = "grey", size = node_size, shape = 21, colour = "black", stroke = node_size/5) +
                scale_edge_color_manual(values = interaction_color) +
                scale_color_manual(values = c(`TRUE` = "#DB7469", `FALSE` = "#557BAA")) +
                scale_x_continuous(limits = c(-1.1, 1.1)) +
                scale_y_continuous(limits = c(-1.1, 1.1)) +
                guides(color = "none", fill = "none") +
                theme_void() +
                theme(legend.position = "none", plot.background = element_rect(fill = "white", color = NA),
                      legend.title = element_blank(), plot.margin = margin(10,30,10,10, unit = "pt"),
                      legend.text = element_text(size = 20)) +
                labs(x = "")
        }
        set.seed(1)
        p1 <- make_example_graph(n_species, c(1,0)) %>%
            activate(nodes) %>%
            mutate(x = c(-1,-1,1,1), y = c(-1,1,-1,1)) %>%
            plot_example_graph(node_size = 10, edge_width = 2) +
            # annotate("text", x = 0, y = 1, label = "Not Emergent", vjust = 0, size = 3) +
            scale_y_continuous(limits = c(-1, 1.3))
        p2 <- make_example_graph(n_species, c(1,1)) %>%
            activate(nodes) %>%
            mutate(x = c(-1,-1,1,1), y = c(-1,1,-1,1)) %>%
            plot_example_graph(node_size = 10, edge_width = 2) +
            scale_y_continuous(limits = c(-1, 1.3))
        p3 <- cowplot::get_legend(p2 + theme(legend.position = "bottom", legend.text = element_text(size = 10), legend.direction = "vertical"))
        pA <- plot_grid(p1, p2, p3, ncol = 1, scale = 1, rel_heights = c(2,2,1), axis = "tb", align = "v") + paint_white_background()
        #ggsave(here::here("plots/Fig1A-two_networks.png"), pA, width = 3, height = 5)
        ggsave(here::here("plots/cartoons/Fig1A-two_networks.pdf"), pA, width = 3, height = 5)

        ## Plot pairs example dynamics
        p_pairs_example_outcomes <- pairs_example_outcomes %>%
            filter(InteractionType != "neutrality") %>%
            left_join(pairs_freq, by = c("Community", "Isolate1", "Isolate2")) %>%
            mutate(Isolate1InitialODFreq = factor(Isolate1InitialODFreq)) %>%
            mutate(InteractionType = factor(InteractionType, c("exclusion", "coexistence"))) %>%
            ggplot(aes(x = Time, y = Isolate1MeasuredFreq, color = Isolate1InitialODFreq, group = Isolate1InitialODFreq)) +
            geom_point(size = 1) +
            geom_line(size = .5) +
            scale_y_continuous(breaks = c(0, .5, 1), limits = c(0,1)) +
            facet_grid(.~InteractionType) +
            theme_bw() +
            theme(panel.spacing = unit(2, "mm"), strip.text.x = element_blank(),
                  panel.border = element_rect(color = 1, fill = NA, size = 1),
                  axis.title = element_text(size = 10), axis.text = element_text(color = 1, size = 8)) +
            guides(color = "none") +
            labs(x = "Time", y = "Frequency")
    }















}

isolates %>%
    ggplot() +
    geom_histogram(aes(x = r_glucose_curver, fill = Fermenter)) +
    theme_classic()


isolates %>%
    ggplot() +
    geom_histogram(aes(x = r_acetate_curver, fill = Fermenter)) +
    theme_classic()


isolates %>%
    group_by(Fermenter) %>%
    summarize(M = mean(leakiness_16hr, na.rm = T), sd = sd(leakiness_16hr, na.rm = T))










