# Figure 3: model
library(tidyverse)
library(tidymodels)
library(tidygraph)
library(cowplot)
library(ggpubr)
library(ggraph)
source(here::here("plotting_scripts/network_functions.R"))

output_dir <- "~/Dropbox/lab/invasion-network/simulation/data/raw9/"
input_independent <- read_csv(paste0(output_dir,"input_independent.csv"), col_types = cols())
input_pairs <- read_csv(paste0(output_dir, "input_pairs.csv"), col_types = cols())
input_row <- input_independent[1,]

# Generate family-species and class-resource matching tibble
sa <- input_independent$sa[1]
ma <- input_independent$ma[1]
sal <- tibble(Family = paste0("F", c(rep(0, sa), rep(1, sa))), Species = paste0("S", 0:(sa * 2 - 1)))
mal <- tibble(Class = paste0("T", c(rep(0, ma), rep(1, ma), rep(2, ma))), Resource = paste0("R", 0:(ma * 3 - 1)))
paint_white_background <- function(x) theme(plot.background = element_rect(color = NA, fill = "white"))

#
df_communities <- read_csv(here::here("data/output/df_communities.csv"), col_types = cols())
df_communities_abundance <- read_csv(here::here("data/output/df_communities_abundance.csv"), col_types = cols()) %>% mutate(Community = factor(Community, df_communities$Community))
df_pairs <- read_csv(here::here("data/output/df_pairs.csv"), col_types = cols())

# Figure 3A. Dynamics of one example community
df_comm_example <- df_communities_abundance %>%
    filter(Community == "W1") %>%
    mutate(Time = ifelse(Time != 0 , (Transfer-1)*10+Time, Time)) # Note that this line changes when the growth cycle time changes
df_comm_example_end <- df_comm_example %>% filter(Time == max(Time)) %>%
    mutate(PlotColor = Family, PlotAlpha = factor(1:n())) %>%
    select(ID, PlotColor, PlotAlpha)

p1 <- df_comm_example %>%
    left_join(df_comm_example_end, by = "ID") %>%
    ggplot() +
    geom_line(aes(x = Time, y = Abundance, color = PlotColor, alpha = PlotAlpha, group = ID), size = 1) +
    scale_x_continuous(breaks = scales::pretty_breaks(n=2), expand = c(0,0)) +
    scale_y_continuous(breaks = scales::pretty_breaks(n=2), expand = c(.1,0)) +
    scale_color_manual(values = category_color, breaks = c("F0", "F1")) +
    scale_alpha_manual(values = c("1" = 1, "2" = 1, "3" = .8)) +
    theme_classic() +
    theme(panel.border = element_rect(color = 1, fill = NA)) +
    guides(alpha = "none", color = "none") +
    labs()

# Figure 3B. End point composition
p2 <- df_communities_abundance %>%
    filter(Transfer == max(Transfer), Time == max(Time)) %>%
    mutate(Community = factor(Community, c("W1", "W0", paste0("W", 2:20)))) %>%
    ggplot() +
    geom_col(aes(x = Community, y = Abundance, fill = Family, group = ID), color = 1, position = "fill") +
    scale_fill_manual(values = category_color, breaks = c("F0", "F1")) +
    scale_x_discrete(expand = c(0,0), label = 1:20) +
    scale_y_continuous(breaks = c(0, .5, 1), expand = c(0,0)) +
    #facet_grid(Time~.) +
    theme_classic() +
    theme(panel.border = element_rect(color = 1, fill = NA)) +
    labs()

p <- plot_grid(p1, p2, nrow = 1, scale = .9, axis = "tb", align = "h", rel_widths = c(1, 1.5)) + paint_white_background()
ggsave(here::here("plots/Fig3AB-community_composition.png"), p, width = 8, height = 3)




# Figure 3C: pairwise outcome of community pairs
## Overall
p1 <- df_pairs %>%
    filter(InteractionType != "no-growth") %>%
    mutate(InteractionType = factor(InteractionType, c("exclusion", "coexistence"))) %>%
    group_by(InteractionType) %>%
    summarize(Count = n()) %>%
    mutate(Fraction = Count / sum(Count), TotalCount = sum(Count)) %>%
    ggplot() +
    geom_col(aes(x = InteractionType, y = Fraction, fill = InteractionType), color = 1) +
    geom_text(aes(x = InteractionType, y = Fraction, label = paste0("n=",Count)), vjust = -1) +
    scale_fill_manual(values = assign_interaction_color()) +
    scale_y_continuous(limits = c(0, 1), breaks = scales::pretty_breaks(n=3)) +
    theme_classic() +
    theme(legend.position = "top", axis.title.x = element_blank()) +
    guides(fill = "none") +
    labs(y = "Count")


## Each community
p2 <- df_pairs %>%
    filter(InteractionType != "no-growth") %>%
    mutate(InteractionType = factor(InteractionType, c("exclusion", "coexistence"))) %>%
    group_by(Community, InteractionType) %>%
    summarize(Count = n()) %>%
    group_by(Community) %>%
    mutate(TotalCount = sum(Count)) %>%
    ggplot() +
    geom_col(aes(x = Community, y = Count, fill = InteractionType), position = "fill", color = 1) +
    geom_text(aes(x = Community, label = TotalCount), y = 1, vjust = 2) +
    scale_fill_manual(values = assign_interaction_color()) +
    scale_y_continuous(breaks = c(0, 0.5, 1)) +
    scale_x_discrete(label = 1:20) +
    theme_classic() +
    theme(legend.position = "top") +
    guides(fill = "none") +
    labs(y = "Fraction")
p <- plot_grid(p1, p2, nrow = 1, axis = "tb", align = "h", labels = c("A", "B"), scale = .9, rel_widths = c(1,1.5)) + paint_white_background()
ggsave(here::here("plots/Fig3C-community_pairs.png"), p, width = 8, height = 3)


# Figure 3D. Hierarchy
df_communities_hierarchy <- read_csv(here::here("data/output/df_communities_hierarchy.csv"), col_types = cols()) %>% mutate(Community = factor(Community, df_communities$Community))
df_communities_hierarchy_randomized <- read_csv(here::here("data/output/df_communities_hierarchy_randomized.csv"), col_types = cols()) %>% mutate(Community = factor(Community, df_communities$Community))

p1 <- df_communities_hierarchy %>%
    drop_na() %>%
    ggplot(aes(x = Metric, y = HierarchyScore)) +
    #geom_hline(yintercept = c(0, 0.5, 1), color = "grey", linetype = 2) +
    geom_boxplot(width = .5, lwd = .8) +
    geom_jitter(shape = 1, size = 2, width = .2, stroke = .8) +
    scale_y_continuous(limits = c(0,1.01), breaks = c(0, .25, .5, .75, 1)) +
    theme_classic() +
    theme(panel.grid.major.y = element_line(color = "grey", linetype = 2),
          axis.text = element_text(size = 10, color = 1),
          panel.border = element_rect(fill = NA, color = 1, size = 1.5)) +
    labs(x = "", y = "Hierarchy Score")

p2 <- df_communities_hierarchy_randomized %>%
    ggplot() +
    geom_histogram(aes(y = HierarchyScore), color = 1, fill = "white") +
    geom_hline(data = df_communities_hierarchy, aes(yintercept = HierarchyScore), color = "red") +
    facet_wrap(Metric~Community) +
    theme_classic()
pD <- p1
ggsave(here::here("plots/Fig3D-community_hierarhcy.png"), pD, width = 3, height = 3)



# Figure 3E. Motif
## Simulation result
df_motif <- read_csv(here::here("data/output/df_motif.csv"))
df_motif_randomized <- read_csv(here::here("data/output/df_motif_randomized.csv"))

## Overall
df_motif_aggregated <- df_motif %>%
    group_by(Motif) %>%
    summarize(Count = sum(Count))
df_motif_randomized_aggregated <- df_motif_randomized %>%
    group_by(Replicate, Motif) %>%
    summarize(Count = sum(Count))

df_motif_percentile_aggregated <- df_motif_randomized_aggregated %>%
    group_by(Motif) %>%
    arrange(desc(Count)) %>%
    slice(b * 0.05, b * 0.95) %>%
    mutate(Percentile = c("p5", "p95")) %>%
    select(Motif, Percentile, Count) %>%
    pivot_wider(names_from = Percentile, values_from = Count)

df_motif_aggregated <- df_motif_aggregated %>%
    left_join(df_motif_percentile_aggregated) %>%
    mutate(Sign = case_when(Count > p95 ~ "top",
                            Count < p5 ~ "bottom",
                            Count < p95 & Count > p5 ~ "n.s."))


## motif diagram
load(here::here("data/output/motif_list.Rdata"))
p_motif_list <- lapply(motif_list, function(x) plot_competitive_network(x, node_size = 3))
p1 <- plot_grid(plotlist = p_motif_list, nrow = 1)

## motif counts
p2 <- df_motif_randomized_aggregated %>%
    ggplot() +
    geom_rect(data = df_motif_aggregated, aes(fill = Sign), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, alpha = .2) +
    geom_vline(xintercept = 0, color = 1) +
    geom_histogram(aes(y = Count), binwidth = 2, color = 1, fill = "white") +
    geom_point(data = df_motif_aggregated, x = 0, aes(y = Count, color = "observed network"), pch = 1, size = 2, stroke = 2, inherit.aes = F) +
    scale_color_manual(values = c("observed network" = "red")) +
    facet_grid(.~Motif, scales = "free_x") +
    scale_fill_manual(values = c("top" = "blue", "bottom" = "red", "n.s." = "grey")) +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 2)) +
    theme_cowplot() +
    theme(panel.background = element_rect(color = 1, size = 1), panel.spacing = unit(0, "mm"),
          strip.background = element_rect(color = NA, fill = NA, size = 1)) +
    guides(color = "none", fill = "none") +
    labs(x = "Probability density", y = "Motif count")

pE <- plot_grid(p1, p2, ncol = 1, rel_heights = c(1,3), axis = "lr", align = "v", labels = LETTERS[1:2]) + theme(plot.background = element_rect(fill = "white", color = NA))
ggsave(here::here("plots/Fig3E-motif_counts.png"), pE, width = 10, height = 5)


## title
p0_1 <- ggdraw() + draw_label("Nontransitive", x = .5, hjust = 0.5) + theme(plot.margin = margin(0, 0, 5, 0))
p0_2 <- ggdraw() + draw_label("Hierarchical", x = .5, hjust = 0.5) + theme(plot.margin = margin(0, 0, 5, 0))
## motif diagram
load(here::here("data/output/motif_list.Rdata"))
p_motif_list <- lapply(motif_list, function(x) plot_competitive_network(x, node_size = 2, edge_width = 1))
p1_1 <- plot_grid(plotlist = p_motif_list[c(1)], nrow = 1, greedy = T)
p1_2 <- plot_grid(plotlist = p_motif_list[c(2,3,5)], nrow = 1, greedy = T)
## Motif count
comm <- "W0"
df_motif_subset <- df_motif %>% filter(Community == comm)
df_motif_randomized_subset <- df_motif_randomized %>% filter(Community == comm)
plot_motif_count <- function(motif_randomized_subset, motif_community_subset,  motif_subset = c(1,2,3)) {
    motif_randomized_subset %>%
        filter(Motif %in% motif_subset) %>%
        group_by(Motif, Replicate) %>%
        mutate(p5 = quantile(Count, 0.05), p95 = quantile(Count, 0.95)) %>%
        ggplot() +
        geom_vline(xintercept = 0, color = 1) +
        geom_hline(yintercept = 0, color = 1, linetype = 2) +
        geom_histogram(aes(y = Fraction, x = after_stat(count / max(count))), alpha = .3, color = 1) +
        geom_point(data = filter(motif_community_subset, Motif %in% motif_subset), aes(x = 0, y = Fraction, color = "observed network"), pch = 1, size = 2, stroke = 2, inherit.aes = F) +
        scale_fill_manual(values = c("head" = "#FF0000A0", "body" = "#A0A0A0A0", "tail" = "#FF0000A0"),
                          labels = c("head" = "top 5%", "body" = "middle", "tail" = "bottom 5%")) +
        scale_color_manual(values = c("observed network" = "red")) +
        facet_grid(.~Motif) +
        scale_x_continuous(breaks = c(0,0.5), labels = c("0", "0.5")) +
        theme_cowplot() +
        theme(panel.background = element_rect(color = 1, size = 1.5, fill = NA),
              panel.spacing = unit(0, "mm"), strip.text = element_blank(),
              legend.position = "top",
              axis.title = element_text(size = 10), axis.text = element_text(size = 8),
              plot.background = element_rect(fill = "white", color = NA)) +
        guides(fill = "none", color = "none") +
        labs(x = "Probability density", y = "Fraction of motif")
}
p2_1 <- plot_motif_count(df_motif_randomized_subset, df_motif_subset, c(1)) + theme(axis.title.x = element_blank())
p2_2 <- plot_motif_count(df_motif_randomized_subset, df_motif_subset, c(2,3,5)) + theme(axis.title = element_blank())
p_left <- plot_grid(p0_1, p1_1, p2_1, ncol = 1, rel_heights = c(.4, 1,2.5), axis = "tblr", align = "v")
p_right <- plot_grid(p0_2, p1_2, p2_2, ncol = 1, rel_heights = c(.4, 1,2.5), axis = "tblr", align = "v")
p_upper <- plot_grid(p_left, p_right, nrow = 1, rel_widths = c(1.5,3), axis = "tblr", align = "h") + paint_white_background()
pF_xlab <- ggdraw() + draw_label("Probability density", x = 0.5, hjust = .5) + theme(plot.margin = margin(0, 0, 0, 0))
pF <- plot_grid(p_upper, pF_xlab, ncol = 1, rel_heights = c(1, .1)) + paint_white_background()
ggsave(here::here("plots/Fig3F-motifs.png"), pF, width = 5, height = 4)





# Figure 3G. Diagonal
df_diag <- read_csv(here::here("data/output/df_diag.csv"))
df_diag_randomized <- read_csv(here::here("data/output/df_diag_randomized.csv"))

df_diag_randomized %>%
    group_by(DistanceToDiagonal, Replicate) %>%
    summarize(Count = sum(CountCoexistence)) %>%
    ggplot(aes(x = DistanceToDiagonal, y = Count, group = DistanceToDiagonal)) +
    geom_boxplot() +
    geom_jitter(height = 0, width = 0.1) +
    geom_point(data = df_diag %>% group_by(DistanceToDiagonal) %>% summarize(Count = sum(CountCoexistence)),
               aes(x = DistanceToDiagonal, y = Count, group = DistanceToDiagonal), color = "red") +
    theme_classic()





#=========================================================================================================

# Figure 3S1: u, l, and D matrices
Dm <- read_csv(paste0(output_dir, "D_seed1.csv"), skip = 1) # D matrix
cm <- read_csv(paste0(output_dir, "c_seed1.csv"), skip = 1) # c matrix
lm <- read_csv(paste0(output_dir, "l_seed1.csv"), skip = 1) # l matrix

## D matrix
Dml <- Dm %>% # D matrix longer
    pivot_longer(cols = starts_with("R"), names_to = "Resource1", values_to = "SecretionFlux") %>%
    rename(Class2 = ...1, Resource2 = ...2) %>%
    left_join(rename_with(mal, ~ paste0(., "1"), everything())) %>%
    mutate(Resource1 = ordered(Resource1, mal$Resource), Resource2 = ordered(Resource2, mal$Resource)) %>%
    select(Class1, Resource1, Class2, Resource2, SecretionFlux)

Dml %>%
    group_by(Class1) %>%
    summarize(TotalFlux = sum(SecretionFlux))

p1 <- Dml %>%
    ggplot() +
    geom_tile(aes(x = Resource1, y = Resource2, fill = SecretionFlux)) +
    geom_segment(aes(color = "sugar"), x = "R0", xend = paste0("R", ma-1), y = "spacer", yend = "spacer", lwd = 2) +
    geom_segment(aes(color = "acid"), x = paste0("R", ma), xend = paste0("R", 2*ma-1), y = "spacer", yend = "spacer", lwd = 2) +
    geom_segment(aes(color = "waste"), x = paste0("R", 2*ma), xend = paste0("R", 3*ma-1), y = "spacer", yend = "spacer", lwd = 2) +
    geom_segment(aes(color = "sugar"), x = "spacer", xend = "spacer", y = "R0", yend = paste0("R", ma-1), lwd = 2) +
    geom_segment(aes(color = "acid"), x = "spacer", xend = "spacer", y = paste0("R", ma), yend = paste0("R", 2*ma-1), lwd = 2) +
    geom_segment(aes(color = "waste"), x = "spacer", xend = "spacer", y = paste0("R", 2*ma), yend = paste0("R", 3*ma-1), lwd = 2) +
    scale_fill_gradient(high = "black", low = "white") +
    scale_color_manual(values = category_color, breaks = c("sugar", "acid", "waste")) +
    scale_x_discrete(limits = c("spacer", mal$Resource)) +
    scale_y_discrete(limits = c("spacer", mal$Resource)) +
    theme_classic() +
    theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank()) +
    guides(fill = guide_colorbar(title = expression(D[i][beta][alpha])),
           color = "none") +
    labs(x = expression(R[alpha]), y = expression(R[beta]))

## c matrix
cml <- cm %>% # c matrix longer
    pivot_longer(cols = starts_with("R"), names_to = "Resource", values_to = "ConsumptionRate") %>%
    rename(Family = ...1, Species = ...2) %>%
    left_join(mal) %>%
    mutate(Species = ordered(Species, sal$Species), Resource = ordered(Resource, mal$Resource)) %>%
    select(Family, Species, Class, Resource, ConsumptionRate)
p2 <- cml %>%
    ggplot() +
    geom_tile(aes(x = Resource, y = Species, fill = ConsumptionRate)) +
    geom_segment(aes(color = "sugar"), x = "R0", xend = paste0("R", ma-1), y = "spacer", yend = "spacer", lwd = 2) +
    geom_segment(aes(color = "acid"), x = paste0("R", ma), xend = paste0("R", 2*ma-1), y = "spacer", yend = "spacer", lwd = 2) +
    geom_segment(aes(color = "waste"), x = paste0("R", 2*ma), xend = paste0("R", 3*ma-1), y = "spacer", yend = "spacer", lwd = 2) +
    geom_segment(aes(color = "fermenter"), x = "spacer", xend = "spacer", y = "S0", yend = paste0("S", sa-1), lwd = 2) +
    geom_segment(aes(color = "respirator"), x = "spacer", xend = "spacer", y = paste0("S", sa), yend = paste0("S", 2*sa-1), lwd = 2) +
    scale_fill_gradient(high = "black", low = "white") +
    scale_color_manual(values = category_color) +
    scale_x_discrete(limits = c("spacer", mal$Resource)) +
    scale_y_discrete(limits = c("spacer", sal$Species)) +
    theme_classic() +
    theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank()) +
    guides(fill = guide_colorbar(title = expression(c[i][alpha])),
           color = "none") +
    labs()

cml %>%
    ggplot() +
    geom_histogram(aes(x = ConsumptionRate, fill = Family)) +
    scale_fill_manual(values = category_color) +
    scale_y_log10() +
    theme_classic() +
    #guides(fill = "none") +
    labs()

## l matrix
lml <- lm %>% # l matrix longer
    pivot_longer(cols = starts_with("R"), names_to = "Resource", values_to = "Leakiness") %>%
    rename(Family = ...1, Species = ...2) %>%
    left_join(mal) %>%
    mutate(Species = ordered(Species, sal$Species), Resource = ordered(Resource, mal$Resource)) %>%
    select(Family, Species, Class, Resource, Leakiness)
p3 <- lml %>%
    ggplot() +
    geom_tile(aes(x = Resource, y = Species, fill = Leakiness)) +
    geom_segment(aes(color = "sugar"), x = "R0", xend = paste0("R", ma-1), y = "spacer", yend = "spacer", lwd = 2) +
    geom_segment(aes(color = "acid"), x = paste0("R", ma), xend = paste0("R", 2*ma-1), y = "spacer", yend = "spacer", lwd = 2) +
    geom_segment(aes(color = "waste"), x = paste0("R", 2*ma), xend = paste0("R", 3*ma-1), y = "spacer", yend = "spacer", lwd = 2) +
    geom_segment(aes(color = "fermenter"), x = "spacer", xend = "spacer", y = "S0", yend = paste0("S", sa-1), lwd = 2) +
    geom_segment(aes(color = "respirator"), x = "spacer", xend = "spacer", y = paste0("S", sa), yend = paste0("S", 2*sa-1), lwd = 2) +
    scale_fill_gradient(high = "black", low = "white") +
    scale_color_manual(values = category_color) +
    scale_x_discrete(limits = c("spacer", mal$Resource)) +
    scale_y_discrete(limits = c("spacer", sal$Species)) +
    theme_classic() +
    theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank()) +
    guides(fill = guide_colorbar(title = expression(l[i][alpha])),
           color = guide_legend(title = "")) +
    labs()

lml %>%
    filter(Class != "T2") %>%
    ggplot() +
    geom_histogram(aes(x = Leakiness, fill = Family)) +
    scale_fill_manual(values = category_color) +
    theme_classic() +
    guides(fill = "none") +
    labs()

p_lower <- plot_grid(p2, p3, labels = c("c matrix", "l matrix"), scale = .9)
pS1 <- plot_grid(p1, p_lower, ncol = 1, labels = c("D matrix", ""), align = "v", axis = "lr", scale = c(0.9, 1)) + paint_white_background()
ggsave(here::here("plots/Fig3S1-matrices.png"), pS1, width = 12, height = 12)






















# =========================================================================================================

# Figure 3A: diagram cartoon
#p_A <- ggdraw() + draw_image(here::here("plots/cartoons/Fig3A.png")) + theme(plot.background = element_rect(fill = "white", color = NA))
p_A <-  ggplot(mtcars, aes(x = wt, y = mpg)) + annotate("text", x = 0 , y = 0, label = "Cartoon for\nmodel") + theme_void() + theme(plot.background = element_rect(fill = "white", color = NA))
ggsave(here::here("plots/Fig3A-functional_groups.png"), p_A, width = 5, height = 5)



# Figure 3C: community composition
## Initial
df_comm_init <- input_independent %>%
    filter(str_detect(init_N0, "selfAssembly")) %>%
    pull(init_N0) %>%
    paste0(output_dir, .) %>%
    read_wide_file() %>%
    full_join(sal, by = c("Family", "Species")) %>%
    replace_na(list(Abundance = 0)) %>%
    mutate(Community = factor(Well, paste0("W", 0:1000)), .keep = "unused") %>%
    mutate(Species = factor(Species, sal$Species)) %>%
    # Remove rare species (relative abundance <0.01)
    group_by(Community) %>%
    mutate(RelativeAbundance = Abundance/sum(Abundance)) %>%
    #filter(RelativeAbundance > 0.01) %>%
    arrange(Community, Species)

## End
df_comm_end <- input_independent %>%
    filter(str_detect(init_N0, "selfAssembly")) %>%
    pull(init_N0) %>% str_replace("_init.csv", "_end.csv") %>%
    paste0(output_dir, .) %>%
    read_wide_file() %>%
    full_join(sal, by = c("Family", "Species")) %>%
    replace_na(list(Abundance = 0)) %>%
    mutate(Community = factor(Well, paste0("W", 0:1000)), .keep = "unused") %>%
    mutate(Species = factor(Species, sal$Species)) %>%
    # Remove rare species (relative abundance <0.01)
    group_by(Community) %>%
    mutate(RelativeAbundance = Abundance/sum(Abundance)) %>%
    filter(RelativeAbundance > 0.01) %>%
    arrange(Community, Species)

##
p_C <- bind_rows(df_comm_init %>% mutate(Time = "init"),
                 df_comm_end %>% mutate(Time = "end")) %>%
    mutate(Time = factor(Time, c("init", "end"))) %>%
    ggplot() +
    geom_col(aes(x = Community, y = Abundance, fill = Family, group = Species), color = 1, position = "fill") +
    scale_fill_manual(values = category_color, breaks = c("F0", "F1")) +
    scale_y_continuous(breaks = c(0, .5, 1)) +
    facet_grid(Time~.) +
    theme_classic() +
    labs()

ggsave(here::here("plots/Fig3C-community_composition.png"), p_C, width = 7, height = 6)


# Figure 3D: Monoculture growth
df_mono <- input_independent %>% filter(str_detect(init_N0, "monoculture")) %>%
    pull(init_N0) %>% str_replace("_init.csv", "_end.csv") %>%
    paste0(output_dir, .) %>%
    read_wide_file() %>%
    full_join(sal, by = c("Family", "Species")) %>%
    replace_na(list(Abundance = 0))

p_D <- df_mono %>%
    mutate(Growth = ifelse(Abundance != 0, "culturable", "unculturable") %>% factor(c("unculturable", "culturable"))) %>%
    group_by(Family, Growth) %>%
    ggplot() +
    geom_bar(aes(x = Family, fill = Family, alpha = Growth), color = 1) +
    scale_fill_manual(values = category_color, breaks = c("F0", "F1")) +
    scale_alpha_manual(values = c("culturable" = 1, "unculturable" = .1)) +
    theme_classic() +
    guides(fill = guide_legend(title = ""),
           alpha = guide_legend(title = "")) +
    labs()

ggsave(here::here("plots/Fig3D-monoculture.png"), p_D, width = 3, height = 3)

# Figure 3E: pairwise outcome of pool pairs
df_pp_init <- input_pairs %>%
    filter(str_detect(init_N0, "poolPairs")) %>%
    pull(init_N0) %>%
    paste0(output_dir, .) %>%
    lapply(function(x) {
        read_wide_file(x) %>%
            #full_join(sal, by = c("Family", "Species")) %>%
            #replace_na(list(Abundance = 0)) %>%
            mutate(Well = factor(Well, paste0("W", 0:1000)), .keep = "unused") %>%
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
    mutate(Time = "Tinit")

df_pp_end <- input_pairs %>%
    filter(str_detect(init_N0, "poolPairs")) %>%
    pull(init_N0) %>% str_replace("_init.csv", "_end.csv") %>%
    paste0(output_dir, .) %>%
    lapply(function(x) {
        read_wide_file(x) %>%
            #full_join(sal, by = c("Family", "Species")) %>%
            #replace_na(list(Abundance = 0)) %>%
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
    mutate(Time = "Tend")

## Determine outcome
pairs_init <- df_pp_init %>% filter(Community == "W0")
pairs_end <- df_pp_end %>% filter(Community == "W0")
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


#
pairs_pool <- determine_interaction(df_pp_init, df_pp_end) %>%
    left_join(rename_with(sal, ~paste0(., 1))) %>%
    left_join(rename_with(sal, ~paste0(., 2))) %>%
    mutate(PairConspecific = with(., case_when(
        (Family1 == Family2) ~ "conspecific",
        (Family1 != Family2) ~ "heterospecific"
    ))) %>%
    mutate(Community = factor(Community, paste0("W", 0:1000))) %>%
    arrange(Community)
temp <- pairs_pool %>%
    filter(InteractionType != "no-growth") %>%
    mutate(InteractionType = factor(InteractionType, c("exclusion", "coexistence"))) %>%
    group_by(Community, InteractionType, PairConspecific) %>%
    summarize(Count = n())

## Overall
p1 <- pairs_pool %>%
    filter(InteractionType != "no-growth") %>%
    mutate(InteractionType = factor(InteractionType, c("exclusion", "coexistence"))) %>%
    group_by(InteractionType) %>%
    summarize(Count = n()) %>%
    mutate(Fraction = Count / sum(Count), TotalCount = sum(Count)) %>%
    ggplot() +
    geom_col(aes(x = InteractionType, y = Fraction, fill = InteractionType), color = 1) +
    geom_text(aes(x = InteractionType, y = Fraction, label = paste0("n=",Count)), vjust = -1) +
    scale_fill_manual(values = assign_interaction_color()) +
    scale_y_continuous(limits = c(0, 1), breaks = scales::pretty_breaks(n=3)) +
    theme_classic() +
    theme(legend.position = "top", axis.title.x = element_blank()) +
    guides(fill = "none") +
    labs(y = "Count")
p1
## Conspecific vs. heterospecific
p2 <- pairs_pool %>%
    filter(InteractionType != "no-growth") %>%
    mutate(InteractionType = factor(InteractionType, c("exclusion", "coexistence"))) %>%
    group_by(InteractionType, PairConspecific) %>%
    summarize(Count = n()) %>%
    group_by(PairConspecific) %>%
    mutate(TotalCount = sum(Count)) %>%
    ggplot() +
    geom_col(aes(x = PairConspecific, y = Count, fill = InteractionType), color = 1, position = "fill") +
    #geom_text(aes(x = PairConspecific, y = Count, label = round(Percentage,2)), vjust = -1) +
    geom_text(aes(x = PairConspecific, label = paste0("n=", TotalCount)), y = 1, vjust = 2) +
    scale_fill_manual(values = assign_interaction_color()) +
    #scale_y_continuous(limits = c(0, 500), breaks = scales::pretty_breaks(n=3)) +
    theme_classic() +
    theme(legend.position = "top", axis.title.x = element_blank()) +
    guides(fill = "none") +
    labs(y = "Coount")


## Each community
p3 <- pairs_pool %>%
    filter(InteractionType != "no-growth") %>%
    mutate(InteractionType = factor(InteractionType, c("exclusion", "coexistence"))) %>%
    group_by(Community, InteractionType) %>%
    summarize(Count = n()) %>%
    group_by(Community) %>%
    mutate(TotalCount = sum(Count)) %>%
    ggplot() +
    geom_col(aes(x = Community, y = Count, fill = InteractionType), position = "fill") +
    geom_text(aes(x = Community, label = TotalCount), y = 1, vjust = 2) +
    scale_fill_manual(values = assign_interaction_color()) +
    scale_y_continuous(breaks = c(0, 0.5, 1)) +
    theme_classic() +
    theme(legend.position = "top") +
    guides(fill = "none") +
    labs(y = "Fraction")
p_upper <- plot_grid(p1, p2, nrow = 1, axis = "tb", align = "h", labels = c("A", "B"), scale = .9)
p_E <- plot_grid(p_upper, p3, nrow = 2, scale = c(1, .9), labels = c("", "C")) + paint_white_background()
ggsave(here::here("plots/Fig3E-pool_pairs.png"), p_E, width = 8, height = 8)


# Figure 3F: u and l versus coexistence in pool pairs
## Subset Dml. Only use the R0 secretion
Dmls <- Dml %>%
    filter(Resource1 == "R0") %>%
    select(Resource = Resource2, SecretionFlux) # secretion from R0
## Subset cml. Only use c on R0
cmls <- cml %>% filter(Resource == "R0")
## Subset lml, Only use the R0 secretion
lmls <- lml %>%
    left_join(Dmls) %>%
    mutate(CrossFeedingPotential = Leakiness * SecretionFlux, .keep = "unused") %>%
    group_by(Species) %>%
    summarize(CrossFeedingPotential = sum(CrossFeedingPotential)) %>% # Sum of R0-to-other flux times leakiness
    select(Species, CrossFeedingPotential)


pairs_pool_meta <- pairs_pool %>%
    left_join(cmls %>% rename_with(~paste0(., 1), everything())) %>%
    left_join(cmls %>% rename_with(~paste0(., 2), everything())) %>%
    left_join(lmls %>% rename_with(~paste0(., 1), everything())) %>%
    left_join(lmls %>% rename_with(~paste0(., 2), everything())) %>%
    mutate(d_ConsumptionRate = ConsumptionRate1 - ConsumptionRate2,
           d_CrossFeedingPotential = CrossFeedingPotential1 - CrossFeedingPotential2,
           .keep = "unused")

## boxplot: r_glu
p1 <- pairs_pool_meta %>%
    #filter(PairConspecific == "conspecific") %>%
    mutate(InteractionType = factor(InteractionType, c("exclusion", "coexistence"))) %>%
    ggplot(aes(x = InteractionType, y = d_ConsumptionRate, fill = InteractionType)) +
    geom_hline(yintercept = 0, linetype = 2, color = "grey") +
    geom_boxplot(width = .5) +
    geom_point(aes(group = InteractionType), shape = 1, size = 1, position = position_jitterdodge(jitter.width = 0.3)) +
    # p value
    ggpubr::stat_compare_means(aes(group = InteractionType), label = "p.format", vjust = 2, hjust = 0, method = "t.test") +
    scale_fill_manual(values = assign_interaction_color()) +
    facet_grid(.~PairConspecific) +
    theme_classic() +
    theme(legend.title = element_blank(), legend.position = "top",
          axis.text = element_text(size = 12, color = 1),
          axis.title = element_text(size = 15, color = 1),
          axis.text.x = element_text(size = 12, color = 1, angle = 30, vjust = 1, hjust = 1),
          panel.border = element_rect(color = 1, fill = NA, size = 1),
          strip.background = element_rect(color = NA, fill = NA)) +
    guides(alpha = "none", fill = "none", color = "none") +
    labs(x = "", y = expression(r[A]-r[B]))

## boxplot: cross-feeding potential
p2 <- pairs_pool_meta %>%
    #filter(PairConspecific == "conspecific") %>%
    mutate(InteractionType = factor(InteractionType, c("exclusion", "coexistence"))) %>%
    ggplot(aes(x = InteractionType, y = d_CrossFeedingPotential, fill = InteractionType)) +
    geom_hline(yintercept = 0, linetype = 2, color = "grey") +
    geom_boxplot(width = 0.5) +
    geom_point(aes(group = InteractionType), shape = 1, size = 1, position = position_jitterdodge(jitter.width = 0.3)) +
    # p value
    ggpubr::stat_compare_means(aes(group = InteractionType), label = "p.format", vjust = 2, hjust = 0, method = "t.test") +
    scale_fill_manual(values = assign_interaction_color()) +
    facet_grid(.~PairConspecific) +
    theme_classic() +
    theme(legend.title = element_blank(), legend.position = "top",
          axis.text = element_text(size = 12, color = 1),
          axis.title = element_text(size = 15, color = 1),
          axis.text.x = element_text(size = 12, color = 1, angle = 30, vjust = 1, hjust = 1),
          panel.border = element_rect(color = 1, fill = NA, size = 1),
          strip.background = element_rect(color = NA, fill = NA)) +
    guides(alpha = "none", fill = "none", color = "none") +
    labs(x = "", y = expression(X[A]-X[B]))

## Scatterplot
p3 <- pairs_pool_meta %>%
    #filter(InteractionType != "no-growth") %>%
    ggplot() +
    geom_vline(xintercept = 0, linetype = 2) +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_point(aes(x = d_ConsumptionRate, y = d_CrossFeedingPotential, color = InteractionType), shape = 21, size = 2, stroke = .5) +
    scale_color_manual(values = c(assign_interaction_color())) +
    facet_grid(.~PairConspecific) +
    theme_classic() +
    theme(legend.position = "top", strip.background = element_blank(), panel.background = element_rect(color = 1)) +
    guides(color = "none") +
    labs()

p_F <- plot_grid(p1, p2, p3, ncol = 1, axis = "rl", align = "v", labels = LETTERS[1:3], scale = .9) + paint_white_background()
ggsave(here::here("plots/Fig3F-pairs_pool_trait.png"), p_F, width = 6, height = 10)

## Stat
### Two sample
pairs_pool_meta %>%
    filter(PairConspecific == "conspecific") %>%
    t_test(d_CrossFeedingPotential ~ InteractionType, order = c("coexistence", "exclusion"))
### glm
pairs_pool_meta %>%
    mutate_if(is.character, as.factor) %>%
    filter(PairConspecific == "conspecific") %>%
    #filter(PairConspecific != "conspecific") %>%
    #filter(!is.na(InteractionType)) %>%
    mutate(InteractionType = ifelse(InteractionType == "coexistence", 1, 0)) %>%
    glm(formula = InteractionType ~  d_ConsumptionRate * d_CrossFeedingPotential, data = ., family = "binomial") %>%
    broom::tidy() %>%
    {.}


# Figure 3G: pool networks
## For pairs from a network
calculate_rank <- function(pairs) {
    pairs %>%
        mutate(Point1 = with(., case_when(
            (InteractionType == "coexistence") ~ 0,
            (InteractionType == "exclusion") ~ 1,
        )),
        Point2 = with(., case_when(
            (InteractionType == "coexistence") ~ 0,
            (InteractionType == "exclusion") ~ -1,
        ))) %>%
        select(Community, starts_with("Species"), starts_with("Point")) %>%
        pivot_longer(cols = starts_with("Species"), names_to = "temp", values_to = "Species") %>%
        mutate(Point = ifelse(temp == "Species1", Point1, Point2)) %>%
        # Competitive score
        group_by(Species) %>%
        summarize(Score = sum(Point, na.rm = T)) %>%
        arrange(desc(Score)) %>%
        # Rank
        mutate(Rank = 1:n(), PlotRank = Rank) %>%
        rename(ID = Species) %>%
        mutate(ID = factor(ID, paste0("S", 0:10000))) %>%
        arrange(ID) %>%
        mutate(Isolate = as.character(1:n())) %>%
        ungroup()
}
pn_names <- pairs_pool %>% pull(Community) %>% unique # pool network names
pn_list <- pairs_pool %>%
    group_by(Community) %>%
    group_split() %>%
    lapply(function(x) {
        tt <- unique(x$Community)
        isolates <- calculate_rank(x)
        pairs <- x %>%
            mutate(across(starts_with("Species"), factor)) %>%
            rename(ID1 = Species1, ID2 = Species2) %>%
            left_join(isolates %>% select(ID, Isolate) %>% rename(ID1 = ID, Isolate1 = Isolate), by = "ID1") %>%
            left_join(isolates %>% select(ID, Isolate) %>% rename(ID2 = ID, Isolate2 = Isolate), by = "ID2") %>%
            mutate(From = Isolate1, To = Isolate2)
        make_network(isolates, pairs)
    }) %>%
    set_names(pn_names)

## Motif
pn_motif <- pn_list %>%
    lapply(function(x) tibble(Motif = 1:7, Count = count_motif(x))) %>%
    bind_rows() %>%
    mutate(Community = rep(pn_names, each = 7))
pn_motif %>%
    group_by(Community) %>%
    summarize(sum(Count))

p_G <- pn_motif %>%
    ggplot() +
    geom_point(aes(x = Motif, y = Count, color = Community), shape = 21, size = 2) +
    geom_line(aes(x = Motif, y = Count, color = Community)) +
    scale_x_continuous(breaks = 1:7) +
    theme_classic() +
    theme() +
    guides(color = "none") +
    labs()
ggsave(here::here("plots/Fig3G-pool_motif.png"), p_G, width = 6, height = 4)


# Figure 3H: individual pool network
## Network of pool networks
## Matrix and graph
p_pn_matrix_list <- lapply(pn_list, function(x) plot_adjacent_matrix(x) + theme(plot.margin = grid::unit(c(5,0,3,0), "mm")))
p_pn_list <- lapply(pn_list, function(x) plot_competitive_network(x, node_size = 2) + theme(plot.background = element_rect(fill = NA)))
p_list <- rep(list(NA), length(p_pn_list))
for (i in 1:length(pn_list)) p_list[[i]] <- ggdraw(p_pn_matrix_list[[i]]) + draw_plot(plot = p_pn_list[[i]], x = -.1, y = -.1, width = 0.7, height = 0.7)
## Motif count
plot_motif_count <- function (x = 1, data = pn_motif, network_names = pn_names) {
    # motif_randomized_subset <- networks_motif_randomized_percentile %>%
    #     filter(Community %in% network_names[x]) %>%
    #     mutate(Community = factor(Community, network_names[x]))
    motif_community_subset <- data %>%
        filter(Community %in% network_names[x]) %>%
        mutate(Community = factor(Community, network_names[x]))

    ggplot() +
        # 5% and 95% percentiles in randomized networks
        # geom_point(data = motif_randomized_subset, aes(x = Motif, y = Count, group = Motif, color = "randomized network")) +
        # geom_segment(data = motif_randomized_subset %>% pivot_wider(id_cols = c(Community, Motif), names_from = Percentile, values_from = Count),
        #              aes(x = Motif, xend = Motif, y = p5, yend = p95, color = "randomized network")) +
        # Observations
        geom_point(data = motif_community_subset, aes(x = Motif, y = Count, color = "observed network")) +
        scale_x_continuous(breaks = 1:7) +
        scale_color_manual(values = c("observed network" = "red", "randomized network" = "black"))+
        #facet_wrap(Community ~., scale = "free_y", nrow = 1)  +
        theme_classic() +
        theme(panel.background = element_rect(color = 1, size = 1), legend.position = "none")
}
p_pn_motif_count_list <- rep(list(NA), length(p_pn_list))
for (i in 1:length(p_pn_list)) {
    if (i %in% c(1, 6, 11)) p_pn_motif_count_list[[i]] <- plot_motif_count(i, pn_motif, pn_names) + theme(axis.title.x = element_blank())
    if (i %in% c(2:5, 7:10, 12:15)) p_pn_motif_count_list[[i]] <- plot_motif_count(i, pn_motif, pn_names) + theme(axis.title = element_blank())
    if (i == 16) p_pn_motif_count_list[[i]] <- plot_motif_count(i, pn_motif, pn_names)
    if (i %in% 17:20) p_pn_motif_count_list[[i]] <- plot_motif_count(i, pn_motif, pn_names) + theme(axis.title.y = element_blank())
}

## Get legend for line
p_temp <- plot_motif_count(1) + theme(legend.position = "right", legend.title = element_blank(), legend.text = element_text(size = 15))
shared_legend_line <- cowplot::get_legend(p_temp)
p_H <- list(p_list[1:5], p_pn_motif_count_list[1:5],
            p_list[6:10], p_pn_motif_count_list[6:10],
            p_list[11:15], p_pn_motif_count_list[11:15],
            p_list[16:20], p_pn_motif_count_list[16:20]) %>%
    unlist(recursive = F) %>%
    plot_grid(plotlist = ., labels = c(pn_names[1:5], rep("", 5), pn_names[6:10], rep("", 5), pn_names[11:15], rep("", 5), pn_names[16:20], rep("", 5)),
              ncol = 5, axis = "tbrl", align = "v") + paint_white_background()
ggsave(here::here("plots/Fig3H-pool_matrix.png"), p_H, width = 10, height = 12)




# Figure 3I: pairwise outcome of community pairs
df_cp_init <- input_pairs %>%
    filter(str_detect(init_N0, "communityPairs")) %>%
    pull(init_N0) %>%
    paste0(output_dir, .) %>%
    lapply(function(x) {
        read_wide_file(x) %>%
            #full_join(sal, by = c("Family", "Species")) %>%
            #replace_na(list(Abundance = 0)) %>%
            mutate(Well = factor(Well, paste0("W", 0:1000)), .keep = "unused") %>%
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
    paste0(output_dir, .) %>%
    lapply(function(x) {
        read_wide_file(x) %>%
            #full_join(sal, by = c("Family", "Species")) %>%
            #replace_na(list(Abundance = 0)) %>%
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

## Determine outcome
pairs_init <- df_cp_init %>% filter(Community == "W0")
pairs_end <- df_cp_end %>% filter(Community == "W0")
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


#
pairs_comm <- determine_interaction(df_cp_init, df_cp_end) %>%
    left_join(rename_with(sal, ~paste0(., 1))) %>%
    left_join(rename_with(sal, ~paste0(., 2))) %>%
    mutate(PairConspecific = with(., case_when(
        (Family1 == Family2) ~ "conspecific",
        (Family1 != Family2) ~ "heterospecific"
    ))) %>%
    mutate(Community = factor(Community, paste0("W", 0:1000))) %>%
    arrange(Community)
temp <- pairs_comm %>%
    filter(InteractionType != "no-growth") %>%
    mutate(InteractionType = factor(InteractionType, c("exclusion", "coexistence"))) %>%
    group_by(Community, InteractionType, PairConspecific) %>%
    summarize(Count = n())

## Overall
p1 <- pairs_comm %>%
    filter(InteractionType != "no-growth") %>%
    mutate(InteractionType = factor(InteractionType, c("exclusion", "coexistence"))) %>%
    group_by(InteractionType) %>%
    summarize(Count = n()) %>%
    mutate(Fraction = Count / sum(Count), TotalCount = sum(Count)) %>%
    ggplot() +
    geom_col(aes(x = InteractionType, y = Fraction, fill = InteractionType), color = 1) +
    geom_text(aes(x = InteractionType, y = Fraction, label = paste0("n=",Count)), vjust = -1) +
    scale_fill_manual(values = assign_interaction_color()) +
    scale_y_continuous(limits = c(0, 1), breaks = scales::pretty_breaks(n=3)) +
    theme_classic() +
    theme(legend.position = "top", axis.title.x = element_blank()) +
    guides(fill = "none") +
    labs(y = "Count")

## Conspecific vs. heterospecific
p2 <- pairs_comm %>%
    filter(InteractionType != "no-growth") %>%
    mutate(InteractionType = factor(InteractionType, c("exclusion", "coexistence"))) %>%
    group_by(InteractionType, PairConspecific) %>%
    summarize(Count = n()) %>%
    group_by(PairConspecific) %>%
    mutate(TotalCount = sum(Count)) %>%
    ggplot() +
    geom_col(aes(x = PairConspecific, y = Count, fill = InteractionType), color = 1, position = "fill") +
    #geom_text(aes(x = PairConspecific, y = Count, label = round(Percentage,2)), vjust = -1) +
    geom_text(aes(x = PairConspecific, label = paste0("n=", TotalCount)), y = 1, vjust = 2) +
    scale_fill_manual(values = assign_interaction_color()) +
    #scale_y_continuous(limits = c(0, 500), breaks = scales::pretty_breaks(n=3)) +
    theme_classic() +
    theme(legend.position = "top", axis.title.x = element_blank()) +
    guides(fill = "none") +
    labs(y = "Coount")


## Each community
p3 <- pairs_comm %>%
    filter(InteractionType != "no-growth") %>%
    mutate(InteractionType = factor(InteractionType, c("exclusion", "coexistence"))) %>%
    group_by(Community, InteractionType) %>%
    summarize(Count = n()) %>%
    group_by(Community) %>%
    mutate(TotalCount = sum(Count)) %>%
    ggplot() +
    geom_col(aes(x = Community, y = Count, fill = InteractionType), position = "fill") +
    geom_text(aes(x = Community, label = TotalCount), y = 1, vjust = 2) +
    scale_fill_manual(values = assign_interaction_color()) +
    scale_y_continuous(breaks = c(0, 0.5, 1)) +
    theme_classic() +
    theme(legend.position = "top") +
    guides(fill = "none") +
    labs(y = "Fraction")
p_upper <- plot_grid(p1, p2, nrow = 1, axis = "tb", align = "h", labels = c("A", "B"), scale = .9)
p_I <- plot_grid(p_upper, p3, nrow = 2, scale = c(1, .9), labels = c("", "C")) + paint_white_background()
ggsave(here::here("plots/Fig3I-community_pairs.png"), p_I, width = 8, height = 8)




# Figure 3J: u and l versus coexistence in community pairs
## Subset Dml. Only use the R0 secretion
Dmls <- Dml %>%
    filter(Resource1 == "R0") %>%
    select(Resource = Resource2, SecretionFlux) # secretion from R0
## Subset cml. Only use c on R0
cmls <- cml %>% filter(Resource == "R0")
## Subset lml, Only use the R0 secretion
lmls <- lml %>%
    left_join(Dmls) %>%
    mutate(CrossFeedingPotential = Leakiness * SecretionFlux, .keep = "unused") %>%
    group_by(Species) %>%
    summarize(CrossFeedingPotential = sum(CrossFeedingPotential)) %>% # Sum of R0-to-other flux times leakiness
    select(Species, CrossFeedingPotential)


pairs_comm_meta <- pairs_comm %>%
    left_join(cmls %>% rename_with(~paste0(., 1), everything())) %>%
    left_join(cmls %>% rename_with(~paste0(., 2), everything())) %>%
    left_join(lmls %>% rename_with(~paste0(., 1), everything())) %>%
    left_join(lmls %>% rename_with(~paste0(., 2), everything())) %>%
    mutate(d_ConsumptionRate = ConsumptionRate1 - ConsumptionRate2,
           d_CrossFeedingPotential = CrossFeedingPotential1 - CrossFeedingPotential2,
           .keep = "unused")

## boxplot: r_glu
p1 <- pairs_comm_meta %>%
    #filter(PairConspecific == "conspecific") %>%
    mutate(InteractionType = factor(InteractionType, c("exclusion", "coexistence"))) %>%
    ggplot(aes(x = InteractionType, y = d_ConsumptionRate, fill = InteractionType)) +
    geom_hline(yintercept = 0, linetype = 2, color = "grey") +
    geom_boxplot(width = .5) +
    geom_point(aes(group = InteractionType), shape = 1, size = 1, position = position_jitterdodge(jitter.width = 0.3)) +
    # p value
    ggpubr::stat_compare_means(aes(group = InteractionType), label = "p.format", vjust = 2, hjust = 0, method = "t.test") +
    scale_fill_manual(values = assign_interaction_color()) +
    facet_grid(.~PairConspecific) +
    theme_classic() +
    theme(legend.title = element_blank(), legend.position = "top",
          axis.text = element_text(size = 12, color = 1),
          axis.title = element_text(size = 15, color = 1),
          axis.text.x = element_text(size = 12, color = 1, angle = 30, vjust = 1, hjust = 1),
          panel.border = element_rect(color = 1, fill = NA, size = 1),
          strip.background = element_rect(color = NA, fill = NA)) +
    guides(alpha = "none", fill = "none", color = "none") +
    labs(x = "", y = expression(r[A]-r[B]))

## boxplot: cross-feeding potential
p2 <- pairs_comm_meta %>%
    #filter(PairConspecific == "conspecific") %>%
    mutate(InteractionType = factor(InteractionType, c("exclusion", "coexistence"))) %>%
    ggplot(aes(x = InteractionType, y = d_CrossFeedingPotential, fill = InteractionType)) +
    geom_hline(yintercept = 0, linetype = 2, color = "grey") +
    geom_boxplot(width = 0.5) +
    geom_point(aes(group = InteractionType), shape = 1, size = 1, position = position_jitterdodge(jitter.width = 0.3)) +
    # p value
    ggpubr::stat_compare_means(aes(group = InteractionType), label = "p.format", vjust = 2, hjust = 0, method = "t.test") +
    scale_fill_manual(values = assign_interaction_color()) +
    facet_grid(.~PairConspecific) +
    theme_classic() +
    theme(legend.title = element_blank(), legend.position = "top",
          axis.text = element_text(size = 12, color = 1),
          axis.title = element_text(size = 15, color = 1),
          axis.text.x = element_text(size = 12, color = 1, angle = 30, vjust = 1, hjust = 1),
          panel.border = element_rect(color = 1, fill = NA, size = 1),
          strip.background = element_rect(color = NA, fill = NA)) +
    guides(alpha = "none", fill = "none", color = "none") +
    labs(x = "", y = expression(X[A]-X[B]))

## Scatterplot
p3 <- pairs_comm_meta %>%
    ggplot() +
    geom_vline(xintercept = 0, linetype = 2) +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_point(aes(x = d_ConsumptionRate, y = d_CrossFeedingPotential, color = InteractionType), shape = 21, size = 2, stroke = .5) +
    scale_color_manual(values = c(assign_interaction_color())) +
    facet_grid(.~PairConspecific) +
    theme_classic() +
    theme(legend.position = "top", strip.background = element_blank(), panel.background = element_rect(color = 1)) +
    guides(color = "none") +
    labs()


p_J <- plot_grid(p1, p2, p3, ncol = 1, axis = "rl", align = "v", labels = LETTERS[1:3], scale = .9) + paint_white_background()
ggsave(here::here("plots/Fig3J-pairs_community_trait.png"), p_J, width = 6, height = 10)

## Stat
### Two sample
pairs_comm_meta %>%
    filter(PairConspecific == "conspecific") %>%
    t_test(d_ConsumptionRate ~ InteractionType, order = c("coexistence", "exclusion"))
pairs_comm_meta %>%
    filter(PairConspecific == "conspecific") %>%
    t_test(d_CrossFeedingPotential ~ InteractionType, order = c("coexistence", "exclusion"))

### glm
pairs_comm_meta %>%
    mutate_if(is.character, as.factor) %>%
    filter(PairConspecific == "conspecific") %>%
    #filter(PairConspecific != "conspecific") %>%
    #filter(!is.na(InteractionType)) %>%
    mutate(InteractionType = ifelse(InteractionType == "coexistence", 1, 0)) %>%
    glm(formula = InteractionType ~  d_ConsumptionRate * d_CrossFeedingPotential, data = ., family = "binomial") %>%
    broom::tidy() %>%
    {.}



# Figure 3K: community networks
## For pairs from a network
cn_names <- pairs_comm %>% pull(Community) %>% unique # community network names
cn_list <- pairs_comm %>%
    group_by(Community) %>%
    group_split() %>%
    lapply(function(x) {
        tt <- unique(x$Community)
        isolates <- calculate_rank(x)
        pairs <- x %>%
            mutate(across(starts_with("Species"), factor)) %>%
            rename(ID1 = Species1, ID2 = Species2) %>%
            left_join(isolates %>% select(ID, Isolate) %>% rename(ID1 = ID, Isolate1 = Isolate), by = "ID1") %>%
            left_join(isolates %>% select(ID, Isolate) %>% rename(ID2 = ID, Isolate2 = Isolate), by = "ID2") %>%
            mutate(From = Isolate1, To = Isolate2)
        make_network(isolates, pairs)
    }) %>%
    set_names(cn_names)

## Motif
cn_motif <- cn_list %>%
    lapply(function(x) tibble(Motif = 1:7, Count = count_motif(x))) %>%
    bind_rows() %>%
    mutate(Community = rep(cn_names, each = 7))
cn_motif %>%
    group_by(Community) %>%
    summarize(sum(Count))

p_K <- cn_motif %>%
    ggplot() +
    geom_point(aes(x = Motif, y = Count, color = Community), shape = 21, size = 2) +
    geom_line(aes(x = Motif, y = Count, color = Community)) +
    scale_x_continuous(breaks = 1:7) +
    theme_classic() +
    theme() +
    guides(color = "none") +
    labs()
ggsave(here::here("plots/Fig3K-community_motif.png"), p_K, width = 6, height = 4)


# Figure 3L: individual community networks
## Matrix and graph
p_cn_matrix_list <- lapply(cn_list, function(x) plot_adjacent_matrix(x) + theme(plot.margin = grid::unit(c(5,0,3,0), "mm")))
p_cn_list <- lapply(cn_list, function(x) plot_competitive_network(x, node_size = 2) + theme(plot.background = element_rect(fill = NA)))
p_list <- rep(list(NA), length(p_cn_list))
for (i in 1:length(cn_list)) p_list[[i]] <- ggdraw(p_cn_matrix_list[[i]]) + draw_plot(plot = p_cn_list[[i]], x = -.1, y = -.1, width = 0.7, height = 0.7)
p_cn_motif_count_list <- rep(list(NA), length(p_cn_list))
for (i in 1:length(p_cn_list)) {
    if (i %in% c(1, 6, 11)) p_cn_motif_count_list[[i]] <- plot_motif_count(i, cn_motif, cn_names) + theme(axis.title.x = element_blank())
    if (i %in% c(2:5, 7:10, 12:15)) p_cn_motif_count_list[[i]] <- plot_motif_count(i, cn_motif, cn_names) + theme(axis.title = element_blank())
    if (i == 16) p_cn_motif_count_list[[i]] <- plot_motif_count(i, cn_motif, cn_names)
    if (i %in% 17:20) p_cn_motif_count_list[[i]] <- plot_motif_count(i, cn_motif, cn_names) + theme(axis.title.y = element_blank())
}

## Get legend for line
p_temp <- plot_motif_count(1) + theme(legend.position = "right", legend.title = element_blank(), legend.text = element_text(size = 15))
shared_legend_line <- cowplot::get_legend(p_temp)
p_L <- list(p_list[1:5], p_cn_motif_count_list[1:5],
            p_list[6:10], p_cn_motif_count_list[6:10],
            p_list[11:15], p_cn_motif_count_list[11:15],
            p_list[16:20], p_cn_motif_count_list[16:20]) %>%
    unlist(recursive = F) %>%
    plot_grid(plotlist = ., labels = c(cn_names[1:5], rep("", 5), cn_names[6:10], rep("", 5), cn_names[11:15], rep("", 5), cn_names[16:20], rep("", 5)),
              ncol = 5, axis = "tbrl", align = "v") + paint_white_background()
ggsave(here::here("plots/Fig3L-community_matrix.png"), p_L, width = 10, height = 12)


# Figure 3M: community vs. pool
## Pairs
p1 <- bind_rows(pairs_pool %>% mutate(Assembly = "pool"), pairs_comm %>% mutate(Assembly = "community")) %>%
    mutate(InteractionType = factor(InteractionType, c("exclusion", "coexistence"))) %>%
    group_by(Assembly, InteractionType) %>%
    summarize(Count = n()) %>%
    mutate(Fraction = Count / sum(Count), TotalCount = sum(Count)) %>%
    ggplot() +
    geom_col(aes(x = Assembly, y = Fraction, fill = InteractionType), color = 1) +
    geom_text(aes(x = Assembly, label = paste0("n=", TotalCount)), y = 1, vjust = 2) +
    scale_fill_manual(values = assign_interaction_color()) +
    theme_classic() +
    theme(axis.title.x = element_blank()) +
    guides(fill = "none") +
    labs(x = "")
p1
## Pairs by family groups
p2 <- bind_rows(pairs_pool %>% mutate(Assembly = "pool"), pairs_comm %>% mutate(Assembly = "community")) %>%
    mutate(PairConspecific = ifelse(Family1 == Family2, "conspecific", "heterospecific")) %>%
    mutate(InteractionType = factor(InteractionType, c("exclusion", "coexistence"))) %>%
    group_by(Assembly, PairConspecific, InteractionType) %>%
    summarize(Count = n(), .groups = "drop_last") %>%
    mutate(Fraction = Count / sum(Count), TotalCount = sum(Count)) %>%
    ggplot() +
    geom_col(aes(x = PairConspecific, y = Fraction, fill = InteractionType), color = 1) +
    geom_text(aes(x = PairConspecific, label = paste0("n=", TotalCount)), y = 1, vjust = 2) +
    scale_fill_manual(values = assign_interaction_color()) +
    scale_y_continuous(limits = c(0, 1), breaks = scales::pretty_breaks(n=3)) +
    facet_grid(.~Assembly) +
    theme_classic() +
    theme(axis.title.x = element_blank(), strip.background = element_rect(fill = NA, color = NA),
          panel.border = element_rect(color = 1, fill = NA, size = 1)) +
    guides(fill = "none") +
    labs(x = "")

## Motifs
load(here::here("data/output/motif_list.Rdata"))
p_motif_list <- lapply(motif_list, function(x) plot_competitive_network(x, node_size = 3))
p_motif <- plot_grid(plotlist = p_motif_list, nrow = 1)

p_temp <- bind_rows(pn_motif %>% mutate(Assembly = "pool"), cn_motif %>% mutate(Assembly = "community")) %>%
    group_by(Assembly, Motif) %>%
    summarize(Count = sum(Count), .groups = "drop_last") %>%
    mutate(Fraction = Count / sum(Count), TotalCount = sum(Count))%>%
    ggplot() +
    geom_point(aes(x = Motif, y = Fraction, color = Assembly), shape = 21, size = 2, stroke = 2) +
    scale_color_manual(values = c("pool" = "grey50", "community" = "red")) +
    scale_x_continuous(breaks = 1:7) +
    scale_y_continuous(limits = c(0, 1), breaks = scales::pretty_breaks(n=3)) +
    facet_grid(.~Motif, scales = "free_x") +
    theme_classic() +
    theme(panel.border = element_rect(color = 1, fill = NA),
          panel.grid.major.x = element_line(linetype = 2, color = "grey10")) +
    guides(color = guide_legend(title = "")) +
    labs()

p3 <- plot_grid(p_motif, p_temp, ncol = 1, axis = "lr", align = "v", rel_heights = c(1, 1.5))
p_upper <- plot_grid(p1, p2, nrow = 1, axis = "tb", align = "h", labels = c("A", "B"), scale = .9)
p_M <- plot_grid(p_upper, p3, ncol = 1, labels = c("", "C"), rel_heights = c(1,1.5), scale = c(1, .9)) + paint_white_background()
ggsave(here::here("plots/Fig3M-community_versus_pool.png"), p_M, width = 10, height = 6)

















