# Figure 3: model
library(tidyverse)
library(tidymodels)
library(tidygraph)
library(cowplot)
library(ggpubr)
library(ggraph)
source(here::here("plotting_scripts/network_functions.R"))


#input_independent <- read_csv(here::here("simulation/input_independent.csv"), col_types = cols())
#input_independent <- read_csv(("~/Dropbox/lab/invasion-network/simulation/data/raw3/input_independent.csv"), col_types = cols())
input_independent <- read_csv(("~/Dropbox/lab/invasion-network/simulation/data/raw4/input_independent.csv"), col_types = cols())
input_pairs <- read_csv(("~/Dropbox/lab/invasion-network/simulation/data/raw4/input_pairs.csv"), col_types = cols())
output_dir <- input_independent$output_dir[1]
category_colors <- c(sugar = "#ED6A5A", acid = "#03CEA4", waste = "#51513D", fermenter = "#8A89C0", respirator = "#FFCB77", F0 = "#8A89C0", F1 = "#FFCB77")
input_row <- input_independent[2,]
# Generate family-species and class-resource matching tibble
sa <- input_independent$sa[1]
ma <- input_independent$ma[1]
sal <- tibble(Family = paste0("F", c(rep(0, sa), rep(1, sa))), Species = paste0("S", 0:(sa * 2 - 1)))
mal <- tibble(Class = paste0("T", c(rep(0, ma), rep(1, ma), rep(2, ma))), Resource = paste0("R", 0:(ma * 3 - 1)))

# Functions
read_wide_file <- function(x) {
    tt <- read_csv(x, col_types = cols()) %>%
        # Remove abundance=0
        mutate_all(~replace(., .==0, NA)) %>%
        pivot_longer(cols = starts_with("W"), names_to = "Well", values_to = "Abundance", values_drop_na = T)
    if ("...1" %in% colnames(tt)) tt <- tt %>% rename(Family = ...1, Species = ...2)
    return(tt)
}
paint_white_background <- function(x) theme(plot.background = element_rect(color = NA, fill = "white"))

# Figure 3A: diagram cartoon
#p_A <- ggdraw() + draw_image(here::here("plots/cartoons/Fig3A.png")) + theme(plot.background = element_rect(fill = "white", color = NA))
p_A <-  ggplot(mtcars, aes(x = wt, y = mpg)) + annotate("text", x = 0 , y = 0, label = "Cartoon for\nmodel") + theme_void() + theme(plot.background = element_rect(fill = "white", color = NA))
ggsave(here::here("plots/Fig3A-functional_groups.png"), p_A, width = 5, height = 5)


# Figure 3B: u, l, and D matrices
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
    scale_color_manual(values = category_colors, breaks = c("sugar", "acid", "waste")) +
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
    scale_color_manual(values = category_colors) +
    scale_x_discrete(limits = c("spacer", mal$Resource)) +
    scale_y_discrete(limits = c("spacer", sal$Species)) +
    theme_classic() +
    theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank()) +
    guides(fill = guide_colorbar(title = expression(c[i][alpha])),
           color = "none") +
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
    scale_color_manual(values = category_colors) +
    scale_x_discrete(limits = c("spacer", mal$Resource)) +
    scale_y_discrete(limits = c("spacer", sal$Species)) +
    theme_classic() +
    theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank()) +
    guides(fill = guide_colorbar(title = expression(l[i][alpha])),
           color = guide_legend(title = "")) +
    labs()
p_lower <- plot_grid(p2, p3, labels = c("c matrix", "l matrix"), scale = .9)
p_B <- plot_grid(p1, p_lower, ncol = 1, labels = c("D matrix", ""), align = "v", axis = "lr", scale = c(0.9, 1)) + paint_white_background()
ggsave(here::here("plots/Fig3B-matrices.png"), p_B, width = 12, height = 12)


# Figure 3C: community composition
df_comm <- input_independent %>% filter(str_detect(init_N0, "selfAssembly")) %>%
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

p_C <- df_comm %>%
    ggplot() +
    geom_col(aes(x = Community, y = Abundance, fill = Family, group = Species), color = 1, position = "fill") +
    scale_fill_manual(values = category_colors, breaks = c("F0", "F1")) +
    scale_y_continuous(breaks = c(0, .5, 1)) +
    theme_classic() +
    labs()

ggsave(here::here("plots/Fig3C-community_composition.png"), p_C, width = 7, height = 3)


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
    scale_fill_manual(values = category_colors, breaks = c("F0", "F1")) +
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
    mutate(Conspecific = with(., case_when(
        (Family1 == Family2) ~ "conspecific",
        (Family1 != Family2) ~ "heterospecific"
    ))) %>%
    mutate(Community = factor(Community, paste0("W", 0:1000))) %>%
    arrange(Community)
temp <- pairs_pool %>%
    filter(InteractionType != "no-growth") %>%
    mutate(InteractionType = factor(InteractionType, c("exclusion", "coexistence"))) %>%
    group_by(Community, InteractionType, Conspecific) %>%
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
    group_by(InteractionType, Conspecific) %>%
    summarize(Count = n()) %>%
    group_by(Conspecific) %>%
    mutate(TotalCount = sum(Count)) %>%
    ggplot() +
    geom_col(aes(x = Conspecific, y = Count, fill = InteractionType), color = 1, position = "fill") +
    #geom_text(aes(x = Conspecific, y = Count, label = round(Percentage,2)), vjust = -1) +
    geom_text(aes(x = Conspecific, label = paste0("n=", TotalCount)), y = 1, vjust = 2) +
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


temp <- pairs_pool %>%
    left_join(cmls %>% rename_with(~paste0(., 1), everything())) %>%
    left_join(cmls %>% rename_with(~paste0(., 2), everything())) %>%
    left_join(lmls %>% rename_with(~paste0(., 1), everything())) %>%
    left_join(lmls %>% rename_with(~paste0(., 2), everything())) %>%
    mutate(d_ConsumptionRate = ConsumptionRate1 - ConsumptionRate2,
           d_CrossFeedingPotential = CrossFeedingPotential1 - CrossFeedingPotential2,
           .keep = "unused")

p_F <- temp %>%
    #filter(InteractionType != "no-growth") %>%
    ggplot() +
    geom_vline(xintercept = 0, linetype = 2) +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_point(aes(x = d_ConsumptionRate, y = d_CrossFeedingPotential, color = InteractionType), shape = 21, size = 2, stroke = .5) +
    scale_color_manual(values = c(assign_interaction_color())) +
    facet_grid(.~Conspecific) +
    theme_classic() +
    theme(legend.position = "top", strip.background = element_blank(), panel.background = element_rect(color = 1)) +
    guides(color = guide_legend(title = "")) +
    labs()

ggsave(here::here("plots/Fig3F-pairs_pool_trait.png"), p_F, width = 6, height = 4)

## Stat
### Two sample
temp %>%
    #filter(InteractionType != "no-growth") %>%
    t_test(d_CrossFeedingPotential ~ InteractionType, order = c("coexistence", "exclusion"))
### glm
temp %>%
    mutate_if(is.character, as.factor) %>%
    #filter(InteractionType != "no-growth") %>%
    filter(Conspecific == "conspecific") %>%
    #filter(Conspecific != "conspecific") %>%
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

network_names <- pairs_pool %>% pull(Community) %>% unique
net_list <- pairs_pool %>%
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
    set_names(network_names)

## Motif
networks_motif <- net_list %>%
    lapply(function(x) tibble(Motif = 1:7, Count = count_motif(x))) %>%
    bind_rows() %>%
    mutate(Community = rep(network_names, each = 7))
networks_motif %>%
    group_by(Community) %>%
    summarize(sum(Count))

p_G <- networks_motif %>%
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
p_net_matrix_list <- lapply(net_list, function(x) plot_adjacent_matrix(x) + theme(plot.margin = grid::unit(c(5,0,3,0), "mm")))
p_net_list <- lapply(net_list, function(x) plot_competitive_network(x, node_size = 2) + theme(plot.background = element_rect(fill = NA)))
p_list <- rep(list(NA), length(p_net_list))
for (i in 1:length(net_list)) p_list[[i]] <- ggdraw(p_net_matrix_list[[i]]) + draw_plot(plot = p_net_list[[i]], x = -.1, y = -.1, width = 0.7, height = 0.7)
## Motif count
plot_motif_count <- function (x = 1) {
    # motif_randomized_subset <- networks_motif_randomized_percentile %>%
    #     filter(Community %in% network_names[x]) %>%
    #     mutate(Community = factor(Community, network_names[x]))
    motif_community_subset <- networks_motif %>%
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
p_motif_count_list <- rep(list(NA), length(p_net_list))
for (i in 1:length(p_net_list)) {
    if (i %in% c(1, 6, 11)) p_motif_count_list[[i]] <- plot_motif_count(i) + theme(axis.title.x = element_blank())
    if (i %in% c(2:5, 7:10, 12:15)) p_motif_count_list[[i]] <- plot_motif_count(i) + theme(axis.title = element_blank())
    if (i == 16) p_motif_count_list[[i]] <- plot_motif_count(i)
    if (i %in% 17:20) p_motif_count_list[[i]] <- plot_motif_count(i) + theme(axis.title.y = element_blank())
}

## Get legend for line
p_temp <- plot_motif_count(1) + theme(legend.position = "right", legend.title = element_blank(), legend.text = element_text(size = 15))
shared_legend_line <- cowplot::get_legend(p_temp)
p_H <- list(p_list[1:5], p_motif_count_list[1:5],
            p_list[6:10], p_motif_count_list[6:10],
            p_list[11:15], p_motif_count_list[11:15],
            p_list[16:20], p_motif_count_list[16:20]) %>%
    unlist(recursive = F) %>%
    plot_grid(plotlist = ., labels = c(network_names[1:5], rep("", 5), network_names[6:10], rep("", 5), network_names[11:15], rep("", 5), network_names[16:20], rep("", 5)),
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
    mutate(Conspecific = with(., case_when(
        (Family1 == Family2) ~ "conspecific",
        (Family1 != Family2) ~ "heterospecific"
    ))) %>%
    mutate(Community = factor(Community, paste0("W", 0:1000))) %>%
    arrange(Community)
temp <- pairs_comm %>%
    filter(InteractionType != "no-growth") %>%
    mutate(InteractionType = factor(InteractionType, c("exclusion", "coexistence"))) %>%
    group_by(Community, InteractionType, Conspecific) %>%
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
    group_by(InteractionType, Conspecific) %>%
    summarize(Count = n()) %>%
    group_by(Conspecific) %>%
    mutate(TotalCount = sum(Count)) %>%
    ggplot() +
    geom_col(aes(x = Conspecific, y = Count, fill = InteractionType), color = 1, position = "fill") +
    #geom_text(aes(x = Conspecific, y = Count, label = round(Percentage,2)), vjust = -1) +
    geom_text(aes(x = Conspecific, label = paste0("n=", TotalCount)), y = 1, vjust = 2) +
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


temp <- pairs_comm %>%
    left_join(cmls %>% rename_with(~paste0(., 1), everything())) %>%
    left_join(cmls %>% rename_with(~paste0(., 2), everything())) %>%
    left_join(lmls %>% rename_with(~paste0(., 1), everything())) %>%
    left_join(lmls %>% rename_with(~paste0(., 2), everything())) %>%
    mutate(d_ConsumptionRate = ConsumptionRate1 - ConsumptionRate2,
           d_CrossFeedingPotential = CrossFeedingPotential1 - CrossFeedingPotential2,
           .keep = "unused")

p_J <- temp %>%
    #filter(InteractionType != "no-growth") %>%
    ggplot() +
    geom_vline(xintercept = 0, linetype = 2) +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_point(aes(x = d_ConsumptionRate, y = d_CrossFeedingPotential, color = InteractionType), shape = 21, size = 2, stroke = .5) +
    scale_color_manual(values = c(assign_interaction_color())) +
    facet_grid(.~Conspecific) +
    theme_classic() +
    theme(legend.position = "top", strip.background = element_blank(), panel.background = element_rect(color = 1)) +
    guides(color = guide_legend(title = "")) +
    labs()

ggsave(here::here("plots/Fig3J-pairs_community_trait.png"), p_J, width = 6, height = 4)

## Stat
### Two sample
temp %>%
    #filter(InteractionType != "no-growth") %>%
    t_test(d_CrossFeedingPotential ~ InteractionType, order = c("coexistence", "exclusion"))
### glm
temp %>%
    mutate_if(is.character, as.factor) %>%
    #filter(InteractionType != "no-growth") %>%
    filter(Conspecific == "conspecific") %>%
    #filter(Conspecific != "conspecific") %>%
    #filter(!is.na(InteractionType)) %>%
    mutate(InteractionType = ifelse(InteractionType == "coexistence", 1, 0)) %>%
    glm(formula = InteractionType ~  d_ConsumptionRate * d_CrossFeedingPotential, data = ., family = "binomial") %>%
    broom::tidy() %>%
    {.}



# Figure 3K: community networks
## For pairs from a network
network_names <- pairs_comm %>% pull(Community) %>% unique
net_list <- pairs_comm %>%
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
    set_names(network_names)

## Motif
networks_motif <- net_list %>%
    lapply(function(x) tibble(Motif = 1:7, Count = count_motif(x))) %>%
    bind_rows() %>%
    mutate(Community = rep(network_names, each = 7))
networks_motif %>%
    group_by(Community) %>%
    summarize(sum(Count))

p_K <- networks_motif %>%
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
p_net_matrix_list <- lapply(net_list, function(x) plot_adjacent_matrix(x) + theme(plot.margin = grid::unit(c(5,0,3,0), "mm")))
p_net_list <- lapply(net_list, function(x) plot_competitive_network(x, node_size = 2) + theme(plot.background = element_rect(fill = NA)))
p_list <- rep(list(NA), length(p_net_list))
for (i in 1:length(net_list)) p_list[[i]] <- ggdraw(p_net_matrix_list[[i]]) + draw_plot(plot = p_net_list[[i]], x = -.1, y = -.1, width = 0.7, height = 0.7)
## Motif count
plot_motif_count <- function (x = 1) {
    # motif_randomized_subset <- networks_motif_randomized_percentile %>%
    #     filter(Community %in% network_names[x]) %>%
    #     mutate(Community = factor(Community, network_names[x]))
    motif_community_subset <- networks_motif %>%
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
p_motif_count_list <- rep(list(NA), length(p_net_list))
for (i in 1:length(p_net_list)) {
    if (i %in% c(1, 6, 11)) p_motif_count_list[[i]] <- plot_motif_count(i) + theme(axis.title.x = element_blank())
    if (i %in% c(2:5, 7:10, 12:15)) p_motif_count_list[[i]] <- plot_motif_count(i) + theme(axis.title = element_blank())
    if (i == 16) p_motif_count_list[[i]] <- plot_motif_count(i)
    if (i %in% 17:20) p_motif_count_list[[i]] <- plot_motif_count(i) + theme(axis.title.y = element_blank())
}

## Get legend for line
p_temp <- plot_motif_count(1) + theme(legend.position = "right", legend.title = element_blank(), legend.text = element_text(size = 15))
shared_legend_line <- cowplot::get_legend(p_temp)
p_L <- list(p_list[1:5], p_motif_count_list[1:5],
            p_list[6:10], p_motif_count_list[6:10],
            p_list[11:15], p_motif_count_list[11:15],
            p_list[16:20], p_motif_count_list[16:20]) %>%
    unlist(recursive = F) %>%
    plot_grid(plotlist = ., labels = c(network_names[1:5], rep("", 5), network_names[6:10], rep("", 5), network_names[11:15], rep("", 5), network_names[16:20], rep("", 5)),
              ncol = 5, axis = "tbrl", align = "v") + paint_white_background()
ggsave(here::here("plots/Fig3L-community_matrix.png"), p_L, width = 10, height = 12)

