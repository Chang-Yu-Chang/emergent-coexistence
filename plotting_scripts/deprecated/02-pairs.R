# Figure 2: pairs alike coexist more than pairs dissimilar

library(tidyverse)
library(tidymodels)
library(cowplot)
library(ggpubr)
library(tidygraph)
library(ggraph)
library(vip)
library(officer)
library(flextable)
source(here::here("plotting_scripts/network_functions.R"))

isolates <- read_csv(here::here("data/output/isolates.csv"), col_types = cols())
communities <- read_csv(here::here("data/output/communities.csv"), col_types = cols())
pairs <- read_csv(here::here("data/output/pairs.csv"), col_types = cols()) %>% mutate(InteractionType = ifelse(InteractionType == "neutrality", "coexistence", InteractionType))
pairs_meta <- read_csv(here::here("data/output/pairs_meta.csv"), col_types = cols()) %>% mutate(InteractionType = ifelse(InteractionType == "neutrality", "coexistence", InteractionType))
load("~/Dropbox/lab/invasion-network/data/output/network.Rdata") # Randomized networks
load("~/Dropbox/lab/invasion-network/data/output/network_randomized.Rdata")

#
pairs_coexistence <- pairs_meta %>%
    filter(!is.na(PairFermenter)) %>%
    filter(Assembly == "self_assembly") %>%
    mutate(PairConspecific = ifelse(PairFermenter == "FF" | PairFermenter == "NN", "conspecific", ifelse(PairFermenter == "FN", "heterospecific", NA))) %>%
    mutate(Rank_d = Rank1 - Rank2, Score_d = Score1 - Score2)
pairs_count <- pairs_coexistence %>%
    group_by(PairConspecific) %>%
    count()

# Figure S10. Lasso regression testing whether pairwise coexistence is predicted by metabolic dissimilarity --------------
# Build a training set
pairs_train <- pairs_meta %>%
    filter(Assembly == "self_assembly") %>%
    mutate(InteractionType = ifelse(InteractionType == "coexistence", 1, 0)) %>%
    select(InteractionType, ends_with("d")) %>%
    mutate(across(where(is.character), as_factor)) %>%
    mutate(Pair = 1:n()) %>%
    # Drop all pairs with na in the features we are interested in
    drop_na()

# Colinearity
#model <- lm(InteractionType ~ ., data = pairs_train)
#car::vif(model)
pairs_train %>%
    select(-InteractionType, Pair) %>%
    cor() %>%
    as_tibble() %>%
    mutate(Feature1 = colnames(.)) %>%
    pivot_longer(-Feature1, names_to = "Feature2", values_to = "r") %>%
    ggplot() +
    geom_histogram(aes(x = r), color = 1, fill = "white") +
    theme_classic()


if(FALSE) {
    ## Prepare recipe
    pairs_rec <- recipe(InteractionType ~ ., data = pairs_train) %>%
        update_role(Pair, new_role = "ID") %>% # Pair identifier
        step_dummy(all_nominal()) %>% # Create dummy variables for categorical variables
        step_zv(all_numeric(), -all_numeric()) %>% # Remove numeric variables that have 0 variance
        step_normalize(all_numeric(), -all_outcomes()) # Normalize (center and rescale) the numeric variables
    pairs_rec <- prepare_recipe(pairs_train)
    summary(pairs_rec)
    pairs_prep <- pairs_rec %>% prep(strings_as_factor = F)

    # Create workflow
    wf <- workflow() %>% add_recipe(pairs_rec)
    ## Specify and fit the model for one penalty value
    lasso_spec <- linear_reg(penalty = 0.1, mixture = 1) %>% set_engine("glmnet")
    lasso_fit <- wf %>% add_model(lasso_spec) %>% fit(data = pairs_train)
    lasso_fit %>%
        extract_fit_parsnip() %>%
        tidy() %>%
        arrange(desc(abs(estimate)))
    ## Specify and fit the model. Tune the lasso parameter: penalty
    set.seed(1)
    pairs_boot <- bootstraps(pairs_train) # Build a set of bootstrap resamples
    tune_spec <- linear_reg(penalty = tune(), mixture = 1) %>% set_engine("glmnet")
    lambda_grid <- grid_regular(penalty(), levels = 50) # a grid of penalty values
    lasso_grid <- tune_grid(wf %>% add_model(tune_spec), resamples = pairs_boot, grid = lambda_grid)
    ## Finalize the workflow. Pick the penalty value with lowest rmse (root mean square error)
    lowest_rmse <- lasso_grid %>% select_best("rmse")
    final_lasso <- finalize_workflow(wf %>% add_model(tune_spec), lowest_rmse)
}

# Correlation plot between all predictors
factor_level <- c(str_subset(names(pairs_train), "X"),
                  str_subset(names(pairs_train), "r_glucose"),
                  str_subset(names(pairs_train), "r_acetate"),
                  str_subset(names(pairs_train), "r_lactate"),
                  str_subset(names(pairs_train), "r_succinate"),
                  str_subset(names(pairs_train), "pH"))
p0_0 <- pairs_train %>%
    select(-Pair, -InteractionType) %>%
    cor() %>%
    as_tibble() %>%
    mutate(Row = names(.)) %>%
    pivot_longer(cols = -Row, names_to = "Column") %>%
    mutate(Row = ordered(Row, factor_level) %>% fct_rev(), Column = ordered(Column, factor_level)) %>%
    ggplot() +
    geom_tile(aes(x = Column, y = Row, fill = value)) +
    scale_fill_gradient2() +
    scale_x_discrete(position = "bottom") +
    scale_y_discrete(position = "right") +
    #scale_y_reverse() +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0), axis.title = element_blank(),
          legend.position = "top") +
    guides(fill = guide_colorbar(title = "Correlation coefficient")) +
    labs()

## Color legend
temp <- tibble(Variable = factor_level) %>%
    mutate(Type = case_when(str_detect(Variable, "^X_") ~ "secretion",
                            str_detect(Variable, "^r_") ~ "growth rate",
                            str_detect(Variable, "^pH") ~ "pH")) %>%
    mutate(temp = n():1)
temp2 <- temp %>%
    group_by(Type) %>%
    summarize(Variable = median(temp), Text = unique(Type))
p0_1 <- temp %>%
    ggplot() +
    geom_tile(aes(x = .5, y = Variable, fill = Type), width = .1) +
    geom_text(data = temp2, aes(x = .75, y = Variable, label = Text, color = Type), angle = 270) +
    scale_x_continuous(limits = c(0,1)) +
    theme_void() +
    theme() +
    guides(fill = "none", color = "none") +
    labs()

p0_2 <- temp %>%
    ggplot() +
    geom_tile(aes(x = Variable, y = .75, fill = Type), height = .1) +
    geom_text(data = temp2, aes(x = Variable, y = .5, label = Text, color = Type)) +
    scale_y_continuous(limits = c(0,1)) +
    theme_void() +
    theme() +
    guides(fill = "none", color = "none") +
    labs()
p0 <- plot_grid(p0_0, p0_1, nrow = 1, rel_widths = c(1, 0.1), axis = "tb", align = "h")

if (FALSE) {
    ## Find penalty parameter with lowest rmse
    p2 <- lasso_grid %>%
        collect_metrics() %>%
        drop_na() %>%
        ggplot(aes(penalty, mean)) +
        geom_vline(xintercept = lowest_rmse$penalty, color = "red") +
        geom_errorbar(aes(ymin = mean - std_err, ymax = mean + std_err), alpha = 0.5) +
        geom_line(size = 1.5) +
        #geom_text(x = lowest_rmse$penalty, y = Inf, label = paste0("penalty=", round(lowest_rmse$penalty, 3)), vjust = 1) +
        facet_wrap(~.metric, scales = "free_y", nrow = 2, labeller = labeller(.metric = toupper)) +
        scale_x_log10() +
        theme_classic() +
        theme(legend.position = "none", strip.background = element_rect(color = NA),
              panel.border = element_rect(fill = NA, color = 1)) +
        labs()

}

# Lasso
# Comprehensive function for implementing lasso
fit_lasso <- function (pairs_train) {
    # Prepare recipe
    pairs_rec <- recipe(InteractionType ~ ., data = pairs_train) %>%
        update_role(Pair, new_role = "ID") %>% # Pair identifier
        step_dummy(all_nominal()) %>% # Create dummy variables for categorical variables
        step_zv(all_numeric(), -all_numeric()) %>% # Remove numeric variables that have 0 variance
        step_normalize(all_numeric(), -all_outcomes()) # Normalize (center and rescale) the numeric variables
    pairs_prep <- pairs_rec %>% prep(strings_as_factor = F)
    # Create workflow
    wf <- workflow() %>% add_recipe(pairs_rec)
    ## Specify and fit the model. Tune the lasso parameter: penalty
    set.seed(1)
    pairs_boot <- bootstraps(pairs_train) # Build a set of bootstrap resamples
    tune_spec <- linear_reg(penalty = tune(), mixture = 1) %>% set_engine("glmnet")
    lambda_grid <- grid_regular(penalty(), levels = 50) # a grid of penalty values
    lasso_grid <- tune_grid(wf %>% add_model(tune_spec), resamples = pairs_boot, grid = lambda_grid)
    ## Finalize the workflow. Pick the penalty value with lowest rmse (root mean square error)
    lowest_rmse <- lasso_grid %>% select_best("rmse")
    final_lasso <- finalize_workflow(wf %>% add_model(tune_spec), lowest_rmse)
    return(list(final_lasso=final_lasso, lowest_rmse=lowest_rmse))
}
plot_lasso <- function(lasso, pairs, lowest_rmse) {
    lasso %>%
        fit(pairs) %>%
        extract_fit_parsnip() %>%
        #tidy() %>%
        #arrange(desc(abs(estimate))) %>%
        vi(lambda = lowest_rmse$penalty) %>%
        mutate(Importance = abs(Importance), Variable = fct_reorder(Variable, Importance)) %>%
        ggplot(aes(x = Importance, y = Variable, fill = Sign)) +
        geom_col(color = 1, width = .8) +
        scale_fill_manual(values = c("NEG" = "black", "POS" = "white")) +
        theme_classic() +
        labs(y = NULL, fill = "")

}
if (FALSE) {
    pairs_train2 <- pairs_meta %>%
        filter(Assembly == "self_assembly") %>%
        mutate(InteractionType = ifelse(InteractionType == "coexistence", 1, 0)) %>%
        select(InteractionType, ends_with("d")) %>%
        select(InteractionType, contains("X_") & !contains("_0hr"), contains("curver")) %>%
        mutate(across(where(is.character), as_factor)) %>%
        mutate(Pair = 1:n()) %>%
        drop_na()
    names(pairs_train2)
    lasso2 <- fit_lasso(pairs_train2)
    p2 <- plot_lasso(lasso2$final_lasso, pairs_train2, lasso2$lowest_rmse) +
        theme(legend.position = "none") +
        ggtitle("All pairs, all d")
    p2


}

## All pairs, all d
pairs_train1 <- pairs_meta %>%
    filter(Assembly == "self_assembly") %>%
    mutate(InteractionType = ifelse(InteractionType == "coexistence", 1, 0)) %>%
    select(InteractionType, ends_with("d")) %>%
    mutate(across(where(is.character), as_factor)) %>%
    mutate(Pair = 1:n()) %>%
    drop_na()
lasso1 <- fit_lasso(pairs_train1)
p1 <- plot_lasso(lasso1$final_lasso, pairs_train1, lasso1$lowest_rmse) +
    theme(legend.position = "none") +
    ggtitle("All pairs, all d")

## FF pairs
pairs_train2 <- pairs_meta %>%
    filter(Assembly == "self_assembly") %>%
    mutate(InteractionType = ifelse(InteractionType == "coexistence", 1, 0)) %>%
    filter(PairFermenter == "FF") %>%
    select(InteractionType, ends_with("d")) %>%
    mutate(across(where(is.character), as_factor)) %>%
    mutate(Pair = 1:n()) %>%
    drop_na()
lasso2 <- fit_lasso(pairs_train2)
p2 <- plot_lasso(lasso2$final_lasso, pairs_train2, lasso2$lowest_rmse) +
    theme(legend.position = "bottom") +
    ggtitle("FF pairs, all d")

## All pairs, only X_d
pairs_train3 <- pairs_meta %>%
    filter(Assembly == "self_assembly") %>%
    mutate(InteractionType = ifelse(InteractionType == "coexistence", 1, 0)) %>%
    dplyr::select(InteractionType, ends_with("d") & !starts_with("r")) %>%
    mutate(across(where(is.character), as_factor)) %>%
    mutate(Pair = 1:n()) %>%
    drop_na()
lasso3 <- fit_lasso(pairs_train3)
p3 <- plot_lasso(lasso3$final_lasso, pairs_train3, lasso3$lowest_rmse) +
    theme(legend.position = "none") +
    ggtitle("All pairs, only X_d")


p_bottom <- plot_grid(p1, p2, p3, nrow = 1, labels = LETTERS[2:4], scale = .9, axis = "tb", align = "h")
p <- plot_grid(p0, p_bottom, ncol = 1, labels = c("A", ""), rel_heights = c(1.5, 1), scale = c(.9, 1)) + paint_white_background()
ggsave(here::here("plots/FigS10-lasso_regression.png"), p, width = 10, height = 14)



# Figure S11. r_mid----
## r
p1 <- isolates %>%
    filter(Assembly == "self_assembly") %>%
    select(Community, Fermenter, r_glucose_midhr) %>%
    drop_na() %>%
    mutate(Fermenter = ifelse(Fermenter, "fermenter", ifelse(!Fermenter, "respirator", NA))) %>%
    #mutate(Group = ifelse(r_glucose_midhr > 0.07, "fermenter", "respirator")) %>%
    ggplot() +
    geom_histogram(aes(x = r_glucose_midhr, fill = Fermenter), color = 1, binwidth = 0.01) +
    geom_vline(xintercept = 0.06, color = "red", linetype = 2) +
    scale_fill_manual(values = fermenter_color) +
    theme_classic() +
    guides(fill = "none") +
    labs()

## r_mid vs. isolate Rank
isolate_ranked <- isolates %>%
    filter(Assembly == "self_assembly") %>%
    select(Community, Fermenter, r_glucose_midhr, Rank) %>%
    group_by(Community) %>%
    mutate(r_ranked = rank(-r_glucose_midhr)) %>% # Put a negative sign such that the highest abs value is rank 1
    mutate(Fermenter = ifelse(Fermenter, "fermenter", ifelse(!Fermenter, "respirator", NA)))
## Stat
cor.test(isolate_ranked$Rank, isolate_ranked$r_ranked) %>% tidy()
## Plot
p2 <- isolate_ranked %>%
    ggplot() +
    geom_smooth(aes(x = r_ranked, y = Rank), formula = y ~ x, method = "lm", color = 1) +
    geom_point(aes(x = r_ranked, y = Rank, color = Fermenter), shape = 21, size = 2, stroke = 1, position = position_jitter(width = .1, height = .1)) +
    scale_x_continuous(breaks = 1:12) +
    scale_y_continuous(breaks = 1:12) +
    scale_color_manual(values = fermenter_color) +
    theme_classic() +
    theme(legend.position = "right", legend.title = element_blank()) +
    labs(x = "r_mid rank", y = "Isolate rank")

p <- plot_grid(p1, p2, nrow = 1, rel_widths = c(1, 1.5), labels = c("A", "B"), axis = "tb", align = "h", scale = .9) + paint_white_background()
ggsave(here::here("plots/FigS11-r_mid.png"), p, width = 6, height = 2.5)


# Figure S12. Leakiness ----
isolates_byproduct <- isolates %>%
    filter(Assembly == "self_assembly") %>%
    select(Community, Isolate, ID, Fermenter, ends_with("hr")) %>%
    pivot_longer(cols = ends_with("hr"), names_pattern = "(.*)_(.*)hr", names_to = c("Measure", "Time"), values_to = "Value") %>%
    filter(!is.na(Value)) %>%
    mutate(Fermenter = ifelse(Fermenter, "fermenter", "respirator"))
## Secretion over time
p1 <- isolates_byproduct %>%
    # Only use acids
    filter(str_detect(Measure, "X_"), !str_detect(Measure, "sum")) %>%
    mutate(Measure = str_replace(Measure, "X_", "")) %>%
    ggplot(aes(x = Time, y = Value, group = ID, color = Fermenter, alpha = ID)) +
    geom_point() +
    geom_line() +
    scale_color_manual(values = category_color, breaks = c("fermenter", "respirator")) +
    #scale_color_npg(labels = c("TRUE" = "Fermenter", "FALSE" = "Respirator"), breaks = c(T, F)) +
    facet_wrap(Measure~., nrow = 1, scales = "free_y") +
    theme_classic() +
    theme(legend.title = element_blank(), legend.position = "top",
          strip.background = element_rect(color = NA, fill = NA),
          panel.border = element_rect(fill = NA, color = 1)) +
    guides(alpha = "none") +
    labs(x = "Time (hr)", y = "Concentration (mM)")
p_top <- get_legend(p1)
p1 <- p1 + theme(legend.position = "none")

## pH over time
p2 <- isolates_byproduct %>%
    filter(Measure == "pH") %>%
    ggplot(aes(x = Time, y = Value, group = ID, color = Fermenter, alpha = ID)) +
    geom_point() +
    geom_line() +
    scale_color_manual(values = category_color, breaks = c("fermenter", "respirator")) +
    facet_wrap(Measure~., nrow = 1, scales = "free_y") +
    theme_classic() +
    theme(legend.title = element_blank(), legend.position = "none",
          strip.background = element_rect(color = NA, fill = NA),
          panel.border = element_rect(fill = NA, color = 1)) +
    guides(alpha = "none") +
    labs(x = "Time (hr)", y = "pH")

## Leakiness at 16 hr
temp <- isolates_byproduct %>%
    filter(Assembly == "self_assembly") %>%
    drop_na(Fermenter, leakiness_16hr)
p3 <- isolates_byproduct %>%
    filter(Time == "16", Measure == "leakiness") %>%
    ggplot() +
    geom_boxplot(aes(x = Fermenter, y = Value, color = Fermenter), outlier.size = 2) +
    geom_jitter(aes(x = Fermenter, y = Value, color = Fermenter), shape = 1, size = 2, width = 0.3) +
    ggpubr::stat_compare_means(aes(group = Fermenter, x = Fermenter, y = Value), method = "t.test", vjust = 1) +
    scale_y_continuous(limits = c(0, 0.65)) +
    scale_color_manual(values = category_color, breaks = c("fermenter", "respirator")) +
    theme_classic() +
    theme(axis.title.x = element_blank(), legend.position = "none", panel.border = element_rect(fill = NA, color = 1)) +
    labs(y = "Leakiness")

p_bottom <- plot_grid(p2, p3, NULL, nrow = 1, axis = "tb", align = "h", labels = c("B", "C"))
p <- plot_grid(p_upper, p1, p_bottom, ncol = 1, rel_heights = c(1, 5, 5), labels = c("", "A", "")) + paint_white_background()
ggsave(here::here("plots/FigS12-leakiness.png"), p, width = 8, height = 6)


# Figure S13. r_mid two group ----
p <- isolates %>%
    filter(Assembly == "self_assembly") %>%
    select(Community, Fermenter, r_glucose_midhr, leakiness_16hr) %>%
    drop_na() %>%
    mutate(Group = ifelse(r_glucose_midhr > 0.09 & leakiness_16hr > 0.1, "fermenter", "respirator")) %>%
    ggplot(aes(x = r_glucose_midhr, y = leakiness_16hr, color = Group)) +
    geom_point(size = 2) +
    stat_ellipse() +
    scale_color_manual(values = category_color, breaks = c("fermenter", "respirator")) +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 3)) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 3)) +
    theme_classic() +
    theme(legend.position = "right") +
    labs(color = "")
ggsave(here::here("plots/FigS13-r_vs_leakiness.png"), p, width = 4, height = 3)



# Table S2. features used in the lasso regression ----
features <- pairs_meta %>%
    filter(Assembly == "self_assembly") %>%
    mutate(InteractionType = ifelse(InteractionType == "coexistence", 1, 0)) %>%
    select(ends_with("d")) %>%
    mutate(across(where(is.character), as_factor))  %>%
    colnames() %>%
    sort()
ft2 <- tibble(Feature = features) %>%
    mutate(Type = str_sub(Feature, 1, 2) %>% str_replace("_", "")) %>%
    mutate(Type = ifelse(Type == "r", "growth rate", ifelse(Type == "X", "secretion", Type))) %>%
    flextable() %>%
    width(j = 2, width = 2)
save_as_image(ft2, here::here("plots/TableS2.png"))









# deprecated ----

if (FALSE) {

# Figure 3A: one example community of crossfeeding networks
## Isolate growth rate and secretion
temp <- isolates %>%
    filter(Assembly == "self_assembly") %>%
    select(ID, Fermenter, starts_with("r_"), starts_with("X_")) %>%
    select(ID, Fermenter, ends_with("midhr"), ends_with("16hr")) %>%
    drop_na(r_glucose_16hr, X_acetate_16hr)
## Community and strain ID
temp2 <- isolates %>%
    mutate(Node = as.character(ID), Fermenter = ifelse(Fermenter, "fermenter", ifelse(Fermenter == FALSE, "respirator", NA))) %>%
    select(Node, Community, Fermenter)
## Fluxes
isolates_in <- temp %>%
    select(ID, Fermenter, ends_with("midhr")) %>%
    pivot_longer(cols = starts_with("r_"), names_to = "CarbonSource", names_pattern = "r_(.+)_", values_to = "GrowthRate") %>%
    filter(!CarbonSource %in% c("ketogluconate", "gluconate")) %>%
    mutate(ID = as.character(ID)) %>%
    mutate(CarbonSource = factor(CarbonSource, c("glucose", "acetate", "lactate", "succinate"))) %>%
    arrange(Fermenter, ID, CarbonSource) %>%
    mutate(FermenterID = rep(1:(n()/4), each = 4))
isolates_out <- temp %>%
    select(ID, Fermenter, ends_with("16hr") & starts_with("X")) %>%
    pivot_longer(cols = starts_with("X_"), names_to = "CarbonSource", names_pattern = "X_(.+)_", values_to = "Secretion") %>%
    filter(CarbonSource != "sum", !CarbonSource %in% c("ketogluconate", "gluconate")) %>%
    mutate(CarbonSource = factor(CarbonSource, c("glucose", "acetate", "lactate", "succinate"))) %>%
    mutate(ID = as.character(ID)) %>%
    arrange(Fermenter, ID, CarbonSource) %>%
    mutate(FermenterID = rep(1:(n()/3), each = 3))

## Make the big bipartite network
nodes <- bind_rows(
    tibble(Node = unique(isolates_in$ID)) %>% mutate(Type = "consumer"),
    tibble(Node = unique(isolates_in$CarbonSource)) %>% mutate(Type = "resource")
) %>%
    left_join(temp2)
edges <- bind_rows(
    isolates_in %>% mutate(from = CarbonSource, to = ID, Strength = GrowthRate, Direction = "consumed", .keep = "unused"),
    isolates_out %>% mutate(from = ID, to = CarbonSource, Strength = Secretion, Direction = "secreted", .keep = "unused")
) %>%
    mutate(BinaryStrength = ifelse(Strength == 0, 0, 1)) %>%
    select(from, to, everything())
## Subset the bipartite network for a community
example_comm <- "C1R2"
nodes_example <- nodes %>% filter(Community == example_comm | Type == "resource")
edges_example <- edges %>% filter(from %in% nodes_example$Node & to %in% nodes_example$Node) %>%
    mutate(from = match(from, nodes_example$Node), to = match(to, nodes_example$Node))

g <- tbl_graph(nodes_example, edges_example) %>%
    activate(nodes) %>%
    mutate(showResourceID = ifelse(Type == "resource", Node, ""),
           showNode = ifelse(Type == "resource", F, T))

node_size = 10
pA <- g %>%
    activate(nodes) %>%
    group_by(Type) %>%
    mutate(x = 1:n(), y = ifelse(Type == "consumer", 1, 0)) %>%
    activate(edges) %>%
    filter(Strength != 0) %>%
    ggraph(layout = "nicely") +
    geom_node_text(aes(label = showResourceID), size = node_size/2) +
    geom_edge_arc(aes(linetype = Direction), strength = 0.03, width = 1, start_cap = circle(node_size/2+1, "mm"), end_cap = circle(node_size/2+1, "mm")) +
    scale_x_continuous(limits = c(.5, 4.5)) +
    scale_y_continuous(limits = c(-.1, 1.1)) +
    theme_void() +
    #theme_bw() +
    theme(legend.title = element_blank(), legend.position = "top", axis.title = element_blank()) +
    #guides(width = "none") +
    labs() +
    draw_image(here::here("plots/cartoons/Fig1B_1.png"), x = 0.5, y = .5, vjust = 0, hjust = 0, clip = "on", scale = .5) +
    draw_image(here::here("plots/cartoons/Fig1B_2.png"), x = 1.5, y = .5, vjust = 0, hjust = 0, clip = "on", scale = .5) +
    draw_image(here::here("plots/cartoons/Fig1B_3.png"), x = 2.5, y = .5, vjust = 0, hjust = 0, clip = "on", scale = .5) +
    draw_image(here::here("plots/cartoons/Fig1B_4.png"), x = 3.5, y = .5, vjust = 0, hjust = 0, clip = "on", scale = .5) +
    paint_white_background()

ggsave(here::here("plots/Fig3A-example_crossfeeding.png"), pA, width = 5, height = 3)


# Figure 3B: d_r explains ranking
pB <- pairs_coexistence %>%
    drop_na(r_glucose_midhr_d) %>%
    ggplot(aes(x = abs(r_glucose_midhr_d), y = abs(Rank_d))) +
    geom_point(aes(color = InteractionType), size = 2, shape = 21, stroke = 1) +
    geom_smooth(method = "lm") +
    scale_color_manual(values = interaction_color) +
    theme_classic() +
    theme(legend.title = element_blank(), legend.position = "top") +
    labs(x = expression(r[1]-r[2]), y = expression(Rank[1]-Rank[2]))

# Figure 3C: d_X explains ranking
pC <- pairs_coexistence %>%
    drop_na(X_sum_16hr_d) %>%
    ggplot(aes(x = X_sum_16hr_d, y = Rank_d)) +
    geom_point(aes(color = InteractionType), size = 2, shape = 21, stroke = 1) +
    geom_smooth(method = "lm") +
    scale_color_manual(values = interaction_color) +
    theme_classic() +
    theme(legend.title = element_blank(), legend.position = "top") +
    labs(x = expression(X[1]-X[2]), y = expression(Rank[1]-Rank[2]))

# Figure 3D: together
if (FALSE) {

    pD <- pairs_coexistence %>%
        drop_na(PairFermenter, r_glucose_midhr_d, X_sum_16hr_d) %>%
        ggplot() +
        geom_vline(xintercept = 0, linetype = 2) +
        geom_hline(yintercept = 0, linetype = 2) +
        geom_point(aes(x = r_glucose_midhr_d, y = X_sum_16hr_d, color = InteractionType), shape = 21, size = 2, stroke = 1) +
        scale_color_manual(values = assign_interaction_color(level = "simple")) +
        theme_classic() +
        theme(legend.title = element_blank(), legend.position = "top",
              axis.text = element_text(color = 1),
              strip.background = element_blank(),
              panel.background = element_rect(color = 1, fill = NA)) +
        labs(x = expression(r[1]-r[2]), y = expression(X[1]-X[2]))
}
## Model prediction
pairs_fit <- pairs_coexistence %>%
    #mutate_if(is.character, as.factor) %>%
    #filter(PairFermenter == "FF") %>%
    filter(!is.na(InteractionType)) %>%
    mutate(InteractionType = ifelse(InteractionType == "coexistence", 1, 0)) %>%
    glm(formula = InteractionType ~  r_glucose_midhr_d + X_sum_16hr_d, data = ., family = "binomial") %>%
    broom::tidy()

pairs_model <- function (glu_d, X_d) {
    pairs_fit$estimate[pairs_fit$term == "(Intercept)"] +
        pairs_fit$estimate[pairs_fit$term == "r_glucose_midhr_d"] * glu_d +
        pairs_fit$estimate[pairs_fit$term == "X_sum_16hr_d"] * X_d
    #pairs_fit$estimate[pairs_fit$term == "r_glucose_midhr_d:X_sum_16hr_d"] * glu_d *X_d
}
x_range <- range(pairs_coexistence$r_glucose_midhr_d, na.rm = T) * 1.1
y_range <- range(pairs_coexistence$X_sum_16hr_d, na.rm = T) * 1.1
pairs_predicted <- tibble(x = seq(x_range[1], x_range[2], length.out = 100), y = seq(y_range[1], y_range[2], length.out = 100)) %>%
    tidyr::expand(x, y) %>%
    mutate(value = pairs_model(x, y))

## Figure
pD <- pairs_predicted %>%
    ggplot() +
    geom_tile(aes(x = x, y = y, fill = value)) +
    geom_vline(xintercept = 0, linetype = 2) +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_point(data = pairs_coexistence %>% drop_na(PairFermenter, r_glucose_midhr_d, X_sum_16hr_d), aes(x = r_glucose_midhr_d, y = X_sum_16hr_d, color = InteractionType), size = 2, shape = 21, stroke = 1) +
    scale_fill_gradient2(low = interaction_color["exclusion"], mid = "white", high = interaction_color["coexistence"]) +
    scale_color_manual(values = interaction_color) +
    scale_x_continuous(expand = c(0,0)) +
    scale_y_continuous(expand = c(0,0)) +
    theme_classic() +
    theme(legend.title = element_blank(), legend.position = "top") +
    guides(fill = "none") +
    labs(x = expression(r[1]-r[2]), y = expression(X[1]-X[2]))


if (FALSE) {
    # Figure 3E: simulation
    pE <- pairs_pool_meta %>%
        #filter(InteractionType != "no-growth") %>%
        ggplot() +
        geom_vline(xintercept = 0, linetype = 2) +
        geom_hline(yintercept = 0, linetype = 2) +
        geom_point(aes(x = d_ConsumptionRate, y = d_CrossFeedingPotential, color = InteractionType), shape = 21, size = 2, stroke = .5) +
        scale_color_manual(values = c(assign_interaction_color())) +
        # facet_grid(.~PairConspecific) +
        theme_classic() +
        theme(legend.position = "top", strip.background = element_blank(), panel.background = element_rect(color = 1)) +
        guides(color = "none") +
        labs()

}

p <- plot_grid(pA, pB, pC, pD, nrow = 1, labels = LETTERS[1:5], scale = c(.8, .9, .9, .9, .9),
               rel_widths = c(1.5,1,1,1), axis = "lrtb", align = "h") + paint_white_background()

ggsave(here::here("plots/Fig3.png"), p, width = 12, height = 3)








#====================================================================================================

# Figure 2S1: fermenter and respirator cartoon
pS1 <- ggdraw() + draw_image(here::here("plots/cartoons/Fig2A.png")) + theme(plot.background = element_rect(fill = "white", color = NA))
ggsave(here::here("plots/Fig2S1-functional_groups.png"), pS1, width = 5, height = 5)

# Figure 2S2: coexistence fraction in conspecific and heterospecific pairs
pS2 <- pairs_coexistence %>%
    group_by(InteractionType, PairConspecific) %>%
    summarize(Count = n()) %>%
    group_by(PairConspecific) %>%
    mutate(Fraction = Count / sum(Count)) %>%
    ggplot() +
    geom_col(aes(x = PairConspecific, y = Fraction, fill = factor(InteractionType, c("exclusion", "coexistence"))), color = 1, width = .7, position = "fill") +
    # sample size
    geom_text(data = pairs_count, aes(x = PairConspecific, y = 1, label = paste0("n=", n)), vjust = 2) +
    scale_fill_manual(values = assign_interaction_color(level = "simple")) +
    scale_y_continuous(breaks = c(0,.5,1), limit = c(0, 1), expand = c(0,0), labels = c("0%", "50%", "100%")) +
    theme_classic() +
    theme(axis.title.x = element_blank(), axis.title.y = element_text(size = 12),
          axis.ticks.x = element_blank(),
          axis.text = element_text(size = 12, color = 1),
          axis.text.x = element_text(size = 12, color = 1, angle = 30, vjust = 1, hjust = 1)
    ) +
    guides(fill = "none") +
    labs(y = "Fracrion of coexistence") +
    draw_image(here::here("plots/cartoons/Fig2B_FF.png"), x = 1, y = 0, scale = .6, vjust = .4, hjust = .5) +
    draw_image(here::here("plots/cartoons/Fig2B_FR.png"), x = 2, y = 0, scale = .6, vjust = .4, hjust = .5) +
    draw_image(here::here("plots/cartoons/Fig2B_RR.png"), x = 1, y = 0, scale = .6, vjust = .25, hjust = .5)

ggsave(here::here("plots/Fig2S2-pairs_alike.png"), pS2, width = 3, height = 3)


# Figure 2S3. Growth curve vs. pairwise outcomes
## Sylvie's growth curve
isolates_curves1 <- read_csv(here::here("data/output/isolates_curves1.csv")) %>%
    filter(CS == "glucose") %>%
    select(ID, Time, OD620) %>%
    # Substract the blank
    group_by(ID) %>%
    mutate(OD620Blank = min(OD620)) %>%
    mutate(OD620 = OD620 - OD620Blank, .keep = "unused")

#$ Jean's growth curves
isolates_curves2 <- read_csv(here::here("data/output/isolates_curves2.csv")) %>%
    filter(CS == "glucose") %>%
    select(ID, Time, OD620)

plot_pairwise_curve <- function(pairs_simple_row) {
    if (nrow(pairs_simple_row$GrowthCurve1[[1]]) != 0 & nrow(pairs_simple_row$GrowthCurve2[[1]]) != 0) {
        temp1 <- pairs_simple_row$GrowthCurve1[[1]] %>%
            mutate(Dominance = ifelse(pairs_simple_row$Dominant == "isolate1", T, F)) %>%
            mutate(Fermenter = pairs_simple_row$Fermenter1)
        temp2 <- pairs_simple_row$GrowthCurve2[[1]] %>%
            mutate(Dominance = ifelse(pairs_simple_row$Dominant == "isolate2", T, F)) %>%
            mutate(Fermenter = pairs_simple_row$Fermenter2)

        bind_rows(temp1, temp2) %>%
            ggplot() +
            geom_line(aes(x = Time, y = OD620, linetype = Dominance, color = Fermenter), size = 1) +
            #scale_y_log10() +
            scale_linetype_manual(values = c("TRUE" = 1, "FALSE" = 2), labels = c("TRUE" = "winner", "FALSE" = "loser")) +
            scale_color_manual(values = c("TRUE" = unname(category_color["fermenter"]), "FALSE" = unname(category_color["respirator"]))) +
            theme_classic() +
            theme(legend.position = c(.1, .5), legend.title = element_blank(),
                  legend.background = element_blank(),
                  axis.title = element_blank(),
                  plot.background = element_rect(size = 2, color = interaction_color[pairs_simple_row$InteractionType])) +
            guides(color = "none", linetype = "none") +
            ggtitle(paste0(pairs_simple_row$Community, "-", pairs_simple_row$ID1, "-", pairs_simple_row$ID2))
    } else {
        return(NULL)
    }
}
plot_pairwise_curve_group <- function(isolates_curves) {
    isolates_curves_ID <- isolates_curves %>% distinct(ID)

    pairs_simple <- pairs %>%
        filter(Assembly == "self_assembly") %>%
        mutate(Dominant = ifelse(From == Isolate1, "isolate1", "isolate2")) %>%
        select(Community, InteractionType, PairFermenter, Isolate1, Isolate2, ID1, ID2, Fermenter1, Fermenter2, Dominant) %>%
        rowwise() %>%
        mutate(GrowthCurve1 = isolates_curves %>% filter(ID == ID1) %>% list()) %>%
        mutate(GrowthCurve2 = isolates_curves %>% filter(ID == ID2) %>% list()) %>%
        drop_na(ID1, ID2, Fermenter1, Fermenter2) %>%
        filter(ID1 %in% isolates_curves_ID$ID, ID2 %in% isolates_curves_ID$ID) %>%
        arrange(InteractionType, PairFermenter)
    p_list <- rep(list(NA), nrow(pairs_simple))
    for (i in 1:nrow(pairs_simple)) p_list[[i]] <- plot_pairwise_curve(pairs_simple[i,]); print(i)
    p <- plot_grid(plotlist = p_list, nrow = 10, scale = .9) + paint_white_background()
    return(p)


}

p1 <- plot_pairwise_curve_group(isolates_curves1)
p2 <- plot_pairwise_curve_group(isolates_curves2)
ggsave(here::here("plots/Fig2S3A-pairwise_growthcurve_sylvie.png"), p1, width = 45, height = 25)
ggsave(here::here("plots/Fig2S3B-pairwise_growthcurve_jean.png"), p2, width = 45, height = 25)


# Figure 2S5. logistic regression on lasso-selected features
## All pairs
feature1 <- lasso1 %>%
    fit(pairs_train1) %>%
    tidy() %>%
    arrange(desc(abs(estimate))) %>%
    filter(term != "(Intercept)", estimate != 0)
p1 <- pairs_train1 %>%
    ggplot(aes(x = r_glucose_midhr_d, y = InteractionType)) +
    geom_point(shape = 21) +
    geom_smooth(method = "glm", method.args = list(family = "binomial")) +
    scale_y_continuous(breaks = c(0,1), labels = c("exclusion", "coexistence")) +
    theme_classic() +
    labs(y = "") +
    ggtitle("All pairs")

## FF pairs
feature2 <- lasso2 %>%
    fit(pairs_train2) %>%
    tidy() %>%
    arrange(desc(abs(estimate))) %>%
    filter(term != "(Intercept)", estimate != 0)

pairs_fit <- pairs_train2 %>%
    glm(formula = InteractionType ~  r_glucose_midhr_d + X_sum_28hr_d, data = ., family = "binomial") %>%
    broom::tidy()
pairs_model <- function (glu_d, X_d) {
    pairs_fit$estimate[pairs_fit$term == "(Intercept)"] +
        pairs_fit$estimate[pairs_fit$term == "r_glucose_midhr_d"] * glu_d +
        pairs_fit$estimate[pairs_fit$term == "X_sum_28hr_d"] * X_d
}
x_range <- range(pairs_coexistence$r_glucose_midhr_d, na.rm = T) * 1.1
y_range <- range(pairs_coexistence$X_sum_28hr_d, na.rm = T) * 1.1
pairs_predicted <- tibble(x = seq(x_range[1], x_range[2], length.out = 100), y = seq(y_range[1], y_range[2], length.out = 100)) %>%
    tidyr::expand(x, y) %>%
    mutate(value = pairs_model(x, y))

p2 <- pairs_predicted %>%
    ggplot() +
    geom_tile(aes(x = x, y = y, fill = value)) +
    geom_vline(xintercept = 0, linetype = 2) +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_point(data = pairs_train2 %>% drop_na(r_glucose_midhr_d, X_sum_28hr_d) %>% mutate(InteractionType = ifelse(InteractionType == 1, "coexistence", "exclusion")) ,
               aes(x = r_glucose_midhr_d, y = X_sum_28hr_d, color = InteractionType),
               size = 2, shape = 21, stroke = 1) +
    scale_fill_gradient2(low = interaction_color["exclusion"], mid = "white", high = interaction_color["coexistence"]) +
    scale_color_manual(values = interaction_color) +
    scale_x_continuous(expand = c(0,0)) +
    scale_y_continuous(expand = c(0,0)) +
    theme_classic() +
    theme(legend.title = element_blank(), legend.position = "right") +
    guides(fill = "none") +
    labs(x = "r_glucose_midhr_d" , y = "X_sum_28hr_d") +
    ggtitle("FF pairs")

pS4 <- plot_grid(p1, p2, nrow = 1, axis = "tbrl", align = "h", labels = c("A", "B"), rel_widths = c(1,1.2))
ggsave(here::here("plots/Fig2S5-regression.png"), pS4 , width = 8, height = 3)



# Figure 2S6. Ranked r vs. ranked
isolate_ranked <- isolates %>%
    filter(Assembly == "self_assembly") %>%
    select(Community, Fermenter, r_glucose_midhr, Rank) %>%
    group_by(Community) %>%
    mutate(r_ranked = rank(-r_glucose_midhr)) %>% # Put a negative sign such that the highest abs value is rank 1
    mutate(Fermenter = ifelse(Fermenter, "fermenter", ifelse(!Fermenter, "respirator", NA)))

isolate_ranked %>%
    lm(Rank ~ r_ranked, data = .) %>%
    tidy()

cor.test(isolate_ranked$Rank, isolate_ranked$r_ranked) %>%
    tidy()

pS6 <- isolate_ranked %>%
    ggplot() +
    geom_smooth(aes(x = r_ranked, y = Rank), formula = y ~ x, method = "lm", color = 1) +
    geom_point(aes(x = r_ranked, y = Rank, color = Fermenter), shape = 21, size = 2, stroke = 1, position = position_jitter(width = .1, height = .1)) +
    scale_x_continuous(breaks = 1:12) +
    scale_y_continuous(breaks = 1:12) +
    scale_color_manual(values = fermenter_color) +
    theme_classic() +
    theme(legend.position = "top", legend.title = element_blank()) +
    labs(x = "r_mid rank", y = "Isolate rank")

ggsave(here::here("plots/Fig2S6-r_rank.png"), pS6, width = 4, height = 4)


# Figure 2S7. r_glu_mid distribution
isolates %>%
    filter(Assembly == "self_assembly") %>%
    drop_na(r_glucose_midhr) %>%
    ggplot() +
    geom_histogram(aes(x = r_glucose_midhr), binwidth = 0.01, color = 1, fill = "white") +
    geom_vline(xintercept = 0.06, linetype = 2) +
    scale_y_continuous(breaks = scales::pretty_breaks(n=3)) +
    theme_classic()



#======================================================================================================





























# Figure 2S3: r_glu_midhr per isolate
pS3 <- pairs_coexistence %>%
    filter(PairConspecific == "conspecific") %>%
    filter(!is.na(r_glucose_midhr_d)) %>%
    mutate(InteractionType = factor(InteractionType, c("exclusion", "coexistence"))) %>%
    ggplot(aes(x = InteractionType, y = r_glucose_midhr_d, fill = InteractionType)) +
    geom_hline(yintercept = 0, linetype = 2, color = "grey") +
    geom_boxplot(width = .5) +
    geom_point(aes(group = InteractionType), shape = 1, size = 1, position = position_jitterdodge(jitter.width = 0.3)) +
    # p value
    ggpubr::stat_compare_means(aes(group = InteractionType), label = "p.format", vjust = 0, method = "t.test") +
    scale_fill_manual(values = assign_interaction_color()) +
    scale_y_continuous(expand = expansion(mult = .1, add = 0), limits = c(-0.1, 0.14), breaks = scales::pretty_breaks(n = 3)) +
    theme_classic() +
    theme(legend.title = element_blank(), legend.position = "top",
          axis.text = element_text(size = 12, color = 1),
          axis.title = element_text(size = 15, color = 1),
          axis.text.x = element_text(size = 12, color = 1, angle = 30, vjust = 1, hjust = 1),
          panel.border = element_rect(color = NA, fill = NA, size = 1)) +
    guides(alpha = "none", fill = "none", color = "none") +
    labs(x = "", y = expression(r[A]-r[B]))

ggsave(here::here("plots/Fig2S3-r_glu.png"), pS3, width = 4, height = 4)

## Stats: among conspecific, whether the dominant has a higher r_glu than subdominant
temp <- pairs_coexistence %>%
    filter(PairConspecific == "conspecific") %>%
    tidyr::drop_na(r_glucose_midhr_d) %>%
    select(InteractionType, r_glucose_midhr_d)
### One sample
temp %>%
    filter(InteractionType == "coexistence") %>%
    t_test(response = r_glucose_midhr_d)
temp %>%
    filter(InteractionType == "exclusion") %>%
    t_test(response = r_glucose_midhr_d)
### CI
x <- filter(temp, InteractionType == "exclusion") %>% pull(r_glucose_midhr_d)
t.test(x)
mean(x) + qt(0.975, length(x) - 1) * sd(x) / sqrt(length(x))
### Two sample
temp %>%
    t_test(r_glucose_midhr_d ~ InteractionType, order = c("coexistence", "exclusion"))



# Figure 2S4: Amount of total acid secretion. X_sum_16hr
pS4 <- pairs_coexistence %>%
    filter(PairConspecific == "conspecific") %>%
    drop_na(X_sum_16hr_d) %>%
    mutate(InteractionType = factor(InteractionType, c("exclusion", "coexistence"))) %>%
    ggplot(aes(x = InteractionType, y = X_sum_16hr_d, fill = InteractionType)) +
    geom_hline(yintercept = 0, linetype = 2, color = "grey") +
    geom_boxplot(width = 0.5) +
    geom_point(aes(group = InteractionType), shape = 1, size = 1, position = position_jitterdodge(jitter.width = 0.3)) +
    # p value
    ggpubr::stat_compare_means(aes(group = InteractionType), label = "p.format", vjust = 0, method = "t.test") +
    scale_fill_manual(values = assign_interaction_color()) +
    scale_y_continuous(expand = expansion(mult = .1, add = 0), limits = c(-20, 20), breaks = scales::pretty_breaks(n = 3)) +
    theme_classic() +
    theme(legend.title = element_blank(), legend.position = "top",
          axis.text = element_text(size = 12, color = 1),
          axis.title = element_text(size = 15, color = 1),
          axis.text.x = element_text(size = 12, color = 1, angle = 30, vjust = 1, hjust = 1),
          panel.border = element_rect(color = NA, fill = NA, size = 1)) +
    guides(alpha = "none", fill = "none", color = "none") +
    labs(x = "", y = expression(X[A]-X[B]))

ggsave(here::here("plots/Fig2S4-secretion_total.png"), pS4, width = 4, height = 4)

## Stats: whether the dominant has a higher X_sum than subdominant
temp <- pairs_coexistence %>%
    filter(PairConspecific == "conspecific") %>%
    tidyr::drop_na(X_sum_16hr_d) %>%
    select(InteractionType, X_sum_16hr_d)
### One sample
temp %>%
    filter(InteractionType == "coexistence") %>%
    t_test(response = X_sum_16hr_d)
temp %>%
    filter(InteractionType == "exclusion") %>%
    t_test(response = X_sum_16hr_d)
### two sample
temp %>% t_test(X_sum_16hr_d ~ InteractionType, order = c("coexistence", "exclusion"))



# Figure 2S5: growth curves, alpha by ranks
isolates_curves <- isolates %>%
    filter(str_detect(Community, "C\\d")) %>%
    mutate(Isolate = factor(Isolate)) %>%
    # Standardize rank
    group_by(Community) %>%
    mutate(Rank = Rank / n()) %>%
    left_join(read_csv(here::here("data/output/isolates_curves.csv"))) %>%
    filter(CS == "glucose") %>%
    group_by(ID, CS, Time) %>%
    arrange(ID, CS, Time)

## Replicate dummy variable
temp <- isolates_curves %>%
    distinct(Community, Isolate, Well) %>%
    group_by(ID, CS, Well) %>%
    summarize() %>%
    group_by(ID, CS) %>%
    mutate(Replicate = factor(1:n()))

pS5 <- isolates_curves %>%
    left_join(temp) %>%
    filter(Replicate == 1) %>%
    mutate(Community = factor(Community, communities$Community)) %>%
    drop_na(Fermenter) %>%
    mutate(Fermenter = ifelse(Fermenter, "fermenter", "respirator")) %>%
    ggplot() +
    geom_line(aes(x = Time, y = OD620, color = Fermenter, group = interaction(Isolate,Replicate)), size = 1) +
    scale_color_manual(values = fermenter_color) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 3)) +
    facet_wrap(Community~., ncol = 4) +
    guides(color = guide_legend(title = "")) +
    theme_classic() +
    theme(panel.border = element_rect(fill = NA, color = 1), legend.position = "top") +
    labs(x = "Time (hr)", y = "OD (620 nm)")

ggsave(here::here("plots/Fig2S5-growth_curves.png"), pS5, width = 8, height = 8)



# Figure 2S6: r_glu_midhr per isolate
## isolate r_glu
p1 <- isolates %>%
    filter(Assembly == "self_assembly") %>%
    drop_na(Fermenter, r_glucose_midhr) %>%
    mutate(Fermenter = ifelse(Fermenter, "fermenter", "respirator")) %>%
    ggplot(aes(x = Fermenter, y = r_glucose_midhr, fill = Fermenter)) +
    geom_boxplot() +
    geom_jitter(shape = 1, size = 1, stroke = 1) +
    ggpubr::stat_compare_means(comparisons = list(c("fermenter", "respirator")), label = "p.signif", method = "t.test", ) +
    scale_fill_manual(values = fermenter_color) +
    scale_y_continuous(limits = c(0, 0.23), breaks = scales::pretty_breaks(n = 3)) +
    theme_classic() +
    theme(axis.title.x = element_blank()) +
    guides(fill = "none") +
    labs(y = expression(r))

## interaction outcome as a function of  difference in r_glu
temp <- pairs_meta %>%
    filter(Assembly == "self_assembly") %>%
    drop_na(InteractionType, PairFermenter, r_glucose_midhr_d) %>%
    mutate(Coexistence = ifelse(InteractionType == "coexistence", 1, ifelse(InteractionType == "exclusion", 0, NA)))
p2 <- temp %>%
    ggplot(aes(x = r_glucose_midhr_d, y = Coexistence)) +
    geom_point(shape = 1, size = 2) +
    geom_smooth(formula = y ~ x, method = "glm") +
    ggpubr::stat_cor(label.x = -.1, label.y = 1.1, method = "pearson") +
    scale_y_continuous(breaks = c(0,1), labels = c("exclusion", "coexistence")) +
    theme_classic() +
    theme(axis.title.y = element_blank()) +
    labs(x = expression(r[A]-r[B]))
## Stats: correlation
cor.test(temp$r_glucose_midhr_d, temp$Coexistence) %>% broom::tidy()


## isolate r_glu in pairFermenter and InteractionType
p3 <- pairs_meta %>%
    filter(Assembly == "self_assembly") %>%
    select(PairFermenter, InteractionType, r_glucose_midhr1, r_glucose_midhr2) %>%
    pivot_longer(cols = starts_with("r_glucose_midhr"), names_pattern = "r_glucose_midhr(.)", names_to = "Isolate", values_to = "r_glucose_midhr") %>%
    drop_na(PairFermenter, r_glucose_midhr) %>%
    mutate(Isolate = ifelse(Isolate == 1, "dominant", "subdominant")) %>%
    ggplot(aes(x = InteractionType, y = r_glucose_midhr, fill = InteractionType, alpha = Isolate)) +
    geom_boxplot() +
    geom_point(aes(group = Isolate), alpha = 1, shape = 1, size = 1, stroke = 1, position = position_jitterdodge(jitter.width = 0.2)) +
    # p value
    ggpubr::stat_compare_means(aes(group = Isolate), label = "p.signif", method = "t.test") +
    scale_alpha_manual(values = c("dominant" = 1, "subdominant" = .2)) +
    scale_fill_manual(values = assign_interaction_color()) +
    scale_y_continuous(expand = expansion(mult = .1, add = 0)) +
    facet_grid(.~PairFermenter, labeller = labeller(PairFermenter = c(FF = "Fermenter-Fermenter", FN = "Fermenter-Respirator", NN = "Respirator-Respirator"))) +
    theme_classic() +
    theme(legend.title = element_blank(), legend.position = "right",
          strip.background = element_blank(),
          panel.border = element_rect(color = 1, fill = NA, size = .5)) +
    guides(fill = "none", alpha = "none") +
    labs(x = "", y = expression(r))

#
p4 <- pairs_meta %>%
    filter(Assembly == "self_assembly") %>%
    drop_na(PairFermenter, r_glucose_midhr_d) %>%
    ggplot(aes(x = InteractionType, y = r_glucose_midhr_d, fill = InteractionType)) +
    geom_hline(yintercept = 0, linetype = 2, color = "grey") +
    geom_boxplot(width = .5) +
    geom_point(aes(group = InteractionType), shape = 1, size = 1, stroke = 1, position = position_jitterdodge(jitter.width = 0.3)) +
    # p value
    ggpubr::stat_compare_means(aes(group = InteractionType), label = "p.signif", method = "t.test") +
    scale_fill_manual(values = assign_interaction_color()) +
    scale_y_continuous(expand = expansion(mult = .1, add = 0), breaks = scales::pretty_breaks(n = 3)) +
    facet_grid(.~PairFermenter, labeller = labeller(PairFermenter = c(FF = "Fermenter-Fermenter", FN = "Fermenter-Respirator", NN = "Respirator-Respirator"))) +
    theme_classic() +
    theme(legend.title = element_blank(), legend.position = "right",
          strip.background = element_blank(),
          panel.border = element_rect(color = 1, fill = NA, size = .5)) +
    guides(alpha = "none", fill = "none", color = "none") +
    labs(x = "", y = expression(r[A]-r[B]))

## Stats
### FF
pairs_meta %>%
    filter(Assembly == "self_assembly") %>%
    drop_na(PairFermenter, r_glucose_midhr_d) %>%
    filter(PairFermenter == "FN") %>%
    select(InteractionType, r_glucose_midhr_d) %>%
    t_test(r_glucose_midhr_d ~ InteractionType, order = c("coexistence", "exclusion"))
### FR
pairs_meta %>%
    filter(Assembly == "self_assembly") %>%
    drop_na(PairFermenter, r_glucose_midhr_d) %>%
    filter(PairFermenter == "FN") %>%
    select(InteractionType, r_glucose_midhr_d) %>%
    t_test(r_glucose_midhr_d ~ InteractionType, order = c("coexistence", "exclusion"))
### RR
pairs_meta %>%
    filter(Assembly == "self_assembly") %>%
    drop_na(PairFermenter, r_glucose_midhr_d) %>%
    filter(PairFermenter == "NN") %>%
    select(InteractionType, r_glucose_midhr_d) %>%
    t_test(r_glucose_midhr_d ~ InteractionType, order = c("coexistence", "exclusion"))


p_upper <- plot_grid(p1, p2, nrow = 1, labels = c(LETTERS[1:2], ""), rel_widths = c(1,1.5), axis = "tb", align = "h", scale = .9)
pS6 <- plot_grid(p_upper, p3, p4, ncol = 1, labels = c("", LETTERS[3:4]),
                 rel_heights = c(1,1,1)) + theme(plot.background = element_rect(fill = "white", color = NA))
ggsave(here::here("plots/Fig2S6-r_glu.png"), pS6, width = 6, height = 8)


## Stats: does difference in r_glu explain pairwise coexistence?
if (FALSE) {
    pairs_meta %>%
        mutate_if(is.character, as.factor) %>%
        filter(!is.na(InteractionType)) %>%
        mutate(InteractionType = ifelse(InteractionType == "coexistence", 1, 0)) %>%
        glm(formula = InteractionType ~  r_glucose_midhr_d * PairFermenter, data = ., family = "binomial") %>%
        broom::tidy() %>%
        {.}

}
temp <- pairs_meta %>%
    filter(Assembly == "self_assembly") %>%
    tidyr::drop_na(r_glucose_midhr_d) %>%
    oelect(PairFermenter, InteractionType, r_glucose_midhr_d)
### one-sample
temp %>%
    filter(PairFermenter == "FF") %>%
    filter(InteractionType == "coexistence") %>%
    t_test(response = r_glucose_midhr_d)
temp %>%
    filter(PairFermenter == "FF") %>%
    filter(InteractionType == "exclusion") %>%
    t_test(response = r_glucose_midhr_d)
temp %>%
    filter(PairFermenter == "FN") %>%
    t_test(response = r_glucose_midhr_d)
temp %>%
    filter(PairFermenter == "NN") %>%
    t_test(response = r_glucose_midhr_d)

### two-sample
temp %>%
    filter(PairFermenter == "FF") %>%
    t_test(r_glucose_midhr_d ~ InteractionType, order = c("coexistence", "exclusion"))
temp %>%
    filter(PairFermenter == "FN") %>%
    t_test(r_glucose_midhr_d ~ InteractionType, order = c("coexistence", "exclusion"))
temp %>%
    filter(PairFermenter == "NN") %>%
    t_test(r_glucose_midhr_d ~ InteractionType, order = c("coexistence", "exclusion"))


# Figure 2S7: rmid_glu vs. number of wins
temp <- isolates %>%
    filter(Assembly == "self_assembly") %>%
    drop_na(r_glucose_midhr)
## raw
p1 <- temp %>%
    ggplot(aes(x = r_glucose_midhr, y = Score)) +
    geom_point(shape = 1, size = 2) +
    geom_smooth(method = "lm") +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 3)) +
    theme_classic() +
    labs(x = expression(r[glu]), y = "Competitive score") +
    ggpubr::stat_cor(label.x = .1, label.y = -10, method = "pearson")
## standardized
p2 <- temp %>%
    ggplot(aes(x = r_glucose_midhr, y = Score/Game)) +
    geom_point(shape = 1, size = 2) +
    geom_smooth(method = "lm") +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 3)) +
    theme_classic() +
    labs(x = expression(r[glu]), y = "Standardized score") +
    ggpubr::stat_cor(label.x = .1, label.y = -1, method = "pearson")
## r1-r2 versus rank1 - rank2
p3 <- pairs_meta %>%
    mutate(Rank_d = Rank1 - Rank2) %>%
    ggplot(aes(x = r_glucose_midhr_d, y = Rank_d)) +
    geom_vline(xintercept = 0, linetype = 2, color = 1) +
    geom_point(aes(color = InteractionType), shape = 21, size = 2, stroke = 1) +
    geom_smooth(method = "lm") +
    scale_color_manual(values = assign_interaction_color()) +
    theme_classic() +
    theme(legend.position = "top", legend.title = element_blank()) +
    labs(x = expression(r[1]-r[2]), y = expression(Rank[1]-Rank[2]))

## rank1 - rank2 vs. r1 - r2
p4 <- pairs_meta %>%
    mutate(Rank_d = Rank1 - Rank2) %>%
    ggplot(aes(x = Rank_d, y = r_glucose_midhr_d)) +
    geom_vline(xintercept = 0, linetype = 2, color = 1) +
    geom_point(aes(color = InteractionType), shape = 21, size = 2, stroke = 1) +
    geom_smooth(method = "lm") +
    scale_color_manual(values = assign_interaction_color()) +
    #scale_x_continuous() +
    theme_classic() +
    theme(legend.position = "top", legend.title = element_blank()) +
    labs(x = expression(Rank[1]-Rank[2]), y = expression(r[1]-r[2]))



pS7 <- plot_grid(p1, p2, p3, p4, labels = LETTERS[1:4], nrow = 2, scale = .9) + theme(plot.background = element_rect(fill = "white", color = NA))
ggsave(here::here("plots/Fig2S7-r_glu_score.png"), pS7, width = 8, height = 8)

## Stats: correlation
cor.test(temp$Score, temp$r_glucose_midhr, alternative = "two.sided", method = "pearson") %>% broom::tidy()
cor.test(temp$Score/temp$Game, temp$r_glucose_midhr, alternative = "two.sided", method = "pearson") %>% broom::tidy()


# Figure 2S8: isolate byproduct acids measure
## Three acids
isolates_byproduct <- isolates %>%
    filter(Assembly == "self_assembly") %>%
    select(Community, Isolate, ID, Fermenter, ends_with("hr")) %>%
    pivot_longer(cols = ends_with("hr"), names_pattern = "(.*)_(.*)hr", names_to = c("Measure", "Time"), values_to = "Value") %>%
    filter(!is.na(Value))
p1 <- isolates_byproduct %>%
    # Only use acids
    filter(str_detect(Measure, "X_"), !str_detect(Measure, "sum")) %>%
    mutate(Measure = str_replace(Measure, "X_", "")) %>%
    ggplot(aes(x = Time, y = Value, group = ID, color = Fermenter, alpha = ID)) +
    geom_point() +
    geom_line() +
    scale_color_npg(labels = c("TRUE" = "Fermenter", "FALSE" = "Respirator"), breaks = c(T, F)) +
    facet_wrap(Measure~., nrow = 1, scales = "free_y") +
    theme_classic() +
    theme(legend.title = element_blank(), legend.position = "top",
          strip.background = element_rect(color = NA, fill = NA),
          panel.border = element_rect(fill = NA, color = 1)) +
    guides(alpha = "none") +
    labs(x = "", y = "Concentration (mM)")
p_upper <- get_legend(p1)
p1 <- p1 + theme(legend.position = "none")

## Total acids production
p2 <- isolates %>%
    select(ID, Fermenter, starts_with("X_sum")) %>%
    pivot_longer(cols = c(-ID, -Fermenter), names_to = c("Measure", "Time"), names_pattern = "(.*)_(.*)hr", values_to = "Value") %>%
    mutate(Measure = ifelse(Measure == "X_sum", "total acids", Measure)) %>%
    ggplot(aes(x = Time, y = Value, group = ID, color = Fermenter, alpha = ID)) +
    geom_point() +
    geom_line() +
    scale_color_npg(labels = c("TRUE" = "Fermenter", "FALSE" = "Respirator"), breaks = c(T, F)) +
    facet_wrap(Measure~., nrow = 1, scales = "free_y") +
    theme_classic() +
    theme(legend.title = element_blank(), legend.position = "none",
          strip.background = element_rect(color = NA, fill = NA),
          panel.border = element_rect(fill = NA, color = 1)) +
    guides(alpha = "none") +
    labs(x = "", y = "Concentration (mM)")

## Leakiness
p3 <- isolates %>%
    select(ID, Fermenter, starts_with("leakiness")) %>%
    pivot_longer(cols = c(-ID, -Fermenter), names_to = c("Measure", "Time"), names_pattern = "(.*)_(.*)hr", values_to = "Value") %>%
    ggplot(aes(x = Time, y = Value, group = ID, color = Fermenter, alpha = ID)) +
    geom_point() +
    geom_line() +
    scale_color_npg(labels = c("TRUE" = "Fermenter", "FALSE" = "Respirator"), breaks = c(T, F)) +
    facet_wrap(Measure~., nrow = 1, scales = "free_y") +
    theme_classic() +
    theme(legend.title = element_blank(), legend.position = "none",
          strip.background = element_rect(color = NA, fill = NA),
          panel.border = element_rect(fill = NA, color = 1)) +
    guides(alpha = "none") +
    labs(x = "Time (hr)", y = "leakiness")

## pH
p4 <- isolates_byproduct %>%
    filter(Measure == "pH") %>%
    ggplot(aes(x = Time, y = Value, group = ID, color = Fermenter, alpha = ID)) +
    geom_point() +
    geom_line() +
    scale_color_npg(labels = c("TRUE" = "Fermenter", "FALSE" = "Respirator"), breaks = c(T, F)) +
    facet_wrap(Measure~., nrow = 1, scales = "free_y") +
    theme_classic() +
    theme(legend.title = element_blank(), legend.position = "none",
          strip.background = element_rect(color = NA, fill = NA),
          panel.border = element_rect(fill = NA, color = 1)) +
    guides(alpha = "none") +
    labs(x = "", y = "pH")

p_lower <- plot_grid(p2, p3, p4, nrow = 1, axis = "tb", align = "h")
pS8 <- plot_grid(p_upper, p1, p_lower, ncol = 1, rel_heights = c(1, 5, 5)) + theme(plot.background = element_rect(fill = "white", color = NA))
ggsave(here::here("plots/Fig2S8-byproduct.png"), pS8, width = 9, height = 6)


# Figure 2S9: leakiness
temp <- isolates %>%
    filter(Assembly == "self_assembly") %>%
    drop_na(Fermenter, leakiness_16hr)
pS9 <- temp %>%
    ggplot() +
    geom_boxplot(aes(x = Fermenter, y = leakiness_16hr, color = Fermenter), outlier.size = 2) +
    geom_jitter(aes(x = Fermenter, y = leakiness_16hr, color = Fermenter), shape = 1, size = 2, width = 0.3) +
    ggpubr::stat_compare_means(aes(group = Fermenter, x = Fermenter, y = leakiness_16hr), method = "t.test") +
    scale_x_discrete(labels = c("TRUE" = "Fermenter", "FALSE" = "Respirator")) +
    scale_y_continuous(limits = c(0, 0.65)) +
    scale_color_npg() +
    theme_classic() +
    theme(axis.title.x = element_blank(), legend.position = "none") +
    labs(y = "Leakiness")
ggsave(here::here("plots/Fig2S9-isolate_leakiness.png"), pS9, width = 3, height = 3)

## Stats: are respirators less leakier than fermenters?
temp %>% t_test(leakiness_16hr ~ Fermenter, order = c("TRUE", "FALSE"))


# Figure 2S10: Amount of total acid secretion. X_sum_16hr
## isolate X_sum
p1 <- isolates %>%
    filter(Assembly == "self_assembly") %>%
    drop_na(Fermenter, X_sum_16hr) %>%
    mutate(Fermenter = ifelse(Fermenter, "fermenter", "respirator")) %>%
    ggplot(aes(x = Fermenter, y = X_sum_16hr, fill = Fermenter)) +
    geom_boxplot() +
    geom_jitter(shape = 1, size = 1, stroke = 1) +
    ggpubr::stat_compare_means(comparisons = list(c("fermenter", "respirator")), label = "p.signif", method = "t.test", ) +
    scale_fill_manual(values = fermenter_color) +
    scale_y_continuous(limits = c(0, 45), breaks = scales::pretty_breaks(n = 3)) +
    theme_classic() +
    theme(axis.title.x = element_blank()) +
    guides(fill = "none") +
    labs(y = expression(X))

## interaction outcome as a function of difference in X_sum
temp <- pairs_meta %>%
    filter(Assembly == "self_assembly") %>%
    drop_na(InteractionType, PairFermenter, X_sum_16hr_d) %>%
    mutate(Coexistence = ifelse(InteractionType == "coexistence", 1, ifelse(InteractionType == "exclusion", 0, NA)))
p2 <- temp %>%
    ggplot(aes(x = X_sum_16hr_d, y = Coexistence)) +
    geom_point(shape = 1, size = 2) +
    geom_smooth(formula = y ~ x, method = "glm") +
    ggpubr::stat_cor(label.x = -20, label.y = 1.1, method = "pearson") +
    scale_y_continuous(breaks = c(0,1), labels = c("exclusion", "coexistence")) +
    theme_classic() +
    theme(axis.title.y = element_blank()) +
    labs(x = expression(X[A]-X[B]))
## Stats: correlation
cor.test(temp$X_sum_16hr_d, temp$Coexistence) %>% broom::tidy()


## isolate r_glu in pairFermenter and InteractionType
p3 <- pairs_meta %>%
    filter(Assembly == "self_assembly") %>%
    select(PairFermenter, InteractionType, X_sum_16hr1, X_sum_16hr2) %>%
    pivot_longer(cols = starts_with("X_sum_16hr"), names_pattern = "X_sum_16hr(.)", names_to = "Isolate", values_to = "X_sum_16hr") %>%
    drop_na(PairFermenter, X_sum_16hr) %>%
    mutate(Isolate = ifelse(Isolate == 1, "dominant", "subdominant")) %>%
    ggplot(aes(x = InteractionType, y = X_sum_16hr, fill = InteractionType, alpha = Isolate)) +
    geom_boxplot() +
    geom_point(aes(group = Isolate), alpha = 1, shape = 1, size = 1, stroke = 1, position = position_jitterdodge(jitter.width = 0.2)) +
    # p value
    ggpubr::stat_compare_means(aes(group = Isolate), label = "p.signif", method = "t.test") +
    scale_alpha_manual(values = c("dominant" = 1, "subdominant" = .2)) +
    scale_fill_manual(values = assign_interaction_color()) +
    scale_y_continuous(expand = expansion(mult = .1, add = 0)) +
    facet_grid(.~PairFermenter, labeller = labeller(PairFermenter = c(FF = "Fermenter-Fermenter", FN = "Fermenter-Respirator", NN = "Respirator-Respirator"))) +
    theme_classic() +
    theme(legend.title = element_blank(), legend.position = "right",
          strip.background = element_blank(),
          panel.border = element_rect(color = 1, fill = NA, size = .5)) +
    guides(fill = "none") +
    labs(x = "", y = expression(X))

#
p4 <- pairs_meta %>%
    filter(Assembly == "self_assembly") %>%
    drop_na(PairFermenter, X_sum_16hr_d) %>%
    ggplot(aes(x = InteractionType, y = X_sum_16hr_d, fill = InteractionType)) +
    geom_hline(yintercept = 0, linetype = 2, color = "grey") +
    geom_boxplot(width = .5) +
    geom_point(aes(group = InteractionType), shape = 1, size = 1, stroke = 1, position = position_jitterdodge(jitter.width = 0.3)) +
    # p value
    ggpubr::stat_compare_means(aes(group = InteractionType), label = "p.signif", method = "t.test") +
    scale_fill_manual(values = assign_interaction_color()) +
    scale_y_continuous(expand = expansion(mult = .1, add = 0), breaks = scales::pretty_breaks(n = 3)) +
    facet_grid(.~PairFermenter, labeller = labeller(PairFermenter = c(FF = "Fermenter-Fermenter", FN = "Fermenter-Respirator", NN = "Respirator-Respirator"))) +
    theme_classic() +
    theme(legend.title = element_blank(), legend.position = "right",
          strip.background = element_blank(),
          panel.border = element_rect(color = 1, fill = NA, size = .5)) +
    guides(alpha = "none", fill = "none", color = "none") +
    labs(x = "", y = expression(X[A]-X[B]))

## Stats
### FF
pairs_meta %>%
    filter(Assembly == "self_assembly") %>%
    drop_na(PairFermenter, X_sum_16hr_d) %>%
    filter(PairFermenter == "FN") %>%
    select(InteractionType, X_sum_16hr_d) %>%
    t_test(X_sum_16hr_d ~ InteractionType, order = c("coexistence", "exclusion"))
### FR
pairs_meta %>%
    filter(Assembly == "self_assembly") %>%
    drop_na(PairFermenter, X_sum_16hr_d) %>%
    filter(PairFermenter == "FN") %>%
    select(InteractionType, X_sum_16hr_d) %>%
    t_test(X_sum_16hr_d ~ InteractionType, order = c("coexistence", "exclusion"))
### RR
pairs_meta %>%
    filter(Assembly == "self_assembly") %>%
    drop_na(PairFermenter, X_sum_16hr_d) %>%
    filter(PairFermenter == "NN") %>%
    select(InteractionType, X_sum_16hr_d) %>%
    t_test(X_sum_16hr_d ~ InteractionType, order = c("coexistence", "exclusion"))


p_upper <- plot_grid(p1, p2, nrow = 1, labels = c(LETTERS[1:2], ""), rel_widths = c(1,1.5), axis = "tb", align = "h", scale = .9)
pS10 <- plot_grid(p_upper, p3, p4, ncol = 1, labels = c("", LETTERS[3:4]),
                  rel_heights = c(1,1,1), axis = "lr", align = "v") + theme(plot.background = element_rect(fill = "white", color = NA))
ggsave(here::here("plots/Fig2S10-secretion_total.png"), pS10, width = 8, height = 8)


# Figure 2S11: Scatter r_glu_d vs. sum_acids_d
pS11 <- pairs_coexistence %>%
    drop_na(PairFermenter, r_glucose_midhr_d, X_sum_16hr_d) %>%
    ggplot() +
    geom_vline(xintercept = 0, linetype = 2) +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_point(aes(x = r_glucose_midhr_d, y = X_sum_16hr_d, color = InteractionType), shape = 21, size = 2, stroke = 1) +
    scale_color_manual(values = assign_interaction_color(level = "simple")) +
    facet_wrap(.~PairConspecific, scale = "free") +
    theme_classic() +
    theme(legend.title = element_blank(), legend.position = "top",
          axis.text = element_text(color = 1),
          strip.background = element_blank(),
          panel.background = element_rect(color = 1, fill = NA)) +
    labs(x = expression(r[A]-r[B]), y = expression(X[A]-X[B]))
ggsave(here::here("plots/Fig2S11-r_glu_secretion_d.png"), pS11, width = 5, height = 3)

# Stats: does difference in r_glu and X_sum explain pairwise coexistence?
pairs_meta %>%
    mutate_if(is.character, as.factor) %>%
    filter(PairFermenter == "FF") %>%
    filter(!is.na(InteractionType)) %>%
    mutate(InteractionType = ifelse(InteractionType == "coexistence", 1, 0)) %>%
    glm(formula = InteractionType ~  r_glucose_midhr_d * X_sum_16hr_d, data = ., family = "binomial") %>%
    broom::tidy() %>%
    {.}
#p_top <- plot_grid(p_A, p_B, nrow = 1, labels = LETTERS[1:2], rel_widths = c(1.5,1), scale = .9)
#p_bottom <- plot_grid(p_C, p_D, nrow = 1, labels = LETTERS[3:5], axis = "tb", align = "h", scale = .9)
#p <- plot_grid(p_top, p_bottom, nrow = 1, rel_widths = c(3,2), axis = "tb", align = "h") + theme(plot.background = element_rect(fill = "white", color = NA))



## Stats: does difference in X_sum explain pairwise coexistence?
pairs_meta %>%
    mutate_if(is.character, as.factor) %>%
    filter(!is.na(InteractionType)) %>%
    #filter(PairFermenter == "FF") %>%
    mutate(InteractionType = ifelse(InteractionType == "coexistence", 1, 0)) %>%
    glm(formula = InteractionType ~  X_sum_16hr_d * PairFermenter, data = ., family = "binomial") %>%
    broom::tidy() %>%
    {.}

## Stats: does difference in r_glu*X_sum explain pairwise coexistence?
pairs_meta %>%
    mutate_if(is.character, as.factor) %>%
    filter(!is.na(InteractionType)) %>%
    filter(PairFermenter == "FF") %>%
    mutate(InteractionType = ifelse(InteractionType == "coexistence", 1, 0)) %>%
    glm(formula = InteractionType ~  r_glucose_midhr_d*X_sum_16hr_d * PairFermenter, data = ., family = "binomial") %>%
    broom::tidy() %>%
    {.}


# Figure 2S12: Scatter r_glu_midhr vs. sum_acids
pS12 <- isolates %>%
    filter(!is.na(Fermenter)) %>%
    ggplot() +
    geom_segment(data = pairs_meta, aes(x = r_glucose_midhr1, xend = r_glucose_midhr2, y = X_sum_16hr1, yend = X_sum_16hr2, linetype = InteractionType)) +
    geom_point(aes(x = r_glucose_midhr, y = X_sum_16hr, color = Fermenter), stroke = 1.5, size = 2, shape = 1) +
    scale_color_npg(labels = c(`TRUE` = "fermenter", `FALSE` = "respirator")) +
    scale_linetype_manual(values = c("coexistence" = 1, "exclusion" = 2)) +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 3)) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 3)) +
    theme_classic() +
    theme(legend.title = element_blank(), legend.position = "right") +
    labs(x = expression(r[glu]), y = expression(2*acetate + 3*lactate + 4*succinate(mM)))
ggsave(here::here("plots/Fig2S12-r_glu_secretion.png"), pS12, width = 6, height = 4)




# Figure S13: r_glu vs. leakiness
pS13 <- isolates %>%
    filter(!is.na(Fermenter), !is.na(leakiness_16hr)) %>%
    ggplot(aes(x = r_glucose_midhr, y = leakiness_16hr, color = Fermenter)) +
    geom_point(shape = 1, size = 2) +
    geom_smooth(method = "lm") +
    scale_color_npg(labels = c("TRUE" = "Fermenter", "FALSE" = "Respirator")) +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 3)) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 3)) +
    theme_classic() +
    theme(legend.title = element_blank(), legend.position = "top") +
    labs(x = expression(r[glu]), y = "Leakiness")
ggsave(here::here("plots/Fig2S13-r_glu_leakiness.png"), pS13, width = 3, height = 3)

## Stats: do the fermenters follow pareto front?
isolates %>%
    filter(!is.na(Fermenter), !is.na(leakiness_16hr)) %>%
    glm(formula = leakiness_16hr ~  r_glucose_midhr * Fermenter, data = .) %>%
    broom::tidy() %>%
    {.}


# Figure 2S14: r_acids
pS14 <- pairs_meta %>%
    filter(!is.na(PairFermenter)) %>%
    mutate(Pair = 1:n()) %>%
    select(Pair, PairFermenter, InteractionType, starts_with("r_") & contains("midhr") & !ends_with("d")) %>%
    pivot_longer(cols = contains("hr"), names_to = c("CarbonSource", "Isolate"), names_pattern = "r_(.*)_midhr(.)", values_drop_na = T, values_to = "r") %>%
    ggplot(aes(x = InteractionType, y = r, fill = Isolate)) +
    geom_boxplot() +
    geom_point(shape = 1, size = 1, position = position_jitterdodge(jitter.width = 0.2)) +
    # p value
    ggpubr::stat_compare_means(aes(group = Isolate), label = "p.signif", vjust = -.5, method = "t.test") +
    scale_fill_npg(labels = c("dominant", "subdominant"), name = "Isolate") +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 4), expand = expansion(mult = .2, add = 0)) +
    facet_grid(CarbonSource~PairFermenter, scales = "free_y", labeller = labeller(PairFermenter = c(FF = "Fermenter-Fermenter", FN = "Fermenter-Respirator", NN = "Respirator-Respirator"))) +
    theme_classic() +
    theme(legend.title = element_blank(), legend.position = "right",
          panel.border = element_rect(color = 1, fill = NA, size = 1)) +
    labs(x = "", y = expression(r))
ggsave(here::here("plots/Fig2S14-r_acids.png"), pS14, width = 8, height = 12)



# Figure 2S15: matrix with species name
append_meta_data <- function(graph, community_name, isolates, pairs) {
    graph %>%
        activate(nodes) %>%
        mutate(Community = community_name) %>%
        left_join(isolates) %>%
        activate(edges) %>%
        mutate(Community = community_name, From = from, To = to) %>%
        left_join(pairs)
}
plot_adjacent_matrix <- function(graph) {
    graph_ranked <- graph %>%
        activate(nodes) %>%
        select(Isolate, PlotRank) %>%
        activate(edges) %>%
        mutate(
            fromRank = .N()$PlotRank[match(from, .N()$Isolate)],
            toRank = .N()$PlotRank[match(to, .N()$Isolate)])
    isolates_comm <- graph %>%
        igraph::vertex.attributes() %>%
        as_tibble()
    n_nodes <- nrow(isolates_comm)

    axis_label <- isolates_comm %>%
        #mutate(Discription = glue("{Genus}\twin = {Win}, draw = {Draw}, lose = {Lose}")) %>%
        #mutate(Discription = paste0(Genus, "\nwin = ", Win,", draw = ", Draw,", lose = ", Lose)) %>%
        mutate(Discription = Genus) %>%
        select(Isolate, PlotRank, Discription) %>%
        arrange(PlotRank) %>%
        pull(Discription)

    # Color lable fermenter and respirator
    color_logic <- isolates_comm %>%
        arrange(PlotRank) %>%
        pull(Fermenter) %>%
        ifelse("#C3423F", "#2978A0")

    graph_ranked %>%
        filter(fromRank <= toRank) %>%
        bind_edges(tibble(from = 1:n_nodes, to = 1:n_nodes, fromRank = 1:n_nodes, toRank = 1:n_nodes, InteractionType = "self")) %>%
        as_tibble() %>%
        ggplot() +
        geom_tile(aes(x = toRank, y = fromRank, fill = InteractionType), width = 0.9, height = 0.9) +
        scale_x_continuous(breaks = 1:n_nodes, labels = axis_label, expand = c(0,0), position = "top") +
        scale_y_reverse(breaks = 1:n_nodes, labels = axis_label, position = "right", expand = c(0,0)) +
        scale_fill_manual(values = assign_interaction_color(level = "matrix")) +
        theme_classic() +
        theme(legend.position = "none",
              axis.title = element_blank(),
              axis.line = element_blank(), axis.ticks = element_blank(),
              axis.text.x.top = element_text(size = 15, angle = 90, face = "bold.italic", color = color_logic),
              axis.text.y = element_text(size = 15, face = "bold.italic", color = color_logic),
              plot.margin = margin(10,10,10,10, "mm")) +
        labs()

}

communities_network <- communities %>%
    filter(str_detect(Community, "C\\d")) %>%
    mutate(graph = net_list[str_detect(names(net_list), "C\\d")]) %>%
    rowwise() %>%
    mutate(graph_meta = list(append_meta_data(graph, Community, isolates, pairs))) %>%
    mutate(p_net_matrix_list = list(plot_adjacent_matrix(graph_meta)))
pS15 <- plot_grid(plotlist = communities_network$p_net_matrix_list, nrow = 4, labels = communities_network$Community) +
    theme(plot.background = element_rect(fill = "white", color = NA))
ggsave(here::here("plots/Fig2S15-matrix_taxa.png"), pS15, width = 20, height = 20)


# Figure 2S16: isolate byproduct acids measure
isolates_byproduct <- isolates %>%
    filter(Assembly == "self_assembly") %>%
    select(Community, Isolate, ID, Fermenter, ends_with("hr")) %>%
    pivot_longer(cols = ends_with("hr"), names_pattern = "(.*)_(.*)hr", names_to = c("Measure", "Time"), values_to = "Value") %>%
    filter(!is.na(Value))
p1 <- isolates_byproduct %>%
    # Only use acids
    filter(str_detect(Measure, "X_"), !str_detect(Measure, "sum")) %>%
    mutate(Measure = str_replace(Measure, "X_", "")) %>%
    ggplot(aes(x = Time, y = Value, group = ID, color = Fermenter, alpha = ID)) +
    geom_point() +
    geom_line() +
    scale_color_npg(labels = c("TRUE" = "Fermenter", "FALSE" = "Respirator"), breaks = c(T, F)) +
    facet_wrap(Measure~., nrow = 1, scales = "free_y") +
    theme_classic() +
    theme(legend.title = element_blank(), legend.position = "top",
          strip.background = element_rect(color = NA, fill = NA),
          panel.border = element_rect(fill = NA, color = 1)) +
    guides(alpha = "none") +
    labs(x = "Time (hr)", y = "Concentration (mM)")
p_upper <- get_legend(p1)
p1 <- p1 + theme(legend.position = "none")

p2 <- isolates_byproduct %>%
    filter(Measure == "pH") %>%
    ggplot(aes(x = Time, y = Value, group = ID, color = Fermenter, alpha = ID)) +
    geom_point() +
    geom_line() +
    scale_color_npg(labels = c("TRUE" = "Fermenter", "FALSE" = "Respirator"), breaks = c(T, F)) +
    facet_wrap(Measure~., nrow = 1, scales = "free_y") +
    theme_classic() +
    theme(legend.title = element_blank(), legend.position = "none",
          strip.background = element_rect(color = NA, fill = NA),
          panel.border = element_rect(fill = NA, color = 1)) +
    guides(alpha = "none") +
    labs(x = "Time (hr)", y = "pH")

p_lower <- plot_grid(p1, p2, nrow = 1, rel_widths = c(3,1), axis = "tb", align = "h")
pS16 <- plot_grid(p_upper, p_lower, ncol = 1, rel_heights = c(1, 5)) + theme(plot.background = element_rect(fill = "white", color = NA))
ggsave(here::here("plots/Fig2S16-byproduct.png"), pS16, width = 9, height = 3)


# Figure 2S17: tree, all
"deprecated. fix it"
load(here::here("data/output/tree.Rdata"))
pS17 <- tree_meta %>%
    #ggtree(branch.length = "none", layout = "circular") +
    ggtree(branch.length = "none") +
    geom_tippoint(aes(colour = Fermenter), size = 3) +
    scale_color_npg(name = "", label = c("TRUE" = "fermenter", "FALSE" = "respirator")) +
    ggnewscale::new_scale_color() +
    geom_tiplab(aes(label = paste0(Community, " ", Genus, " sp.")), offset = 1) +
    scale_x_continuous(limits = c(0, 30)) +
    theme_void() +
    theme(plot.background = element_rect(fill = "white", color = NA)) +
    labs()
ggsave(here::here("plots/Fig2S17-tree.png"), pS17, width = 10, height = 10)

# Figure 2S18: tree for each community
plot_tree <- function(tree) {
    tree %>%
        #ggtree(branch.length = "none") +
        ggtree() +
        geom_tippoint(aes(colour = Fermenter), size = 5) +
        scale_color_npg(name = "", label = c("TRUE" = "Fermenter", "FALSE" = "Respirator")) +
        new_scale_color() +
        geom_tiplab(aes(label = paste0(Genus, " sp.")), offset = 0.05) +
        scale_x_continuous(expand = expansion(.5, .1)) +
        theme_void() +
        theme(legend.position = "none")
}
here::here("data/output/tree.Rdata")
## Get legend
l <- get_legend(p_L  + theme(legend.text = element_text(size = 20)))
pS18 <- plot_grid(plotlist = c(communities_tree$tree_plot, list(l)), labels = communities_tree$Community, ncol = 2) +
    theme(plot.background = element_rect(fill = "white", color = NA))
ggsave(here::here("plots/Fig2S18-tree_community.png"), pS18, width = 15, height = 15)


# Figure 2S19: growth trait differences in a community
isolates_plot <- isolates %>% filter(Assembly == "self_assembly") %>% filter(!is.na(Fermenter))
pairs_plot <- pairs_meta %>% filter(Assembly == "self_assembly")
pS19 <- isolates_plot %>%
    ggplot(aes(x = Rank, y = r_glucose_midhr)) +
    geom_segment(data = pairs_plot, aes(x = Rank1, xend = Rank2, y = r_glucose_midhr1, yend = r_glucose_midhr2, linetype = InteractionType)) +
    geom_point(aes(color = Fermenter), shape = 1, size = 3, stroke = 2) +
    scale_color_npg(name = "", label = c("TRUE" = "fermenter", "FALSE" = "respirator")) +
    scale_linetype_manual(name = "", values = c("coexistence" = 1, "exclusion" = 0)) +
    scale_x_continuous(breaks = 1:12) +
    facet_wrap(Community ~., scales = "free") +
    theme_classic() +
    theme(panel.border = element_rect(color = 1, fill = NA, size = 1)) +
    labs(x = "Rank", y = expression(r[glu]))
ggsave(here::here("plots/Fig2S19-r_glu_rank.png"), pS19, width = 10, height = 8)


# Figure 2S20: r_glu versus rank for C11R2
isolates_plot <- isolates %>% filter(Assembly == "self_assembly", Community == "C11R2") %>% filter(!is.na(Fermenter))
pairs_plot <- pairs_meta %>% filter(Assembly == "self_assembly", Community == "C11R2")
pS20 <- isolates_plot %>%
    ggplot(aes(x = Rank, y = r_glucose_midhr)) +
    geom_segment(data = pairs_plot, aes(x = Rank1, xend = Rank2, y = r_glucose_midhr1, yend = r_glucose_midhr2, linetype = InteractionType)) +
    geom_point(aes(color = Fermenter), shape = 1, size = 3, stroke = 2) +
    scale_color_npg(name = "", label = c("TRUE" = "fermenter", "FALSE" = "respirator")) +
    scale_linetype_manual(name = "", values = c("coexistence" = 1, "exclusion" = 0)) +
    scale_x_continuous(breaks = 1:12) +
    theme_classic() +
    theme(panel.border = element_rect(color = 1, fill = NA, size = 1)) +
    labs(x = "Rank", y = expression(r[glu]))
ggsave(here::here("plots/Fig2S20-r_glu_rank_C11R2.png"), pS20, width = 4, height = 3)


# Figure 2S21. Genus level metabolic traits
temp <- isolates %>%
    drop_na(Genus) %>%
    filter(Family %in% c("Enterobacteriaceae", "Pseudomonadaceae")) %>%
    arrange(Family) %>%
    mutate(Genus = factor(Genus, levels = unique(Genus)))

p1 <- temp %>%
    ggplot(aes(x = Genus, y = r_glucose_midhr, color = Family)) +
    geom_boxplot() +
    geom_jitter() +
    theme_classic() +
    labs()

p2 <- temp %>%
    ggplot(aes(x = Genus, y = X_sum_16hr, color = Family)) +
    geom_boxplot() +
    geom_jitter() +
    theme_classic() +
    labs()

pS21 <- plot_grid(p1, p2, ncol = 1, axis = "lr", align = "v", labels = c("A", "B"))
ggsave(here::here("plots/Fig2S21-isolate_trait.png"), pS21, width = 8, height = 6)


# Figure 2S22. Pairwise coexistence at genus resolution
temp <- pairs_meta %>%
    filter(Assembly == "self_assembly") %>%
    filter(PairFermenter == "FF") %>%
    filter(Family1 == Family2) %>%
    #select(InteractionType, Family1, Family2, PairFermenter, Genus1, Genus2) %>%
    drop_na(Genus1, Genus2) %>%
    mutate(InteractionType = factor(InteractionType, c("exclusion", "coexistence"))) %>%
    mutate(PairGenusConspecific = ifelse(Genus1 == Genus2, "conspecific", "heterospecific")) %>%
    unite(col = "PairGenus", Genus1, Genus2)

pS22 <- temp %>%
    group_by(PairGenusConspecific, InteractionType) %>%
    summarize(Count = n()) %>%
    mutate(TotalCount = sum(Count)) %>%
    ggplot() +
    #geom_col(aes(x = PairGenusConspecific, y = Count, fill = InteractionType)) +
    geom_col(aes(x = PairGenusConspecific, y = Count, fill = InteractionType), position = "fill") +
    geom_text(aes(x = PairGenusConspecific, label = paste0("n=", TotalCount)), y = 1, vjust = 2) +
    scale_fill_manual(values = assign_interaction_color()) +
    theme_classic() +
    theme(legend.position = "top", legend.title = element_blank()) +
    labs(x = "Within-fermenter pairs", y = "Fraction")

ggsave(here::here("plots/Fig2S22-isolate_trait.png"), pS22, width = 3, height = 3)


# Figure 2S23: cross-feeding networks
## Using the secretion profile from Sylvie
temp <- isolates %>%
    filter(Assembly == "self_assembly") %>%
    select(ID, Fermenter, starts_with("r_"), starts_with("X_")) %>%
    select(ID, Fermenter, ends_with("midhr"), ends_with("16hr")) %>%
    drop_na(r_glucose_16hr, X_acetate_16hr)
temp2 <- isolates %>%
    mutate(Node = as.character(ID), Fermenter = ifelse(Fermenter, "fermenter", ifelse(Fermenter == FALSE, "respirator", NA)), .keep = "used") %>%
    select(Node, Fermenter)

isolates_in <- temp %>%
    select(ID, Fermenter, ends_with("midhr")) %>%
    pivot_longer(cols = starts_with("r_"), names_to = "CarbonSource", names_pattern = "r_(.+)_", values_to = "GrowthRate") %>%
    filter(!CarbonSource %in% c("ketogluconate", "gluconate")) %>%
    mutate(ID = as.character(ID)) %>%
    mutate(CarbonSource = factor(CarbonSource, c("glucose", "acetate", "lactate", "succinate"))) %>%
    arrange(Fermenter, CarbonSource)
isolates_out <- temp %>%
    select(ID, Fermenter, ends_with("16hr") & starts_with("X")) %>%
    pivot_longer(cols = starts_with("X_"), names_to = "CarbonSource", names_pattern = "X_(.+)_", values_to = "Secretion") %>%
    filter(CarbonSource != "sum", !CarbonSource %in% c("ketogluconate", "gluconate")) %>%
    mutate(CarbonSource = factor(CarbonSource, c("glucose", "acetate", "lactate", "succinate"))) %>%
    mutate(ID = as.character(ID)) %>%
    arrange(Fermenter, CarbonSource)
nodes <- bind_rows(tibble(Node = unique(isolates_in$ID)) %>% mutate(Type = "consumer", y = 1, x = 1:n()),
                   tibble(Node = unique(isolates_in$CarbonSource)) %>% mutate(Type = "resource", y = 0, x = seq(10, 50, length = n()))) %>%
    left_join(temp2)
edges <- bind_rows(isolates_in %>% mutate(from = CarbonSource, to = ID, Strength = GrowthRate, Direction = "consumed", .keep = "unused"),
                   isolates_out %>% mutate(from = ID, to = CarbonSource, Strength = Secretion, Direction = "secreted", .keep = "unused")) %>%
    mutate(BinaryStrength = ifelse(Strength == 0, 0, 1))

p1 <- edges %>%
    ggplot() +
    geom_histogram(aes(x = Strength, fill = Direction), alpha = .5, position = "dodge") +
    facet_grid(.~Direction, scales = "free_x") +
    theme_classic() +
    theme(panel.border = element_rect(color = 1, fill = NA)) +
    guides(fill = "none") +
    labs()

g <- tbl_graph(nodes, edges) %>%
    activate(nodes) %>%
    mutate(showResourceID = ifelse(Type == "resource", Node, "o"))
node_size = 2
p2 <- g %>%
    activate(edges) %>%
    filter(Direction == "consumed") %>%
    filter(Strength != 0) %>%
    ggraph(layout = "nicely") +
    geom_node_text(aes(label = showResourceID, color = Fermenter), size = 5) +
    geom_edge_link(aes(color = Direction), width = .2, start_cap = circle(node_size/2+1, "mm"), end_cap = circle(node_size/2+1, "mm")) +
    annotate(geom = "text", label = "respirator", x = 10, y = 1, color = "#FFCB77", vjust = -1, size = 5) +
    annotate(geom = "text", label = "fermenter", x = 40, y = 1, color = "#8A89C0", vjust = -1, size = 5) +
    scale_color_manual(values = fermenter_color) +
    scale_fill_manual(values = c("consumer" = "grey50", "resource" = "grey90")) +
    scale_edge_color_manual(values = c("consumed" = "#557BAA", "secreted" = "#DB7469")) +
    scale_y_continuous(limits = c(0, 1.1)) +
    theme_void() +
    theme(legend.title = element_blank(), legend.position = "none") +
    guides(width = "none") +
    labs() +
    ggtitle("consumption network")

p3 <- g %>%
    activate(edges) %>%
    filter(Direction == "secreted") %>%
    filter(Strength != 0) %>%
    ggraph(layout = "nicely") +
    geom_node_text(aes(label = showResourceID, color = Fermenter), size = 5) +
    geom_edge_link(aes(color = Direction), width = .2, start_cap = circle(node_size/2+1, "mm"), end_cap = circle(node_size/2+1, "mm")) +
    annotate(geom = "text", label = "respirator", x = 10, y = 1, color = "#FFCB77", vjust = -1, size = 5) +
    annotate(geom = "text", label = "fermenter", x = 40, y = 1, color = "#8A89C0", vjust = -1, size = 5) +
    scale_color_manual(values = fermenter_color) +
    scale_fill_manual(values = c("consumer" = "grey50", "resource" = "grey90")) +
    scale_edge_color_manual(values = c("consumed" = "#557BAA", "secreted" = "#DB7469")) +
    scale_y_continuous(limits = c(0, 1.1)) +
    theme_void() +
    theme(legend.title = element_blank(), legend.position = "none") +
    guides(width = "none") +
    labs() +
    ggtitle("secretion network")


pS23 <- plot_grid(p1, p2, p3, ncol = 1, scale = .9, labels = c("A", "", "")) + paint_white_background()
ggsave(here::here("plots/Fig2S23-crossfeeding_network.png"), pS23, width = 8, height = 10)


# Figure 2S24: cross-feeding. Using the cross-feeding assay from Sylvie
gr <- read_csv("~/Dropbox/lab/invasion-network/data/raw/growth_rate/growth_od_processed_SE.csv")
cs <- read_csv("~/Dropbox/lab/invasion-network/data/raw/growth_rate/crossfeeding_processed_SE.csv")

gr %>%
    filter(gluConc == 0.2) %>%
    filter(str_detect(CommunityReplicate, "C\\d")) %>%
    filter(replicate == "R1") %>%
    filter(strain_type == "F") %>%
    ggplot(aes(x = timepoint, y = abs.mean, group = SangerID)) +
    geom_point() +
    geom_line() +
    facet_wrap(.~gluConc) +
    theme_classic() +
    labs()

gr_subset <- gr %>%
    filter(gluConc == 0.2) %>%
    filter(str_detect(CommunityReplicate, "C\\d")) %>%
    filter(replicate == "R1") %>%
    filter(strain_type == "F") %>%
    filter(timepoint == "T3") %>%
    select(CommunityReplicate, strain_type, SangerID, Genus, replicate, timepoint, abs.mean)

cs_subset <- cs %>%
    filter(gluConc == 0.2) %>%
    filter(str_detect(CommunityReplicate, "C\\d")) %>%
    filter(replicate == "R1") %>%
    select(CommunityReplicate, strain_type, SangerID, Genus, replicate, timepoint, abs.mean)

temp <- bind_rows(gr_subset, cs_subset) %>%
    pivot_wider(id_cols = CommunityReplicate, names_from = strain_type, values_from = abs.mean) %>%
    mutate(Dummy = "haha") %>%
    mutate(Community = factor(CommunityReplicate, communities$Community), .keep = "unused") %>%
    drop_na(Community)

pS24 <- temp %>%
    ggplot() +
    geom_tile(aes(x = Community, y = Dummy, fill = R)) +
    scale_fill_gradient(low = "black", high = "grey90", limits = c(0, 0.15)) +
    scale_y_discrete(expand = c(0,0)) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, vjust = .5),
          axis.text.y = element_blank(),
          axis.title = element_blank(),
          axis.line.y = element_blank(),
          axis.ticks.y = element_blank()) +
    guides(fill = guide_colourbar(title = "OD at 16hr")) +
    labs()

ggsave(here::here("plots/Fig2S24-crossfeeding_assay.png"), pS24, width = 4, height = 2)










# Figure 2S24: cross-feeding networks in heatmaps
temp <- isolates %>%
    filter(Assembly == "self_assembly") %>%
    select(ID, Fermenter, starts_with("r_"), starts_with("X_")) %>%
    select(ID, Fermenter, ends_with("midhr"), ends_with("16hr")) %>%
    drop_na(r_glucose_16hr, X_acetate_16hr)
temp2 <- isolates %>%
    mutate(Node = as.character(ID), Fermenter = ifelse(Fermenter, "fermenter", ifelse(Fermenter == FALSE, "respirator", NA)), .keep = "used") %>%
    select(Node, Fermenter)

isolates_in <- temp %>%
    select(ID, Fermenter, ends_with("midhr")) %>%
    pivot_longer(cols = starts_with("r_"), names_to = "CarbonSource", names_pattern = "r_(.+)_", values_to = "GrowthRate") %>%
    filter(!CarbonSource %in% c("ketogluconate", "gluconate")) %>%
    mutate(ID = as.character(ID)) %>%
    mutate(CarbonSource = factor(CarbonSource, c("glucose", "acetate", "lactate", "succinate"))) %>%
    arrange(Fermenter, ID, CarbonSource) %>%
    mutate(FermenterID = rep(1:(n()/4), each = 4))
isolates_out <- temp %>%
    select(ID, Fermenter, ends_with("16hr") & starts_with("X")) %>%
    pivot_longer(cols = starts_with("X_"), names_to = "CarbonSource", names_pattern = "X_(.+)_", values_to = "Secretion") %>%
    filter(CarbonSource != "sum", !CarbonSource %in% c("ketogluconate", "gluconate")) %>%
    mutate(CarbonSource = factor(CarbonSource, c("glucose", "acetate", "lactate", "succinate"))) %>%
    mutate(ID = as.character(ID)) %>%
    arrange(Fermenter, ID, CarbonSource) %>%
    mutate(FermenterID = rep(1:(n()/3), each = 3))

# Growth rate
p1 <- isolates_in %>%
    ggplot() +
    geom_tile(aes(x = CarbonSource, y = FermenterID, fill = GrowthRate), color = 1) +
    scale_fill_gradient2() +
    scale_y_continuous(limits = c(0, 60), breaks = 1:59) +
    theme_classic() +
    theme(legend.position = "top") +
    labs()

# Secretion
p2 <- isolates_out %>%
    ggplot() +
    geom_tile(aes(x = CarbonSource, y = FermenterID, fill = Secretion), color = 1) +
    scale_fill_gradient2() +
    scale_y_continuous(limits = c(0, 60), breaks = 1:59) +
    theme_classic() +
    theme(legend.position = "top", axis.title.y = element_blank()) +
    labs()

pS24 <- plot_grid(p1, p2, ncol = 2, scale = .9, labels = c("", "")) + paint_white_background()
ggsave(here::here("plots/Fig2S24-cross_feeding.png"), pS24, width = 8, height = 10)





}
