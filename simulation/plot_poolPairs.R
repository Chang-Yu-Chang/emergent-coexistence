library(tidyverse)

n_treatment = 3*4
n_rep = 20
# Generate input csv
input_poolPairs <- tibble(
    output_dir = "~/Dropbox/lab/invasion-network/simulation/data/raw/",
    exp_id = 1:(n_rep*n_treatment),
    seed = rep(1:n_rep, each = n_treatment),
    q = 0.9,
    q2 = c(0.1, 0.5, 0.9) %>% rep(4) %>% rep(n_rep),
    vamp = c(0.5, 1, 1.5, 2) %>% rep(each = 3) %>% rep(n_rep),
    S = 100,
    n_pairs = 100
)
write_csv(input_poolPairs, file = "input_poolPairs.csv")


# Read data
input_poolPairs <- read_csv("input_poolPairs.csv")
target_exp_id <- input_poolPairs %>% filter(seed == 1) %>% pull(exp_id)
pattern_exp_id <- paste(target_exp_id, collapse = "|") %>% paste0("(", . , ")")
#pattern_exp_id <- "\\d+"

# Subset the target exp_id
df_init <- read_simulation_data("~/Dropbox/lab/invasion-network/simulation/data/raw", paste0("poolPair_", pattern_exp_id, "_init")) %>% mutate(Time = 0)
df_end <- read_simulation_data("~/Dropbox/lab/invasion-network/simulation/data/raw", paste0("poolPair_", pattern_exp_id, "_end")) %>% mutate(Time = 1)

# Determine interactions
df_interaction <- bind_rows(df_init, df_end) %>%
    pivot_wider(names_from = Time, names_prefix = "T", values_from = Abundance) %>%
    replace_na(list(T1 = 0)) %>%
    group_by(Experiment, exp_id, Well) %>%
    summarize(InteractionType = ifelse(all(T1 != 0), "coexistence", ifelse(all(T1 == 0), NA, "exclusion")),
              PairFermenter = ifelse(all(Family == "F0"), "FF", ifelse(all(Family == "F1"), "NN", "FN"))) %>%
    filter(!is.na(InteractionType)) %>%
    left_join(select(input_poolPairs, exp_id, seed, vamp, q2))

# Check sample size
df_interaction %>%
    group_by(Experiment, exp_id, seed, vamp, q2) %>%
    summarize(Count = n())

# Plot
# Plot by vamp
p <- df_interaction %>%
    group_by(Experiment, seed, vamp, q2, PairFermenter, InteractionType) %>%
    summarize(Count = n()) %>%
    filter(seed == 1) %>%
    ggplot() +
    geom_col(aes(x = PairFermenter, y= Count, fill = InteractionType), color = 1) +
    #geom_col(aes(x = PairFermenter, y= Count, fill = InteractionType), position = "fill") +
    facet_grid(vamp~q2, scales = "free_y") +
    scale_fill_manual(values = interaction_color, breaks = c("coexistence", "exclusion")) +
    theme_classic() +
    labs(x = "", fill = "")
p
ggsave("plots/poolPairs_seed1.png", p, width = 10, height = 10)

df_summary <- df_interaction %>%
    left_join(select(input_csv, exp_id, vamp, seed)) %>%
    group_by(Experiment, seed, vamp, PairFermenter, InteractionType) %>%
    summarize(Count = n()) %>%
    group_by(Experiment, seed, vamp, PairFermenter) %>%
    arrange(Experiment, seed, vamp, PairFermenter) %>%
    mutate(Fraction = Count / sum(Count)) %>%
    filter(InteractionType == "coexistence") %>%
    group_by(Experiment, vamp, PairFermenter) %>%
    summarize(MeanFraction = mean(Fraction), SdFraction = sd(Fraction), SampleSize = n())
    {.}

df_summary %>%
    ggplot() +
    geom_col(aes(x = PairFermenter, y = MeanFraction)) +
    #geom_col(aes(x = PairFermenter, y= Count, fill = InteractionType), position = "fill") +
    facet_wrap(.~vamp) +
    scale_y_continuous(limits = c(0, 1)) +
    scale_fill_manual(values = interaction_color, breaks = c("coexistence", "exclusion")) +
    theme_classic() +
    labs(x = "", fill = "")





# Plot by q
df_interaction %>%
    left_join(select(input_csv, exp_id, q, q2, seed)) %>%
    group_by(Experiment, seed, q, q2, PairFermenter, InteractionType) %>%
    summarize(Count = n()) %>%
    #filter(seed == ) %>%
    ggplot() +
    geom_col(aes(x = PairFermenter, y= Count, fill = InteractionType)) +
    #geom_col(aes(x = PairFermenter, y= Count, fill = InteractionType), position = "fill") +
    facet_grid(q~q2, scales = "free_y") +
    scale_fill_manual(values = interaction_color, breaks = c("coexistence", "exclusion")) +
    theme_classic() +
    labs(x = "", fill = "")



#
df_interaction %>%
    filter(Seed %in% 1:5) %>%
    ggplot() +
    geom_bar(aes(x = PairFermenter, fill = InteractionType), position = "fill") +
    facet_grid(Seed~., scales = "free_y") +
    scale_fill_manual(values = interaction_color, breaks = c("coexistence", "exclusion")) +
    theme_classic() +
    labs(x = "", fill = "")




if (FALSE) {
read_csv("data/raw/1_end.csv", col_types = cols()) %>%
    select(Species = `...2`, starts_with("W")) %>%
    na_if(0) %>%
    pivot_longer(cols = -Species, names_to = "Well", values_to = "Abundance", values_drop_na = T)

}

