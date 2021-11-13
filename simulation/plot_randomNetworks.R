library(tidyverse)
source("misc.R")

# Generate input csv
n_treatment = 3*3*5
n_rep = 1
input_randomNetworks <- tibble(
    output_dir = "~/Dropbox/lab/invasion-network/simulation/data/raw2/",
    exp_id = 1:(n_treatment*n_rep),
    seed = 1,
    q = 0.9,
    q2 = c(0.1, 0.5, 0.9) %>% rep(5) %>% rep(3) %>% rep(n_rep),
    l1 = c(0, 0.1, 0.5, 0.9, 1) %>% rep(3) %>% rep(each = 3) %>% rep(n_rep),
    l2 = 0.1,
    vamp = c(0.5, 1, 1.5) %>% rep(each = 3) %>% rep(each = 5) %>% rep(n_rep),
    S = 100,
    n_pairs = 100
)

write_csv(input_randomNetworks, file = "input_randomNetworks.csv")



# Read data
input_poolPairs <- read_csv("input_poolPairs.csv")
target_exp_id <- input_poolPairs %>%
    filter(seed == 1) %>%
    filter(q2 %in% c(0.1, 0.5, 0.9), vamp %in% c(0.5, 1, 1.5)) %>%
    pull(exp_id)
pattern_exp_id <- paste(target_exp_id, collapse = "|") %>% paste0("poolPair_(", . , ")")
#pattern_exp_id <- "poolPair_\\d+"
#pattern_exp_id <- "poolPair_1"

# Subset the target exp_id
df_init <- read_simulation_data(input_poolPairs$output_dir[1], paste0(pattern_exp_id, "_init")) %>% mutate(Time = "init")
df_end <- read_simulation_data(input_poolPairs$output_dir[1], paste0(pattern_exp_id, "_end")) %>% mutate(Time = "end")

# Isolate list
isolates <- df_init %>% distinct(Experiment, exp_id, Family, Species) %>% arrange(Experiment, exp_id, Family, Species)

# Compute frequency changes
df_freq <- bind_rows(df_init, df_end) %>%
    pivot_wider(names_from = Time, names_prefix = "T", values_from = Abundance) %>%
    arrange(Experiment, exp_id, Well) %>%
    group_by(Experiment, exp_id, Well) %>%
    mutate(Tend = Tend/sum(Tend, na.rm = T)) %>%
    replace_na(list(Tend = 0)) %>%
    mutate(Sign = ifelse(Tend > Tinit, 1, 0))
df_freq_w <- df_freq %>%
    select(Experiment, exp_id, Species, Well, Tinit, Tend) %>%
    group_by(Experiment, exp_id, Well) %>%
    mutate(temp = c("1", "2")) %>%
    pivot_wider(names_from = temp, values_from = c(Species, Tinit, Tend), names_sep = "")
df_sign <- df_freq %>%
    mutate(temp = c("1", "2")) %>%
    select(-Tinit, -Tend) %>%
    pivot_wider(names_from = temp, values_from = c(Species, Family, Sign), names_sep = "") %>%
    select(Experiment, exp_id, Species1, Species2, Well, Sign1)

# Determine interactions
df_interaction <- df_sign %>%
    group_by(Experiment, exp_id, Species1, Species2) %>%
    mutate(temp = c("Freq1", "Freq2")) %>%
    pivot_wider(id_cols = -Well, names_from = temp, values_from = Sign1) %>%
    mutate(InteractionType = ifelse(Freq1 == 1 & Freq2 == 1, "exclusion",
                                    ifelse(Freq1 == 0 & Freq2 == 0, "exclusion",
                                           ifelse(Freq1 == 1 & Freq2 == 0, "coexistence", NA))),
           From = ifelse(Freq1 == 1 & Freq2 == 1, "species1",
                         ifelse(Freq1 == 0 & Freq2 == 0, "species2", NA))) %>%
    select(-Freq1, -Freq2) %>%
    left_join(rename_with(isolates, ~ paste0(., "1"), !contains("Experiment") & !contains("exp_id"))) %>%
    left_join(rename_with(isolates, ~ paste0(., "2"), !contains("Experiment") & !contains("exp_id"))) %>%
    mutate(PairFermenter = ifelse(Family1 == "F0" & Family2 == "F0", "FF", ifelse(Family1 == "F1" & Family2 == "F1", "NN", "FN"))) %>%
    left_join(select(input_poolPairs, exp_id, seed, vamp, q2, l1))
df_interaction

"
check the sign of changes and if F0 is always the winner
"

# Check sample size
df_interaction %>%
    group_by(Experiment, exp_id, seed, vamp, q2, l1) %>%
    summarize(Count = n())

# Plot
# Plot by vamp and q2
l = c(0, 0.1, 0.5, 0.9, 1)
for (i in 1:5) {
    p <- df_interaction %>%
        group_by(Experiment, seed, vamp, q2, l1, PairFermenter, InteractionType) %>%
        summarize(Count = n()) %>%
        filter(l1 == l[i]) %>%
        ggplot() +
        geom_col(aes(x = PairFermenter, y= Count, fill = InteractionType), color = 1) +
        #geom_col(aes(x = PairFermenter, y= Count, fill = InteractionType), position = "fill") +
        facet_grid(vamp~q2, scales = "free_y", labeller = label_both) +
        scale_fill_manual(values = interaction_color, breaks = c("coexistence", "exclusion")) +
        theme_classic() +
        theme(legend.position = "top") +
        labs(x = "", fill = "") +
        ggtitle(paste0("l1 = ", l[i]))
    ggsave(paste0("plots/poolPairs_seed1_l1_", i, ".png"), p, width = 10, height = 10)
}
p
#ggsave("plots/poolPairs_seed1.png", p, width = 10, height = 10)



if (FALSE) {

df_summary <- df_interaction %>%
    left_join(select(input_poolPairs, exp_id, vamp, seed)) %>%
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

}











