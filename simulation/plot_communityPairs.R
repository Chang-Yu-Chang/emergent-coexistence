library(tidyverse)
source("misc.R")

# Generate input csv
input_selfAssembly <- read_csv("input_selfAssembly.csv")
input_communityPairs <- input_selfAssembly #%>% filter(q2 %in% c(0.1, 0.5, 0.9), vamp %in% c(0.5, 1, 1.5))
write_csv(input_communityPairs, file = "input_communityPairs.csv")

# Read data
input_communityPairs <- read_csv("input_communityPairs.csv")
target_exp_id <- input_communityPairs %>%
    filter(seed == 1) %>%
    pull(exp_id)
pattern_exp_id <- paste(target_exp_id, collapse = "|") %>% paste0("communityPairs_(", . , ")")
#pattern_exp_id <- "communityPairs_\\d+"

# Subset the target exp_id
df_init <- read_simulation_data(input_communityPairs$output_dir[1], paste0(pattern_exp_id, "_init"), cs_output_format = F) %>% mutate(Time = "init")
df_end <- read_simulation_data(input_communityPairs$output_dir[1], paste0(pattern_exp_id, "_end")) %>% mutate(Time = "end")

"
plot and find the decent q and q2 (fixed) and a range of vamp
plot the community pairs versus pool pairs
"
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

# Check sample size
df_interaction %>%
    group_by(Experiment, exp_id, seed, vamp, q2, l1) %>%
    summarize(Count = n())

# Plot the community pairs
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
    ggsave(paste0("plots/communityPairs_seed1_l1_", i, ".png"), p, width = 10, height = 10)
}
p
ggsave("plots/communityPairs_seed1.png", p, width = 10, height = 10)





# Match pairs to communities
target_exp_id
pattern_exp_id <- paste(target_exp_id, collapse = "|") %>% paste0("selfAssembly_(", . , ")")
df_end_comm <- read_simulation_data("~/Dropbox/lab/invasion-network/simulation/data/raw", paste0(pattern_exp_id, "_end")) %>% mutate(Time = 1)

df_end_comm %>%
    left_join(input_communityPairs) %>%
    ggplot() +
    geom_col(aes(x = Well, y = Abundance, fill = Family, alpha = Species), color = 1, position = "fill") +
    facet_grid(vamp~q2) +
    theme_classic() +
    guides(alpha = "none")













