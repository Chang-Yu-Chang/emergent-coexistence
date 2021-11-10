library(tidyverse)
source("misc.R")

# Generate input csv
input_selfAssembly <- read_csv("input_selfAssembly.csv")
input_communityPairs <- input_selfAssembly %>%
    filter(q2 %in% c(0.1, 0.5, 0.9), vamp %in% c(0.5, 1, 1.5))
write_csv(input_communityPairs, file = "input_communityPairs.csv")

# Read data
input_communityPairs <- read_csv("input_communityPairs.csv")
target_exp_id <- input_communityPairs %>% pull(exp_id)
#pattern_exp_id <- paste(target_exp_id, collapse = "|") %>% paste0("communityPairs_(", . , ")")
pattern_exp_id <- "communityPairs_\\d+"

# Subset the target exp_id
df_init <- read_simulation_data("~/Dropbox/lab/invasion-network/simulation/data/raw", paste0(pattern_exp_id, "_init"), cs_output_format = F) %>% mutate(Time = 0)
df_end <- read_simulation_data("~/Dropbox/lab/invasion-network/simulation/data/raw", paste0(pattern_exp_id, "_end")) %>% mutate(Time = 1)


"
plot and find the decent q and q2 (fixed) and a range of vamp
plot the community pairs versus pool pairs
"


# Determine interactions
df_interaction <- bind_rows(df_init, df_end) %>%
    pivot_wider(names_from = Time, names_prefix = "T", values_from = Abundance) %>%
    replace_na(list(T1 = 0)) %>%
    group_by(Experiment, exp_id, Well) %>%
    summarize(InteractionType = ifelse(all(T1 != 0), "coexistence", ifelse(all(T1 == 0), NA, "exclusion")),
              PairFermenter = ifelse(all(Family == "F0"), "FF", ifelse(all(Family == "F1"), "NN", "FN"))) %>%
    filter(!is.na(InteractionType)) %>%
    left_join(select(input_communityPairs, exp_id, seed, vamp, q2))

# Check sample size
df_interaction %>%
    group_by(Experiment, exp_id, seed, vamp, q2) %>%
    summarize(Count = n())

#
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
ggsave("plots/communityPairs_seed1.png", p, width = 10, height = 10)














