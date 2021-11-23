library(tidyverse)
source("misc.R")

# Read data
input_randomNetworks <- read_csv("input_randomNetworks.csv")
target_exp_id <- input_randomNetworks %>%
    filter(seed == 1) %>%
    #filter(q2 %in% c(0.5, 0), vamp %in% c(0.5, 1, 1.5)) %>%
    #filter(vamp == 1, l1 == 0.5) %>%
    pull(exp_id)
pattern_exp_id <- paste(target_exp_id, collapse = "|") %>% paste0("randomNetworks_(", . , ")")

# Subset the target exp_id
df_init <- read_simulation_data(input_randomNetworks$output_dir[1], paste0(pattern_exp_id, "_init"), cs_output_format = F) %>% mutate(Time = "init")
df_end <- read_simulation_data(input_randomNetworks$output_dir[1], paste0(pattern_exp_id, "_end")) %>% mutate(Time = "end")
isolates <- df_init %>% distinct(Experiment, exp_id, Family, Species) %>% arrange(Experiment, exp_id, Family, Species)

# Compute frequency changes
df_freq <- bind_rows(df_init, df_end) %>%
    pivot_wider(names_from = Time, names_prefix = "T", values_from = Abundance) %>%
    mutate(Well = factor(Well, levels = paste0("W", 0:1000))) %>%
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
    select(Experiment, exp_id, Species1, Species2, Family1, Family2, Well, Sign1) %>%
    ungroup() %>%
    # Generate a unique pair id to sidestep the issue of duplicated pairs
    mutate(Pair = rep(1:(nrow(.)/2), each = 2))

# Determine interactions
df_interaction <- df_sign %>%
    group_by(Experiment, exp_id, Pair, Species1, Species2, Family1, Family2) %>%
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
    left_join(select(input_randomNetworks, exp_id, seed, vamp, q2, l1))
df_interaction

write_csv(df_interaction, file = "~/Dropbox/lab/invasion-network/simulation/data/temp/pairs_randomNetworks.csv")


if (FALSE) {
    df_init %>%
        mutate(Species = factor(Species, levels = paste0("S", 0:1000))) %>%
        mutate(Well = factor(Well, levels = paste0("W", 0:1000))) %>%
        arrange(exp_id, Well)

    # Isolate list
    isolates <- df_init %>% distinct(Experiment, exp_id, Family, Species) %>% arrange(Experiment, exp_id, Family, Species)

    # Compute frequency changes
    df_freq <- bind_rows(df_init, df_end) %>%
        pivot_wider(names_from = Time, names_prefix = "T", values_from = Abundance) %>%
        mutate(Well = factor(Well, levels = paste0("W", 0:1000))) %>%
        arrange(Experiment, exp_id, Well) %>%
        group_by(Experiment, exp_id, Well) %>%
        mutate(Tend = Tend/sum(Tend, na.rm = T)) %>%
        replace_na(list(Tend = 0)) %>%
        ungroup() %>%
        mutate(Pair = rep(1:(nrow(.)/4), each = 4)) %>%
        group_by(Experiment, exp_id, Pair) %>%
        # some pairs have 0.1/0.9 in well 1 and some have 0.9/0.1
        #mutate(temp = ifelse(Tinit == c(0.9, 0.1, 0.1, 0.9), c(2,1,2,1), ifelse(Tinit == c(0.1, 0.9, 0.9, 0.1), c(1,2,1,2)))) %>%
        mutate(Freq = ifelse(Tinit == c(0.9, 0.1, 0.1, 0.9), c("Freq2", "Freq2", "Freq1", "Freq1"), ifelse(Tinit == c(0.1, 0.9, 0.9, 0.1), c("Freq1", "Freq1", "Freq2", "Freq2")))) %>%
        mutate(Sign = ifelse(Tend > Tinit, 1, 0))
    # df_freq_w <- df_freq %>%
    #     select(Experiment, exp_id, Species, Well, Freq, Tinit, Tend) %>%
    #     pivot_wider(names_from = Freq, values_from = c(Species, Tinit, Tend), names_sep = "")
    df_sign <- df_freq %>%
        select(-Tinit, -Tend) %>%
        arrange(Pair, Freq) %>%
        group_by(Pair, Freq) %>%
        mutate(temp = c("1", "2")) %>%
        select(-Well) %>%
        pivot_wider(names_from = temp, values_from = c(Species, Family, Sign), names_sep = "") %>%
        select(Experiment, exp_id, Species1, Species2, Freq, Sign1)

    # Determine interactions
    df_interaction <- df_sign %>%
        group_by(Experiment, exp_id, Pair, Species1, Species2) %>%
        pivot_wider(names_from = Freq, values_from = Sign1) %>%
        mutate(InteractionType = ifelse(Freq1 == 1 & Freq2 == 1, "exclusion",
                                        ifelse(Freq1 == 0 & Freq2 == 0, "exclusion",
                                               ifelse(Freq1 == 1 & Freq2 == 0, "coexistence", NA))),
               From = ifelse(Freq1 == 1 & Freq2 == 1, "species1",
                             ifelse(Freq1 == 0 & Freq2 == 0, "species2", NA))) %>%
        select(-Freq1, -Freq2) %>%
        left_join(rename_with(isolates, ~ paste0(., "1"), !contains("Experiment") & !contains("exp_id"))) %>%
        left_join(rename_with(isolates, ~ paste0(., "2"), !contains("Experiment") & !contains("exp_id"))) %>%
        mutate(PairFermenter = ifelse(Family1 == "F0" & Family2 == "F0", "FF", ifelse(Family1 == "F1" & Family2 == "F1", "NN", "FN"))) %>%
        left_join(select(input_randomNetworks, exp_id, seed, vamp, q2, l1))
    df_interaction

    write_csv(df_interaction, file = "~/Dropbox/lab/invasion-network/simulation/data/temp/pairs_randomNetworks.csv")

}

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
        facet_grid(vamp~q2, scales = "free_y", labeller = label_both) +
        scale_fill_manual(values = interaction_color, breaks = c("coexistence", "exclusion")) +
        theme_classic() +
        theme(legend.position = "top") +
        labs(x = "", fill = "") +
        ggtitle(paste0("l1 = ", l[i]))
    ggsave(paste0("plots/randomNetworks_seed1_l1_", i, ".png"), p, width = 10, height = 10)
}



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

}











