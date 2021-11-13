# Plot the difference in growth rates on glucose and acetate
library(tidyverse)
library(cowplot)

# Data ----
isolates <- read_csv(here::here("data/output/isolates.csv"))
pairs_meta <- read_csv(here::here("data/output/pairs_meta.csv"))

if (FALSE) {
    isolates_preference <- read_csv(here::here("data/temp/isolates_preference.csv"))
    pairs <- read_csv(here::here("data/output/pairs.csv")); pairs$InteractionType[pairs$InteractionType == "neutrality"] <- "coexistence"
    pairs_freq <- read_csv(here::here("data/output/pairs_freq.csv"))
    # Growth rate data from Jean
    isolates_growth <- read_csv(here::here("data/raw/growth_rate/Growthcurver.csv"))
    # Growth rate data from Sylvie
    isolates_growth_syl <- read_csv(here::here("data/raw/growth_rate/Estrela_2021_isolates_grmax.csv"))
    # Byproduct measurement on glucose. Data from Sylvie
    isolates_byproduct <- read_csv(here::here("data/raw/growth_rate/By_Products_Glucose.csv")) %>%
        select(OD620_16h = OD620, ID = SangerID, Glucose_perc, acetate_mM, succinate_mM, lactate_mM, gluconate_mM, ketogluconate_mM)
    isolates_byproduct_time <- read_csv(here::here("data/raw/growth_rate/Estrela_2021_isolates_ph_OAs.csv")) %>%
        select(ID = SangerID, Time = time_hours, OD620, pH, Glucose_perc, acetate_mM, succinate_mM, lactate_mM)
    # OD data from Jean. Filter the 16hr data only
    isolates_OD_DW <- read_csv(here::here("data/raw/growth_rate/OD_Data_MH2.csv")) %>%
        filter(Time == 16) %>% select(ID = Sequence, CS = Carbon_Source, OD) %>%
        mutate(CS = sub("[LD]-", "", CS) %>% tolower() %>% paste0("OD620_16h_", .))


    #' 50-50 frequency changes. To determine the dominance in coexistence pairs.
    #' The isolate that increases frequency in the 50-50 competition is the dominant strain
    pairs_coexist_dominant <- pairs %>%
        filter(InteractionType == "coexistence") %>%
        select(Community, Isolate1, Isolate2) %>%
        left_join(pairs_freq) %>%
        filter(Isolate1InitialODFreq == 50) %>%
        select(Community, Isolate1, Isolate2, Time, Isolate1MeasuredFreq) %>%
        pivot_wider(names_from = "Time", values_from = "Isolate1MeasuredFreq") %>%
        mutate(Isolate1Dominant = T8 > T0) %>%
        select(Community, Isolate1, Isolate2, Isolate1Dominant, T0, T8)
    pairs <- left_join(pairs, pairs_coexist_dominant)

    # Swap so that fermenter is always isolate1
    for (i in 1:nrow(pairs)) {
        if (is.na(pairs$Fermenter1[i]) | is.na(pairs$Fermenter2[i])) next
        if (!pairs$Fermenter1[i] & pairs$Fermenter2[i]) {
            temp <- pairs$Isolate1[i]
            pairs$Isolate1[i] <- pairs$Isolate2[i]
            pairs$Isolate2[i] <- temp
            temp <- pairs$Fermenter1[i]
            pairs$Fermenter1[i] <- pairs$Fermenter2[i]
            pairs$Fermenter2[i] <- temp
            temp <- pairs$Family1[i]
            pairs$Family1[i] <- pairs$Family2[i]
            pairs$Family2[i] <- temp
            temp <- pairs$Genus1[i]
            pairs$Genus1[i] <- pairs$Genus2[i]
            pairs$Genus2[i] <- temp
        }
    }
    # Swap so that winner is also isolate1
    for (i in 1:nrow(pairs)) {
        #if (is.na(pairs$Fermenter1[i]) | is.na(pairs$Fermenter2[i])) next
        if (pairs$InteractionType[i] == "exclusion" & pairs$From[i] != pairs$Isolate1[i]) {
            temp <- pairs$Isolate1[i]
            pairs$Isolate1[i] <- pairs$Isolate2[i]
            pairs$Isolate2[i] <- temp
            temp <- pairs$Fermenter1[i]
            pairs$Fermenter1[i] <- pairs$Fermenter2[i]
            pairs$Fermenter2[i] <- temp
            temp <- pairs$Family1[i]
            pairs$Family1[i] <- pairs$Family2[i]
            pairs$Family2[i] <- temp
            temp <- pairs$Genus1[i]
            pairs$Genus1[i] <- pairs$Genus2[i]
            pairs$Genus2[i] <- temp
        }
    }
    # Swap so that in the coexistence pairs, isolate1 is the dominant strain
    for (i in 1:nrow(pairs)) {
        if (is.na(pairs$Isolate1Dominant[i])) next
        if (pairs$Isolate1Dominant[i] == FALSE) {
            temp <- pairs$Isolate1[i]
            pairs$Isolate1[i] <- pairs$Isolate2[i]
            pairs$Isolate2[i] <- temp
            temp <- pairs$Fermenter1[i]
            pairs$Fermenter1[i] <- pairs$Fermenter2[i]
            pairs$Fermenter2[i] <- temp
            temp <- pairs$Family1[i]
            pairs$Family1[i] <- pairs$Family2[i]
            pairs$Family2[i] <- temp
            temp <- pairs$Genus1[i]
            pairs$Genus1[i] <- pairs$Genus2[i]
            pairs$Genus2[i] <- temp
            temp <- pairs$From[i]
            pairs$From[i] <- pairs$To[i]
            pairs$To[i] <- temp
        }
    }

    # Jean's growth curve data. Use the fitted Rmid
    isolates_growth_w <- isolates_growth %>%
        separate(col = SID, sep = "_", into  = c("ID", "CS"), convert = T) %>%
        select(-SangerID, -Family) %>%
        select(ID, CS, RMid) %>%
        filter(CS %in% c("D-Glucose", "Acetate", "D-Lactate", "Succinate", "Gluconate", "2-Ketogluconate")) %>%
        pivot_wider(names_from = CS, values_from = RMid) %>%
        right_join(select(isolates, Community, Isolate, ID))

    pairs_meta <- pairs %>%
        left_join(rename_with(isolates_growth_w, ~ paste0(., "1"), !contains("Community"))) %>%
        left_join(rename_with(isolates_growth_w, ~ paste0(., "2"), !contains("Community"))) %>%
        # Difference in glucose, acetate, lactate, succinate growth rate
        mutate(r_glu_d = `D-Glucose1` - `D-Glucose2`, r_ace_d = Acetate1 - Acetate2,
               r_lac_d = `D-Lactate1` - `D-Lactate2`, r_suc_d = Succinate1 - Succinate2)

    # Isolate OD and glucose byproduct at 16hr. Sylvie's data
    pairs_meta <- pairs_meta %>%
        left_join(rename_with(isolates_byproduct, ~ paste0(., "1"), !contains("Community"))) %>%
        left_join(rename_with(isolates_byproduct, ~ paste0(., "2"), !contains("Community")))

    # Isolate's OD at 16 hr in DW96 plate. Jean's data
    isolates_OD_DW_w <- isolates_OD_DW %>%
        filter(!is.na(ID), ID != "N/A") %>%
        mutate(ID = as.numeric(ID)) %>%
        pivot_wider(names_from = CS, values_from = OD)
    pairs_meta <- pairs_meta %>%
        left_join(rename_with(isolates_OD_DW_w, ~ paste0(., "1"), !contains("Community"))) %>%
        left_join(rename_with(isolates_OD_DW_w, ~ paste0(., "2"), !contains("Community")))
    # Match the preferred CS
    pairs_meta <- pairs_meta %>%
        left_join(rename_with(isolates_preference, ~ paste0(., "1"))) %>%
        left_join(rename_with(isolates_preference, ~ paste0(., "2"))) %>%
        mutate(MatchedCS = PreferredCS1 == PreferredCS2)
    # isolate pH
    isolates_pH <- isolates_byproduct_time %>%
        select(ID, Time, pH) %>%
        pivot_wider(names_from = Time, values_from = pH, names_glue = "{Time}hr")
    pairs_meta <- pairs_meta %>%
        left_join(rename_with(isolates_pH, ~ paste0(., "1"))) %>%
        left_join(rename_with(isolates_pH, ~ paste0(., "2"))) %>%
        mutate(pH0hr_d = `0hr1` - `0hr2`, pH16hr_d = `16hr1` - `16hr2`,
               pH28hr_d = `28hr1` - `28hr2`, pH48hr_d = `48hr1` - `48hr2`)
    # Total amount of secretion
    isolates_byproduct_time_sum <- isolates_byproduct_time %>%
        group_by(ID, Time) %>%
        mutate(ByproductSum = sum(2 * acetate_mM, 4 * succinate_mM, 3* lactate_mM)) %>%
        select(ID, Time, ByproductSum) %>%
        pivot_wider(names_from = Time, values_from = ByproductSum, names_glue = "ByproductSum_{Time}hr")
    pairs_meta <- pairs_meta %>%
        left_join(rename_with(isolates_byproduct_time_sum, ~ paste0(., "1"), !contains("Time"))) %>%
        left_join(rename_with(isolates_byproduct_time_sum, ~ paste0(., "2"), !contains("Time"))) %>%
        mutate(byproduct_sum_16hr_d = ByproductSum_16hr1 - ByproductSum_16hr2,
               byproduct_sum_28hr_d = ByproductSum_28hr1 - ByproductSum_28hr2,
               byproduct_sum_48hr_d = ByproductSum_48hr1 - ByproductSum_48hr2)


}

# Plot ----
interaction_type <- c("exclusion", "coexistence")
interaction_color = c("#DB7469", "#557BAA")
names(interaction_color) <- interaction_type

# scatterplot r_glu vs. X_sum
isolates %>%
    ggplot(aes(x = r_glucose, y = X_sum_28hr, color = Fermenter)) +
    geom_point() +
    theme_classic()


# Scatterplot r_glu_d ----
## r_glu_d vs. r_ace_d
p1 <- pairs_meta %>%
    filter(!is.na(PairFermenter)) %>%
    filter(!is.na(r_glu_d)) %>%
    ggplot() +
    geom_abline(slope = -1, intercept = 0) +
    geom_point(aes(x = r_glu_d, y = r_ace_d, shape = InteractionType, color = PairFermenter), size = 2) +
    scale_shape_manual(values = c("coexistence" = 16, "exclusion" = 1)) +
    scale_color_discrete(labels = c(FF = "Fermenter-Fermenter", FN = "Fermenter-Respirator", NN = "Respirator-Respirator")) +
    theme_classic() +
    theme(legend.title = element_blank(), legend.position = "top", legend.direction = "horizontal") +
    guides(shape = "none") +
    labs(x = expression(r[A_glu] - r[B_glu]), y = expression(r[A_ace] - r[B_ace]))

## r_glu_d vs. r_lac_d
p1_2 <- pairs_meta %>%
    filter(!is.na(PairFermenter)) %>%
    filter(!is.na(r_glu_d)) %>%
    ggplot() +
    geom_abline(slope = -1, intercept = 0) +
    geom_point(aes(x = r_glu_d, y = r_lac_d, shape = InteractionType, color = PairFermenter), size = 2) +
    scale_shape_manual(values = c("coexistence" = 16, "exclusion" = 1)) +
    scale_color_discrete(labels = c(FF = "Fermenter-Fermenter", FN = "Fermenter-Respirator", NN = "Respirator-Respirator")) +
    theme_classic() +
    theme(legend.title = element_blank(), legend.position = "top", legend.direction = "horizontal") +
    guides(color = "none") +
    labs(x = expression(r[A_glu] - r[B_glu]), y = expression(r[A_lac] - r[B_lac]))

## r_glu_d vs. r_suc_d
p1_3 <- pairs_meta %>%
    filter(!is.na(PairFermenter)) %>%
    filter(!is.na(r_glu_d)) %>%
    ggplot() +
    geom_abline(slope = -1, intercept = 0) +
    geom_point(aes(x = r_glu_d, y = r_suc_d, shape = InteractionType, color = PairFermenter), size = 2) +
    scale_shape_manual(values = c("coexistence" = 16, "exclusion" = 1)) +
    scale_color_discrete(labels = c(FF = "Fermenter-Fermenter", FN = "Fermenter-Respirator", NN = "Respirator-Respirator")) +
    theme_classic() +
    theme(legend.title = element_blank(), legend.position = "none", legend.direction = "horizontal") +
    labs(x = expression(r[A_glu] - r[B_glu]), y = expression(r[A_suc] - r[B_suc]))

p <- plot_grid(p1, p1_2, p1_3, nrow = 1, align = "hv", axis = "tblr")
ggsave(here::here("plots/Fig_pairs-growth_rate_scatter.png"), p, width = 15, height = 5)

# Boxplot ----
## r_glu_d
p2 <- pairs_meta %>%
    filter(!is.na(PairFermenter)) %>%
    filter(!is.na(r_glu_d)) %>%
    ggplot(aes(x = InteractionType, y = r_glu_d, fill = InteractionType)) +
    geom_boxplot() +
    geom_jitter(shape = 1, size = 2) +
    geom_hline(yintercept = 0, color = "red") +
    scale_fill_manual(values = c("coexistence"=grey(0.5), "exclusion"=grey(0.9))) +
    facet_grid(.~PairFermenter, labeller = labeller(PairFermenter = c(FF = "Fermenter-Fermenter", FN = "Fermenter-Respirator", NN = "Respirator-Respirator"))) +
    theme_classic() +
    theme(legend.title = element_blank()) +
    labs(x = "", y = expression(r[A_glu] - r[B_glu]))

## r_ace_d
p2_2 <- pairs_meta %>%
    filter(!is.na(PairFermenter)) %>%
    filter(!is.na(r_glu_d)) %>%
    ggplot(aes(x = InteractionType, y = r_ace_d, fill = InteractionType)) +
    geom_boxplot() +
    geom_jitter(shape = 1, size = 2) +
    geom_hline(yintercept = 0, color = "red") +
    scale_fill_manual(values = c("coexistence"=grey(0.5), "exclusion"=grey(0.9))) +
    facet_grid(.~PairFermenter, labeller = labeller(PairFermenter = c(FF = "Fermenter-Fermenter", FN = "Fermenter-Respirator", NN = "Respirator-Respirator"))) +
    theme_classic() +
    theme(legend.title = element_blank()) +
    labs(x = "", y = expression(r[A_ace] - r[B_ace]))

## r_lac_d
p2_3 <- pairs_meta %>%
    filter(!is.na(PairFermenter)) %>%
    filter(!is.na(r_glu_d)) %>%
    ggplot(aes(x = InteractionType, y = r_lac_d, fill = InteractionType)) +
    geom_boxplot() +
    geom_jitter(shape = 1, size = 2) +
    geom_hline(yintercept = 0, color = "red") +
    scale_fill_manual(values = c("coexistence"=grey(0.5), "exclusion"=grey(0.9))) +
    facet_grid(.~PairFermenter, labeller = labeller(PairFermenter = c(FF = "Fermenter-Fermenter", FN = "Fermenter-Respirator", NN = "Respirator-Respirator"))) +
    theme_classic() +
    theme(legend.title = element_blank()) +
    labs(x = "", y = expression(r[A_lac] - r[B_lac]))

## r_suc_d
p2_4 <- pairs_meta %>%
    filter(!is.na(PairFermenter)) %>%
    filter(!is.na(r_glu_d)) %>%
    ggplot(aes(x = InteractionType, y = r_suc_d, fill = InteractionType)) +
    geom_boxplot() +
    geom_jitter(shape = 1, size = 2) +
    geom_hline(yintercept = 0, color = "red") +
    scale_fill_manual(values = c("coexistence"=grey(0.5), "exclusion"=grey(0.9))) +
    facet_grid(.~PairFermenter, labeller = labeller(PairFermenter = c(FF = "Fermenter-Fermenter", FN = "Fermenter-Respirator", NN = "Respirator-Respirator"))) +
    theme_classic() +
    theme(legend.title = element_blank()) +
    labs(x = "", y = expression(r[A_suc] - r[B_suc]))

p <- plot_grid(p2, p2_2, p2_3, p2_4, ncol = 1)
ggsave(here::here("plots/Fig_pairs-growth_rate_boxplot.png"), p, width = 8, height = 15)

# r_glu
p <- pairs_meta %>%
    filter(!is.na(PairFermenter)) %>%
    select(PairFermenter, InteractionType, `D-Glucose1`, `D-Glucose2`) %>%
    pivot_longer(cols = starts_with("D-"), names_to = "Isolate", values_to = "r_glu") %>%
    ggplot(aes(x = InteractionType, y = r_glu, color = Isolate)) +
    geom_boxplot() +
    geom_point(aes(color = Isolate), shape = 1, size = 2, position = position_jitterdodge(jitter.width = 0.2)) +
    scale_color_discrete(labels = c("A", "B"), name = "Isolate") +
    #scale_fill_manual(values = c("coexistence"=grey(0.5), "exclusion"=grey(0.9))) +
    facet_grid(.~PairFermenter, scales = "free_y",
               labeller = labeller(PairFermenter = c(FF = "Fermenter-Fermenter", FN = "Fermenter-Respirator", NN = "Respirator-Respirator"))) +
    theme_classic() +
    theme(legend.title = element_blank()) +
    labs(x = "", y = expression(r[glu]))
p
ggsave(here::here("plots/Fig_pairs-growth_rate_boxplot2.png"), width = 8, height = 3)


## OD620_16h_d, on glucose ----
# From Sylvie's deep well data
pairs_meta %>%
    filter(!is.na(PairFermenter)) %>%
    mutate(OD620_16h_d = OD620_16h1 - OD620_16h2) %>%
    ggplot(aes(x = InteractionType, y = OD620_16h_d, fill = InteractionType)) +
    geom_boxplot() +
    geom_jitter(shape = 1, size = 2) +
    geom_hline(yintercept = 0, color = "red") +
    scale_fill_manual(values = c("coexistence"=grey(0.5), "exclusion"=grey(0.9))) +
    facet_grid(.~PairFermenter, labeller = labeller(PairFermenter = c(FF = "Fermenter-Fermenter", FN = "Fermenter-Respirator", NN = "Respirator-Respirator"))) +
    theme_classic() +
    theme(legend.title = element_blank()) +
    labs(x = "", y = expression(OD[A_16h] - OD[B_16h]))

## OD620_16_d on glucose. From Jean's data
pairs_meta %>%
    filter(!is.na(PairFermenter)) %>%
    mutate(OD620_16h_glu_d = OD620_16h_glucose1 - OD620_16h_glucose2) %>%
    ggplot(aes(x = InteractionType, y = OD620_16h_glu_d, fill = InteractionType)) +
    geom_boxplot() +
    geom_jitter(shape = 1, size = 2) +
    geom_hline(yintercept = 0, color = "red") +
    scale_fill_manual(values = c("coexistence"=grey(0.5), "exclusion"=grey(0.9))) +
    facet_grid(.~PairFermenter, labeller = labeller(PairFermenter = c(FF = "Fermenter-Fermenter", FN = "Fermenter-Respirator", NN = "Respirator-Respirator"))) +
    theme_classic() +
    theme(legend.title = element_blank()) +
    labs(x = "", y = expression(OD[A_glu] - OD[B_glu]))

## OD620_16_d on succinate
pairs_meta %>%
    filter(!is.na(PairFermenter)) %>%
    mutate(OD620_16h_suc_d = OD620_16h_succinate1 - OD620_16h_succinate2) %>%
    ggplot(aes(x = InteractionType, y = OD620_16h_suc_d, fill = InteractionType)) +
    geom_boxplot() +
    geom_jitter(shape = 1, size = 2) +
    geom_hline(yintercept = 0, color = "red") +
    scale_fill_manual(values = c("coexistence"=grey(0.5), "exclusion"=grey(0.9))) +
    facet_grid(.~PairFermenter, labeller = labeller(PairFermenter = c(FF = "Fermenter-Fermenter", FN = "Fermenter-Respirator", NN = "Respirator-Respirator"))) +
    theme_classic() +
    theme(legend.title = element_blank()) +
    labs(x = "", y = expression(OD[A_suc] - OD[B_suc]))


# Barplot r_glu_d ----
## Whether the winner is the faster grower overall.
p3 <- pairs_meta %>%
    filter(!is.na(PairFermenter)) %>%
    filter(!is.na(r_glu_d)) %>%
    filter(InteractionType == "exclusion") %>%
    group_by(InteractionType, FastWinner = r_glu_d > 0) %>%
    ggplot() +
    geom_bar(aes(x = InteractionType, fill = FastWinner), color = 1) +
    scale_fill_discrete(labels = c("TRUE"="winner is faster", "FALSE" = "winner is slower")) +
    scale_y_continuous(limits = c(0, 100)) +
    theme_classic() +
    theme(legend.position = "none") +
    labs(x = "", fill = "growth rate on glucose")
## Whether winner is the faster grower, in FF, FN, NN groups
p4 <- pairs_meta %>%
    filter(!is.na(PairFermenter)) %>%
    filter(!is.na(r_glu_d)) %>%
    filter(InteractionType == "exclusion") %>%
    group_by(InteractionType, PairFermenter, FastWinner = r_glu_d > 0) %>%
    ggplot() +
    geom_bar(aes(x = PairFermenter, fill = FastWinner), color = 1) +
    scale_fill_discrete(labels = c("TRUE"="winner is faster", "FALSE" = "winner is slower")) +
    scale_x_discrete(labels = c(FF = "Fermenter-Fermenter", FN = "Fermenter-Respirator", NN = "Respirator-Respirator")) +
    scale_y_continuous(limits = c(0, 100)) +
    theme_classic() +
    theme(legend.position = "top") +
    labs(x = "", fill = "growth rate on glucose")

p <- plot_grid(p3, p4, align = "hv", axis = "tblr", nrow = 1, rel_widths = c(1, 2))
ggsave(here::here("plots/Fig_pairs-growth_rate_barplot.png"), p, width = 8, height = 5)



# Barplot matched CS ----
## Matched preferred CS
p5 <- pairs_meta %>%
    filter(!is.na(PairFermenter)) %>%
    ungroup() %>%
    group_by(InteractionType, PairFermenter, MatchedCS) %>%
    summarize(Count = n()) %>%
    ggplot() +
    geom_bar(aes(x = InteractionType, y = Count , fill = MatchedCS), color = 1, stat = "identity") +
    facet_grid(.~PairFermenter, labeller = labeller(PairFermenter = c(FF = "Fermenter-Fermenter", FN = "Fermenter-Respirator", NN = "Respirator-Respirator"))) +
    theme_classic() +
    labs(x = "")
p5_2 <- pairs_meta %>%
    filter(!is.na(PairFermenter)) %>%
    ungroup() %>%
    group_by(InteractionType, PairFermenter, MatchedCS) %>%
    summarize(Count = n()) %>%
    mutate(Fration = Count / sum(Count)) %>%
    ggplot() +
    geom_bar(aes(x = InteractionType, y = Fration , fill = MatchedCS), color = 1, stat = "identity") +
    facet_grid(.~PairFermenter, labeller = labeller(PairFermenter = c(FF = "Fermenter-Fermenter", FN = "Fermenter-Respirator", NN = "Respirator-Respirator"))) +
    theme_classic() +
    labs(x = "")

x <- pairs_meta %>%
    ungroup() %>%
    group_by(InteractionType, MatchedCS, PairFermenter) %>%
    summarize(Count = n()) %>%
    filter(PairFermenter == "FF")
#chisq.test(matrix(c(21, 19, 64, 39), 2, 2)) %>% broom::tidy()
chisq.test(matrix(c(14, 9, 11, 29), 2, 2)) %>% broom::tidy()

ggsave(here::here("plots/Fig_pairs-matches_barplot.png"), p5, width = 8, height = 5)
ggsave(here::here("plots/Fig_pairs-matches_barplot_fraction.png"), p5_2, width = 8, height = 5)


# Dissimilarity using the growth rates on different acids ----
isolates_growth_w <- isolates_growth %>%
    separate(col = SID, sep = "_", into  = c("ID", "CS"), convert = T) %>%
    select(-SangerID, -Family) %>%
    select(ID, CS, RMid) %>%
    #filter(CS %in% c("D-Glucose", "Acetate", "D-Lactate", "Succinate", "Gluconate", "2-Ketogluconate")) %>%
    #filter(CS %in% c("D-Glucose", "D-Lactate", "Succinate")) %>%
    filter(CS %in% c("D-Glucose", "Acetate", "D-Lactate", "Succinate")) %>%
    pivot_wider(names_from = CS, values_from = RMid)
ID <- isolates_growth_w$ID

## Sign determined by the overall sum of the growth rates on acids
isolates_growth_total <- isolates_growth_w %>%
    # Exclude glucose
    select(-`D-Glucose`) %>%
    pivot_longer(-ID) %>%
    group_by(ID) %>%
    summarize(r_total = sum(value))

pairs_temp <- combn(ID, 2) %>% t() %>% as_tibble(.name_repair = "minimal") %>% setNames(c("ID1", "ID2"))
pairs_r_total <- pairs_temp %>%
    bind_rows(setNames(pairs_temp, c("ID2", "ID1"))) %>%
    left_join(rename_with(isolates_growth_total, ~ paste0(., "1"))) %>%
    left_join(rename_with(isolates_growth_total, ~ paste0(., "2"))) %>%
    mutate(ID1_has_higher_r_total = r_total1 > r_total2, r_total_d = r_total1 - r_total2) %>%
    {.}

## Distance between the two IDs
isolates_growth_dis <- isolates_growth_w %>%
    # Exclude glucose
    select(-ID, -`D-Glucose`) %>%
    cluster::daisy(metric="euclidean") %>%
    as.matrix() %>% as_tibble() %>%
    mutate(ID1 = ID) %>%
    pivot_longer(cols = -ID1, names_to = "ID2", names_transform = list(ID2 = as.integer), values_to = "r_dissim") %>%
    mutate(ID2 = ID[ID2]) %>%
    # Add signs
    left_join(pairs_r_total) %>%
    filter(!is.na(r_total1)) %>%
    mutate(r_dissim = ifelse(ID1_has_higher_r_total, r_dissim, -r_dissim))

## Scatterplot of r_glu_d against +-|r_dis|, for all isolates
isolates_growth_dis %>%
    left_join(rename_with(select(isolates_growth_w, ID, `D-Glucose`), ~ paste0(., "1"))) %>%
    left_join(rename_with(select(isolates_growth_w, ID, `D-Glucose`), ~ paste0(., "2"))) %>%
    mutate(r_glu_d = `D-Glucose1` - `D-Glucose2`) %>%
    ggplot(aes(x = r_glu_d, y = r_dissim)) +
    geom_point(size = 2, shape = 1) +
    #geom_smooth(method = "lm") +
    theme_classic() +
    theme(legend.title = element_blank(), legend.position = "top", legend.direction = "horizontal") +
    guides(shape = "none") +
    labs(x = expression(r[A_glu] - r[B_glu]), y = expression(r[dis_acids]))

## Scatterplot of r_glu_d against +-|r_dis|, for only my isolates
p6 <- pairs_meta %>%
    left_join(isolates_growth_dis) %>%
    filter(!is.na(PairFermenter)) %>%
    ggplot(aes(x = r_glu_d, y = r_dissim, shape = InteractionType, color = PairFermenter)) +
    geom_point(size = 2) +
    scale_shape_manual(values = c("coexistence" = 16, "exclusion" = 1)) +
    scale_color_discrete(labels = c(FF = "Fermenter-Fermenter", FN = "Fermenter-Respirator", NN = "Respirator-Respirator")) +
    theme_classic() +
    theme(legend.title = element_blank(), legend.position = "top", legend.direction = "horizontal") +
    guides(shape = "none") +
    labs(x = expression(r[A_glu] - r[B_glu]), y = expression(r[dis_acids]))

p6_2 <- pairs_meta %>%
    left_join(isolates_growth_dis) %>%
    filter(!is.na(PairFermenter)) %>%
    ggplot(aes(x = r_glu_d, y = abs(r_dissim), shape = InteractionType, color = PairFermenter)) +
    geom_vline(xintercept = 0, linetype = 2) +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_point(size = 2) +
    scale_shape_manual(values = c("coexistence" = 16, "exclusion" = 1)) +
    scale_color_discrete(labels = c(FF = "Fermenter-Fermenter", FN = "Fermenter-Respirator", NN = "Respirator-Respirator")) +
    theme_classic() +
    theme(legend.title = element_blank(), legend.position = "right", legend.direction = "vertical") +
    #guides(shape = "none") +
    labs(x = expression(r[A_glu] - r[B_glu]), y = expression("pairwise distance in r"[acids]))

p <- plot_grid(p6, p6_2, ncol = 2)

ggsave(here::here("plots/Fig_pairs-r_dissim_scatter.png"), p6_2, width = 10, height = 5)


## Scatterplot of r_glu_d against r_total_d. The difference in the total r_acids sum
pairs_meta %>%
    left_join(isolates_growth_dis) %>%
    filter(!is.na(PairFermenter)) %>%
    ggplot() +
    #geom_abline(slope = -1, intercept = 0) +
    geom_point(aes(x = r_glu_d, y = r_total_d, shape = InteractionType, color = PairFermenter), size = 2) +
    scale_shape_manual(values = c("coexistence" = 16, "exclusion" = 1)) +
    scale_color_discrete(labels = c(FF = "Fermenter-Fermenter", FN = "Fermenter-Respirator", NN = "Respirator-Respirator")) +
    theme_classic() +
    theme(legend.title = element_blank(), legend.position = "top", legend.direction = "horizontal") +
    guides(shape = "none") +
    labs(x = expression(r[A_glu] - r[B_glu]), y = expression(r[dis_acids]))

## boxplot
pairs_meta %>%
    left_join(isolates_growth_dis) %>%
    filter(!is.na(PairFermenter)) %>%
    ggplot(aes(x = InteractionType, y = r_dissim, fill = InteractionType)) +
    geom_boxplot() +
    geom_jitter(shape = 1, size = 2) +
    geom_hline(yintercept = 0, color = "red") +
    scale_fill_manual(values = c("coexistence"=grey(0.5), "exclusion"=grey(0.9))) +
    facet_grid(.~PairFermenter, labeller = labeller(PairFermenter = c(FF = "Fermenter-Fermenter", FN = "Fermenter-Respirator", NN = "Respirator-Respirator"))) +
    theme_classic() +
    theme(legend.title = element_blank()) +
    labs(x = "", y = expression(r[dis_acids]))


# pH ----
# Boxplot pH
p <- pairs_meta %>%
    select(InteractionType, PairFermenter, ends_with("hr_d")) %>%
    pivot_longer(cols = starts_with("pH"), names_to = "Time", values_to = "pH_d") %>%
    filter(!is.na(PairFermenter)) %>%
    ggplot(aes(x = InteractionType, y = pH_d, fill = InteractionType)) +
    geom_boxplot() +
    geom_jitter(shape = 1, size = 2) +
    geom_hline(yintercept = 0, color = "red") +
    scale_fill_manual(values = c("coexistence"=grey(0.5), "exclusion"=grey(0.9))) +
    facet_grid(Time~PairFermenter, scales = "free_y",
               labeller = labeller(PairFermenter = c(FF = "Fermenter-Fermenter", FN = "Fermenter-Respirator", NN = "Respirator-Respirator"))) +
    theme_classic() +
    theme(legend.title = element_blank()) +
    labs(x = "", y = expression(pH[A] - pH[B]))

ggsave(here::here("plots/Fig_pairs-pH_boxplot.png"), p, width = 8, height = 8)

# Scatterplot, pH vs. growth rate for isolates
isolates_growth %>%
    separate(SID, into = c("ID", "CS"), sep = "_", convert = T) %>%
    filter(CS == "D-Glucose") %>%
    select(ID, CS, RMid) %>%
    left_join(isolates_byproduct_time) %>%
    select(ID, Time, RMid, pH, OD620) %>%
    filter(!is.na(Time)) %>%
    pivot_wider(id_cols = c(ID, RMid), names_from = Time, values_from = c(pH, OD620)) %>%
    right_join(isolates) %>%
    ggplot() +
    geom_point(aes(x = RMid, y = pH_48, color = Fermenter), shape = 1, size = 3) +
    theme_classic()

# Amount of total acid secretion ----
## Boxplot difference in the total amount of secretion
p <- pairs_meta %>%
    select(Community, ID1, ID2, r_glu_d, InteractionType, PairFermenter, contains("_sum_")) %>%
    pivot_longer(cols = ends_with("hr_d"), names_to = "Time", values_to = "byproduct_sum_d") %>%
    mutate(Time = sub("byproduct_sum_", "", sub("hr_d", "", Time))) %>%
    filter(!is.na(Time), !is.na(PairFermenter)) %>%
    ggplot(aes(x = InteractionType, y = byproduct_sum_d, fill = InteractionType)) +
    geom_boxplot() +
    #geom_jitter(shape = 1, size = 2) +
    geom_hline(yintercept = 0, color = "red") +
    scale_fill_manual(values = c("coexistence"=grey(0.5), "exclusion"=grey(0.9))) +
    facet_grid(Time~PairFermenter, scales = "free_y",
               labeller = labeller(PairFermenter = c(FF = "Fermenter-Fermenter", FN = "Fermenter-Respirator", NN = "Respirator-Respirator"))) +
    theme_classic() +
    theme(legend.title = element_blank()) +
    labs(x = "", y = expression(sum[A] - sum[B]))
ggsave(here::here("plots/Fig_pairs-total_acids_boxplot.png"), p, width = 8, height = 8)

## Boxplot, r_glu_d vs. sum_A and sum_B
pairs_meta_sum <- pairs_meta %>%
    select(Community, ID1, ID2, r_glu_d, InteractionType, PairFermenter, contains("ByproductSum")) %>%
    pivot_longer(cols = contains("Sum"), names_to = "Time", values_to = "byproduct_sum") %>%
    separate(col = Time, into = c("Time", "Isolate"), sep = "hr") %>%
    mutate(Time = sub("ByproductSum_", "", sub("", "", Time))) %>%
    filter(!is.na(PairFermenter), Time != 0)
p <- pairs_meta_sum %>%
    ggplot(aes(x = InteractionType, y = byproduct_sum, fill = Isolate)) +
    geom_boxplot() +
    #geom_jitter(aes(fill = Isolate), shape = 1, size = 2, height = 0, width = 0.1) +
    scale_fill_discrete(labels = c("A", "B"), name = "Isolate") +
    #scale_fill_manual(values = c("coexistence"=grey(0.5), "exclusion"=grey(0.9))) +
    facet_grid(Time~PairFermenter, scales = "free_y",
               labeller = labeller(PairFermenter = c(FF = "Fermenter-Fermenter", FN = "Fermenter-Respirator", NN = "Respirator-Respirator"))) +
    theme_classic() +
    theme(legend.title = element_blank()) +
    labs(x = "", y = expression(sum[acids]))
ggsave(here::here("plots/Fig_pairs-total_acids_boxplot2.png"), p, width = 8, height = 8)
## Scatterplot, r_glu_d vs. sum_d
p <- pairs_meta %>%
    filter(!is.na(PairFermenter)) %>%
    select(Community, ID1, ID2, r_glu_d, InteractionType, PairFermenter, contains("_sum_")) %>%
    pivot_longer(cols = ends_with("hr_d"), names_to = "Time", values_to = "byproduct_sum_d") %>%
    mutate(Time = sub("byproduct_sum_", "", sub("hr_d", "", Time))) %>%
    ggplot() +
    geom_vline(xintercept = 0, linetype = 2) +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_point(aes(x = r_glu_d, y = byproduct_sum_d, shape = InteractionType, color = PairFermenter), size = 2) +
    scale_shape_manual(values = c("coexistence" = 16, "exclusion" = 1)) +
    scale_color_discrete(labels = c(FF = "Fermenter-Fermenter", FN = "Fermenter-Respirator", NN = "Respirator-Respirator")) +
    facet_grid(.~Time, labeller = labeller(Time = c(`16` = "16hr", `28` = "28hr", `48` = "48hr"))) +
    theme_classic() +
    theme(legend.title = element_blank(), legend.position = "top", legend.direction = "horizontal") +
    labs(x = expression(r[A_glu] - r[B_glu]), y = expression(sum[A]-sum[B]))
ggsave(here::here("plots/Fig_pairs-total_acids_scatter.png"), p, width = 10, height = 4)


# Scatter r_glu vs. sum_acids ----
pairs_glu <- pairs_meta %>%
    select(Community, ID1, ID2, PairFermenter, InteractionType, `D-Glucose1`, `D-Glucose2`, Fermenter1, Fermenter2) %>%
    pivot_longer(cols = contains("Glucose"), names_to = "Isolate", values_to = "r_glu") %>%
    mutate(Isolate = sub("D-Glucose", "", Isolate))

pairs_sum <- pairs_meta_sum %>%  filter(Time == 16) %>% select(-r_glu_d)

pairs_glu_sum <- pairs_glu %>% left_join(pairs_sum)
pairs_glu_sum_w <- pairs_glu_sum %>% pivot_wider(names_from = Isolate, values_from = c(r_glu, byproduct_sum))
pairs_glu_sum %>%
    ggplot() +
    geom_segment(data = pairs_glu_sum_w, aes(x = r_glu_1, xend = r_glu_2, y = byproduct_sum_1, yend = byproduct_sum_2, color = InteractionType)) +
    geom_point(aes(x = r_glu, y = byproduct_sum, shape = Fermenter), shape = 16, size = 2) +
    scale_color_manual(values = interaction_color) +
    #scale_y_log10() +
    theme_classic()




# Statistics ---
## r_glu + byproduct amount
pairs_meta %>%
    filter(PairFermenter == "FF") %>%
    filter(!is.na(InteractionType), !is.na(r_glu_d)) %>%
    mutate(InteractionType = ifelse(InteractionType == "coexistence", 1, 0)) %>%
    glm(formula = InteractionType ~ r_glu_d * byproduct_sum_16hr_d, data = ., family = "binomial") %>%
    broom::tidy()

pairs_meta %>%
    filter(PairFermenter == "FF") %>%
    filter(!is.na(InteractionType), !is.na(r_glu_d)) %>%
    mutate(InteractionType = ifelse(InteractionType == "coexistence", 1, 0)) %>%
    glm(formula = InteractionType ~ r_glu_d * pH48hr_d, data = ., family = "binomial") %>%
    broom::tidy()





pairs_meta %>%
    filter(PairFermenter == "FF") %>%
    filter(!is.na(InteractionType), !is.na(r_glu_d)) %>%
    mutate(InteractionType = ifelse(InteractionType == "coexistence", 1, 0)) %>%
    glm(formula = InteractionType ~ r_glu_d * r_lac_d * r_suc_d, data =.,  family = "binomial") %>%
    broom::tidy()


x <- pairs_meta %>% filter(PairFermenter == "FF", InteractionType == "coexistence") %>% pull(r_glu_d)
y <- pairs_meta %>% filter(PairFermenter == "FF", InteractionType == "exclusion") %>% pull(r_glu_d)

t.test(x,y) %>% broom::tidy()
t.test(x) %>% broom::tidy()
t.test(y) %>% broom::tidy()

x <- pairs_meta %>% filter(PairFermenter == "FN", InteractionType == "coexistence") %>% pull(r_glu_d)
y <- pairs_meta %>% filter(PairFermenter == "FN", InteractionType == "exclusion", r_glu_d > -0.1) %>% pull(r_glu_d)
t.test(x,y)
t.test(x)
t.test(y)


x <- pairs_meta %>% filter(PairFermenter == "NN", InteractionType == "coexistence") %>% pull(r_glu_d)
y <- pairs_meta %>% filter(PairFermenter == "NN", InteractionType == "exclusion") %>% pull(r_glu_d)
t.test(x,y)
t.test(x)
t.test(y)



