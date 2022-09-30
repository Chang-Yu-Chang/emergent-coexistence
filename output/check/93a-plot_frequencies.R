#' This script reads the aggregated csv from93-cleanup_frequencies and for each pair
library(tidyverse)
library(cowplot)

folder_script <- "~/Desktop/Lab/emergent-coexistence/output/check/"
folder_main <- "~/Dropbox/lab/emergent-coexistence/data/raw/plate_scan/emergent_coexistence_plate_scan_check/"
pairs_freq_ID <- read_csv(paste0(folder_main, "meta/00-pairs_freq_ID.csv"), show_col_types = F)
pairs_T0_boots <- read_csv(paste0(folder_main, "meta/93-pairs_T0_boots.csv"), show_col_types = F)
pairs_T8_boots <- read_csv(paste0(folder_main, "meta/93-pairs_T8_boots.csv"), show_col_types = F)
pairs_T8 <- read_csv(paste0(folder_main, "meta/93-pairs_T8.csv"), show_col_types = F)
plates_no_colony <- c(
    "B2_T8_C11R1_5-95_2_8",
    "B2_T8_C11R1_5-95_2_9",
    "B2_T8_C11R1_5-95_8_2",
    "B2_T8_C11R1_5-95_9_8",
    "B2_T8_C11R1_50-50_2_8",
    "B2_T8_C11R1_50-50_2_9",
    "C2_T8_C11R2_50-50_2_10",
    "C2_T8_C11R2_50-50_9_13"
)

pairs_no_colony <- c(
    "C11R1_2_8",
    "C11R1_2_9",
    "C11R1_8_9",
    "C11R2_2_10"
)


pairs_ID <- distinct(pairs_freq_ID, Batch, Community, Isolate1, Isolate2)

i=1
if (FALSE) {
# 1. T8 bootstraps ----
for (i in 1:nrow(pairs_ID)) {
    community <- pairs_ID$Community[i]
    isolate1 <- pairs_ID$Isolate1[i]
    isolate2 <- pairs_ID$Isolate2[i]
    pair_name <- paste0(community, "_", isolate1, "_", isolate2)
    if (pair_name %in% pairs_no_colony) {cat("\nT8 has no colony, skip pair\t", pair_name); next}

    # Increase or decrease sign
    pair_boots <- bind_rows(
        filter(pairs_T0_boots, Community == community, Isolate1 == isolate1, Isolate2 == isolate2),
        filter(pairs_T8_boots, Community == community, Isolate1 == isolate1, Isolate2 == isolate2),
    ) %>%
        select(Community, Isolate1, Isolate2, Isolate1InitialODFreq, Time, BootstrapID, Isolate1CFUFreq) %>%
        pivot_wider(names_from = Time, values_from = Isolate1CFUFreq) %>%
        mutate(FreqChange = T8-T0) %>%
        mutate(FreqChangeSign = case_when(
            FreqChange > 0 ~ "increase",
            FreqChange == 0 ~ "same",
            FreqChange < 0 ~ "decrease",
        )) %>%
        pivot_longer(cols = c(T0, T8), names_to = "Time", values_to = "Isolate1CFUFreq")

    # Table of increase and decrease
    pair_boots_table <- pair_boots %>%
        filter(Time == "T0") %>%
        mutate(FreqChangeSign = factor(FreqChangeSign, c("increase", "same", "decrease"))) %>%
        group_by(Community, Isolate1, Isolate2, Isolate1InitialODFreq, FreqChangeSign, .drop = F) %>%
        count(name = "Count") %>%
        group_by(Community, Isolate1, Isolate2, Isolate1InitialODFreq) %>%
        mutate(Fraction = Count / sum(Count))

    #
    color_names <- c("increase" = "#EF798A", "same" = grey(0.8), "decrease" = "#7D82B8")
    fill_names <- c("increase" = "#EF798A", "same" = grey(0.8), "decrease" = "#7D82B8")

    # T0 vs. T8 interaction plots
    p1 <- pair_boots %>%
        ggplot() +
        #geom_hline(yintercept = c(0, 1), linetype = 1) +
        geom_point(aes(x = Time, y = Isolate1CFUFreq, color = FreqChangeSign), shape = 21) +
        ## Add color on increase or decrease
        geom_line(aes(x = Time, y = Isolate1CFUFreq, color = FreqChangeSign, group = BootstrapID), size = .1) +
        scale_color_manual(values = color_names) +
        scale_y_continuous(limits = c(0,1), expand = c(0, 0)) +
        facet_grid(.~Isolate1InitialODFreq) +
        theme_classic() +
        theme(panel.border = element_rect(color = 1, fill = NA)) +
        ggtitle(pair_name)

    # Barplot
    p2 <- pair_boots_table %>%
        ggplot() +
        geom_hline(yintercept = c(0, 1), linetype = 1) +
        geom_col(aes(x = FreqChangeSign, y = Fraction, fill = FreqChangeSign), color = 1) +
        geom_hline(yintercept = c(0.05, 0.95), linetype = 2, color = "red") +
        scale_y_continuous(expand = c(0, 0)) +
        scale_fill_manual(values = color_names) +
        facet_grid(.~Isolate1InitialODFreq) +
        theme_classic() +
        theme(panel.border = element_rect(color = 1, fill = NA))


    p <- plot_grid(p1, p2, nrow = 2, axis = "lr", align = "v")
    ggsave(paste0(folder_main, "check/meta-93-T0_T8_frequencies/bootstrapped/", pair_name, ".png"), p, width = 10, height = 5)
    write_csv(pair_boots_table, paste0(folder_main, "check/meta-93-T0_T8_frequencies/bootstrapped/", pair_name, ".csv"))
    cat("\n", pair_name, "\t", i, "/", nrow(pairs_ID))
}


}
# 2. T8 random forest prediction -----
i=75
for (i in 1:nrow(pairs_ID)) {
#for (i in 74:nrow(pairs_ID)) {
    community <- pairs_ID$Community[i]
    isolate1 <- pairs_ID$Isolate1[i]
    isolate2 <- pairs_ID$Isolate2[i]
    pair_name <- paste0(community, "_", isolate1, "_", isolate2)
    if (pair_name %in% pairs_no_colony) {cat("\nT8 has no colony, skip pair\t", pair_name); next}

    # Increase or decrease sign
    pairs_T0_boot <- filter(pairs_T0_boots, Community == community, Isolate1 == isolate1, Isolate2 == isolate2) %>%
        select(Community, Isolate1, Isolate2, Isolate1InitialODFreq, Time, BootstrapID, Isolate1CFUFreq)
    pairs_T8_boot <- filter(pairs_T8, Community == community, Isolate1 == isolate1, Isolate2 == isolate2) %>%
        select(Community, Isolate1, Isolate2, Isolate1InitialODFreq, Time, Isolate1CFUFreq) %>%
        slice(rep(1:n(), each = 1000)) %>%
        mutate(BootstrapID = rep(1:1000, n()/1000)) %>%
        select(Community, Isolate1, Isolate2, Isolate1InitialODFreq, Time, BootstrapID, Isolate1CFUFreq)

    pair_boots <- bind_rows(pairs_T0_boot, pairs_T8_boot) %>%
        select(Community, Isolate1, Isolate2, Isolate1InitialODFreq, Time, BootstrapID, Isolate1CFUFreq) %>%
        pivot_wider(names_from = Time, values_from = Isolate1CFUFreq) %>%
        mutate(FreqChange = T8-T0) %>%
        mutate(FreqChangeSign = case_when(
            FreqChange > 0 ~ "increase",
            FreqChange == 0 ~ "same",
            FreqChange < 0 ~ "decrease",
        )) %>%
        pivot_longer(cols = c(T0, T8), names_to = "Time", values_to = "Isolate1CFUFreq")

    # Table of increase and decrease
    pair_boots_table <- pair_boots %>%
        filter(Time == "T0") %>%
        mutate(FreqChangeSign = factor(FreqChangeSign, c("increase", "same", "decrease"))) %>%
        group_by(Community, Isolate1, Isolate2, Isolate1InitialODFreq, FreqChangeSign, .drop = F) %>%
        count(name = "Count") %>%
        group_by(Community, Isolate1, Isolate2, Isolate1InitialODFreq) %>%
        mutate(Fraction = Count / sum(Count))

    #
    color_names <- c("increase" = "#EF798A", "same" = grey(0.8), "decrease" = "#7D82B8")
    fill_names <- c("increase" = "#EF798A", "same" = grey(0.8), "decrease" = "#7D82B8")

    # T0 vs. T8 interaction plots
    p1 <- pair_boots %>%
        ggplot() +
        #geom_hline(yintercept = c(0, 1), linetype = 1) +
        geom_point(aes(x = Time, y = Isolate1CFUFreq, color = FreqChangeSign), shape = 21) +
        ## Add color on increase or decrease
        geom_line(aes(x = Time, y = Isolate1CFUFreq, color = FreqChangeSign, group = BootstrapID), size = .1) +
        scale_color_manual(values = color_names) +
        scale_y_continuous(limits = c(0,1), expand = c(0, 0)) +
        facet_grid(.~Isolate1InitialODFreq) +
        theme_classic() +
        theme(panel.border = element_rect(color = 1, fill = NA)) +
        ggtitle(pair_name)

    # Barplot
    p2 <- pair_boots_table %>%
        ggplot() +
        geom_hline(yintercept = c(0, 1), linetype = 1) +
        geom_col(aes(x = FreqChangeSign, y = Fraction, fill = FreqChangeSign), color = 1) +
        geom_hline(yintercept = c(0.05, 0.95), linetype = 2, color = "red") +
        scale_y_continuous(expand = c(0, 0)) +
        scale_fill_manual(values = color_names) +
        facet_grid(.~Isolate1InitialODFreq) +
        theme_classic() +
        theme(panel.border = element_rect(color = 1, fill = NA))


    p <- plot_grid(p1, p2, nrow = 2, axis = "lr", align = "v")
    ggsave(paste0(folder_main, "check/meta-93-T0_T8_frequencies/classified/", pair_name, ".png"), p, width = 10, height = 5)
    write_csv(pair_boots_table, paste0(folder_main, "check/meta-93-T0_T8_frequencies/classified/", pair_name, ".csv"))
    cat("\n", pair_name, "\t", i, "/", nrow(pairs_ID))
}
