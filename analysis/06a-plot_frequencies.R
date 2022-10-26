#' This script reads the aggregated csv from 06-cleanup_frequencies and for each pair
library(tidyverse)
library(cowplot)
source(here::here("analysis/00-metadata.R"))

pairs_ID <- read_csv(paste0(folder_data, "temp/00c-pairs_ID.csv"), show_col_types = F)
pairs_T0_boots <- read_csv(paste0(folder_data, "temp/06-pairs_T0_boots.csv"), show_col_types = F) # bootstraps using T0 mean and sd
pairs_T8_boots <- read_csv(paste0(folder_data, "temp/06-pairs_T8_boots.csv"), show_col_types = F) # bootstraps using random forest object probabilities=


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
    ggsave(paste0(folder_pipeline, "images/", pairs_ID$Batch[i], "-09-bootstrap/", pair_name, ".png"), p, width = 10, height = 5)
    write_csv(pair_boots_table, paste0(folder_pipeline, "images/", pairs_ID$Batch[i], "-09-bootstrap/", pair_name, ".csv"))
    cat("\n", pair_name, "\t", i, "/", nrow(pairs_ID))
}


