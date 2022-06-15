# CFU and CASEU comparison
library(tidyverse)
library(cowplot)
library(gridExtra)
source(here::here("plotting_scripts/misc.R"))

isolates <- read_csv("~/Dropbox/lab/emergent-coexistence/data/output/isolates.csv", col_types = cols())
pairs <- read_csv("~/Dropbox/lab/emergent-coexistence/data/output/pairs.csv", col_types = cols()) %>% mutate(InteractionType = ifelse(InteractionType == "neutrality", "coexistence", InteractionType))
pairs_freq <- read_csv("~/Dropbox/lab/emergent-coexistence/data/output/pairs_freq.csv", col_types = cols()) %>% mutate(Time = factor(Time, c("Tini", "Tend")))
pairs_example_outcomes_finer <- read_csv("~/Dropbox/lab/emergent-coexistence/data/output/pairs_example_outcomes_finer.csv", col_types = cols())
communities <- read_csv("~/Dropbox/lab/emergent-coexistence/data/output/communities.csv", col_types = cols())
community_factor <- communities$Community
communities_size <- communities$CommunitySize
communities_hierarchy <- read_csv("~/Dropbox/lab/emergent-coexistence/data/output/communities_hierarchy.csv", col_types = cols())
load("~/Dropbox/lab/emergent-coexistence/data/output/communities_network.Rdata") # communities_network


pairs_interaction <- read_csv("~/Dropbox/lab/emergent-coexistence/data/output/pairs_interaction.csv", col_types = cols()) %>%
    mutate(InteractionType = ifelse(InteractionType == "neutrality", "coexistence", InteractionType))
pairs_ID <- read_csv("~/Dropbox/lab/emergent-coexistence/data/temp/pairs_ID.csv", col_types = cols())


pairs_interaction %>%
    mutate(Set = factor(Set, c("CFUonly", "CASEUonly", "CFUandCASEU"))) %>%
    select(Set, PairID, InteractionType) %>%
    pivot_wider(names_from = Set, values_from = InteractionType) %>%
    mutate(ResultMatch = CASEUonly == CFUandCASEU) %>%
    group_by(ResultMatch) %>%
    summarize(Count = n())


pairs %>%
    filter(Community == "C7R1") %>%
    select(PairID, Community, Isolate1, Isolate2, InteractionType)


pairs_example_freq <- pairs %>%
    mutate(InteractionType = factor(InteractionType, c("exclusion", "coexistence"))) %>%
    mutate(InteractionTypeFiner = ordered(InteractionTypeFiner, c("competitive exclusion", "mutual exclusion", "stable coexistence", "frequency-dependent coexistence", "neutrality"))) %>%
    arrange(InteractionType, InteractionTypeFiner) %>%
    # Add the coordinate in the grid
    bind_rows(tibble(InteractionType = rep(NA, 4), InteractionTypeFiner = rep(NA, 4))) %>%
    filter(!is.na(tibble(InteractionType))) %>%
    # Join the frequency data
    select(Set, PairID, InteractionType, InteractionTypeFiner) %>%
    left_join(filter(pairs_freq, Set == "CFUandCASEU"), by = c("Set", "PairID")) %>%
    select(PairID, InteractionType, InteractionTypeFiner, Isolate1InitialODFreq, Time, Isolate1MeasuredFreq) %>%
    mutate(Isolate1InitialODFreq = factor(Isolate1InitialODFreq),
           Time = factor(Time, c("Tini", "Tend")))


pairs_example_freq %>%
    filter(PairID == 29) %>%
    view




pairs_freq %>%
    filter(Set == "CFUandCASEU") %>%
    group_by(PairID) %>%
    filter(Time == "Tend") %>%
    summarize(Mix = length(unique(RawDataType))) %>%
    filter(Mix == 2)
    group_by(Mix) %>%
    summarize(Count = n())
    filter(is.na(Isolate1MeasuredFreq))






pairs_freq %>%
    filter(PairID %in% c(120, 144)) %>%
    filter(Set == "CFUandCASEU") %>%
    filter(Time == "Tend")  %>% view


pairs_freq %>%
    filter(Community == "C8R4") %>%
    filter(Set == "CFUandCASEU") %>%
    #filter(Time == "Tend")  %>%
    view


pairs_example_freq <- pairs_interaction %>%
    left_join(pairs_ID) %>%
    select(Set, InteractionType, PairID) %>%
    mutate(InteractionType = factor(InteractionType, c("exclusion", "coexistence"))) %>%
    arrange(Set, InteractionType) %>%
    right_join(pairs_freq) %>%
    select(Set, PairID, InteractionType, Isolate1InitialODFreq, Time, Isolate1MeasuredFreq) %>%
    mutate(Set = factor(Set, c("CFUonly", "CASEUonly", "CFUandCASEU"))) %>%
    arrange(PairID, Set)

pairs_ID_grid <- pairs %>%
    mutate(InteractionType = factor(InteractionType, c("exclusion", "coexistence"))) %>%
    mutate(InteractionTypeFiner = ordered(InteractionTypeFiner, c("competitive exclusion", "mutual exclusion", "stable coexistence", "frequency-dependent coexistence", "neutrality"))) %>%
    arrange(InteractionType, InteractionTypeFiner) %>%
    # Add the coordinate in the grid
    bind_rows(tibble(InteractionType = rep(NA, 4), InteractionTypeFiner = rep(NA, 4))) %>%
    filter(!is.na(tibble(InteractionType))) %>%
    # Join the frequency data
    select(Set, PairID, InteractionType, InteractionTypeFiner) %>%
    left_join(filter(pairs_freq, Set == "CFUandCASEU"), by = c("Set", "PairID")) %>%
    select(PairID, InteractionType, InteractionTypeFiner, Isolate1InitialODFreq, Time, Isolate1MeasuredFreq) %>%
    mutate(Isolate1InitialODFreq = factor(Isolate1InitialODFreq),
           Time = factor(Time, c("Tini", "Tend"))) %>%
    arrange(InteractionTypeFiner, PairID) %>%
    distinct(PairID)

plot_sidebyside_freq <- function(tb) {

    temp <- distinct(tb, Set, InteractionType)
    if (nrow(temp) != 3) cat("The outcome of one data type is not unique")

    tb %>%
        mutate(Isolate1InitialODFreq = factor(Isolate1InitialODFreq)) %>%
        mutate(Time = factor(Time, c("Tini", "Tend"))) %>%
        mutate(Set = factor(Set, c("CFUonly", "CASEUonly", "CFUandCASEU"))) %>%
        ggplot() +
        geom_rect(data = temp, size = .5, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, aes(fill = InteractionType, linetype = Set), color = 1, alpha = .4) +
        geom_hline(size = .2, yintercept = c(0,1), linetype = 1, color = "grey90") +
        geom_point(size = .4, aes(x = Time, y = Isolate1MeasuredFreq, color = Isolate1InitialODFreq, group = Isolate1InitialODFreq)) +
        geom_line(size = .4, aes(x = Time, y = Isolate1MeasuredFreq, color = Isolate1InitialODFreq, group = Isolate1InitialODFreq)) +
        scale_y_continuous(breaks = c(0, .5, 1), limits = c(-.1, 1.1)) +
        scale_color_manual(values = frequency_color, label = c("95%", "50%", "5%")) +
        scale_fill_manual(values = c(assign_interaction_color(), "NA" = "grey90")) +
        scale_linetype_manual(values = c("CFUonly" = 1, "CASEUonly" = 2, "CFUandCASEU" = 3), name = "") +
        facet_grid(.~Set) +
        theme_bw() +
        theme(panel.spacing = unit(0, "mm"),
              strip.background = element_blank(),
              strip.text = element_blank(),
              #panel.border = element_rect(color = 1, fill = NA, size = 1),
              panel.border = element_blank(),
              panel.grid = element_blank(),
              panel.grid.minor.y = element_blank(),
              axis.title = element_blank(), axis.text = element_blank(),
              axis.ticks = element_blank(),
              plot.background = element_blank(),
              plot.title = element_text(size = 6, margin = margin(0,0,0,0)),
              plot.margin = margin(0,0,0,0, "mm")) +
        guides(color = "none", fill = "none", linetype = "none") +
        labs(x = "Time", y = "Frequency") +
        ggtitle(unique(tb$PairID))
}
temp_list <- pairs_example_freq %>%
    filter(PairID %in% 4:6) %>%
    mutate(PairID = factor(PairID, pairs_ID_grid$PairID)) %>%
    arrange(Set, PairID) %>%
    group_split(PairID) %>%
    lapply(plot_sidebyside_freq)













