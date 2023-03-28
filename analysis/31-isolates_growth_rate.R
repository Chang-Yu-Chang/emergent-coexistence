library(tidyverse)
library(cowplot)
library(broom)
source(here::here("analysis/00-metadata.R"))

isolates <- read_csv(paste0(folder_data, "output/isolates.csv"), show_col_types = F)
communities <- read_csv(paste0(folder_data, "temp/00c-communities.csv"), show_col_types = F)
isolates_growth_syl <- read_csv(paste0(folder_data, "raw/growth_rate/Estrela_2021_isolates_grmax.csv"), col_types = cols()) %>%
    filter(cs == "glucose") %>%
    rename(ID = SangerID)

if (FALSE) {

# Data from Jean
isolates_curves <- read_csv(paste0(folder_data, "raw/growth_rate/20GC_Data.csv"), col_types = cols()) %>%
    select(ID = SangerID, CS, Time, OD620) %>%
    filter(ID %in% isolates$ID) %>%
    # Clean the Carbon Source names
    mutate(CS = tolower(CS)) %>%
    mutate(CS = str_replace(CS, "d-", "") %>% str_replace("l-", "") %>% str_replace("2-", "")) %>%
    # Clean the Time points
    mutate(Time = cut_width(Time, 1, boundary = 0, labels = F))  %>%
    distinct(ID, CS, Time, .keep_all = T) %>%
    filter(CS != "ddh20")

# Growth rates using the time points 12, 16, 28 hr
calculate_r <- function(N0, N1, T0, T1) (log10(N1)-log10(N0)) / (T1 - T0)
isolates_curves_T0 <- isolates_curves %>%
    group_by(ID, CS) %>%
    filter(Time == min(Time)) %>%
    mutate(T0 = ifelse(Time == min(Time), 0, Time), N0 = OD620) %>%
    # Assign a minimum OD value to prevent error in log(0)
    mutate(N0 = ifelse(N0 == 0, 0.001, N0)) %>%
    select(-Time, -OD620)
isolates_growth_jean <- isolates_curves %>%
    filter(Time %in% c(12, 16, 28)) %>%
    group_by(ID, CS) %>%
    arrange(ID, CS, Time) %>%
    mutate(T1 = Time, N1 = OD620) %>%
    select(-Time, -OD620) %>%
    left_join(isolates_curves_T0) %>%
    # Remove negative OD reads
    filter(N1 > 0) %>%
    # Calculate r
    mutate(r = calculate_r(N0, N1, T0, T1)) %>%
    select(ID, CS, Time = T1, r) %>%
    # Remove contamination
    filter(r>0) %>%
    # Average
    group_by(ID, CS, Time) %>%
    summarize(r = mean(r)) %>%
    pivot_wider(names_from = c(Time, CS), values_from = r, names_glue = "r_{CS}_{Time}hr") %>%
    ungroup()
}

# Ranked glucose growth rate
isolates <- isolates %>%
    left_join(isolates_growth_syl) %>%
    #left_join(isolates_growth_jean) %>%
    group_by(Community) %>%
    # Rank glucose maximum growth rate
    drop_na(gr_max) %>%
    mutate(rank_gr_max = rank(-gr_max)) # Top growth rate is rank 1
    # # Rank glucose growth rate at 16hr
    # drop_na(r_glucose_16hr) %>%
    # mutate(rank_r_glucose_16hr = rank(-r_glucose_16hr)) # Top growth rate is rank 1


# Growth vs rank
p <- isolates %>%
    ggplot(aes(x = rank_gr_max, y = Rank, group = ExpID)) +
    geom_jitter(shape = 21, size = 2, stroke = 1, width = 0.1, height = 0.1) +
    scale_x_continuous(breaks = 1:12) +
    scale_y_continuous(breaks = 1:12) +
    theme_classic() +
    guides(fill = guide_legend(title = "")) +
    labs(x = expression(ranked~growth~rate), y = "competitive rank")

ggsave(paste0(folder_data, "temp/31-01-growth_vs_rank.png"), p, width = 4, height = 4)

cor.test(isolates$Rank, isolates$rank_gr_max, method = "spearman", alternative = "two.sided", exact = FALSE) %>%
    tidy()


#  Growth vs. rank, colored by family
p <- isolates %>%
    ggplot(aes(x = rank_gr_max, y = Rank, group = ExpID, color = Family)) +
    geom_jitter(shape = 21, size = 2, stroke = 1, width = 0.1, height = 0.1) +
    scale_color_manual(values = RColorBrewer::brewer.pal(6, "Set2")) +
    scale_x_continuous(breaks = 1:12) +
    scale_y_continuous(breaks = 1:12) +
    theme_classic() +
    guides(fill = guide_legend(title = "")) +
    labs(x = expression(ranked~growth~rate), y = "competitive rank")
ggsave(paste0(folder_data, "temp/31-02-growth_vs_rank_family.png"), p, width = 5, height = 4)


# Growth vs. rank, facet by community
p <- isolates %>%
    mutate(Community = factor(Community, communities$Community)) %>%
    ggplot(aes(x = rank_gr_max, y = Rank, group = ExpID, color = Family)) +
    geom_point(shape = 21, size = 2, stroke = 1) +
    scale_color_manual(values = RColorBrewer::brewer.pal(6, "Set2")) +
    scale_x_continuous(breaks = 1:12) +
    scale_y_continuous(breaks = 1:12) +
    facet_wrap(~Community, ncol = 3, scales = "free") +
    theme_bw() +
    theme(
        panel.grid.minor = element_blank()
    ) +
    guides() +
    labs()

ggsave(paste0(folder_data, "temp/31-03-growth_vs_rank_comm.png"), p, width = 10, height = 8)

cor_test_long <- function (long_data) cor.test(long_data$Rank, long_data$rank_gr_max, method = "spearman")
isolates %>%
    mutate(Community = factor(Community, communities$Community)) %>%
    nest(data = c(-Community)) %>%
    mutate(fit = map(data, ~ cor_test_long(.x)),
           tidied = map(fit, tidy)) %>%
    unnest(tidied) %>%
    filter(p.value < 0.05)


# Growth vs. rank, normalized within community
isolates_norm <- isolates %>%
    group_by(Community) %>%
    left_join(communities) %>%
    mutate(RankNorm = Rank / CommunitySize, rank_gr_max_norm = rank_gr_max / CommunitySize)

p <- isolates_norm %>%
    ggplot(aes(x = rank_gr_max_norm, y = RankNorm)) +
    geom_abline(slope = 1, intercept = 0, color = "red", linetype = 2) +
    geom_point(shape = 21, size = 2, stroke = 1) +
    theme_classic() +
    theme() +
    guides() +
    labs()

ggsave(paste0(folder_data, "temp/31-04-growth_vs_rank_norm.png"), p, width = 4, height = 4)

cor.test(isolates_norm$RankNorm, isolates_norm$rank_gr_max_norm,
         method = "pearson", alternative = "two.sided", exact = FALSE) %>%
    tidy()

# Growth vs. rank, normalized within family and colored by family
p <- isolates_norm %>%
    mutate(Community = factor(Community, communities$Community)) %>%
    ggplot(aes(x = rank_gr_max_norm, y = RankNorm, group = ExpID, color = Family)) +
    geom_abline(slope = 1, intercept = 0, color = "red", linetype = 2) +
    geom_point(shape = 21, size = 2, stroke = 1) +
    scale_color_manual(values = RColorBrewer::brewer.pal(6, "Set2")) +
    theme_bw() +
    theme(
        panel.grid.minor = element_blank()
    ) +
    guides(fill = guide_legend(title = "")) +
    labs()

ggsave(paste0(folder_data, "temp/31-05-growth_vs_rank_norm_family.png"), p, width = 5, height = 4)

# Growth vs. rank, normalized within family and colored by family
p <- isolates_norm %>%
    mutate(Community = factor(Community, communities$Community)) %>%
    ggplot(aes(x = rank_gr_max_norm, y = RankNorm, group = ExpID, color = Family)) +
    geom_abline(slope = 1, intercept = 0, color = "red", linetype = 2) +
    geom_point(shape = 21, size = 2, stroke = 1) +
    scale_color_manual(values = RColorBrewer::brewer.pal(6, "Set2")) +
    facet_wrap(~Community, ncol = 3, scales = "free") +
    theme_bw() +
    theme(
        panel.grid.minor = element_blank()
    ) +
    guides(fill = guide_legend(title = "")) +
    labs(x = expression(ranked~growth~rate), y = "competitive rank")

ggsave(paste0(folder_data, "temp/31-06-growth_vs_rank_norm_family_comm.png"), p, width = 10, height = 8)

isolates_norm %>%
    mutate(Community = factor(Community, communities$Community)) %>%
    nest(data = c(-Community)) %>%
    mutate(fit = map(data, ~ cor_test_long(.x)),
           tidied = map(fit, tidy)) %>%
    unnest(tidied)


# Growth rate vs. ranked growth rate per community
p <- isolates_norm %>%
    mutate(Community = factor(Community, communities$Community)) %>%
    ggplot(aes(x = rank_gr_max, y = gr_max, group = ExpID, color = Family)) +
    geom_point(shape = 21, size = 2, stroke = 1) +
    scale_color_manual(values = RColorBrewer::brewer.pal(6, "Set2")) +
    scale_x_continuous(breaks = 1:12) +
    facet_wrap(~Community, ncol = 3, scales = "free") +
    theme_bw() +
    theme(
        panel.grid.minor = element_blank()
    ) +
    guides(fill = guide_legend(title = "")) +
    labs()

ggsave(paste0(folder_data, "temp/31-07-growth_vs_rankgrowth_comm.png"), p, width = 10, height = 8)










