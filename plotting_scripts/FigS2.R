library(tidyverse)
library(cowplot)
library(broom)
source(here::here("analysis/00-metadata.R"))
#source(here::here("plotting_scripts/FigS10.R")) # use the function defined in Fig.S10

communities <- read_csv(paste0(folder_data, "temp/00c-communities.csv"), show_col_types = F) %>%
    mutate(Community = factor(Community, Community))
communities_abundance <- read_csv(paste0(folder_data, "raw/community_ESV/Emergent_Comunity_Data.csv"), show_col_types = F) %>%
    filter(Carbon_Source == "Glucose" | Carbon_Source == "Original") %>%
    mutate(Community = factor(paste0("C", Inoculum, "R", Replicate), paste0("C", rep(1:12, each = 8), "R", rep(1:8, 12))))%>%
    # bin_ESV_names() %>%
    # clean_ESV_names() %>%
    arrange(Community, Family, Transfer, ESV)
communities_abundance_temporal <- communities_abundance %>%
    filter(Transfer != 0) %>%
    # Filter for those that has temporal data
    filter(Inoculum %in% c(2,6) | Replicate == 4) %>%
    select(Community, Transfer, ESV_ID, Relative_Abundance) %>%
    arrange(Community, Transfer, ESV_ID)

communities_abundance_temporal_complete <- communities_abundance_temporal %>%
    distinct(Community, ESV_ID) %>%
    slice(rep(1:n(), each = 12)) %>%
    mutate(Transfer = rep(1:12, n()/12))

ESV_stable <- communities_abundance_temporal %>% # ESVs that data from T8-12
    filter(Transfer %in% c(8:12)) %>%
    pivot_wider(id_cols = c(Community, ESV_ID), names_from = Transfer, names_prefix = "T", values_from = Relative_Abundance) %>%
    filter(!is.na(T8), !is.na(T9), !is.na(T10), !is.na(T11), !is.na(T12)) %>%
    #filter(Transfer %in% 11:12) %>%
    distinct(Community, ESV_ID) %>%
    mutate(CommunityESV = paste0(Community, ESV_ID))

# Calculate fitness
communities_abundance_fitness <- communities_abundance_temporal %>%
    right_join(communities_abundance_temporal_complete) %>%
    group_by(Community, ESV_ID) %>%
    arrange(Community, ESV_ID, Transfer) %>%
    mutate(Fitness = log(lead(Relative_Abundance) / Relative_Abundance))


#
communities_abundance_fitness_all <- communities_abundance_fitness %>%
    drop_na() %>%
    mutate(CommunityESV = paste0(Community, ESV_ID)) %>%
    filter(CommunityESV %in% ESV_stable$CommunityESV)

# Linear regression
tb_lm <- communities_abundance_fitness_all %>%
    nest(data = c(-Community, -ESV_ID)) %>%
    mutate(fit = map(data, ~ lm(Fitness ~ Relative_Abundance, data = .x)),
           tidied = map(fit, tidy),
           r.squared = map(fit, function(x) summary(x)$adj.r.squared)) %>%
    unnest(r.squared) %>%
    unnest(tidied)

tb_lm %>%
    filter(term == "Relative_Abundance") %>%
    select(Community, ESV_ID, estimate, std.error, p.value, r.squared)

ESV_sig <- tb_lm %>%
    filter(term == "Relative_Abundance") %>%
    select(Community, ESV_ID, estimate, std.error, p.value, r.squared) %>%
    filter(p.value < 0.05, estimate < 0) %>%
    distinct(Community, ESV_ID) %>%
    mutate(CommunityESV = paste0(Community, ESV_ID))


tb_lm %>%
    filter(term == "Relative_Abundance") %>%
    select(Community, ESV_ID, estimate, std.error, p.value, r.squared) %>%
    filter(p.value < 0.05, estimate < 0) %>%
    pull(r.squared) %>%
    range

# Abundance vs. fitness for ESVs in final communities, T1-12 ----
p <- communities_abundance_fitness_all %>%
    ggplot() +
    geom_smooth(data = filter(communities_abundance_fitness_all, CommunityESV %in% ESV_sig$CommunityESV),
                aes(x = Relative_Abundance, y = Fitness), method = "lm", formula = y~x, se = F) +
    geom_point(aes(x = Relative_Abundance, y = Fitness), shape = 21) +
    geom_hline(yintercept = 0, linetype = 2) +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 3)) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 3)) +
    facet_wrap(Community~ESV_ID, scales = "free", ncol = 9) +
    theme_classic() +
    theme(axis.text = element_text(size = 8, angle = 30, hjust = 1),
          axis.title = element_text(size = 15),
          strip.text = element_text(size = 8),
          panel.border = element_rect(color = 1, fill = NA)) +
    labs(x = expression(x[i]), y = expression(log(x[i+1]/x[i])))

ggsave(here::here("plots/FigS2-species_fitness_all_transfers.png"), p, width = 12, height = 15)


# ESVs that are ephemeral
communities_abundance_ephemeral <- communities_abundance_fitness %>%
    drop_na() %>%
    mutate(CommunityESV = paste0(Community, ESV_ID)) %>%
    filter(!(CommunityESV %in% ESV_stable$CommunityESV)) %>%
    nest(data = c(-Community, -ESV_ID)) %>%
    mutate(fit = map(data, ~ lm(Fitness ~ Relative_Abundance, data = .x)),
           tidied = map(fit, tidy),
           r.squared = map(fit, function(x) summary(x)$adj.r.squared)) %>%
    unnest(r.squared) %>%
    unnest(tidied) %>%
    filter(term == "Relative_Abundance") %>%
    select(Community, ESV_ID, estimate, std.error, p.value, r.squared) %>%
    filter(p.value < 0.05, estimate < 0)


communities_abundance_fitness_ephemeral <- communities_abundance_fitness %>%
    drop_na() %>%
    mutate(CommunityESV = paste0(Community, ESV_ID)) %>%
    filter(CommunityESV %in% paste0(communities_abundance_ephemeral$Community, communities_abundance_ephemeral$ESV_ID))

communities_abundance_fitness_ephemeral %>%
    ggplot() +
    geom_smooth(data = communities_abundance_fitness_ephemeral,
                aes(x = Relative_Abundance, y = Fitness), method = "lm", formula = y~x, se = F) +
    geom_point(aes(x = Relative_Abundance, y = Fitness), shape = 21) +
    geom_hline(yintercept = 0, linetype = 2) +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 3)) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 3)) +
    facet_wrap(Community~ESV_ID, scales = "free", ncol = 6) +
    theme_classic() +
    theme(axis.text = element_text(size = 8, angle = 30, hjust = 1),
          axis.title = element_text(size = 15),
          strip.text = element_text(size = 8),
          panel.border = element_rect(color = 1, fill = NA)) +
    labs(x = expression(x[i]), y = expression(log(x[i+1]/x[i])))



if (FALSE) {
communities <- read_csv(paste0(folder_data, "temp/00c-communities.csv"), show_col_types = F) %>%
    mutate(Community = factor(Community, Community))
communities_abundance <- read_csv(paste0(folder_data, "raw/community_ESV/Emergent_Comunity_Data.csv"), show_col_types = F) %>%
    filter(Carbon_Source == "Glucose" | Carbon_Source == "Original") %>%
    mutate(Community = factor(paste0("C", Inoculum, "R", Replicate), paste0("C", rep(1:12, each = 8), "R", rep(1:8, 12))))%>%
    # bin_ESV_names() %>%
    # clean_ESV_names() %>%
    arrange(Community, Family, Transfer, ESV)

communities_abundance_temporal <- communities_abundance %>%
    filter(Transfer != 0) %>%
    # Filter for those that has temporal data
    filter(Inoculum %in% c(2,6) | Replicate == 4) %>%
    select(Community, Transfer, ESV_ID, Relative_Abundance) %>%
    arrange(Community, Transfer, ESV_ID)

communities_abundance_temporal_complete <- communities_abundance_temporal %>%
    distinct(Community, ESV_ID) %>%
    slice(rep(1:n(), each = 12)) %>%
    mutate(Transfer = rep(1:12, n()/12))

# ESV that are present at T12
ESV_stable <- communities_abundance_temporal %>%
    filter(Transfer == 12) %>%
    select(Community, ESV_ID) %>%
    mutate(StableESV = "")


communities_abundance_temporal %>%
    #right_join(communities_abundance_temporal_complete) %>%
    filter(Transfer %in% c(8:12)) %>%
    pivot_wider(id_cols = c(Community, ESV_ID), names_from = Transfer, names_prefix = "T", values_from = Relative_Abundance) %>%
    right_join(ESV_stable, by = join_by(Community, ESV_ID))


communities_abundance_fitness <- communities_abundance_temporal %>%
    right_join(communities_abundance_temporal_complete) %>%
    group_by(Community, ESV_ID) %>%
    arrange(Community, ESV_ID, Transfer) %>%
    #mutate(lead(Relative_Abundance))
    mutate(Fitness = log(lead(Relative_Abundance) / Relative_Abundance))


# communities_abundance_fitness %>%
#     filter(Transfer %in% 8:12) %>%
#     group_by(Community, ESV_ID) %>%
#     drop_na() %>%
#     #filter(Community %in% c("C1R4")) %>%
#     ggplot() +
#     geom_point(aes(x = Relative_Abundance, y = Fitness), shape = 21) +
#     geom_hline(yintercept = 0, linetype = 2) +
#     facet_wrap(Community~ESV_ID) +
#     theme_classic() +
#     theme(axis.text = element_text(size = 15),
#           axis.title = element_text(size = 15),
#           panel.border = element_rect(color = 1, fill = NA)) +
#     labs(x = expression(x[i]), y = expression(log(x[i+1]/x[i])))



# Calculate Malthusian fitness
communities_abundance_fitness <- communities_abundance_sp %>%
    bin_ESV_names() %>%
    clean_ESV_names() %>%
    filter(Community %in% c("C1R4", "C2R6", "C2R8", "C8R4")) %>%
    filter(Transfer %in% c(1, 9:12)) %>%
    mutate(Time = case_when(
        Transfer == 1 ~ "init",
        Transfer %in% 9:12 ~ "end"
    )) %>%
    group_by(Community, ESV_ID, Time) %>%
    summarize(Relative_Abundance = mean(Relative_Abundance, na.rm = T)) %>%
    arrange(Community, ESV_ID, Time) %>%
    pivot_wider(names_from = Time, names_prefix = "T", values_from = Relative_Abundance) %>%
    drop_na() %>%
    mutate(Fitness = log(Tend/Tinit)) %>%
    select(Community, ESV_ID, Fitness)

# x_T1 vs. log(x_T12/x_T1)
temp <- communities_abundance_T1 %>%
    select(Community, ESV_ID, Relative_Abundance) %>%
    filter(Community %in% c("C1R4", "C2R6", "C2R8", "C8R4")) %>%
    group_by(Community, ESV_ID) %>%
    summarize(Relative_Abundance = sum(Relative_Abundance)) %>%
    ungroup()

p <- communities_abundance_fitness %>%
    left_join(temp) %>%
    ggplot(aes(x = Relative_Abundance, y = Fitness)) +
    geom_point(shape = 21, size = 3, stroke = 1) +
    geom_hline(yintercept = 0, linetype = 2) +
    theme_classic() +
    theme(axis.text = element_text(size = 15),
          axis.title = element_text(size = 15)) +
    labs(x = expression(x[init]), y = expression(log(x[end]/x[init])))

ggsave(here::here("plots/FigS11-species_abundance.png"), p, width = 4, height = 4)
    # log(x[Ti+1]/x[i]) over transfers
    communities_abundance_time <- communities_abundance_sp %>%
        select(ESV_ID, Transfer, Relative_Abundance) %>%
        arrange(ESV_ID, Transfer) %>%
        group_by(ESV_ID) %>%
        # time step change
        mutate(Fitness = log(lead(Relative_Abundance) / Relative_Abundance))

    communities_abundance_time %>%
        ggplot(aes(x = Relative_Abundance, y = Fitness)) +
        geom_point(shape = 21, size = 3, stroke = 1) +
        geom_hline(yintercept = 0, linetype = 2) +
        #facet_wrap(~ESV_ID) +
        #facet_wrap(~Transfer, scales = "free") +
        theme_classic() +
        theme(panel.border = element_rect(color = 1, fill = NA)) +
        labs(x = "x_T", y = "log(x_{Ti+1}/x_{Ti})")

    p3 <- communities_abundance_time %>%
        filter(Transfer != 12) %>%
        ggplot() +
        geom_boxplot(aes(x = Transfer, y = Fitness, group = Transfer), outlier.color = NA) +
        geom_hline(yintercept = 0, linetype = 2) +
        geom_point(aes(x = Transfer, y = Fitness), position = position_jitter(width = 0.1, height = 0),
                   size = 2, shape = 21, stroke = 1) +
        scale_x_continuous(breaks = 1:12) +
        theme_classic() +
        theme() +
        labs(x = "Transfer", y = expression(log(x[Ti+1]/x[i])))

}
