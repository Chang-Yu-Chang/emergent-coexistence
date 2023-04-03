library(tidyverse)
library(cowplot)
library(broom)
source(here::here("analysis/00-metadata.R"))

communities <- read_csv(paste0(folder_data, "temp/00c-communities.csv"), show_col_types = F) %>%
    mutate(Community = factor(Community, Community))
communities_abundance <- read_csv(paste0(folder_data, "raw/community_ESV/Emergent_Comunity_Data.csv"), show_col_types = F) %>%
    # Remove Leucine and Citrate communities
    filter(Carbon_Source == "Glucose" | Carbon_Source == "Original") %>%
    mutate(Community = factor(paste0("C", Inoculum, "R", Replicate), paste0("C", rep(1:12, each = 8), "R", rep(1:8, 12))))%>%
    arrange(Community, Family, Transfer, ESV)

communities_temporal <- communities_abundance %>%
    distinct(Community, Transfer) %>%
    arrange(Community, Transfer) %>%
    filter(Transfer == 4)

curate_abundant_genus <- function () {
    abundant_genus <- c(
        "Enterobacteriaceae","Klebsiella", "Enterobacter", "Raoultella", "Citrobacter", "Salmonella", "Pantoea", "Yersinia",
        "Pseudomonadaceae", "Pseudomonas", "Azomonas", "Aeromonas", "Acinetobacter", "Enterococcus",
        "Stenotrophomonas", "Delftia", "Serratia"
    )
    abundant_genus <- paste0(rep(abundant_genus, each = 101),
                             rep(c("", paste0(".", 1:100)), length(abundant_genus)))
    return(abundant_genus)
}
abundant_genus <- curate_abundant_genus()
bin_ESV_names <- function (comm_abundance) {
    comm_abundance %>%
        #mutate(ESV_ID = str_replace(ESV_ID, "Enterobacteriaceae", "Enterobacter")) %>%
        mutate(ESV_ID = case_when(
            str_detect(ESV_ID, "Enterobacteriaceae") ~ str_replace(ESV_ID, "\\.\\d+", ""),
            str_detect(ESV_ID, "Pantoea") ~ str_replace(ESV_ID, "\\.\\d+", ""),
            str_detect(ESV_ID, "Raoultella") ~ str_replace(ESV_ID, "\\.\\d+", ""),
            str_detect(ESV_ID, "Salmonella") ~ str_replace(ESV_ID, "\\.\\d+", ""),
            str_detect(ESV_ID, "Citrobacter") ~ str_replace(ESV_ID, "\\.\\d+", ""),
            str_detect(ESV_ID, "Enterobacter") ~ str_replace(ESV_ID, "\\.\\d+", ""),
            str_detect(ESV_ID, "Klebsiella") ~ str_replace(ESV_ID, "\\.\\d+", ""),
            str_detect(ESV_ID, "Acinetobacter") ~ str_replace(ESV_ID, ".\\d+", ""),
            #str_detect(ESV_ID, "Pseudomonadaceae") ~ str_replace(ESV_ID, "Pseudomonadaceae", "Pseudomonadas"),
            str_detect(ESV_ID, "Pseudomonadaceae") ~ str_replace(ESV_ID, "\\.\\d+", ""),
            str_detect(ESV_ID, "Pseudomonas.1\\d+") ~ str_replace(ESV_ID, "\\.1\\d+", "\\.1"),
            str_detect(ESV_ID, "Pseudomonas.2\\d+") ~ str_replace(ESV_ID, "\\.2\\d+", "\\.2"),
            str_detect(ESV_ID, "Pseudomonas.3\\d+") ~ str_replace(ESV_ID, "\\.3\\d+", "\\.3"),
            str_detect(ESV_ID, "Pseudomonas.4\\d+") ~ str_replace(ESV_ID, "\\.4\\d+", "\\.4"),
            str_detect(ESV_ID, "Pseudomonas.5\\d+") ~ str_replace(ESV_ID, "\\.5\\d+", "\\.5"),
            str_detect(ESV_ID, "Pseudomonas.5\\d+") ~ str_replace(ESV_ID, "\\.5\\d+", "\\.5"),
            str_detect(ESV_ID, "Azomonas") ~ str_replace(ESV_ID, ".\\d+", ""),
            str_detect(ESV_ID, "Stenotrophomonas") ~ str_replace(ESV_ID, ".\\d+", ""),
            str_detect(ESV_ID, "Aeromonas") ~ str_replace(ESV_ID, ".\\d+", ""),
            str_detect(ESV_ID, "Delftia") ~ str_replace(ESV_ID, ".\\d+", ""),
            str_detect(ESV_ID, "Enterococcus") ~ str_replace(ESV_ID, ".\\d+", ""),
            T ~ ESV_ID
        ))
}
clean_ESV_names <- function (comm_abundance) {
    comm_abundance %>%
        mutate(ESV_ID = ifelse(!(ESV_ID %in% abundant_genus), "Other", ESV_ID)) %>%
        return()
}
get_ESV_colors <- function (comm_abundance) {
    #comm_abundance <- temp
    communities_abundance_ESV_ID <- comm_abundance %>%
        #filter(Transfer == 12) %>%
        distinct(Family, ESV_ID) %>%
        mutate(ESV_ID = factor(ESV_ID, abundant_genus)) %>%
        arrange(Family, ESV_ID)

    # Entero
    temp1 <- communities_abundance_ESV_ID %>%
        filter(Family %in% c("Enterobacteriaceae")) %>%
        drop_na() %>%
        #mutate(ESV_color = viridis::viridis_pal(option = "viridis")(n())) %>%
        #mutate(ESV_color = scales::div_gradient_pal(low = "#313797", mid = "#ffffbf", high = "#a50026", space = "Lab")(seq(0, 1, length.out = n()))) %>%
        filter(ESV_ID != "Other")

    if (nrow(temp1) >= 12) {
        temp1 <- temp1 %>% mutate(ESV_color = c(rev(RColorBrewer::brewer.pal(11, "RdYlBu")), rev(RColorBrewer::brewer.pal(9, "OrRd")))[1:n()])
    } else if (nrow(temp1) >= 3) {
        temp1 <- temp1 %>% mutate(ESV_color = rev(RColorBrewer::brewer.pal(n(), "RdYlBu")))
    } else if (nrow(temp1) < 3) {
        temp1 <- temp1 %>% mutate(ESV_color = rev(RColorBrewer::brewer.pal(3, "RdYlBu"))[1:n()])
    }

    # Pseudo
    temp2 <- communities_abundance_ESV_ID %>%
        filter(Family == c("Pseudomonadaceae")) %>%
        drop_na() %>%
        #mutate(ESV_color = scales::seq_gradient_pal(low = "#d73027", high = "#fef1e7", space = "Lab")(seq(0, 1, length.out = n()))) %>%
        #mutate(ESV_color = scales::div_gradient_pal(low = "#00451a", mid = "#f6f6f7", high = "#41004a", space = "Lab")(seq(0, 1, length.out = n()))) %>%
        filter(ESV_ID != "Other")

    if (nrow(temp2) >= 12) {
        temp2 <- temp2 %>% mutate(ESV_color = c(rev(RColorBrewer::brewer.pal(11, "PRGn")), rev(RColorBrewer::brewer.pal(9, "YlGn")))[1:n()])
    } else if (nrow(temp2) >= 3) {
        temp2 <- temp2 %>% mutate(ESV_color = rev(RColorBrewer::brewer.pal(n(), "PRGn")))
    } else if (nrow(temp2) < 3) {
        temp2 <- temp2 %>% mutate(ESV_color =rev(RColorBrewer::brewer.pal(3, "PRGn"))[1:n()])
    }

    # Other
    temp3 <- communities_abundance_ESV_ID %>%
        filter(!(Family %in% c("Enterobacteriaceae", "Pseudomonadaceae")) & ESV_ID != "Other") %>%
        mutate(ESV_color = RColorBrewer::brewer.pal(9, "Set1")[1:n()]) %>%
        filter(ESV_ID != "Other")

    ESV_colors <- c(temp1$ESV_color %>% setNames(temp1$ESV_ID),
                    temp2$ESV_color %>% setNames(temp2$ESV_ID),
                    temp3$ESV_color %>% setNames(temp3$ESV_ID),
                    Other = "#999999")


    return(ESV_colors)

}


# T12 composition. Replicate reflect the actual order. Fig. S6 in Goldford2018 ----
communities_abundance_T12 <- communities_abundance %>%
    filter(Transfer == 12) %>%
    bin_ESV_names() %>%
    clean_ESV_names()

ESV_colors <- communities_abundance_T12 %>% get_ESV_colors

p <- communities_abundance_T12 %>%
    mutate(ESV_ID = ifelse(ESV_ID %in% names(ESV_colors[-length(ESV_colors)]), ESV_ID, "Other")) %>%
    mutate(ESV_ID = factor(ESV_ID, rev(names(ESV_colors)))) %>%
    ggplot() +
    geom_col(aes(x = Replicate, y = Relative_Abundance, fill = ESV_ID), width = .8, color = NA) +
    scale_x_continuous(breaks = 1:8, limits = c(0.5,8.5), expand = c(0,0)) +
    scale_y_continuous(breaks = seq(0, 1, 0.5), expand = c(0,0)) +
    scale_fill_manual(values = ESV_colors, breaks = names(ESV_colors)) +
    coord_flip() +
    facet_wrap(~Inoculum, labeller = labeller(Inoculum = label_both), ncol = 3, scales = "free_y") +
    theme_classic() +
    theme(strip.background = element_rect(fill = NA, color = NA),
          panel.border = element_rect(fill = NA, color = 1, linewidth = 1)) +
    guides(color = "none") +
    labs(x = "replicate", y = "relative abundance")

ggsave(paste0(folder_data, "temp/35-01-Goldford2018_FigS6_actual_replicate.png"), p, width = 10, height = 8)

# T12 composition. Replicate reordered according to the ESV abundance Fig.S6 in Goldford2018 ----
communities_abundance_reordered <- communities_abundance %>%
    filter(Transfer == 12, !is.na(Community)) %>%
    select(ESV_ID, Community, Inoculum, Replicate, Relative_Abundance) %>%
    bin_ESV_names() %>%
    clean_ESV_names() %>%
    group_by(ESV_ID, Community, Inoculum, Replicate) %>%
    summarize(Relative_Abundance = sum(Relative_Abundance)) %>%
    ungroup() %>%
    pivot_wider(names_from = ESV_ID, values_from = Relative_Abundance, values_fill = 0) %>%
    # Replicate ordered according to ESV abundance
    arrange(Inoculum, Klebsiella, Enterobacteriaceae, Raoultella, Citrobacter, Pseudomonas) %>%
    group_by(Inoculum) %>%
    mutate(ReplicateReordered = 1:n()) %>%
    select(Community, Inoculum, Replicate, ReplicateReordered)


communities_abundance_T12_ordered <- communities_abundance %>%
    left_join(communities_abundance_reordered) %>%
    filter(Transfer == 12) %>%
    bin_ESV_names() %>%
    clean_ESV_names()

ESV_colors <- communities_abundance_T12_ordered %>% get_ESV_colors

p <- communities_abundance_T12_ordered %>%
    mutate(ESV_ID = ifelse(ESV_ID %in% names(ESV_colors[-length(ESV_colors)]), ESV_ID, "Other")) %>%
    mutate(ESV_ID = factor(ESV_ID, rev(names(ESV_colors)))) %>%
    ggplot() +
    geom_col(aes(x = ReplicateReordered, y = Relative_Abundance, fill = ESV_ID), width = .8, color = NA) +
    scale_x_continuous(breaks = 1:8, limits = c(0.5,8.5), expand = c(0,0)) +
    scale_y_continuous(breaks = seq(0, 1, 0.5), expand = c(0,0)) +
    scale_fill_manual(values = ESV_colors, breaks = names(ESV_colors)) +
    coord_flip() +
    facet_wrap(~Inoculum, labeller = labeller(Inoculum = label_both), ncol = 3, scales = "free_y") +
    theme_classic() +
    theme(strip.background = element_rect(fill = NA, color = NA),
          panel.border = element_rect(fill = NA, color = 1, linewidth = 1)) +
    guides(color = "none") +
    labs(x = "replicate (reordered)", y = "relative abundance")

ggsave(paste0(folder_data, "temp/35-02-Goldford2018_FigS6.png"), p, width = 10, height = 8)



# Temporal dynamics of all CXXR4 communities. One replicate per inoculum. Fig.S2 in Goldford2018 ----
communities_abundance_R4 <- communities_abundance %>%
    filter(Community %in% communities_temporal$Community | Transfer == 0) %>%
    filter(Replicate == 4 | Transfer == 0) %>%
    bin_ESV_names() %>%
    clean_ESV_names()

ESV_colors <- communities_abundance_R4 %>% get_ESV_colors

p <- communities_abundance_R4 %>%
    mutate(Community = factor(Community, paste0("C", rep(1:12, each = 8), "R", rep(1:8, 12)))) %>%
    mutate(ESV_ID = ifelse(ESV_ID %in% names(ESV_colors[-length(ESV_colors)]), ESV_ID, "Other")) %>%
    mutate(ESV_ID = factor(ESV_ID, rev(names(ESV_colors)))) %>%
    ggplot() +
    geom_col(aes(x = Transfer, y = Relative_Abundance, fill = ESV_ID), width = 1, linewidth = .1) +
    scale_x_continuous(breaks = seq(0, 12, 2), limits = c(-0.5, 12.5), expand = c(0,0)) +
    scale_y_continuous(breaks = seq(0, 1, 0.2), expand = c(0,0)) +
    scale_fill_manual(values = ESV_colors, breaks = names(ESV_colors)) +
    #scale_color_manual(values = ESV_colors, breaks = names(ESV_colors)) +
    facet_wrap(~Inoculum, ncol = 3, scales = "free", labeller = labeller(Inoculum = label_both)) +
    theme_classic() +
    theme(strip.background = element_rect(fill = NA, color = NA),
          panel.border = element_rect(fill = NA, color = 1, linewidth = 1)) +
    guides(color = "none") +
    labs(x = "transfer", y = "relative abundance")
ggsave(paste0(folder_data, "temp/35-03-Goldford2018_FigS2.png"), p, width = 10, height = 8)


# Temporal dynamics of all C2RXX communities.  ----
communities_abundance_T0 <- communities_abundance %>%
    filter(Transfer == 0, Inoculum == 2) %>%
    slice(rep(1:n(), each = 8)) %>%
    mutate(Community = paste0("C2R", rep(1:8, n()/8)))

communities_abundance_C2 <- communities_abundance %>%
    filter(Community %in% communities_temporal$Community) %>%
    filter(Inoculum == 2) %>%
    bin_ESV_names() %>%
    clean_ESV_names()

ESV_colors <- communities_abundance_C2 %>% get_ESV_colors

p <- communities_abundance_C2 %>%
    bind_rows(communities_abundance_T0) %>%
    mutate(Community = factor(Community, paste0("C", rep(1:12, each = 8), "R", rep(1:8, 12)))) %>%
    mutate(ESV_ID = ifelse(ESV_ID %in% names(ESV_colors[-length(ESV_colors)]), ESV_ID, "Other")) %>%
    mutate(ESV_ID = factor(ESV_ID, rev(names(ESV_colors)))) %>%
    ggplot() +
    geom_col(aes(x = Transfer, y = Relative_Abundance, fill = ESV_ID), width = 1) +
    scale_x_continuous(breaks = seq(0, 12, 2), limits = c(-0.5, 12.5), expand = c(0,0)) +
    scale_y_continuous(breaks = seq(0, 1, 0.2), expand = c(0,0)) +
    scale_fill_manual(values = ESV_colors, breaks = names(ESV_colors)) +
    facet_wrap(~Community, ncol = 2, scales = "free") +
    theme_classic() +
    theme(strip.background = element_rect(fill = NA, color = NA),
          panel.border = element_rect(fill = NA, color = 1, linewidth = 1)) +
    guides(color = "none") +
    labs(x = "transfer", y = "relative abundance")
ggsave(paste0(folder_data, "temp/35-04-Goldford2018_Inoculum2.png"), p, width = 6, height = 8)


# Temporal dynamics of all C6RXX communities. Fig S7 in Goldford2018 ----
communities_abundance_T0 <- communities_abundance %>%
    filter(Transfer == 0, Inoculum == 6) %>%
    slice(rep(1:n(), each = 8)) %>%
    mutate(Community = paste0("C6R", rep(1:8, n()/8)))

communities_abundance_C6 <- communities_abundance %>%
    filter(Community %in% communities_temporal$Community) %>%
    filter(Inoculum == 6) %>%
    bin_ESV_names() %>%
    clean_ESV_names()

ESV_colors <- communities_abundance_C6 %>% get_ESV_colors

p <- communities_abundance_C6 %>%
    bind_rows(communities_abundance_T0) %>%
    mutate(Community = factor(Community, paste0("C", rep(1:12, each = 8), "R", rep(1:8, 12)))) %>%
    mutate(ESV_ID = ifelse(ESV_ID %in% names(ESV_colors[-length(ESV_colors)]), ESV_ID, "Other")) %>%
    mutate(ESV_ID = factor(ESV_ID, rev(names(ESV_colors)))) %>%
    ggplot() +
    geom_col(aes(x = Transfer, y = Relative_Abundance, fill = ESV_ID), width = 1) +
    scale_x_continuous(breaks = seq(0, 12, 2), limits = c(-0.5, 12.5), expand = c(0,0)) +
    scale_y_continuous(breaks = seq(0, 1, 0.2), expand = c(0,0)) +
    scale_fill_manual(values = ESV_colors, breaks = names(ESV_colors)) +
    facet_wrap(~Community, ncol = 2, scales = "free") +
    theme_classic() +
    theme(strip.background = element_rect(fill = NA, color = NA),
          panel.border = element_rect(fill = NA, color = 1, linewidth = 1)) +
    guides(color = "none") +
    labs(x = "transfer", y = "relative abundance")
ggsave(paste0(folder_data, "temp/35-05-Goldford2018_FigS7.png"), p, width = 6, height = 8)




# Temporal dynamics of only the four communities used for pairwise competition experiments ----
communities_abundance_T0 <- bind_rows(
    communities_abundance %>% filter(Transfer == 0, Inoculum == 1) %>% mutate(Community = "C1R4"),
    communities_abundance %>% filter(Transfer == 0, Inoculum == 2) %>% mutate(Community = "C2R6"),
    communities_abundance %>% filter(Transfer == 0, Inoculum == 2) %>% mutate(Community = "C2R8"),
    communities_abundance %>% filter(Transfer == 0, Inoculum == 8) %>% mutate(Community = "C8R4"),
)

temp <- communities_abundance %>%
    filter(Community %in% communities_temporal$Community) %>%
    filter(Community %in% communities$Community) %>%
    bin_ESV_names() %>%
    clean_ESV_names()

ESV_colors <- temp %>% get_ESV_colors

p <- temp %>%
    bind_rows(communities_abundance_T0) %>%
    mutate(Community = factor(Community, paste0("C", rep(1:12, each = 8), "R", rep(1:8, 12)))) %>%
    mutate(ESV_ID = ifelse(ESV_ID %in% names(ESV_colors[-length(ESV_colors)]), ESV_ID, "Other")) %>%
    mutate(ESV_ID = factor(ESV_ID, rev(names(ESV_colors)))) %>%
    ggplot() +
    geom_col(aes(x = Transfer, y = Relative_Abundance, fill = ESV_ID)) +
    scale_x_continuous(breaks = 1:12, limits = c(0.5, 12.5), expand = c(0,0)) +
    scale_y_continuous(breaks = seq(0, 1, 0.2), expand = c(0,0)) +
    scale_fill_manual(values = ESV_colors) +
    facet_wrap(~Community, ncol = 1, scales = "free") +
    theme_classic() +
    theme(strip.background = element_rect(fill = NA, color = NA),
          panel.border = element_rect(fill = NA, color = 1, linewidth = 1)) +
    guides(color = "none", fill = guide_legend(ncol = 1)) +
    labs(x = "", y = "replicate")
ggsave(paste0(folder_data, "temp/35-07-four_communities_temporal.png"), p, width = 5, height = 8)


# Test plotting the fitness ----
communities <- read_csv(paste0(folder_data, "output/communities_remained.csv"), show_col_types = F) %>%
    mutate(Community = factor(Community, Community))
communities_abundance <- read_csv(paste0(folder_data, "raw/community_ESV/Emergent_Comunity_Data.csv"), show_col_types = F) %>%
    filter(Carbon_Source == "Glucose" | Carbon_Source == "Original") %>%
    mutate(Community = factor(paste0("C", Inoculum, "R", Replicate), paste0("C", rep(1:12, each = 8), "R", rep(1:8, 12))))%>%
    arrange(Community, Family, Transfer, ESV)
communities_abundance_temporal <- communities_abundance %>%
    filter(Transfer != 0) %>%
    # Filter for those that has temporal data
    filter(Inoculum %in% c(2,6) | Replicate == 4) %>%
    select(Community, Transfer, ESV_ID, Relative_Abundance) %>%
    arrange(Community, Transfer, ESV_ID)

# A complete list of ESV-Transfer
communities_abundance_temporal_complete <- communities_abundance_temporal %>%
    distinct(Community, ESV_ID) %>%
    slice(rep(1:n(), each = 12)) %>%
    mutate(Transfer = rep(1:12, n()/12))
nrow(distinct(communities_abundance_temporal_complete, Community, ESV_ID)) # 755 unique ESVs
nrow(communities_abundance_temporal_complete) # 755 ESVs * 12 transfers = 9060 rows

# Calculate fitness
communities_abundance_fitness <- communities_abundance_temporal %>%
    right_join(communities_abundance_temporal_complete) %>%
    group_by(Community, ESV_ID) %>%
    arrange(Community, ESV_ID, Transfer) %>%
    mutate(Fitness = log(lead(Relative_Abundance) / Relative_Abundance))


# Check abundance vs. invasion fitness, ESVs present in stable communities, linear fit ----
ESV_stable <- communities_abundance_temporal %>%
    filter(Transfer %in% c(9:12)) %>%
    pivot_wider(id_cols = c(Community, ESV_ID), names_from = Transfer, names_prefix = "T", values_from = Relative_Abundance) %>%
    filter(
        # Species with complete presence at T8-12
        (!is.na(T9) & !is.na(T10) & !is.na(T11) & !is.na(T12) ) |
            # Community C4R4, C9R4 have missing time points at T11. For these communities, add back the species
            (Community %in% c("C4R4", "C9R4") & (!is.na(T9) & !is.na(T10) & !is.na(T12))) |
            # Community C5R4 have missing time points at T9. For these communities, add back the species
            (Community %in% c("C5R4") & (!is.na(T10) & !is.na(T11) & !is.na(T12)))
    ) %>%
    distinct(Community, ESV_ID) %>%
    mutate(CommunityESV = paste0(Community, ESV_ID))

communities_abundance_fitness_stable <- communities_abundance_fitness %>%
    drop_na() %>%
    mutate(CommunityESV = paste0(Community, ESV_ID)) %>%
    filter(CommunityESV %in% ESV_stable$CommunityESV)

# Correlation
cor_test_long <- function (long_data) cor.test(long_data$Relative_Abundance, long_data$Fitness, method = "spearman", exact = F, alternative = c("less"))
tb_cor_stable <- communities_abundance_fitness_stable %>%
    nest(data = c(-Community, -ESV_ID)) %>%
    mutate(fit = map(data, ~ cor_test_long(.x)),
           tidied = map(fit, tidy)) %>%
    unnest(tidied)
ESV_sig_stable <- tb_cor_stable %>%
    select(Community, ESV_ID, estimate, p.value) %>%
    filter(p.value < 0.05, estimate < 0) %>%
    mutate(CommunityESV = paste0(Community, ESV_ID)) %>%
    mutate(Significance = "p<0.05")

nrow(tb_cor_stable) # 99 stable ESVs
nrow(ESV_sig_stable) # 52 shows significant negative correlation
range(ESV_sig_stable$estimate) # range of correlation [-0.95, -0.53]
tb_cor_stable %>% filter(Community %in% c("C4R4", "C5R4", "C9R4")) %>% nrow() # 9 ESVs from C4R4, C5R4, C9R4 have only three time points

# Empirical equilibrium frequency
eq_freq_stable <- communities_abundance_temporal %>%
    mutate(CommunityESV = paste0(Community, ESV_ID)) %>%
    filter(Transfer %in% 9:12) %>%
    filter(CommunityESV %in% ESV_stable$CommunityESV) %>%
    arrange(Community, ESV_ID, CommunityESV, Transfer) %>%  # 387 rows. 99 * 4 - 9 = 387. 9 ESVs from C4R4, C5R4, C9R4 have only three time points
    group_by(Community, ESV_ID, CommunityESV) %>%
    summarize(EquilibriumAbundance = mean(Relative_Abundance)) # 99 ESVs


# Linear model predicted equilibrium frequency
xintercept_stable <- communities_abundance_fitness_stable %>%
    #filter(CommunityESV %in% ESV_sig_stable$CommunityESV) %>%
    nest(data = c(-Community, -ESV_ID, -CommunityESV)) %>%
    mutate(fit = map(data, ~ lm(Fitness ~ Relative_Abundance, data = .x)),
           tidied = map(fit, tidy)) %>%
    unnest(tidied) %>%
    select(Community, ESV_ID, CommunityESV, term, estimate) %>%
    pivot_wider(names_from = term, values_from = estimate) %>%
    # X intercept
    mutate(Xintercept = -`(Intercept)` / Relative_Abundance) %>%
    select(Community, ESV_ID, CommunityESV, Xintercept, Slope = Relative_Abundance, Yintercept = `(Intercept)`) %>%
    mutate(ESVType = "stable") %>%
    # Mark significant correlation
    left_join(select(tb_cor_stable, Community, ESV_ID, estimate, p.value))

p <- communities_abundance_fitness_stable %>%
    ggplot() +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_point(aes(x = Relative_Abundance, y = Fitness), shape = 21) +
    geom_smooth(data = communities_abundance_fitness_stable %>% left_join(ESV_sig_stable) %>% replace_na(list(Significance = "p>=0.05")),
                aes(x = Relative_Abundance, y = Fitness, color = Significance), method = stats::lm, formula = y ~ x, se = F) +
    # Linear model predicted
    geom_vline(data = xintercept_stable, aes(xintercept = Xintercept), color = "green", linetype = 2) +
    # Mean of T9-12
    geom_vline(data = communities_eq_freq_stable, aes(xintercept = EquilibriumAbundance), color = "navyblue", linetype = 2) +
    scale_color_manual(values = c("p<0.05" = "pink", "p>=0.05" = grey(0.8))) +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 3)) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 3)) +
    facet_wrap(Community~ESV_ID, scales = "free", ncol = 9) +
    theme_classic() +
    theme(
        axis.text = element_text(size = 8, angle = 30, hjust = 1),
        axis.title = element_text(size = 15),
        strip.text = element_text(size = 8),
        panel.border = element_rect(color = 1, fill = NA),
        legend.position = "top"
    ) +
    guides(color = "none") +
    labs(x = expression(x[i]), y = expression(log(x[i+1]/x[i])))

ggsave(paste0(folder_data, "temp/35-08-species_fitness_stable.png"), p, width = 12, height = 15)

# Check abundance vs. invasion fitness, ESVs present in stable communities, Polynomial fit, T1-12 ----

p <- communities_abundance_fitness_stable %>%
    ggplot() +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_point(aes(x = Relative_Abundance, y = Fitness), shape = 21) +
    geom_smooth(data = communities_abundance_fitness_stable %>% left_join(ESV_sig_stable) %>% replace_na(list(Significance = "p>=0.05")),
                aes(x = Relative_Abundance, y = Fitness, color = Significance),
                method = stats::loess, span = 0.9, formula = y ~ x, se = F) +
    # Linear model predicted
    geom_vline(data = xintercept_stable, aes(xintercept = Xintercept), color = "green", linetype = 2) +
    # Mean of T9-12
    geom_vline(data = communities_eq_freq_stable, aes(xintercept = EquilibriumAbundance), color = "navyblue", linetype = 2) +
    scale_color_manual(values = c("p<0.05" = "pink", "p>=0.05" = grey(0.8))) +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 3)) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 3)) +
    facet_wrap(Community~ESV_ID, scales = "free", ncol = 9) +
    theme_classic() +
    theme(
        axis.text = element_text(size = 8, angle = 30, hjust = 1),
        axis.title = element_text(size = 15),
        strip.text = element_text(size = 8),
        panel.border = element_rect(color = 1, fill = NA),
        legend.position = "top"
    ) +
    guides(color = "none") +
    labs(x = expression(x[i]), y = expression(log(x[i+1]/x[i])))

ggsave(paste0(folder_data, "temp/35-08a-species_fitness_stable_logfit.png"), p, width = 12, height = 15)


# Check abundance vs. invasion fitness, extinct ESVs, linear fit T1-12 ----
ESV_ephemeral <- communities_abundance_fitness %>%
    drop_na() %>%
    mutate(CommunityESV = paste0(Community, ESV_ID)) %>%
    filter(!(CommunityESV %in% ESV_stable$CommunityESV)) %>%
    group_by(CommunityESV, Community, ESV_ID) %>%
    count() %>%
    filter(n>=5) %>%
    # Remove the artifact C10R4 Stenotrophomonas
    filter(CommunityESV != "C10R4Stenotrophomonas")

communities_abundance_fitness_ephemeral <- communities_abundance_fitness %>%
    drop_na() %>%
    mutate(CommunityESV = paste0(Community, ESV_ID)) %>%
    filter(CommunityESV %in% ESV_ephemeral$CommunityESV)

# Correlation
tb_cor_ephemeral <- communities_abundance_fitness_ephemeral %>%
    nest(data = c(-Community, -ESV_ID)) %>%
    mutate(fit = map(data, ~ cor_test_long(.x)),
           tidied = map(fit, tidy)) %>%
    unnest(tidied)

ESV_sig_ephemeral <- tb_cor_ephemeral %>%
    select(Community, ESV_ID, estimate, p.value) %>%
    filter(p.value < 0.05, estimate < 0) %>%
    mutate(CommunityESV = paste0(Community, ESV_ID)) %>%
    mutate(Significance = "p<0.05")

nrow(tb_cor_ephemeral) # 46 ephemeral ESVs
nrow(ESV_sig_ephemeral) # 10 ephemeral ESVs shows significant negative correlation
range(ESV_sig_ephemeral$estimate) # range of correlation [-1, -0.714]

# Linear model predicted equilibrium frequency
xintercept_ephemeral <- communities_abundance_fitness_ephemeral %>%
    #filter(CommunityESV %in% ESV_sig_ephemeral$CommunityESV) %>%
    nest(data = c(-Community, -ESV_ID, -CommunityESV)) %>%
    mutate(fit = map(data, ~ lm(Fitness ~ Relative_Abundance, data = .x)),
           tidied = map(fit, tidy)) %>%
    unnest(tidied) %>%
    select(Community, ESV_ID, CommunityESV, term, estimate) %>%
    pivot_wider(names_from = term, values_from = estimate) %>%
    # X intercept
    mutate(Xintercept = -`(Intercept)` / Relative_Abundance) %>%
    select(Community, ESV_ID, CommunityESV, Xintercept, Slope = Relative_Abundance, Yintercept = `(Intercept)`) %>%
    mutate(ESVType = "transient") %>%
    # Mark significant correlation
    left_join(select(tb_cor_ephemeral, Community, ESV_ID, estimate, p.value))


p <- communities_abundance_fitness_ephemeral %>%
    ggplot() +
    geom_smooth(data = communities_abundance_fitness_ephemeral %>% left_join(ESV_sig_ephemeral) %>% replace_na(list(Significance = "p>=0.05")),
                aes(x = Relative_Abundance, y = Fitness, color = Significance),
                method = stats::lm, formula = y ~ x, se = F) +
    geom_point(aes(x = Relative_Abundance, y = Fitness), shape = 21) +
    geom_hline(yintercept = 0, linetype = 2) +
    scale_color_manual(values = c("p<0.05" = "pink", "p>=0.05" = grey(0.8))) +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 3)) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 3)) +
    facet_wrap(Community~ESV_ID, scales = "free", ncol = 9) +
    theme_classic() +
    theme(
        axis.text = element_text(size = 8, angle = 30, hjust = 1),
        axis.title = element_text(size = 15),
        strip.text = element_text(size = 8),
        panel.border = element_rect(color = 1, fill = NA),
        legend.position = "top"
    ) +
    guides(color = "none") +
    labs(x = expression(x[i]), y = expression(log(x[i+1]/x[i])))

ggsave(paste0(folder_data, "temp/35-09-species_fitness_extinct.png"), p, width = 12, height = 9)


# Check one ESV's fitness ----

p <- communities_abundance_fitness %>%
    drop_na() %>%
    mutate(CommunityESV = paste0(Community, ESV_ID)) %>%
    #filter(Community == "C8R4", ESV_ID == "Raoultella") %>%
    filter(Community == "C8R4", ESV_ID %in% c("Raoultella", "Klebsiella", "Pseudomonas.10")) %>%
    ggplot() +
    geom_smooth(aes(x = Relative_Abundance, y = Fitness), method = "lm", formula = y~x, se = T, color = "black") +
    geom_point(aes(x = Relative_Abundance, y = Fitness), shape = 21) +
    geom_hline(yintercept = 0, linetype = 2) +
    scale_x_continuous(breaks = c(0.001, 0.01, 0.1, 1), trans = "log10") +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 3)) +
    facet_wrap(~ESV_ID, nrow = 1) +
    theme_classic() +
    theme(
        axis.text = element_text(size = 8, angle = 30, hjust = 1),
        axis.title = element_text(size = 15),
        axis.title.x = element_blank(),
        strip.text = element_text(size = 6),
        panel.border = element_rect(color = 1, fill = NA),
        legend.position = "top"
    ) +
    labs(x = expression(x[i]), y = expression(log(x[i+1]/x[i])))

ggsave(paste0(folder_data, "temp/35-10-check_one_ESV.png"), p, width = 6, height = 3)


# Check the ESVs in the 13 communities -----
tb_cor <- communities_abundance_fitness_stable %>%
    nest(data = c(-Community, -ESV_ID)) %>%
    mutate(fit = map(data, ~ cor_test_long(.x)),
           tidied = map(fit, tidy)) %>%
    unnest(tidied)

ESV_sig <- tb_cor %>%
    select(Community, ESV_ID, estimate, p.value) %>%
    filter(p.value < 0.05, estimate < 0) %>%
    mutate(CommunityESV = paste0(Community, ESV_ID)) %>%
    mutate(Significance = "p<0.05")
plot_comm_ESV <- function (communities_abundance_fitness_stable, comm) {
    n_facets = 9
    communities_abundance_fitness_stable_comm <- communities_abundance_fitness_stable %>%
        left_join(ESV_sig) %>%
        replace_na(list(Significance = "p>=0.05")) %>%
        filter(Community %in% comm)

    n_ESVs <- length(unique(communities_abundance_fitness_stable_comm$ESV_ID))

    p1 <- communities_abundance_fitness_stable %>%
        filter(Community %in% comm) %>%
        ggplot() +
        geom_smooth(data = communities_abundance_fitness_stable_comm,
                    aes(x = Relative_Abundance, y = Fitness, color = Significance),
                    method = "lm", formula = y ~ x, se = F) +
        geom_point(aes(x = Relative_Abundance, y = Fitness), shape = 21) +
        geom_hline(yintercept = 0, linetype = 2) +
        scale_color_manual(values = c("p<0.05" = "maroon", "p>=0.05" = grey(0.8))) +
        #scale_x_continuous(breaks = c(0.001, 0.1, 1), trans = "log10") +
        scale_x_continuous(breaks = c(0.001, 0.01, 0.1, 1), trans = "log10") +
        scale_y_continuous(breaks = scales::pretty_breaks(n = 3)) +
        facet_wrap(~ESV_ID, nrow = 1) +
        theme_classic() +
        theme(
            axis.text = element_text(size = 8, angle = 30, hjust = 1),
            axis.title = element_text(size = 15),
            axis.title.x = element_blank(),
            strip.text = element_text(size = 6),
            panel.border = element_rect(color = 1, fill = NA),
            legend.position = "top"
        ) +
        guides(color = "none") +
        labs(x = expression(x[i]), y = expression(log(x[i+1]/x[i]))) +
        ggtitle(comm)

    plot_grid(p1, NULL, rel_widths = c(n_ESVs, n_facets - n_ESVs))

}

p1 <- plot_comm_ESV(communities_abundance_fitness_stable, "C1R4")
p2 <- plot_comm_ESV(communities_abundance_fitness_stable, "C2R6")
p3 <- plot_comm_ESV(communities_abundance_fitness_stable, "C2R8")
p4 <- plot_comm_ESV(communities_abundance_fitness_stable, "C8R4")

p <- plot_grid(p1,p2,p3,p4 + theme(axis.title.x = element_text(size = 5)), ncol = 1) + paint_white_background()

ggsave(paste0(folder_data, "temp/35-12-species_fitness_stable_comm_ordered.png"), p, width = 10, height = 6)


# Histogram: linear model predicted equilibrium frequency of stable vs. transient ESVs ----
xintercept <- bind_rows(xintercept_stable, xintercept_ephemeral)
table(xintercept$ESVType) # 99 stable ESVs and 46 transient ESVs

xintercept %>%
    filter(p.value < 0.05, estimate < 0) %>%
    pull(ESVType) %>%
    table() # For siginficant negaative correlation, 52 stable ESVs, 10 transient ESVs

p <- xintercept %>%
    filter(p.value < 0.05, estimate < 0) %>%
    ggplot() +
    geom_histogram(aes(x = Xintercept, fill = ESVType), color = 1, binwidth = 0.05) +
    scale_fill_manual(values = RColorBrewer::brewer.pal(3, "Set3"), label = c(stable = "stable ESV", transient = "transient ESV")) +
    scale_x_continuous(limits = c(-0.1,1.1), breaks = seq(0,1,.1), expand = c(0,0)) +
    scale_y_continuous(expand = c(0,0)) +
    theme_classic() +
    theme(
        legend.position = c(0.8, 0.8)
    ) +
    guides(fill = guide_legend(title = "")) +
    labs(x = "equilibirum frequency")
ggsave(paste0(folder_data, "temp/35-13-lm_eq_stable_vs_transient.png"), p, width = 4, height = 3)
#ggsave(here::here("plots/FigS3a-ESV_eq_freq.png"), p, width = 4, height = 3)


# Correlation plot: for all stable ESVs, mean of T9-12 vs. lm predicted eq freq. color by siginficance ----
# x: the equilibrium frequency calculated as mean of T9-T12 vs.
# y: the equilibrium frequency predicted by the negative linear model

eq_freq_stable <- communities_eq_freq_stable %>%
    left_join(xintercept_stable) %>%
    drop_na(Xintercept)

p <- eq_freq_stable %>%
    mutate(ShowNegFreqDep = case_when(
        estimate < 0 & p.value < 0.05 ~ "neg freq dep",
        T ~ "no evidence of neg freq dep"
    )) %>%
    ggplot() +
    geom_abline(intercept = 0, slope = 1, linetype = 2, color = "red") +
    geom_point(aes(x = EquilibriumAbundance, y = Xintercept, color = ShowNegFreqDep), shape = 21, size = 2, stroke = 1) +
    scale_color_manual(values = c("neg freq dep" = "black", "no evidence of neg freq dep" = grey(0.8))) +
    scale_x_continuous(limits = c(0,1)) +
    scale_y_continuous(limits = c(0,1)) +
    theme_classic() +
    theme(
        panel.border = element_rect(color = 1, fill = NA, linewidth = 0.5),
        #legend.position = c(0.7, 0.2),
        legend.background = element_rect(color = NA, fill = NA)
    ) +
    guides() +
    labs(x = "mean of T9-T12", y = "predicetd by linear model")

ggsave(paste0(folder_data, "temp/35-14-ESV_eq_freq_predicted.png"), p, width = 6, height = 4)

cor.test(eq_freq_stable$EquilibriumAbundance, eq_freq_stable$Xintercept, method = "pearson") %>%
    tidy()


# Correlation plot: for all stable ESVs, mean of T9-12 vs. lm predicted eq freq. color by communities used in this study ----
# x: the equilibrium frequency calculated as mean of T9-T12 vs.
# y: the equilibrium frequency predicted by the negative linear model

eq_freq_stable_comm <- eq_freq_stable %>%
    mutate(InThisStudy = case_when(
        Community %in% communities$Community ~ "four communities in current study",
        T ~ "other 22 communities"
    )) %>%
    arrange(desc(InThisStudy))

eq_freq_stable_comm %>%
    group_by(InThisStudy, Community) %>%
    count() %>%
    pull(InThisStudy) %>%
    table() # 4 communities in current study, 22 other communities

eq_freq_stable_comm %>%
    group_by(InThisStudy, Community) %>%
    count() # 18 red points. Or 18 ESVs in the current study


p <- eq_freq_stable_comm %>%
    ggplot() +
    geom_abline(intercept = 0, slope = 1, linetype = 2, color = "black") +
    geom_point(aes(x = EquilibriumAbundance, y = Xintercept, color = InThisStudy), shape = 21, size = 2, stroke = 1) +
    scale_color_manual(values = c("four communities in current study" = "red", "other 22 communities" = grey(0.8))) +
    scale_x_continuous(limits = c(0,1)) +
    scale_y_continuous(limits = c(0,1)) +
    theme_classic() +
    theme(
        panel.border = element_rect(color = 1, fill = NA, linewidth = 0.5),
        #legend.position = c(0.7, 0.2),
        legend.background = element_rect(color = NA, fill = NA)
    ) +
    guides() +
    labs(x = "mean of T9-T12", y = "predicetd by linear model")
# There are 16 red points because two points have negative predicted ESV eq freq from the linear model
ggsave(paste0(folder_data, "temp/35-15-ESV_eq_freq_predicted_comm.png"), p, width = 6, height = 4)

cor.test(eq_freq_stable$EquilibriumAbundance, eq_freq_stable$Xintercept, method = "pearson") %>%
    tidy()

table(eq_freq_stable_comm$Slope < 0) # 95 ESVs have slope <0, 4 ESVs have slope > 0
table(eq_freq_stable_comm$Xintercept > 0) # 94 ESVs have the predicted x intercept >0 , 5 ESVs have x intercept < 0

eq_freq_stable_comm %>%
    mutate(NegativeSlope = Slope < 0) %>%
    mutate(PositiveXintercept = Xintercept > 0) %>%
    group_by(NegativeSlope, PositiveXintercept) %>%
    count()

# For both stable and ephemeral ESVs, compare average fitness value vs. equilibirum freq ----
# Fitness
communities_abundance_fitness_mean <- bind_rows(
    communities_abundance_fitness_stable %>% mutate(ESVType = "stable"),
    communities_abundance_fitness_ephemeral %>% mutate(ESVType = "transient")
) %>%
    group_by(ESVType, Community, ESV_ID) %>%
    summarize(MeanFitness = mean(Fitness), NumberPoint = n(), SdFitness = sd(Fitness))
table(communities_abundance_fitness_mean$ESVType) # 99 stable ESVs and 46 transient ESVs

# Linear model predicted equilibrium frequency
xintercept <- bind_rows(xintercept_stable, xintercept_ephemeral)
table(xintercept$ESVType) # 99 stable ESVs and 46 transient ESVs

fitness_eq_freq <- left_join(communities_abundance_fitness_mean, xintercept)
nrow(fitness_eq_freq) # 99+46 = 145 rows

p <- fitness_eq_freq %>%
    ggplot() +
    geom_point(aes(x = MeanFitness, y = Xintercept, color = ESVType), shape = 21, size = 2, stroke = 1) +
    #geom_segment(aes(x = MeanFitness - 2*SdFitness, xend = MeanFitness + 2*SdFitness, y = Xintercept, yend = Xintercept, color = ESVType)) +
    geom_vline(xintercept = 0, linewidth = 0.1) +
    geom_hline(yintercept = 0, linewidth = 0.1) +
    scale_color_manual(values = RColorBrewer::brewer.pal(3, "Set1"), label = c(stable = "stable ESV", transient = "transient ESV")) +
    theme_classic() +
    theme(
        panel.border = element_rect(color = 1, fill = NA, linewidth = .5)
    ) +
    guides() +
    labs(x = "average fitness value", y = "predicted equilibrium frequency from linear model")

ggsave(paste0(folder_data, "temp/35-16-stable_vs_transient.png"), p, width = 5, height = 4)


# The number of data points in each ESV panel
p <- fitness_eq_freq %>%
    mutate(NumberPoint = factor(NumberPoint, 5:11)) %>%
    mutate(ESVType = factor(ESVType)) %>%
    group_by(ESVType, NumberPoint, .drop = F) %>%
    summarize(Count = n()) %>%
    ggplot() +
    geom_col(aes(x = NumberPoint, y = Count, fill = ESVType), position = position_dodge()) +
    scale_fill_manual(values = RColorBrewer::brewer.pal(3, "Set1"), label = c(stable = "stable ESV", transient = "transient ESV")) +
    theme_classic() +
    theme() +
    guides() +
    labs()
ggsave(paste0(folder_data, "temp/35-17-stable_vs_transient_samplesize.png"), p, width = 5, height = 4)

# Distribution of fitness value without error bar
p1 <- fitness_eq_freq %>%
    arrange(ESVType) %>%
    mutate(CommunityESV = factor(CommunityESV)) %>%
    ggplot() +
    geom_point(aes(x = CommunityESV, y = MeanFitness, color = ESVType), shape = 21, size = 2, stroke = 1) +
    #geom_segment(aes(x = CommunityESV, xend = CommunityESV, y = MeanFitness - SdFitness, yend = MeanFitness + SdFitness, color = ESVType)) +
    scale_color_manual(values = RColorBrewer::brewer.pal(3, "Set1"), label = c(stable = "stable ESV", transient = "transient ESV")) +
    coord_flip() +
    scale_x_discrete(position = "top") +
    theme_classic() +
    theme(
        panel.grid.major.y = element_line(color = grey(0.9)),
        legend.position = "top",
        axis.text.y = element_text(hjust = 0)
        #axis.text.x = element_text(size = 6, angle = 45, hjust = 0)
    ) +
    guides() +
    labs()

# Distribution of fitness value with error bar
p2 <- fitness_eq_freq %>%
    arrange(ESVType) %>%
    mutate(CommunityESV = factor(CommunityESV)) %>%
    ggplot() +
    geom_point(aes(x = CommunityESV, y = MeanFitness, color = ESVType), shape = 21, size = 2, stroke = 1) +
    geom_segment(aes(x = CommunityESV, xend = CommunityESV, y = MeanFitness - SdFitness, yend = MeanFitness + SdFitness, color = ESVType)) +
    scale_color_manual(values = RColorBrewer::brewer.pal(3, "Set1"), label = c(stable = "stable ESV", transient = "transient ESV")) +
    coord_flip() +
    scale_x_discrete(position = "top") +
    theme_classic() +
    theme(
        panel.grid.major.y = element_line(color = grey(0.9)),
        legend.position = "top",
        axis.text.y = element_text(hjust = 0)
        #axis.text.x = element_text(size = 6, angle = 45, hjust = 0)
    ) +
    guides() +
    labs()

p <- plot_grid(p1, p2, nrow = 1, align = "h", axis = "tb")
ggsave(paste0(folder_data, "temp/35-18-stable_vs_transient_fitness_value.png"), p, width = 10, height = 15)
























































