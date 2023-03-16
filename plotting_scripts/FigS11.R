library(tidyverse)
library(cowplot)
source(here::here("analysis/00-metadata.R"))
#source(here::here("plotting_scripts/FigS10.R"))

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
    # Set rare taxa into Other
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


communities_abundance_fitness %>%
    filter(Transfer %in% 8:12) %>%
    group_by(Community, ESV_ID) %>%
    drop_na() %>%
    #filter(Community %in% c("C1R4")) %>%
    ggplot() +
    geom_point(aes(x = Relative_Abundance, y = Fitness), shape = 21) +
    geom_hline(yintercept = 0, linetype = 2) +
    facet_wrap(Community~ESV_ID) +
    theme_classic() +
    theme(axis.text = element_text(size = 15),
          axis.title = element_text(size = 15),
          panel.border = element_rect(color = 1, fill = NA)) +
    labs(x = expression(x[i]), y = expression(log(x[i+1]/x[i])))



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

if (FALSE) {
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
