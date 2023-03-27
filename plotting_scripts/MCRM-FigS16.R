library(tidyverse)
library(cowplot)
source(here::here("analysis/00-metadata.R"))
source(here::here("simulation/01-generate_input.R"))

# Parameters ----
input_parameters <- read_csv(here::here("simulation/01-input_parameters.csv"), col_types = cols())
input_monocultures <- read_csv(here::here("simulation/02a-input_monocultures.csv"), col_types = cols())
input_communities <- read_csv(here::here("simulation/02b-input_communities.csv"), col_types = cols())
input_communitiesWithoutCrossfeeding <- read_csv(here::here("simulation/02c-input_communitiesWithoutCrossfeeding.csv"), col_types = cols())

#mcrm_family_colors <- RColorBrewer::brewer.pal(10, "Set3") %>% setNames(paste0("F", 0:9))
mcrm_family_colors <- c(RColorBrewer::brewer.pal(9, "Set1"), RColorBrewer::brewer.pal(3, "Set2")[1] ) %>% setNames(paste0("F", 0:9))

family_names <- paste0("family ", 1:10) %>% setNames(paste0("F", 0:9))
resource_names <- LETTERS[1:10] %>% setNames(paste0("R", 0:9))
n_timesteps <- input_communities$n_timesteps[1]
n_timepoints <- input_communities$n_timepoints[1]


# Panel A for simulaton procedure



# 1. simulation species-level communities
communities_abundance <- read_csv(paste0(folder_simulation, "aggregated/12-communities_abundance.csv"), col_types = cols()) %>%
    mutate(Community = ordered(Community, paste0("W", 0:(input_communities$n_wells[1]-1)))) %>%
    mutate(Time = ordered(Time, c("init", paste0("T", 1:20), "end"))) %>%
    arrange(Community, Time)

# Barplot final time point
communities_richness <- read_csv(paste0(folder_simulation, "aggregated/12-communities_richness.csv"), col_types = cols()) %>%
    #mutate(Community = factor(Community, paste0("W", 0:(nrow(input_withinCommunityPairs)-1)))) %>%
    mutate(PairSize = choose(Richness, 2)) %>%
    arrange(desc(Richness)) %>%
    slice(1:20) %>%
    mutate(CommunityLabel = 1:20)

p0 <- communities_abundance %>%
    filter(Time == max(Time), Abundance > 0) %>%
    group_by(Community) %>%
    summarize(Richness = n()) %>%
    ggplot() +
    geom_histogram(aes(x = Richness), binwidth = 1, color = 1, fill = NA) +
    scale_x_continuous(breaks = 1:10) +
    scale_y_continuous(breaks = seq(0, 50, 5)) +
    theme_classic() +
    theme() +
    labs()

#
p1 <- communities_abundance %>%
    filter(Time == max(Time), Abundance > 0) %>%
    filter(Community %in% communities_richness$Community) %>%
    group_by(Community) %>%
    mutate(RelativeAbundance = Abundance / sum(Abundance)) %>%
    group_by(Community, Family) %>%
    # Species alpha
    mutate(CommunitySpecies = n():1) %>%
    left_join(communities_richness) %>%
    ggplot() +
    geom_col(aes(x = CommunityLabel, y = RelativeAbundance, fill = Family, alpha = CommunitySpecies)) +
    #geom_text(data = communities_richness, aes(x = CommunityLabel, label = Richness), y = 1.1) +
    #annotate("text", x = 1:20, y = 1.15, label = communities_richness$Richness, size = 4) +
    #annotate("text", x = 21, y = 1.1, label = c("n. of species"), size = 4, hjust = 0) +
    scale_fill_manual(values = mcrm_family_colors) +
    scale_alpha_continuous(range = c(1, 0.3)) +
    scale_x_continuous(breaks = seq(2, 20, 2), expand = c(0,0.1)) +
    scale_y_continuous(breaks = seq(0,1, 0.5), expand = c(0,0), limits = c(0,1.1)) +
    coord_cartesian(xlim = c(0.5, 20.5), ylim = c(0, 1), clip = "off") +
    theme_classic() +
    theme(plot.margin = unit(c(10,10,5,5), "mm"),
          panel.border = element_rect(color = 1, fill = NA)) +
    guides(color = "none", fill = "none", alpha = "none") +
    labs(x = "simulated community", y = "relative abundance") +
    ggtitle('"species"-level composition')

# 2. Simulation family level
p2 <- communities_abundance %>%
    filter(Time == max(Time), Abundance > 0) %>%
    filter(Community %in% communities_richness$Community) %>%
    group_by(Community) %>%
    mutate(RelativeAbundance = Abundance / sum(Abundance)) %>%
    # Species alpha
    group_by(Community, Family) %>%
    mutate(CommunitySpecies = n():1) %>%
    left_join(communities_richness) %>%
    ggplot() +
    geom_col(aes(x = CommunityLabel, y = RelativeAbundance, fill = Family)) +
    scale_fill_manual(values = mcrm_family_colors) +
    scale_alpha_continuous(range = c(1, 0.3)) +
    scale_x_continuous(breaks = seq(2, 20, 2), expand = c(0,0.1)) +
    scale_y_continuous(breaks = seq(0,1, 0.5), expand = c(0,0), limits = c(0,1.1)) +
    coord_cartesian(xlim = c(0.5, 20.5), ylim = c(0, 1), clip = "off") +
    theme_classic() +
    theme(plot.margin = unit(c(10,10,5,5), "mm"),
          panel.border = element_rect(color = 1, fill = NA)) +
    guides(color = "none", fill = guide_legend(title = "Family", ncol = 1),
           alpha = "none") +
    labs(x = "simulated community", y = "relative abundance") +
    ggtitle('"family"-level composition')

# 3. Actual data species level ----
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
    bin_ESV_names() %>%
    clean_ESV_names() %>%
    arrange(Community, Family, Transfer, ESV)

temp <- communities_abundance %>%
    filter(Transfer == 12) %>%
    filter(Community %in% communities$Community) %>%
    filter(Relative_Abundance > 0.01) %>%
    bin_ESV_names() %>%
    clean_ESV_names()

ESV_colors <- temp %>% get_ESV_colors

family_colors <- RColorBrewer::brewer.pal(length(unique(temp$Family)), "Set1") %>% setNames(unique(temp$Family))

p3 <- temp %>%
    mutate(ESV_ID = ifelse(ESV_ID %in% names(ESV_colors[-length(ESV_colors)]), ESV_ID, "Other")) %>%
    mutate(ESV_ID = factor(ESV_ID, rev(names(ESV_colors)))) %>%
    mutate(Family = factor(Family, names(family_colors))) %>%
    group_by(Community, Family) %>%
    mutate(CommunitySpecies = n():1) %>%
    left_join(communities) %>%
    ggplot() +
    geom_col(aes(x = CommunityLabel, y = Relative_Abundance, fill = Family, alpha = CommunitySpecies), color = NA) +
    scale_x_continuous(breaks = seq(2, 20, 2), expand = c(0,0.1)) +
    scale_y_continuous(limits = c(0,1), breaks = seq(0, 1, 0.5), expand = c(0,0)) +
    scale_fill_manual(values = family_colors) +
    scale_alpha_continuous(range = c(1, 0.3)) +
    theme_classic() +
    theme(strip.background = element_rect(fill = NA, color = NA),
          panel.border = element_rect(fill = NA, color = 1, linewidth = 1)) +
    guides(color = "none", alpha = "none") +
    labs(x = "empirical community", y = "relative abundance") +
    ggtitle('ESV-level composition')

p4 <- temp %>%
    mutate(ESV_ID = ifelse(ESV_ID %in% names(ESV_colors[-length(ESV_colors)]), ESV_ID, "Other")) %>%
    mutate(ESV_ID = factor(ESV_ID, rev(names(ESV_colors)))) %>%
    mutate(Family = factor(Family, names(family_colors))) %>%
    group_by(Community, Family) %>%
    mutate(CommunitySpecies = n():1) %>%
    left_join(communities) %>%
    ggplot() +
    geom_col(aes(x = CommunityLabel, y = Relative_Abundance, fill = Family), color = NA) +
    scale_x_continuous(breaks = seq(2, 20, 2), expand = c(0,0.1)) +
    scale_y_continuous(limits = c(0,1), breaks = seq(0, 1, 0.5), expand = c(0,0)) +
    scale_fill_manual(values = family_colors) +
    scale_alpha_continuous(range = c(1, 0.3)) +
    theme_classic() +
    theme(strip.background = element_rect(fill = NA, color = NA),
          panel.border = element_rect(fill = NA, color = 1, linewidth = 1)) +
    guides(color = "none", alpha = "none") +
    labs(x = "empirical community", y = "relative abundance") +
    ggtitle('family-level composition')

# Assemble panels ----
p <- plot_grid(
    plot_grid(NULL, NULL, ncol = 2, rel_widths = c(1,1), labels = LETTERS[c(1,2)], scale = .9),
    plot_grid(
        plot_grid(
            p1 + guides(fill = "none"),
            p3 + guides(fill = "none"),
            ncol = 1, axis = "tblr", align = "vh",
            labels = LETTERS[c(3,5)], scale = .9
        ),
        plot_grid(
            p2,
            p4,
            ncol = 1, axis = "tblr", align = "vh",
            labels = LETTERS[c(4,6)], scale = .9
        ),
        ncol = 2, rel_heights = c(1,1), rel_widths = c(1,1.3)
    ),
    ncol = 1, rel_heights = c(1, 3), scale = c(1, 1)
) +
    paint_white_background()
ggsave(here::here("plots/FigS16-in_silico_assembly.png"), p, width = 10, height = 8)










if (FALSE) {
    # 3. Communities without crossfeeding -----
    communitiesWithoutCrossfeeding_abundance <- read_csv(paste0(folder_simulation, "aggregated/12-communitiesWithoutCrossfeeding_abundance.csv"), col_types = cols()) %>%
        mutate(Community = ordered(Community, paste0("W", 0:(input_communitiesWithoutCrossfeeding$n_wells[1]-1)))) %>%
        mutate(Time = ordered(Time, c("init", paste0("T", 1:n_timepoints), "end"))) %>%
        arrange(Community, Time)

    # Barplot final time point
    communitiesWithoutCrossfeeding_abundance_abundant <- communitiesWithoutCrossfeeding_abundance %>%
        filter(Time == max(Time)) %>%
        filter(Abundance != 0) %>%
        group_by(Community) %>%
        mutate(TotalAbundance = sum(Abundance)) %>%
        filter(Abundance > 0.01 * sum(Abundance))
    communities_names <- tibble(Community = paste0("W", 0:19),
                                Community_new = factor(1:20, 1:20))

    #communitiesWithoutCrossfeeding_richness <- read_csv(paste0(folder_simulation, "11-aggregated/communitiesWithoutCrossfeeding_richness.csv"), col_types = cols())
    communitiesWithoutCrossfeeding_abundance_richness <- communitiesWithoutCrossfeeding_abundance_abundant %>%
        filter(Time == max(Time)) %>%
        group_by(Community, .drop = F) %>%
        summarize(Richness = n()) %>%
        left_join(communities_names)

    p3 <- communitiesWithoutCrossfeeding_abundance_abundant %>%
        filter(Time == max(Time)) %>%
        group_by(Community) %>%
        mutate(RelativeAbundance = Abundance / sum(Abundance)) %>%
        left_join(communities_names) %>%
        ggplot() +
        geom_col(aes(x = Community_new, y = RelativeAbundance, fill = Family), color = 1) +
        geom_text(data = communitiesWithoutCrossfeeding_abundance_richness, aes(x = Community_new, label = Richness), y = 1.1) +
        annotate("text", x = 1:20, y = 1.15, label = communities_richness$Richness, size = 4) +
        annotate("text", x = 21, y = 1.1, label = c("n. of species"), size = 4, hjust = 0) +
        scale_fill_manual(values = mcrm_family_colors) +
        scale_x_discrete(expand = c(0,.1)) +
        scale_y_continuous(breaks = seq(0,1,.5), expand = c(0,0), limits = c(0,1.1)) +
        coord_cartesian(xlim = c(0.5, 20.5), ylim = c(0, 1), clip = "off") +
        theme_classic() +
        theme(plot.margin = unit(c(1,.5,.5,.5), "cm"),
              panel.border = element_rect(color = 1, fill = NA)) +
        guides(color = "none", fill = guide_legend(title = "Family", ncol = 2)) +
        labs(x = "community", y = "relative abundance")


}

