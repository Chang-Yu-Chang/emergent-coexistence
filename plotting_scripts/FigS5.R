library(tidyverse)
library(cowplot)
library(broom)
source(here::here("analysis/00-metadata.R"))

communities <- read_csv(paste0(folder_data, "temp/00c-communities.csv"), show_col_types = F)
isolates <- read_csv(paste0(folder_data, "output/isolates.csv"), show_col_types = F)
sequences_alignment <- read_csv(paste0(folder_data, "temp/32-sequences_alignment.csv"), show_col_types = F)
communities_abundance_T12 <- read_csv(paste0(folder_data, "temp/32-communities_abundance_T12.csv"), show_col_types = F)
isolates_abundance <- read_csv(paste0(folder_data, "temp/32-isolates_abundance.csv"), show_col_types = F)


# Number of matchens
sequences_alignment %>%
    filter(ConsensusLength > 200, BasePairMismatch <= 4) %>%
    distinct(Community, Isolate)

algn_Sanger_ESV <- sequences_alignment %>%
    filter(ConsensusLength > 200, BasePairMismatch <= 4) %>%
    # For each Sanger, find the ESV which it has least base pair mismatch with
    group_by(ExpID) %>%
    filter(BasePairMismatch == min(BasePairMismatch)) %>%
    # When a Sanger has two ESV that both has zero, pick the first one
    slice(1) %>%
    ungroup()

nrow(algn_Sanger_ESV)

algn_Sanger_ESV %>% # Number of mismatches for each isolate Sanger
    group_by(BasePairMismatch) %>%
    count()

algn_Sanger_ESV1 <- algn_Sanger_ESV %>%
    group_by(Community, CommunityESVID) %>%
    arrange(BasePairMismatch) %>%
    slice(1)

table(algn_Sanger_ESV1$BasePairMismatch) # table of base pair mismatch
range(algn_Sanger_ESV1$ConsensusLength) # range of consensus length
nrow(algn_Sanger_ESV1) # Numer of unique matches

clean_isolate_names <- function (x) {
    x %>%
        group_by(Community) %>%
        mutate(ExpIDGenus = paste0(ExpID, "-", Genus)) %>%
        mutate(IsolateGenus = paste0(Isolate, "-", Genus))
}

isolates_algn_bp_mismatch <- isolates %>%
    select(Community, Isolate, ExpID, Genus) %>%
    left_join(algn_Sanger_ESV) %>%
    left_join(communities) %>%
    mutate(Community = factor(Community, communities$Community)) %>%
    clean_isolate_names

isolates_algn_bp_mismatch1 <- isolates %>%
    select(Community, Isolate, ExpID, Genus) %>%
    left_join(algn_Sanger_ESV1) %>%
    left_join(communities) %>%
    mutate(Community = factor(Community, communities$Community)) %>%
    clean_isolate_names

# Heatmap, ESV and sanger alignment, binned by mismatch
p <- sequences_alignment %>%
    left_join(communities) %>%
    #mutate(Community = factor(Community, communities$Community)) %>%
    clean_isolate_names %>%
    mutate(BasePairMismatchCategory = case_when(
        BasePairMismatch == 0 ~ "0",
        BasePairMismatch == 1 ~ "1",
        BasePairMismatch == 2 ~ "2",
        BasePairMismatch == 3 ~ "3",
        BasePairMismatch == 4 ~ "4",
        BasePairMismatch >= 5 ~ ">=5"
    ) %>% ordered(c(0:4, ">=5"))
    ) %>%
    ggplot() +
    geom_tile(aes(x = IsolateGenus, y = CommunityESVID, fill = BasePairMismatchCategory)) +
    # Each Sanger has one ESV match
    geom_tile(data = isolates_algn_bp_mismatch, aes(x = IsolateGenus, y = CommunityESVID), fill = NA, color = "gold", linewidth = 1) +
    # # Each ESV has one Sanger
    geom_tile(data = isolates_algn_bp_mismatch1, aes(x = IsolateGenus, y = CommunityESVID), fill = NA, color = "red", linewidth = 1) +
    #facet_wrap(.~Community, scales = "free", ncol = 4) +
    scale_fill_manual(values = rev(RColorBrewer::brewer.pal(n = 6, "Blues"))) +
    scale_x_discrete(expand = c(0,0)) +
    scale_y_discrete(expand = c(0,0)) +
    facet_wrap(.~CommunityLabel, scales = "free", ncol = 3) +
    theme_classic() +
    theme(
        panel.border = element_rect(fill = NA, color = 1, linewidth = 0.5),
        axis.text.x = element_text(angle = 45, vjust = 1.05, hjust = 1, size = 6),
        axis.text.y = element_text(size = 6),
        legend.position = c(.6, .05),
        strip.background = element_blank(),
        plot.title = element_text(hjust = 0.5)
    ) +
    guides(fill = guide_legend(nrow = 1, override.aes = list(color = "black"))) +
    labs(x = "isolate Sanger sequence", y = "ESV") +
    ggtitle("community")


ggsave(here::here("plots/FigS4-Sanger_ESV_alignment.png"), p, width = 8, height = 10)


if (FALSE) {


    # ESV enus level ESV  ----
    isolates_abundance <- read_csv(paste0(folder_data, "temp/32-isolates_abundance.csv"), col_types = cols())

    curate_abundant_genus <- function () {
        abundant_genus <- c(
            "Enterobacteriaceae","Klebsiella", "Enterobacter", "Raoultella", "Citrobacter", "Salmonella", "Pantoea", "Yersinia",
            "Pseudomonadaceae", "Pseudomonas", "Azomonas", "Aeromonas", "Acinetobacter", "Enterococcus", "Bordetella",
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
                #str_detect(ESV_ID, "Salmonella") ~ str_replace(ESV_ID, "\\.\\d+", ""),
                str_detect(ESV_ID, "Citrobacter") ~ str_replace(ESV_ID, "\\.\\d+", ""),
                str_detect(ESV_ID, "Enterobacter") ~ str_replace(ESV_ID, "\\.\\d+", ""),
                str_detect(ESV_ID, "Klebsiella") ~ str_replace(ESV_ID, "\\.\\d+", ""),
                #str_detect(ESV_ID, "Acinetobacter") ~ str_replace(ESV_ID, ".\\d+", ""),
                #str_detect(ESV_ID, "Pseudomonadaceae") ~ str_replace(ESV_ID, "Pseudomonadaceae", "Pseudomonadas"),
                str_detect(ESV_ID, "Pseudomonadaceae") ~ str_replace(ESV_ID, "\\.\\d+", ""),
                #str_detect(ESV_ID, "Pseudomonas.1\\d+") ~ str_replace(ESV_ID, "\\.1\\d+", "\\.1"),
                str_detect(ESV_ID, "Pseudomonas.2\\d+") ~ str_replace(ESV_ID, "\\.2\\d+", "\\.2"),
                str_detect(ESV_ID, "Pseudomonas.3\\d+") ~ str_replace(ESV_ID, "\\.3\\d+", "\\.3"),
                str_detect(ESV_ID, "Pseudomonas.4\\d+") ~ str_replace(ESV_ID, "\\.4\\d+", "\\.4"),
                str_detect(ESV_ID, "Pseudomonas.5\\d+") ~ str_replace(ESV_ID, "\\.5\\d+", "\\.5"),
                str_detect(ESV_ID, "Pseudomonas.5\\d+") ~ str_replace(ESV_ID, "\\.5\\d+", "\\.5"),
                str_detect(ESV_ID, "Azomonas") ~ str_replace(ESV_ID, ".\\d+", ""),
                #str_detect(ESV_ID, "Stenotrophomonas") ~ str_replace(ESV_ID, ".\\d+", ""),
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
        communities_abundance_ESV_ID <- comm_abundance %>%
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
            mutate(ESV_color = c(RColorBrewer::brewer.pal(8, "Set2"), RColorBrewer::brewer.pal(8, "Set3"))[1:n()]) %>%
            filter(ESV_ID != "Other")

        ESV_colors <- c(temp1$ESV_color %>% setNames(temp1$ESV_ID),
                        temp2$ESV_color %>% setNames(temp2$ESV_ID),
                        temp3$ESV_color %>% setNames(temp3$ESV_ID),
                        Other = "#999998")


        return(ESV_colors)

    }

    isolates_abundance_factor <- isolates_abundance %>%
        drop_na() %>%
        rename(ESV_ID = CommunityESVID)
    ESV_colors <- isolates_abundance_factor %>% get_ESV_colors()

    p2 <- isolates_abundance_factor %>%
        #bind_rows(isolates_abundance_other) %>%
        mutate(ESV_ID = factor(ESV_ID, names(ESV_colors))) %>%
        left_join(communities, by = "Community") %>%
        ggplot() +
        geom_col(aes(x = CommunityLabel, y = RelativeAbundance, fill = ESV_ID), color = "#000000", position = position_stack(reverse = T)) +
        theme_bw() +
        #scale_fill_manual(values = ESV_colors) +
        scale_fill_manual(values = c(RColorBrewer::brewer.pal(8, "Set1"), RColorBrewer::brewer.pal(7, "Set2"), RColorBrewer::brewer.pal(7, "Set3"))) +
        scale_x_continuous(breaks = 1:13, expand = c(0,.3)) +
        scale_y_continuous(breaks = c(0, .5, 1), expand = c(0,0), limits = c(0, 1)) +
        coord_cartesian(ylim = c(0, 1), clip = "off") +
        theme(panel.grid = element_blank(),
              axis.title = element_text(size = 17),
              axis.text = element_text(color = 1, size = 17),
              panel.border = element_rect(color = 1, linewidth = .5),
              plot.background = element_blank(),
              panel.background = element_blank()) +
        guides(alpha = "none", fill = guide_legend(title = "base pair mismatch")) +
        labs(x = "community", y = "relative abundance")

    # Family level ESV
    family_names <- unique(isolates_abundance$Family)
    family_names <- family_names[!is.na(family_names)]
    family_colors <- c(
        "Enterobacteriaceae" = rev(RColorBrewer::brewer.pal(11, "RdYlBu"))[2],
        "Pseudomonadaceae" = rev(RColorBrewer::brewer.pal(11, "PRGn"))[2],
        RColorBrewer::brewer.pal(length(family_names)-2, "Set2") %>% setNames(family_names[!(family_names %in% c("Enterobacteriaceae", "Pseudomonadaceae"))]),
        "Other" = "#999999"
    )

    p3 <- isolates_abundance %>%
        left_join(communities) %>%
        drop_na() %>%
        #bind_rows(mutate(isolates_abundance_other, Family = "Other") %>% left_join(communities)) %>%
        mutate(Family = factor(Family, names(family_colors))) %>%
        ggplot() +
        geom_col(aes(x = CommunityLabel, y = RelativeAbundance, fill = Family), color = "#000000", position = position_stack(reverse = T)) +
        theme_bw() +
        scale_fill_manual(values = family_colors) +
        scale_x_continuous(breaks = 1:13, expand = c(0,.3)) +
        scale_y_continuous(breaks = c(0, .5, 1), expand = c(0,0), limits = c(0, 1)) +
        coord_cartesian(ylim = c(0, 1), clip = "off") +
        theme(panel.grid = element_blank(),
              axis.title = element_text(size = 17),
              axis.text = element_text(color = 1, size = 17),
              panel.border = element_rect(color = 1, linewidth = 0.5),
              plot.background = element_blank(),
              panel.background = element_blank()) +
        guides(alpha = "none") +
        labs(x = "community", y = "relative abundance")

    p5 <- isolates %>%
        left_join(communities) %>%
        group_by(Community) %>%
        summarize(Total = sum(RelativeAbundance, na.rm = T)) %>%
        mutate(Xaxis = 1) %>%
        ggplot(aes(Xaxis, y = Total)) +
        geom_boxplot(width = 0.5, outlier.color = NA) +
        geom_jitter(width = 0.2, shape = 21, stroke = 1, size = 1.5) +
        scale_x_continuous(breaks = 1, limits = c(.5,1.5)) +
        scale_y_continuous(limits = c(0,1), breaks = seq(0,1,.5)) +
        theme_classic() +
        theme(
            panel.grid.major.x = element_line(color = grey(0.9)),
            panel.grid.major.y = element_line(color = grey(0.9)),
            panel.grid.minor.y = element_line(color = grey(0.9)),
            axis.text.x = element_blank(),
            panel.border = element_rect(color = 1, fill = NA),
            axis.text.y = element_text(color = 1, size = 10),
            axis.title = element_text(size = 10)
        ) +
        labs(x = "", y = "total abundance")

    p <- plot_grid(
        p1,
        plot_grid(p2 + guides(fill = guide_legend(ncol = 3, title = "ESV")) + theme(
            legend.key.size = unit(5, "mm"),
            legend.title = element_text(size = 10),
            axis.text = element_text(color = 1, size = 10),
            axis.title = element_text(size = 10)),
            NULL,
            nrow = 1, rel_widths = c(40,1)
        ),
        plot_grid(p3 + guides(fill = guide_legend(ncol = 2)) + theme(
            legend.key.size = unit(5, "mm"),
            legend.title = element_text(size = 10),
            axis.text = element_text(color = 1, size = 10),
            axis.title = element_text(size = 10)),
            p5,
            nrow = 1, rel_widths = c(3.9,1), axis = "tb", align = "h", labels = c("", "D"), label_x = -.2, label_y = 1.1
        ),
        ncol = 1, axis = "lr", scale = c(.95, .9, .9), labels = LETTERS[1:3], rel_heights = c(4,1,1)
    ) + paint_white_background()

}


# ESV richness per community
communities_abundance <- read_csv(paste0(folder_data, "raw/community_ESV/Emergent_Comunity_Data.csv"), show_col_types = F) %>%
    filter(Carbon_Source == "Glucose" | Carbon_Source == "Original") %>%
    mutate(Community = factor(paste0("C", Inoculum, "R", Replicate), paste0("C", rep(1:12, each = 8), "R", rep(1:8, 12))))%>%
    arrange(Community, Family, Transfer, ESV)

communities_abundance %>%
    filter(Community %in% communities$Community) %>%
    filter(Transfer == 12) %>%
    #filter(Relative_Abundance > 0.01) %>%
    group_by(Community) %>%
    count() %>%
    left_join(communities) %>%
    arrange(CommunityLabel)

# Matched ESV richness per community
isolates_abundance %>%
    drop_na() %>%
    group_by(Community) %>%
    count()

# Total abundance
isolates %>%
    group_by(Community) %>%
    summarize(Total = sum(RelativeAbundance, na.rm = T)) %>%
    summarize(Mean = mean(Total))



















