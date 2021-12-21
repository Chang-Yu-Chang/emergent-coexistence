# Figure S1

library(tidygraph)
library(ggraph)
library(tidyverse)
library(data.table)
library(cowplot)

sequences_abundance <- fread(here::here("data/temp/sequences_abundance.csv"))
communities <- fread("../data/output/communities.csv")
communities_name <- communities %>% pull(Community)


family_name <- c("Pseudomonadaceae", "Enterobacteriaceae", "Aeromonadaceae",  "Sphingobacteriaceae", "Xanthomonadaceae", "Moraxellaceae", "Alcaligenaceae", "Comamonadaceae", "Other") #        "Enterococcaceae", "Oxalobacteraceae", "Comamonadaceae", "Porphyromonadaceae", "Flavobacteriaceae", "Nocardiaceae", "Sphingomonadaceae", "Brucellaceae")
family_color <- c("#E21F27", "#397EB8", "#4fb049", "#984f9f", "#a94624", "#8fd1c6", "orange2", "firebrick", "#989898")
names(family_color) <- family_name
genus_name <- c("Pseudomonas", "Klebsiella", "Citrobacter", "Enterobacter", "Raoultella", "Stenotrophomonas", "Aeromonas") #, "Yersinia", "Buttiauxella", "Enterococcus", "Erwinia", "Sphingobacterium", "Bordetella", "Acinetobacter", "Salmonella", "Xanthomonas", "Pedobacter", "Pantoea", "Herbaspirillum", "Comamonas", "Dysgonomonas", "Delftia", "Flavobacterium", "Rhodococcus", "Novosphingobium", "Ochrobactrum", "Achromobacter", "Providencia"
genus_color <- c("#7d1416", "#3eb7c4", "#fdf8ce", "#225fa9", "#a2d6b3", "#f57e2d", "#8fd1c6")
names(genus_color) <- genus_name


plot_abundance <- function(df, label_x = "Community", label_y = "RelativeAbundance", fill = "CommunityESVID") {
    ggplot(df) +
        geom_bar(aes_string(x = label_x, y = label_y, fill = fill),
                 position = "stack", stat = "identity", col = 1) +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 90))
}

p <- sequences_abundance %>%
    filter(AlignmentType == "local") %>%
    filter(AllowMismatch == Inf) %>%
    filter(BasePairMismatch <= 4) %>%
    mutate(Community = ordered(Community,  communities_name)) %>%
    plot_abundance(label_x = "Community", label_y = " RelativeAbundance", fill = "Family") +
    labs(x = "", y = "Relative abundance") +
    scale_fill_manual(values = family_color) +
    scale_x_discrete(expand=c(0,0)) +
    scale_y_continuous(expand=c(0,0), limits = c(0,1), breaks = seq(0,1, .25)) +
    NULL

ggsave("../plots/FigS1.png", plot = p, width = 6, height = 4)

if (FALSE) {
    isolates <- fread("../data/output/isolates.csv")
    isolates_abundance <- fread(here::here("data/temp/isolates_abundance.csv"))
    pairs <- fread("../data/output/pairs.csv")
    pairs_melted <- fread("../data/output/pairs_melted.csv")
    isolates_random <- fread("../data/output/isolates_random.csv")
    pairs_random <- fread("../data/output/pairs_random.csv")
    pairs_random_melted <- fread("../data/output/pairs_random_melted.csv")
    communities <- fread("../data/output/communities.csv")
    community_names_ordered_by_size <- communities %>% arrange(CommunitySize) %>% pull(Community)
    community_abundance <- fread("../data/temp/communities_abundance.csv")


    # Panel: relative abundance of communities
    community_abundance_subset <- community_abundance %>%
        mutate(CommunityESVID = factor(CommunityESVID)) %>%
        filter(Community %in% community_names_ordered_by_size, Transfer == 12) %>%
        filter(!grepl("Glu-T12-X11-Rep1", SampleID), !grepl("Glu-T12-X1-Rep2-AA", SampleID),
               !grepl("Glu-T12-X1-Rep4-AA", SampleID), !grepl("Glu-T12-X1-Rep6-AA", SampleID),
               !grepl("Glu-T12-X4-Rep1", SampleID), !grepl("Glu-T12-X7-Rep1", SampleID)) %>%
        as_tibble()

    library(rRDPData)
    predict_RDP_taxonomy <- function(sequence_set) {
        pred <- predict(rdp(), sequence_set, confidence = 0) #  Predict 16s
        conf_score <- attr(pred, "confidence") %>% as.data.frame()  # Confidence score
        colnames(pred) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus")
        colnames(conf_score) <- paste0(c("Domain", "Phylum", "Class", "Order", "Family", "Genus"), "Score")
        pred$ID <- rownames(pred) %>% as.numeric()
        conf_score$ID <- rownames(pred) %>% as.numeric()
        return(pred)
    }
    esv_set <- DNAStringSet(community_abundance_subset$ESV)
    names(esv_set) <- community_abundance_subset$CommunityESVID # Rename the sequence
    esv_RDP <- predict_RDP_taxonomy(esv_set) %>%
        mutate(CommunityESVID = factor(ID)) %>%
        select(CommunityESVID, Family, Genus) %>%
        left_join(select(community_abundance_subset, CommunityESVID, Community))

    p_family <- left_join(community_abundance_subset, esv_RDP) %>%
        filter(RelativeAbundance >= 0.01) %>%
        #mutate(Family = as.character(Family)) %>%
        select(Family, Genus, Community, Transfer, CommunityESVID, RelativeAbundance) %>%
        mutate(Family = ifelse(Family %in% family_name, Family, "Other") %>% ordered(family_name)) %>%
        mutate(Community = ordered(Community, community_names_ordered_by_size)) %>%
        ggplot() +
        geom_bar(aes(x = Community, y = RelativeAbundance, fill = Family), stat = "identity", color = "grey20") +
        scale_x_discrete(expand = c(0,0)) +
        scale_y_continuous(expand = c(0,0), breaks = seq(0, 1, by = 0.2)) +
        scale_fill_manual(values = family_color) +
        theme_cowplot() +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
        labs(x = "community", y = "relative abundance")

    p_genus <- left_join(community_abundance_subset, esv_RDP) %>%
        filter(RelativeAbundance >= 0.01) %>%
        mutate(Genus = as.character(Genus)) %>%
        select(Family, Genus, Community, Transfer, CommunityESVID, RelativeAbundance) %>%
        mutate(Genus = ifelse(Genus == "Salmonella", "Klebsiella", Genus)) %>%
        mutate(Genus = ifelse(Genus %in% genus_name, Genus, "Other") %>% ordered(genus_name)) %>%
        mutate(Community = ordered(Community, community_names_ordered_by_size)) %>%
        ggplot() +
        geom_bar(aes(x = Community, y = RelativeAbundance, fill = Family, alpha = Genus), position = position_stack(reverse = T), stat = "identity", color = "grey20") +
        scale_x_discrete(expand = c(0,0)) +
        scale_y_continuous(expand = c(0,0), breaks = seq(0, 1, by = 0.2)) +
        scale_fill_manual(values = family_color) +
        theme_cowplot() +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
        labs(x = "community", y = "relative abundance")

    ps0 <- plot_grid(p_family, p_genus, ncol = 1, axis = "lr", align = "hv")
    ggsave("../plots/FigS0.png", plot = ps0, width = 6, height = 8)

}
