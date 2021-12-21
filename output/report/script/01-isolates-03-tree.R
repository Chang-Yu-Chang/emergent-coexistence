library(tidyverse)
library(treedataverse)
library(ggsci)
library(ggnewscale)
library(cowplot)
load(here::here("data/temp/isolates_sanger_seq.Rdata"))
communities <- read_csv(here::here("data/output/communities.csv")) %>% filter(str_detect(Community, "C\\d"))

# Isolate attributes
isolates <- read_csv(here::here("data/output/isolates.csv")) %>%
    filter(!is.na(ID)) %>%
    mutate(ID = as.character(ID)) %>%
    distinct(ID, .keep_all = T) %>%
    filter(Assembly == "self_assembly")

# Remove tips that are not in the isolate list
tip_to_remove <- which(!(tree$tip.label %in% isolates$ID))
tree1 <- tree %>% drop.tip(tip_to_remove)

# Append meta data
tree_meta <- tree1 %>%
    as_tibble() %>%
    mutate(ID = label) %>%
    left_join(isolates) %>%
    as.treedata()

trim_by_community <- function (tree, comm) {
    tips_to_drop <- tree %>%
        as_tibble() %>%
        filter(Community != comm) %>%
        pull(node)
    drop.tip(tree, tips_to_drop)
}
plot_tree <- function(tree) {
    tree %>%
        #ggtree(branch.length = "none") +
        ggtree() +
        geom_tippoint(aes(colour = Fermenter), size = 5) +
        scale_color_npg(name = "", label = c("TRUE" = "Fermenter", "FALSE" = "Respirator")) +
        new_scale_color() +
        geom_tiplab(aes(label = paste0(Genus, " sp.")), offset = 0.05) +
        scale_x_continuous(expand = expansion(.5, .1)) +
        theme_void() +
        theme(legend.position = "none")
}

communities_tree <- communities %>%
    select(comm = Community, everything()) %>%
    mutate(tree_meta = list(tree_meta)) %>%
    rowwise() %>%
    mutate(tree_comm = tree_meta %>% trim_by_community(comm) %>% list()) %>%
    mutate(tree_plot = plot_tree(tree_comm) %>% list) %>%
    select(Community = comm, everything())

p <- plot_grid(plotlist = communities_tree$tree_plot, labels = communities$Community, ncol = 2) +
    theme(plot.background = element_rect(fill = "white", color = NA))
ggsave(here::here("output/report/figure/01-tree.png"), p, width = 15, height = 15)

save(communities_tree, tree_meta, file = here::here("data/output/tree.Rdata"))








if (FALSE) {
tree_meta %>%
    #ggtree(branch.length = "none", layout = "circular") +
    ggtree(branch.length = "none") +
    geom_tippoint(aes(colour = Fermenter)) +
    scale_color_npg(name = "", label = c("TRUE" = "Fermenter", "FALSE" = "Respirator")) +
    new_scale_color() +
    geom_tiplab(aes(label = paste0(Community, " ", Genus, " sp."))) +
    #scale_alpha_manual(values = seq(0, 1, length.out = 13)) +
    #scale_x_continuous(limits = c(0, 30)) +
    #scale_x_continuous(expand = expansion(mult = .1, add = 0)) +
    theme_bw() +
    labs()


p







}
