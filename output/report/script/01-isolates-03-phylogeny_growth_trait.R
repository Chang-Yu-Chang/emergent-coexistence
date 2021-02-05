#' Plot the isolates phylogenetic tree, competitive rank and growth traits

# Read data
isolates <- fread(root$find_file("data/output/isolates.csv")) %>% as_tibble()
load(root$find_file("data/temp/isolates_sanger_seq.Rdata"))

# Mak tree ----
## Make tree for all 68 isolates used in competitive assays
tree

## Make tree for only large communities
isolates_seq_large <- isolates_seq[names(isolates_seq) %in% unlist(isolates[isolates$Community %in% c("C1R7", "C11R1", "C11R2"), "ExpID"])]
aln <- msa::msa(isolates_seq_large)
aln2 <- msaConvert(aln, type = "seqinr::alignment")
d <- seqinr::dist.alignment(aln2, "identity")
tree_large <- ape::nj(d) # Neighbor-joining tree estimation

# Selected traits to plot
selected_traits <- c("OD620", "Rank", "glucose_LagTime", "glucose_MaxGrowthRate", "glucose_MaxOD")
selected_traits_short <- c("ODe", "rank", "lag", "rate", "ODg")

# Plot tree ----
## Plot all isolates
isolates_subset <- scale(isolates[,selected_traits]) %>% as.data.frame() %>% setNames(selected_traits_short)
rownames(isolates_subset) <- isolates$ExpID

p_tree <- ggtree(tree) %<+% isolates +
  geom_tiplab(aes(label = paste0(Community, " ", Isolate, " ", Genus), color = Family), align = T, size = 3)  +
  theme(legend.position = "top") +
  geom_hilight(node = 109, fill = "pink", alpha = 0.5) +
  NULL

p_heat <- gheatmap(p_tree, isolates_subset, width=0.5, font.size=3, offset = .2,
  colnames_position = "top", colnames_offset_y = .1, colnames = T,
  low = "white", high = "red") +
  NULL

# Plot isolates from large communities
isolates_large <- filter(isolates, Community %in% c("C1R7", "C11R1", "C11R2"))
isolates_large_subset <- scale(isolates_large[,selected_traits]) %>% as.data.frame() %>% setNames(selected_traits_short)
rownames(isolates_large_subset) <- isolates_large$ExpID

p_tree_large <- ggtree(tree_large) %<+% isolates_large +
  geom_tiplab(aes(label = paste0(Community, " ", Isolate, " ", Genus), color = Family), align = T, size = 3)  +
  theme(legend.position = "top") +
  geom_hilight(node = 109, fill = "pink", alpha = 0.5) +
  NULL

p_heat_large <- gheatmap(p_tree_large, isolates_large_subset, width=0.5, font.size=3, offset = .2,
  colnames_position = "top", colnames_offset_y = .1, colnames = T,
  low = "white", high = "red") +
  NULL












