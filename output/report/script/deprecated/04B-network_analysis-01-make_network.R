#' Competitive ranks and network plot

# Read data
df_pair_interaction <- fread(root$find_file("output/report/data/temp/df_pair_interaction.txt"))

# Prepare pairs and isolates list for making network ----
## Isolate list for making network
df_isolates <- df_pair_interaction %>%
  select(Community, Isolate1, Isolate2) %>%
  pivot_longer(cols = starts_with("Isolate"), names_to = "temp", values_to = "Isolate") %>%
  select(Community, Isolate) %>%
  distinct(Community, Isolate) %>%
  mutate(ID = ordered(Isolate))


# Function for computing the isolate competitive ranks and plotting matrix ----
#' In the experimental data, isolates are called 1, 2, 3, to whatever the size of the community
#' However, in the simulated communities, the species are called by their exact ID in the species pool
#' To adapt for this issue, I modified two network functions from `misc/network_functions` by
#' replacing the variable names "Isolate" by "ID"

## Calculate isolate ranks
tournament_rank <- function(net) {
  # Extract edges from igraph object
  temp_pairs_interaction <-
    igraph::as_data_frame(net) %>%
    filter(!(InteractionType == "coexistence" & (from > to))) %>%
    select(from, to, InteractionType)

  # Community name and size
  comm <- unique(V(net)$Community)
  comm_size <- length(V(net))
  isolate_list <- ordered(V(net)$name, 0:10000)

  # Isolates' ranks in the tournament
  tour_rank <- data.frame(
    Community = comm,
    ID = isolate_list,
    # Win
    Win = filter(temp_pairs_interaction, InteractionType == "exclusion") %>%
      select(from) %>% unlist() %>% factor(isolate_list) %>%table() %>% as.vector(),
    # Lose
    Lose = filter(temp_pairs_interaction, InteractionType == "exclusion") %>%
      select(to) %>% unlist() %>% factor(isolate_list) %>%table() %>% as.vector(),
    # Draw; Note that I consider neturality and bistability as draw in the tournament
    Draw = filter(temp_pairs_interaction, InteractionType %in% c("coexistence", "neutrality", "bistability")) %>%
      select(from, to) %>% unlist() %>% factor(isolate_list) %>%table()%>% as.vector())

  # Arrange the df by score
  tour_rank <- tour_rank %>%
    mutate(Score = Win - Lose + 0 * Draw, Game = Win + Lose + Draw) %>%
    arrange(desc(Score))

  # Calcluate rank by score; same scores means the same ranks
  temp_score <- ordered(tour_rank$Score, levels = sort(unique(tour_rank$Score), decreasing = T))
  temp_score_table <- table(temp_score)
  temp <- NULL; temp_counter = 1
  for (i in 1:length(temp_score_table)) {
    temp <- c(temp, rep(temp_counter, temp_score_table[i]))
    temp_counter <- temp_counter + temp_score_table[i]
  }

  tour_rank$Rank <- temp
  tour_rank$PlotRank <- 1:comm_size
  return(tour_rank)
}

## Plot pairwise matrix from igraph object, matrix can be ordered in chosen isolate traits
adj_matrix_plot <- function (
  net, # igraph object
  isolate_label = c("name", "ExpID", "Family", "Fermenter", "GramPositive"),  # Isolate labels to be shown on the pairwise matrix plot
  axis_order_by = "competitive_rank", # How to order the axis. Other option are "ID", "OD620", or other growth traits
  show_title_legend = TRUE,
  show_axis_label = TRUE,
  show_half_matrix = FALSE,
  standardize_tile_size = FALSE
) {
  # Make pairwise matrix ----
  # Community name and size
  comm <- unique(V(net)$Community)

  # Make pairwise matrix from igraph object
  net_m <- adj_from_net(net)

  ## Pairwise matrix axis. Now the axis is ordered by isolate ID
  #' To order the axis by isolate traits, first I convert the logical variables (e.g., Fermenter) into binary characters
  if (length(isolate_label) != 0 & any(is.logical(vertex_attr(net)))) {

    temp_index_logical <- which(lapply(vertex_attr(net)[isolate_label], is.logical) %>% unlist)
    for (i in 1:length(temp_index_logical)) {
      vertex_attr(net)[isolate_label[temp_index_logical]][[i]] <-
        ifelse(vertex_attr(net)[isolate_label[temp_index_logical]][[i]],
          names(temp_index_logical)[i],
          paste0("Non-", names(temp_index_logical)[i]))
    }

    ### Make axis label by pasting isolate variables
    isolates_label_axis <-
      vertex_attr(net)[isolate_label] %>% unlist() %>%
      matrix(ncol = length(isolate_label)) %>%
      apply(1, paste0, collapse = "_") # Paste the isolate variables into a vector of single characters

  } else isolates_label_axis <- V(net)$name
  # Sort matrix and assign matrix axis labels ----
  ## Matrix axis labels
  colnames(net_m) <- isolates_label_axis
  rownames(net_m) <- isolates_label_axis

  ## Compute isolate ranks based on growth traits or competitive scores
  ## Isolate ranks; the order for plotting isolates is when `PlotRank` is sorted
  if (axis_order_by == "ID") {
    isolates_rank_axis <- isolates_label_axis
  } else if (axis_order_by == "competitive_rank") {
    temp_tournament <- tournament_rank(net) %>% arrange(Community, ID)
    isolates_rank_axis <- isolates_label_axis[order(temp_tournament$PlotRank)]
  } else {
    temp_traits <- isolate_trait_rank(net, axis_order_by) %>% arrange(Community, ID)
    isolates_rank_axis <- isolates_label_axis[order(temp_traits$TraitPlotRank)] # Isolate rank
  }

  ## Update the matix axis labels by ranks
  temp_order <- order(isolates_rank_axis, isolates_label_axis)
  net_m <- net_m[temp_order, temp_order]

  # Plot pairwise matrix by `geom_tile` ----
  ## Fill colors for tiles by interaction types
  interaction_type <- c("win", "coexistence", "lose", "bistability", "neutrality", "self", "undefined")
  #myColor = c("#73C966", "#557BAA", "#DB7469", "#EECF6D", "#8650C4", "black", "grey80")
  myColor = c("#DB7469", "#557BAA", "#73C966", "#EECF6D", "#8650C4", "black", "grey80")
  names(myColor) <- interaction_type

  ## Reformat the matrix into a df for to be plot
  net_m_df <- t(net_m) %>% melt(.) %>%
    # Ordered the matrix for filter
    mutate(Var1 = ordered(Var1, level = isolates_rank_axis),
      Var2 = ordered(Var2, level = isolates_rank_axis)) %>%
    # Conditonal filter; filter when only showing the upper half of the matrix
    purrr::when (
      show_half_matrix == T ~ filter(., Var1 >= Var2),
      ~ .
    ) %>%
    # Order the axis labels
    mutate(
      Var1 = ordered(Var1, level = isolates_rank_axis),
      # Reverse Var2 for axis direction
      Var2 = ordered(Var2, level = rev(isolates_rank_axis)),
      value = ordered(value, level = interaction_type)
    ) %>%
    as_tibble()

  ## Add empty tiles to standardize the tile sizes across communities
  if (standardize_tile_size) {
    temp_size_diff <- 13 - length(isolates_label_axis)
    net_m_df <- net_m_df %>%
      mutate(Var1 = ordered(Var1, level = c(isolates_rank_axis, 1:temp_size_diff)),
        Var2 = ordered(Var2, level = rev(c(isolates_rank_axis, 1:temp_size_diff)))) %>%
      rbind(data.frame(Var1 = rep(1:temp_size_diff), Var2 = rep(1:temp_size_diff), value = NA))
  }

  ## Plotting
  pp <- net_m_df %>%
    ggplot() +
    geom_tile(aes(x = Var1, y = Var2, fill = value), color = "white", lwd = 1) +
    scale_fill_manual(values = myColor) +
    # Move x axis tick to top
    scale_x_discrete(position = "top", ) +
    theme(
      axis.title = element_blank(),
      axis.text = element_blank(),
      panel.grid.major = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank(),
      axis.ticks = element_blank(),
      legend.position = "none",
      legend.direction = "none",
      legend.title = element_blank()) +
    # Conditional theme; show the axis labels
    {if (show_axis_label) theme(axis.text.y = element_text(hjust = 1))} +
    # Conditonal theme; show title and labels
    {if (show_title_legend) {
      theme(legend.position = "top",legend.direction = "horizontal")
      ggtitle(paste0(comm, ", axis order by ", axis_order_by))
    }} +
    NULL
  return(pp)
}


# Make networks ----
communities_name <- unique(df_isolates$Community)
net_list <- rep(list(NA), length(communities_name))
names(net_list) <- communities_name

for (i in 1:length(net_list)) net_list[[i]] <- network_make(df_isolates, df_pair_interaction, comm = communities_name[i])

# Plot the networks ----
p_net_matrix_list <- rep(list(NA), length(net_list))
names(p_net_matrix_list) <- communities_name

for (i in 1:length(net_list)) p_net_matrix_list[[i]] <- adj_matrix_plot(net_list[[i]], isolate_label = c("ID"), show_half_matrix = T, axis_order_by = "competitive_rank", standardize_tile_size = F)






