# Miscellaneous R functions


# Color palatte ----
# Assign interaction link colors
assign_interaction_color <- function (level = "simple") {
    if (level == "simple") {
        interaction_type <- c("exclusion", "coexistence")
        interaction_color <- c("#DB7469", "#557BAA")
        names(interaction_color) <- interaction_type
        return(interaction_color)
    }
    if (level == "matrix") {
        interaction_type <- c("exclusion", "coexistence", "exclusion violating rank", "bistability", "neutrality", "self", "undefined")
        interaction_color <- c("#DB7469", "#557BAA", "#8CB369", "#EECF6D", "#8650C4", "black", "grey80")
        names(interaction_color) <- interaction_type
        return(interaction_color)
    }
    if (level == "finer") {
        interaction_type <- c("competitive exclusion", "stable coexistence", "mutual exclusion", "frequency-dependent coexistence", "neutrality", "exclusion violating rank")
        interaction_color <- c("#DB7469", "#557BAA", "#FFBC42", "#B9FAF8", "#8650C4", "#8CB369")
        names(interaction_color) <- interaction_type
        return(interaction_color)
    }
}

# Assign motif colors
assign_motif_color <- function () {
    motif_type <- c("others", "0-scored")
    motif_color = c("#DB7469", "#557BAA")
    names(motif_color) <- motif_type
    return(motif_color)
}

# Assign category colors
assign_category_color <- function() {
    c(sugar = "#ED6A5A", acid = "#03CEA4", waste = "#51513D", fermenter = "#8A89C0", respirator = "#FFCB77", F0 = "#8A89C0", F1 = "#FFCB77")
}

motif_color <- assign_motif_color()
interaction_color <- assign_interaction_color()
category_color <- assign_category_color()
fermenter_color <- c("fermenter" = "#8A89C0", "respirator" = "#FFCB77")
dominant_color <- c("dominant" = "grey20", "subdominant" = "grey90")
#frequency_color <- c( "95"="#1B1B3A", "50"="#693668", "5"="#D9DBF1")
frequency_color <- c( "95"="#292F36", "50"="#9F87AF", "5"="#CCCCCC")

# Paint white backgorund for plot_grid
paint_white_background <- function(x) theme(plot.background = element_rect(color = NA, fill = "white"))


# Plot ----
prepare_plate_draw <- function(plate) {
    plate %>%
        mutate(Row = match(str_sub(Well, 1, 1), LETTERS[1:8]),
               Col = str_sub(Well, 2, 3) %>% as.numeric(),
               TextLabel = ifelse(MixIsolate, paste(Isolate1, Isolate2, sep = "_"), Isolate1),
               FillLabel = paste0(Community, ifelse(MixIsolate, "", "_single"))) %>%
        unite(col = "PlateLabel", Batch, PlateLayout, Plate) %>%
        select(PlateLabel, Community, Row, Col, TextLabel, FillLabel)
}
draw_plate_empty <- function  (well_size = 6, text_size = 3) {
    tibble(Row = rep(1:8, each = 12), Col = rep(1:12, 8)) %>%
        ggplot() +
        geom_point(aes(x = Col, y = Row), shape = 1, size = well_size, color = "black", fill = NA) +
        scale_x_continuous(name = "", breaks = 1:12, labels = 1:12, limits = c(1, 12), position = "top") +
        scale_y_reverse(name = "", lim = c(8, 1), breaks = 1:8, labels = LETTERS[1:8]) +
        theme_bw() +
        theme(panel.grid.minor = element_blank())
}
draw_plate_from_df <- function(plate, well_size = 6, text_size = 3) {
    # Check that the input plate tibble has to be a whole plate
    stopifnot(is.tibble(plate), nrow(plate) == 96)

    # Remove blank
    plate <- plate %>%
        filter(!str_detect(TextLabel, "blank")) %>%
        filter(!str_detect(FillLabel, "blank"))

    empty_well <- tibble(Row = rep(1:8, each = 12), Col = rep(1:12, 8))
    # Plot segments
    g_empty_well <- geom_point(data = empty_well, aes(x = Col, y = Row),
                               shape = 1, size = well_size, color = "black", fill = NA)

    g_plate <- geom_point(data = plate, aes(x = Col, y = Row, fill = FillLabel),
                          pch = 21, size = well_size)


    # Text for fill
    plate_fill_text <- plate %>%
        group_by(FillLabel) %>%
        summarize(Col = mean(Col), Row = mean(Row)) %>%
        mutate(Row = ifelse((Row %% 1) == 0, Row - 0.5, Row))
    g_text_label <- geom_text(data = plate, aes(x = Col, y = Row, label = TextLabel), size = text_size)
    g_text_fill <- geom_text(data = plate_fill_text, aes(x = Col, y = Row, label = FillLabel))


    # Plot plate
    p <- ggplot() + g_plate + g_empty_well + g_text_label + g_text_fill +
        scale_x_continuous(name = "", breaks = 1:12, labels = 1:12, limits = c(1, 12), position = "top") +
        scale_y_reverse(name = "", lim = c(8, 1), breaks = 1:8, labels = LETTERS[1:8]) +
        theme_bw() +
        theme() +
        guides(fill = "none", color = "none")

    return(p)

}



# Network ----
# Make network from pairs and isolates data
make_network <- function(isolates, pairs) {
    # Nodes
    nodes <- isolates %>% select(Isolate, ID, Rank, PlotRank)

    # Edges
    ## Remove no-growth
    edges <- pairs %>%
        filter(InteractionType %in% c("coexistence", "exclusion")) %>%
        mutate(from=From, to=To) %>% select(from, to, InteractionType)
    edges_coext <- edges[edges$InteractionType == "coexistence",]
    edges_coext[,c("from", "to")] <- edges_coext[,c("to", "from")] # Add the mutual edges for coexistence links
    edges <- rbind(edges, edges_coext)

    # Network
    graph <- tbl_graph(nodes = nodes, edges = edges, directed = T)

    return(graph)
}

# Summarize the pairwise
summarize_network_pairs <- function(graph){
    number_nodes <- graph %>% activate(nodes) %>%
        igraph::V() %>% length()

    number_pairs <- graph %>% activate(edges) %>%
        igraph::E() %>% length()

    number_coexistence <- graph %>% activate(edges) %>%
        filter(InteractionType == "coexistence") %>%
        igraph::E() %>% length()

    summary_stat <- tibble(
        NumberNodes = number_nodes,
        FractionCoexistence = number_coexistence/number_pairs
    )

    return(summary_stat)

}

# Summarize motifs
summarize_network_motif <- function(graph) {
    number_nodes <- graph %>% activate(nodes) %>%
        igraph::V() %>% length()

    motif_counts <- count_motif(graph)

    tibble(
        CommunitySize = number_nodes,
        Motif = factor(1:7),
        Count = motif_counts,
        TotalMotifCount = rep(sum(motif_counts)),
        RelativeMotifCount = Count/ TotalMotifCount)
}

# Plot the competitive network. Take output from make_network()
plot_competitive_network <- function(g, node_size = 10, edge_width = 1, g_layout = "circle") {
    # g_layout
    if (g_layout == "hierarchy") {
        graph_layout <- create_layout(g, g_layout = "")
    } else if (g_layout == "example_motif") {
        g <- g
    } else {
        graph_layout <- create_layout(g, g_layout)
        mean_x_coord <- mean(graph_layout$x)
        mean_y_coord <- mean(graph_layout$x)
        g <- g %>% activate(nodes) %>% mutate(x = graph_layout$x - mean_x_coord, y = graph_layout$y - mean_y_coord)
    }


    # Nodes
    nodex_axis_x <- activate(g, nodes) %>% pull(x) %>% range()
    nodex_axis_y <- activate(g, nodes) %>% pull(y) %>% range()

    # Edges
    interaction_type <- c("exclusion", "coexistence", "lose", "bistability", "neutrality", "self", "undefined")
    interaction_color = c("#DB7469", "#557BAA", "#73C966", "#EECF6D", "#8650C4", "black", "grey80")
    names(interaction_color) <- interaction_type

    if (g_layout == "linear" & length( activate(g, nodes) %>% pull(x)) == 2) {
        g %>%
            mutate(Isolate = factor(Isolate)) %>%
            ggraph(layout = "nicely") +
            geom_node_point(aes(fill = Isolate), size = node_size, shape = 21, colour = "black", stroke = node_size/5) +
            geom_edge_link(aes(color = InteractionType), width = edge_width,
                           arrow = arrow(length = unit(edge_width, "mm"), type = "closed", angle = 30, ends = "last"),
                           start_cap = circle(node_size/2, "mm"),
                           end_cap = circle(node_size/2, "mm")) +
            scale_edge_color_manual(values = interaction_color) +
            #scale_fill_manual(values = c("white", "grey40")) +
            scale_x_continuous(limits = nodex_axis_x*1.2) +
            theme_graph() +
            theme(
                legend.position = "none",
                legend.direction = "none",
                legend.title = element_blank(),
                panel.background = element_blank(),
                strip.text = element_blank(),
                plot.margin=unit(c(3,3,3,3),"mm")
            )
    } else if (g_layout == "linear"){
        g %>%
            ggraph(layout = "nicely") +
            geom_node_point(fill = "grey", size = node_size, shape = 21, colour = "black", stroke = node_size/5) +
            geom_edge_arc(aes(color = InteractionType), width = edge_width,
                          arrow = arrow(length = unit(edge_width, "mm"), type = "closed", angle = 30, ends = "last"),
                          start_cap = circle(node_size/2, "mm"),
                          end_cap = circle(node_size/2, "mm")) +
            scale_edge_color_manual(values = interaction_color) +
            scale_x_continuous(limits = nodex_axis_x*1.2) +
            #scale_y_continuous(limits = nodex_axis_y*1.2) +
            theme_graph() +
            theme(
                legend.position = "none",
                legend.direction = "none",
                legend.title = element_blank(),
                panel.background = element_blank(),
                strip.text = element_blank(),
                plot.margin=unit(c(3,3,3,3),"mm")
            )
    } else {
        g %>%
            ggraph(layout = "nicely") +
            geom_node_point(fill = "grey", size = node_size, shape = 21, colour = "black", stroke = node_size/5) +
            geom_edge_link(aes(color = InteractionType), width = edge_width,
                           arrow = arrow(length = unit(edge_width, "mm"), type = "closed", angle = 30, ends = "last"),
                           start_cap = circle(node_size/2, "mm"),
                           end_cap = circle(node_size/2, "mm")) +
            scale_edge_color_manual(values = interaction_color) +
            scale_x_continuous(limits = nodex_axis_x*1.3) +
            scale_y_continuous(limits = nodex_axis_y*1.3) +
            theme_graph() +
            theme(
                legend.position = "none",
                legend.direction = "none",
                legend.title = element_blank(),
                panel.background = element_blank(),
                strip.text = element_blank(),
                plot.margin=unit(c(3,3,3,3),"mm")
            )

    }
}

# Randomize the network
randomize_network <- function(graph){
    # Step1: remove the bidirection of coexistence
    graph1 <- graph %>%
        activate(edges) %>%
        reroute(from = to, to = from, subset = from > to) %>%
        distinct()

    # Step2: shuffle interaction types
    graph2 <- graph1 %>%
        activate(edges) %>%
        mutate(InteractionType = InteractionType[order(runif(n=n()))])

    # Step3: for pairs that coexist, add back the reverse direction links
    pairs_coexistence_reversed <- graph2 %>%
        activate(edges) %>%
        filter(InteractionType == "coexistence") %>%
        reroute(from = to, to = from) %>%
        activate(nodes) %>%
        select(Isolate)

    graph_coexistence <- graph_join(graph2, pairs_coexistence_reversed, by = "Isolate") %>%
        filter(InteractionType == "coexistence")

    # Step4: for exclusionary pairs, 50% of those has a flip direction
    n_exclusion_pairs <- graph2 %>% filter(InteractionType == "exclusion") %>% pull(InteractionType) %>% length()
    graph_exclusion <- graph2 %>%
        activate(edges) %>%
        filter(InteractionType == "exclusion") %>%
        reroute(from = to, to = from, subset = sample(c(T,F), n_exclusion_pairs, replace = T, prob = c(0.5,0.5))) %>%
        activate(nodes) %>%
        select(Isolate)

    graph3 <- graph_join(graph_coexistence, graph_exclusion, by = "Isolate") %>%
        arrange(from, to)

    return(graph3)
}



# Make a random network
make_random_network <- function (
    n, # Number of nodes
    p # p denotes the probability of pairwise coexistence. 1-p denotes the probability of exclusion
) {
    edges <- as_tibble(t(combn(n,2))) %>%
        setnames(c("from", "to")) %>%
        mutate(temp = sample(1:3, size = choose(n, 2), replace = T, prob = c(p, (1-p)/2, (1-p)/2)))

    # 1: coexistence, 2: keep the direction, 3: flip the direction
    edges_coexistence <- filter(edges, temp == 1)
    edges_coexistence_another_direction <- tibble(from = edges_coexistence$to, to = edges_coexistence$from, interaction = 1)
    edges_coexistence <- tibble(from = edges_coexistence$from, to = edges_coexistence$to, interaction = 1)

    edges_exclusion_fixed <- filter(edges, temp == 2) %>% select(-temp) %>% mutate(interaction = 1)

    edges_exclusion_to_be_flipped <- filter(edges, temp == 3)
    edges_exclusion_flipped <- tibble(from = edges_exclusion_to_be_flipped$to, to = edges_exclusion_to_be_flipped$from, interaction = 1)

    edges <- bind_rows(edges_coexistence, edges_coexistence_another_direction, edges_exclusion_fixed, edges_exclusion_flipped)

    graph <- as_tbl_graph(edges, directed = T)

    return(graph)
}

# Count motif
count_motif <- function(net) igraph::triad_census(net)[c(10, 9, 12, 14, 13, 15, 16)]

# Count node degree. only coexistence links considered
count_degree <- function(net) {
    net %>%
        activate(edges) %>%
        filter(InteractionType == "coexistence") %>%
        filter(from < to) %>%
        igraph::degree()
}

# Count number of connected components. One network has one value
count_component <- function(net) {
    component_result <- net %>%
        activate(edges) %>%
        filter(InteractionType == "coexistence") %>%
        filter(from < to) %>%
        igraph::components()

    # tibble(NumberCluster = component_result$no, # number of components
    #        SizeCluster= component_result$csize) %>% # Size of components
    return(component_result$no)
}

# Count the number of pairwise coexistence as a function of rank difference
count_diag_coexistence <- function (net, isolates_comm = NULL, observation = T) {
    # Extract link attributes from the network
    pairs_comm <- as_edgelist(net) %>%
        magrittr::set_colnames(c("Isolate1", "Isolate2")) %>%
        as_tibble() %>%
        mutate(InteractionType = edge.attributes(net)$InteractionType) %>%
        filter(!(InteractionType == "coexistence" & Isolate1 > Isolate2))

    # Isolate ranks. Kept the same as the observation
    if (observation) isolates_comm <- vertex.attributes(net) %>% as_tibble()
    if (observation == F) stopifnot(!is.null(isolates_comm))

    #
    pairs_comm %>%
        # Join the isolate information
        left_join(rename_with(isolates_comm, ~paste0(., "1"), everything()), by = "Isolate1") %>%
        left_join(rename_with(isolates_comm, ~paste0(., "2"), everything()), by = "Isolate2") %>%
        # Calculate rank difference
        mutate(RankDifference = abs(Rank1 - Rank2)) %>%
        # Calculate the fraction of coexistence at each distance
        mutate(InteractionType = factor(InteractionType)) %>%
        group_by(RankDifference, InteractionType, .drop = FALSE) %>%
        summarize(Count = n(), .groups = "keep") %>%
        group_by(RankDifference) %>%
        mutate(TotalCount = sum(Count)) %>%
        mutate(Fraction = Count/TotalCount) %>%
        filter(InteractionType == "coexistence")
}

# Hierarchy metric h1
compute_hierarchy1 <- function(pairs_mock) {
    pairs_temp <- pairs_mock %>%
        select(Isolate1, Isolate2, InteractionType, From, To)
    isolates_tournament <- tournament_rank(pairs_temp) %>% select(Isolate, Score)

    pairs_temp %>%
        left_join(rename_with(isolates_tournament, ~ paste0(., "1")), by = "Isolate1") %>%
        left_join(rename_with(isolates_tournament, ~ paste0(., "2")), by = "Isolate2") %>%
        filter(InteractionType == "exclusion") %>%
        mutate(WinnerScore = ifelse(From == Isolate1, Score1, Score2),
               LoserScore = ifelse(From == Isolate1, Score2, Score1)) %>%
        mutate(FollowRank = (WinnerScore > LoserScore) %>% factor(c(T,F))) %>%
        count(FollowRank, .drop = F, name = "Count") %>%
        mutate(FractionFollowRank = Count / sum(Count)) %>%
        filter(FollowRank == T) %>%
        pull(FractionFollowRank) %>%
        return()
}
randomize_pairs1 <- function(x) {
    # Shuffle pairs
    rng <- order(runif(nrow(x),0,1))
    x$InteractionType <- x$InteractionType[rng]
    # Shuffle dominance
    rng <- sample(1:nrow(x), size = nrow(x)/2, replace = F)
    temp <- x$From[rng]
    x$From[rng] <- x$To[rng]
    x$To[rng] <- temp

    return(x)
}



# Hierarchy metric h2
compute_comp_score2 <- function(pairs_comm, target_isolate) {
    pairs_comm %>%
        filter(Isolate1 == target_isolate | Isolate2 == target_isolate) %>%
        mutate(Freq = ifelse(Isolate1 == target_isolate, Isolate1MeasuredFreq, 1-Isolate1MeasuredFreq)) %>%
        pull(Freq) %>%
        mean()
}
compute_hierarchy2 <- function(isolates_mock, pairs_mock) {
    ranking <- isolates_mock %>%
        rowwise() %>%
        mutate(CompetitiveScore = compute_comp_score2(pairs_mock, Isolate)) %>%
        arrange(desc(CompetitiveScore)) %>%
        # Ranking by competitive score
        pull(Isolate)

    pairs_mock %>%
        mutate(Freq1 = Isolate1MeasuredFreq, Freq2 = 1-Freq1) %>%
        select(Isolate1, Isolate2, Freq1, Freq2) %>%
        mutate(Pair = 1:n()) %>%
        pivot_longer(cols = c(-Pair), names_to = ".value", names_pattern = "(.+)[12]") %>%
        mutate(Isolate = ordered(Isolate, ranking)) %>%
        group_by(Pair) %>%
        arrange(Pair, Isolate) %>%
        slice(1) %>%
        pull(Freq) %>%
        mean() %>%
        return()

}
randomize_pairs2 <- function(x) {
    # Shuffle pairs
    rng <- order(runif(nrow(x),0,1))
    x$Isolate1MeasuredFreq <- x$Isolate1MeasuredFreq[rng]
    # Shuffle dominance
    rng <- sample(1:nrow(x), size = nrow(x)/2, replace = F)
    x$Isolate1MeasuredFreq[rng] <- 1 - x$Isolate1MeasuredFreq[rng]

    return(x)
}



# Plot adjacent matrix
plot_adjacent_matrix <- function(graph, show.legend = F, show.axis = F, show_label = "ID") {
    graph_ranked <- graph %>%
        activate(nodes) %>%
        select(Isolate, PlotRank) %>%
        activate(edges) %>%
        mutate(fromRank = .N()$PlotRank[match(from, .N()$Isolate)],
               toRank = .N()$PlotRank[match(to, .N()$Isolate)])

    n_nodes <- igraph::vcount(graph_ranked)
    n_exclusion_violation <- graph_ranked %>%
        activate(edges) %>%
        filter(fromRank > toRank, InteractionType == "exclusion") %>%
        igraph::ecount()

    if (n_exclusion_violation != 0) {
        pairs_exclusion_violation <- graph_ranked %>%
            activate(edges) %>%
            filter(fromRank > toRank, InteractionType == "exclusion") %>%
            mutate(temp = fromRank, fromRank = toRank, toRank = temp) %>%
            select(-temp) %>%
            mutate(InteractionType = "exclusion violating rank") %>%
            as_tibble

        graph_ranked <- graph_ranked %>%
            activate(edges) %>%
            filter(fromRank <= toRank) %>%
            bind_edges(pairs_exclusion_violation)
    }

    interaction_color <- assign_interaction_color("matrix")
    graph_ranked %>%
        filter(fromRank <= toRank) %>%
        bind_edges(tibble(from = 1:n_nodes, to = 1:n_nodes, fromRank = 1:n_nodes, toRank = 1:n_nodes, InteractionType = "self")) %>%
        as_tibble() %>%
        ggplot() +
        geom_tile(aes(x = toRank, y = fromRank, fill = InteractionType), width = 0.9, height = 0.9) +
        scale_y_reverse(breaks = 1:n_nodes) +
        scale_fill_manual(values = interaction_color) +
        {if (show.axis) { theme_minimal()} else theme_void()} +
        {if (show.legend) {theme(legend.position = "top")} else theme(legend.position = "none")} +
        NULL

}


# Compute the rank of an isolate based on its number of wins and lose in pairwise competition
tournament_rank <- function(pairs) {
    if (igraph::is.igraph(pairs)) {
        net <- pairs
        pairs <- as_tibble(get.edgelist(net)) %>%
            setNames(c("From", "To")) %>%
            mutate(InteractionType = get.edge.attribute(net)$InteractionType) %>%
            rowwise() %>%
            mutate(Isolate1 = min(From, To), Isolate2 = max(From, To))
    }
    isolate_name <- pairs %>% select(Isolate1, Isolate2) %>% unlist %>% unique %>% sort()
    # Isolates' ranks in the tournament
    tour_rank <- data.frame(
        Isolate = isolate_name,
        # Win
        Win = filter(pairs, InteractionType == "exclusion") %>%
            select(From) %>% unlist() %>% factor(isolate_name) %>% table() %>% as.vector(),
        # Lose
        Lose = filter(pairs, InteractionType == "exclusion") %>%
            select(To) %>% unlist() %>% factor(isolate_name) %>% table() %>% as.vector(),
        # Draw; Note that I consider neturality and bistability as draw in the tournament
        Draw = filter(pairs, InteractionType %in% c("coexistence", "neutrality", "bistability")) %>%
            select(From, To) %>% unlist() %>% factor(isolate_name) %>% table() %>% as.vector())

    # Arrange the df by score
    tour_rank <- tour_rank %>%
        mutate(Score = Win - Lose + 0 * Draw, Game = Win + Lose + Draw) %>%
        arrange(desc(Score))

    # Calculate rank by score; same scores means the same ranks
    temp_score <- ordered(tour_rank$Score, levels = sort(unique(tour_rank$Score), decreasing = T))
    temp_score_table <- table(temp_score)
    temp <- NULL; temp_counter = 1
    for (i in 1:length(temp_score_table)) {
        temp <- c(temp, rep(temp_counter, temp_score_table[i]))
        temp_counter <- temp_counter + temp_score_table[i]
    }

    tour_rank$Rank <- temp
    tour_rank$PlotRank <- 1:nrow(tour_rank)
    return(tour_rank)
}


