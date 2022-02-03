# R functions making, randomizing, plotting network

# Assign interaction link colors
assign_interaction_color <- function (level = "simple") {
    if (level == "simple") {
        interaction_type <- c("exclusion", "coexistence")
        interaction_color <- c("#DB7469", "#557BAA")
        names(interaction_color) <- interaction_type
        return(interaction_color)
    }
    if (level == "matrix") {
        interaction_type <- c("exclusion", "coexistence", "lose", "bistability", "neutrality", "self", "undefined")
        interaction_color <- c("#DB7469", "#557BAA", "#73C966", "#EECF6D", "#8650C4", "black", "grey80")
        names(interaction_color) <- interaction_type
        return(interaction_color)
    }
    if (level == "finer") {
        interaction_type <- c("competitive exclusion", "stable coexistence", "neutrality", "mutual exclusion", "frequency-dependent coexistence")
        interaction_color <- c("#DB7469", "#557BAA", "#8650C4", "#F7AEF8", "#72DDF7")
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


# Make network from pairs and isolates data
make_network <- function(isolates, pairs) {
    # Nodes
    nodes <- isolates %>% select(Isolate, Rank, PlotRank)

    # Edges
    edges <- pairs %>% mutate(from=From, to=To) %>% select(from, to, InteractionType) #Isolate1, Isolate2, Isolate1Freq)
    edges_coext <- edges[edges$InteractionType == "coexistence",]
    edges_coext[,c("from", "to")] <- edges_coext[,c("to", "from")] # Add the mutual edges for coexistence links
    edges <- rbind(edges, edges_coext)

    # Network
    graph <- tbl_graph(nodes = nodes , edges = edges, directed = T)

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
plot_competitive_network <- function(g, node_size = 10, g_layout = "circle") {
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
            geom_edge_link(aes(color = InteractionType, fill = InteractionType), width = node_size/10,
                           arrow = arrow(length = unit(node_size/2, "mm"), type = "closed", angle = 30, ends = "last"),
                           start_cap = circle(node_size/2+1, "mm"),
                           end_cap = circle(node_size/2+1, "mm")) +
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
            geom_edge_arc(aes(color = InteractionType, fill = InteractionType), width = node_size/10,
                          arrow = arrow(length = unit(node_size/2, "mm"), type = "closed", angle = 30, ends = "last"),
                          start_cap = circle(node_size/2+1, "mm"),
                          end_cap = circle(node_size/2+1, "mm")) +
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
            geom_edge_link(aes(color = InteractionType, fill = InteractionType), width = node_size/10,
                           arrow = arrow(length = unit(node_size/2, "mm"), type = "closed", angle = 30, ends = "last"),
                           start_cap = circle(node_size/2+1, "mm"),
                           end_cap = circle(node_size/2+1, "mm")) +
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

# Plot adjacent matrix
plot_adjacent_matrix <- function(graph, show.legend = F, show.axis = F, show_label = "ID") {
    graph_ranked <- graph %>%
        activate(nodes) %>%
        select(Isolate, PlotRank) %>%
        activate(edges) %>%
        mutate(
            fromRank = .N()$PlotRank[match(from, .N()$Isolate)],
            toRank = .N()$PlotRank[match(to, .N()$Isolate)])

    # fromLabel = .N()[match(from, .N()$Isolate), show_label],
    # toLabel = .N()[match(to, .N()$Isolate), show_label])

    #label_ID <- igraph::get.vertex.attribute(graph)$ID %>% setNames(igraph::get.vertex.attribute(graph)$PlotRank)

    n_nodes <- igraph::vcount(graph_ranked)
    interaction_type <- c("exclusion", "coexistence", "lose", "bistability", "neutrality", "self", "undefined")
    interaction_color = c("#DB7469", "#557BAA", "#73C966", "#EECF6D", "#8650C4", "black", "grey80")
    names(interaction_color) <- interaction_type

    graph_ranked %>%
        filter(fromRank <= toRank) %>%
        #mutate(fromLabel = factor(fromLabel), toLabel = factor(toLabel)) %>%
        bind_edges(tibble(from = 1:n_nodes, to = 1:n_nodes, fromRank = 1:n_nodes, toRank = 1:n_nodes, InteractionType = "self")) %>%
        as_tibble() %>%
        ggplot() +
        geom_tile(aes(x = toRank, y = fromRank, fill = InteractionType), width = 0.9, height = 0.9) +
        #geom_tile(aes(x = toLabel, y = fromLabel, fill = InteractionType), width = 0.9, height = 0.9) +
        #scale_x_discrete(labels= label_ID) %>%
        #scale_x_continuous(position = "top", breaks = 1:n_nodes) +
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


