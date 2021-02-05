#' pairwise-competition_04-network-function.R
#' Description: this script tidy up the pairwise and isolate data for later analyses and visualizations
#' Input objects: community_name, community_size

# Take one community for example.
comm = "C1R7"
size_temp <- community_size[community_name==comm]
m <- matrix(NA, size_temp, size_temp)
pairs2_temp <- pairs2 %>% filter(Community == comm)
tour <- as.data.frame(matrix(NA, size_temp, 4)) %>% setNames(c("Isolate", "Win", "Lose", "Draw"))
tour$Isolate <- 1:size_temp

# Tournament and calculate the competitive score. The code is the same as in `tournament()` ----
tour <- data.frame(Community=comm, Isolate=1:size_temp,
                   # Win
                   Win=filter(pairs2_temp, InteractionType == "exclusion") %>%
                     select(From) %>% unlist() %>% factor(1:7) %>%table() %>% as.vector(),
                   # Lose
                   Lose=filter(pairs2_temp, InteractionType == "exclusion") %>%
                     select(To) %>% unlist() %>% factor(1:7) %>%table() %>% as.vector(), 
                   # Draw
                   Draw=filter(pairs2_temp, InteractionType == "coexistence") %>%
                     select(From, To) %>% unlist() %>% factor(1:7) %>%table()%>% as.vector())

tour <- mutate(tour, Score=Win-Lose, Game=Win+Lose+Draw) %>%
  arrange(desc(Score))
rank_temp <- tour$Isolate

# Make network. The code is the same as `network_make()` ----
# Nodes
nodes <- filter(isolates2, Community==comm) %>% select(Isolate, everything())

# Edges
edges <- mutate(pairs2_temp, from=From, to=To) %>% select(from, to, InteractionType)
edges_coext <- edges[edges$InteractionType == "coexistence",]; edges_coext[,c("from", "to")]  <- edges_coext[,c("to", "from")] # Add the mutual edges for coexistence links
edges <- rbind(edges, edges_coext)

# Make network
net <- graph_from_data_frame(vertices=nodes, d=edges, directed = T)

# Edges attributes
edge_colrs <- c("coral", NA, "grey90")
E(net)$arrow.size <- .2
E(net)$edge.color[E(net)$InteractionType == "exclusion"] <- 1
E(net)$edge.color[E(net)$InteractionType == "coexistence"] <- 2
E(net)$edge.color[is.na(E(net)$InteractionType)] <- 3


# Plot network. The code is the same as `network_plot()` ----
E(net)$edge.color <- ifelse(which_mutual(net), "blue", "red")
plot.igraph(net, vertex.size = 30, vertex.color = "grey80",
            edge.arrow.size = .5, edge.arrow.width = 1.5, 
            edge.width = 2, edge.color = E(net)$edge.color, 
            layout = do.call("layout_in_circle", list(net)))

# Pairwise matrix. The code is the same as `adj_from_net()` ----
# Get adjacent matrix
net_m <- get.adjacency(net, attr = "InteractionType", sparse = F)

# Win
temp <- which(net_m=="exclusion", arr.ind = T) %>% as.data.frame()
for(i in 1:nrow(temp)) net_m[temp$col[i], temp$row[i]] <- "lose"
net_m[net_m=="exclusion"] <- "win"

# Diagonal 
diag(net_m) <- "self"
net_m[net_m == ""] <- NA

# Plot pairwise matrix. The code is the same as `adj_matrix_plot()` ----
myColor = c("#A6C966", "#557BAA", "#DD6E74", "black")
names(myColor) <- c("win", "coexistence", "lose", "self")

colnames(net_m) <- paste0(V(net)$name)#, " ", V(net)$Family)
rownames(net_m) <- paste0(V(net)$name)#, " ", V(net)$Family)
isolate_rank <- paste0(as.character(rank_temp))#, " ", V(net)$Family[match(as.character(rank_temp), V(net)$name)])

t(net_m) %>%
  melt(.) %>%
  mutate(Var1 = ordered(Var1, level = isolate_rank),
         Var2 = ordered(Var2, level = rev(isolate_rank)),
         value = ordered(value, level = c("win", "coexistence", "lose", "self"))) %>%
  ggplot(aes(x = Var1, y = Var2, fill = value)) +
  geom_tile() +
  scale_fill_manual(values =myColor) +
  # Add grid
  geom_vline(xintercept = seq(.5, length(V(net))-.5, 1), color = "white") +
  geom_hline(yintercept = seq(.5, length(V(net))-.5, 1), color = "white") +
  scale_x_discrete(position = "top") + # Move x axis tick to top
  theme(
    axis.text = element_blank(),
    axis.title.x = element_blank(),
    #    axis.text.x = element_text(angle = 90, hjust = 1),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "top",
    legend.direction = "horizontal",
    legend.title = element_blank()) +
  ggtitle(comm)



