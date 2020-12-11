library(rvest)
library(tidygraph)
library(ggraph)
library(tidyverse)

pub_url <- "http://www.sanchezlaboratory.com/pubs"

pub_xpath <- '//*[@id="block-yui_3_17_2_1_1567798471690_21270"]/div/p'
html_content <- read_html(pub_url)
pub_list <- html_content %>%
    html_nodes(xpath = pub_xpath) %>%
    html_text()

coauthor_list <- pub_list %>%
    gsub('\\#', "", .) %>%
    gsub('\\*', "", .) %>%
    strsplit(split = "\\(\\d{4}\\)") %>% sapply(function(x) x[2]) %>%
    gsub("^\\.\\s", "", .) %>%
    gsub("^\\s", "", .) %>%
    strsplit(split = "\\(\\d+\\)") %>% sapply(function(x) x[1]) %>%
     strsplit(split = "\\. mSystem") %>% sapply(function(x) x[1]) %>%
     strsplit(split = "\\. mSystem") %>% sapply(function(x) x[1]) %>%
    strsplit(split = "\\. PNAS") %>% sapply(function(x) x[1]) %>%
    strsplit(split = "\\. Cell System") %>% sapply(function(x) x[1]) %>%
    strsplit(split = "\\. Lecture Notes in Computer Science") %>% sapply(function(x) x[1]) %>%
    strsplit(split = "\\. Biophysical Journal") %>% sapply(function(x) x[1]) %>%
    strsplit(split = "\\. Natural Computing") %>% sapply(function(x) x[1]) %>%
    strsplit(split = "\\. Science") %>% sapply(function(x) x[1]) %>%
    strsplit(split = "\\. ISME") %>% sapply(function(x) x[1]) %>%
    strsplit(split = "\\. PLOS") %>% sapply(function(x) x[1]) %>%
    strsplit(split = "\\.\\s+PLoS") %>% sapply(function(x) x[1]) %>%
    strsplit(split = "\\. In review") %>% sapply(function(x) x[1]) %>%
    strsplit(split = "\\. BioRxiv") %>% sapply(function(x) x[1]) %>%
    strsplit(split = " BioRxiv") %>% sapply(function(x) x[1]) %>%
    strsplit(split = "\\. arXiv") %>% sapply(function(x) x[1]) %>%
    strsplit(split = "\\. Ecoevo") %>% sapply(function(x) x[1]) %>%
    strsplit(split = "\\. Evolution") %>% sapply(function(x) x[1]) %>%
    strsplit(split = "\\. Journal") %>% sapply(function(x) x[1]) %>%
    strsplit(split = "\\. Current Opinion") %>% sapply(function(x) x[1]) %>%
    {.}


sanchez_lab <- c(
    "Chang-Yu Chang",
    "Jean Vila",
    "Maria Rebolleda-Gomez",
    "Nanxi Lu",
    "Sylvie Estrela",
    "Alicia Sanchez-Gorostiaga",
    "Djordje Bajic",
    "Alvaro Sanchez",
    "Madeline Bender",
    "Juan Diaz-Colunga",
    "Nora Pyenson",
    "Jackie Folmar"
)

temp_list <- coauthor_list %>% strsplit(split = ",")
for (i in 1:length(temp_list)) {
    x = temp_list[[i]]
    if (i == 4) x <- c("Chang-Yu Chang", "Jean C.C. Vila", x[-1])
    if (any(grepl("&", x))) temp_list[[i]] <- strsplit(x, "&") %>% unlist()
    if (any(grepl("and", x))) temp_list[[i]] <- strsplit(x, "and") %>% unlist()

}

paper_author_list <-  temp_list %>%
    lapply(function(x) tibble(Author = x)) %>%
    bind_rows(.id = "Paper") %>%
    mutate(Author = gsub("^\\s+", "", Author)) %>%
    mutate(Author = gsub("\\s+$", "", Author))

paper_author_list_sanchez <- paper_author_list %>%
    mutate(Author = ifelse(Author %in% c("Jean CC Vila", "Jean C.C. Vila"), "Jean Vila", Author)) %>%
    mutate(Author = ifelse(Author %in% c("María Rebolleda-Gomez"), "Maria Rebolleda-Gomez", Author)) %>%
    mutate(Author = ifelse(Author %in% c("Djordje Bajić", "Djordje Bajic"), "Djordje Bajic", Author)) %>%
    filter(Author %in% sanchez_lab)

df_temp <- paper_author_list_sanchez %>%
    group_by(Paper) %>%
    mutate(AuthorNumber = n()) %>%
    filter(AuthorNumber >= 2)

df_edges <- df_temp %>%
    split.data.frame(f = .$Paper) %>%
    lapply(function(x) {
        temp <- t(combn(x$Author, 2))
        tibble(from = temp[,1], to = temp[,2])
        }) %>%
    bind_rows(.id = "Paper") %>%
    arrange(from, to) %>%
    group_by(from, to) %>%
    summarise(CoauthoredPaper = n())


graph <- tbl_graph(nodes = tibble(name = sanchez_lab), edges = df_edges)

node_size = 15
p1 <- graph %>%
    ggraph(layout = "circle") +
    geom_edge_link(aes(edge_width = CoauthoredPaper), color = "grey30", alpha = 0.5,
        start_cap = circle(node_size/2+1, "mm"),
        end_cap = circle(node_size/2+1, "mm")) +
    geom_node_point(size = node_size, fill = "white", color = "grey30", shape = 21) +
    geom_node_text(aes(label = name), size = 5, color = "black") +
    theme_graph() +
    scale_x_continuous(limits = c(-2, 2)) +
    scale_y_continuous(limits = c(-2, 2)) +
    theme(legend.position = "none")
#        plot.margin = unit(rep(10,4), units = "mm"))

ggsave("../plots/lab_network.png", plot = p1, width = 10, height = 10)

#nodes <- tibble(Name <- c("Alvaro", "Djordje", "Alicia", "Sylvie", "Nanxi", "Chang-Yu", "Jean", "Nora", "Juan"))

if (FALSE) {
    #gs_id = "vwkZIIMAAAAJ"
    gs_id = "dt8_0JAAAAAJ"
    get_profile(id = gs_id)
    coauthor_network <- get_coauthors(gs_id, n_coauthors = 20)

    sanchez_lab <- c("Chang-Yu Chang",
        "Nanxi Lu",
        "Jean Vila",
        "Sylvie Estrela",
        "Alicia Sanchez-Gorostiaga",
        "Djordje Bajić",
        "Alvaro Sanchez",
        "Maria Rebolleda-Gomez"
    )

    coauthor_network_cleaned <- coauthor_network %>%
        filter(!(coauthors %in% c("Sort By Title", "Sort By Citations", " Sort By Year")),
            !(author %in% c("Sort By Title", "Sort By Citations", " Sort By Year"))) %>%
        filter(author %in% sanchez_lab, coauthors %in% sanchez_lab)


    graph <- tbl_graph(edges = coauthor_network_cleaned)

    graph %>%
        ggraph(layout = "circle") +
        geom_node_point(aes(label = name), size = 10, color = "white") +
        geom_node_text(aes(label = name), size = 5) +
        geom_edge_link() +
        theme_graph()
}
