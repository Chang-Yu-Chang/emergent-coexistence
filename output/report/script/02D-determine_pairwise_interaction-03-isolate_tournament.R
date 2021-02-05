# Tournament ranks

library(tidyverse)
library(data.table)
tournament_rank <- function(pairs) {
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

#
isolates_ID_match <- fread(here::here("data/temp/isolates_ID_match.csv"))
pairs_interaction <- fread(here::here("data/temp/pairs_interaction.csv"))
communities <- fread(here::here("data/output/communities.csv"))
communities_name <- communities$Community
tournaments <- rep(list(NA), length(communities_name)) %>% setNames(communities_name)

for (i in 1:length(communities_name)) {
    tournaments[[i]] <- tournament_rank(pairs = filter(pairs_interaction, Community == communities_name[i]))
}
isolates_tournament <- tournaments %>%
    rbindlist(idcol = "Community") %>%
    arrange(Community, Isolate) %>%
    as_tibble()

fwrite(isolates_tournament, file = here::here("data/temp/isolates_tournament.csv"))
