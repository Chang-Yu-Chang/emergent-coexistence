library(tidyverse)

# Process plates
plates <- read_csv(here::here("data/output/plates.csv")) %>%
    separate(col = Well, into = c("Row", "Column"), sep = 1) %>%
    mutate(Row = match(Row, LETTERS), Column = as.numeric(Column)) %>%
    mutate(Pair = paste0(Isolate1, "_", Isolate2))

plates$Community[plates$Isolate1 == "blank"] <- NA
plates$Pair[plates$Isolate1 == "blank"] <- NA


if (FALSE) {

plate1 <- plates %>%
    filter(PlateLayout == 933, Plate == "P1")

community_label <- plate1 %>%
    filter(MixIsolate == T) %>%
    group_by(Community) %>%
    summarize(Community_x = mean(Column), Community_y = mean(Row)) %>%
    mutate(Community_y = ifelse(Community_y %% 1 == 0, Community_y + 0.5, Community_y))


well_size = 5
p <- plate1 %>%
    ggplot() +
    geom_point(aes(x = Column, y = Row, fill = Community), shape = 21, size = well_size) +
    geom_text(data = community_label, aes(x = Community_x, y = Community_y, label = Community, color = Community)) +
    scale_fill_discrete(na.value = "white") +
    scale_x_continuous(breaks = 1:12, position = "top") +
    scale_y_reverse(breaks = 1:8, labels = setNames(LETTERS[1:8], 1:8)) +
    theme_minimal() +
    theme(panel.grid.minor = element_blank(), panel.border = element_rect(color = "black", fill = NA),
          plot.background = element_rect(color = NA, fill = "white"), legend.position = "none",
          axis.ticks = element_line(color = "black"),
          axis.title = element_blank()) +
    ggtitle("")
p

ggsave(here::here("output/report/figure/plate_layout.png"), p, width = 4, height = 3)


}


# Extract the label position
community_label <- plates %>%
    #filter(PlateLayout %in% c(933, 444), Plate == "P1") %>%
    filter(Plate == "P1") %>%
    filter(MixIsolate == T) %>%
    group_by(PlateLayout, Community) %>%
    mutate(PlateLayout = ordered(PlateLayout, levels = c("933", "444", "13A", "13B", "75", "5543", "C11R1"))) %>%
    summarize(Community_x = mean(Column), Community_y = mean(Row)) %>%
    mutate(Community_y = ifelse(Community_y %% 1 == 0, Community_y + 0.5, Community_y))

# Plot
well_size = 5
p <- plates %>%
    #filter(PlateLayout %in% c(933, 444), Plate == "P1") %>%
    filter(Plate == "P1") %>%
    mutate(PlateLayout = ordered(PlateLayout, levels = c("933", "444", "13A", "13B", "75", "5543", "C11R1"))) %>%
    ggplot() +
    geom_point(aes(x = Column, y = Row, fill = Community), shape = 21, size = well_size) +
    geom_text(data = community_label, aes(x = Community_x, y = Community_y, label = Community, color = Community)) +
    scale_fill_discrete(na.value = "white") +
    scale_x_continuous(breaks = 1:12, position = "bottom") +
    scale_y_reverse(breaks = 1:8, labels = setNames(LETTERS[1:8], 1:8)) +
    facet_wrap(~PlateLayout, scales = "free", ncol = 2) +
    theme_minimal() +
    theme(panel.grid.minor = element_blank(), panel.border = element_rect(color = "black", fill = NA),
          plot.background = element_rect(color = NA, fill = "white"), legend.position = "none",
          strip.text = element_text(size = 15, face = "bold"),
          axis.ticks = element_line(color = "black"),
          axis.title = element_blank()) +
    ggtitle("")

ggsave(here::here("output/report/figure/plate_layout.png"), p, width = 8, height = 12)




