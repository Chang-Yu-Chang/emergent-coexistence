# FigS3

pairs <- fread("../data/output/pairs.csv")
p <- pairs %>%
    mutate(PairFermenter = ifelse(PairFermenter == "FN", "Fermenter-Respirator",
                                  ifelse(PairFermenter == "FF", "Fermenter-Fermenter",
                                         ifelse(PairFermenter == "NN", "Respirator-Respirator", PairFermenter)))) %>%
    filter(PairFermenter != "") %>%
    filter(InteractionType != "neutrality") %>%
    group_by(InteractionType, PairFermenter) %>%
    summarize(Count = n()) %>%
    group_by(PairFermenter) %>%
    mutate(RelativeCount = Count / sum(Count) * 100, Count = sum(Count)) %>%
    ggplot() +
    geom_col(aes(x = PairFermenter, y = RelativeCount, fill = InteractionType)) +
    geom_text(aes(x = PairFermenter, label = paste0("n=", Count)), y = 100, vjust = 1.5, color = "white") +
    scale_fill_manual(values = c("coexistence" = "#557BAA", "exclusion" = "#DB7469")) +
    scale_y_continuous(expand = c(0,0)) +
    theme_cowplot() +
    theme(legend.title = element_blank(), axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    labs(x = "", y = "Fraction (%)")

ggsave("../plots/FigS3.png", p, width = 5, height = 4)

