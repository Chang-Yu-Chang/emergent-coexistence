#' Plot isolates distributions

# Fitted growth rates distribution
p_isolates_growth_rate <- isolates_melted %>%
  filter(Family %in% c("Enterobacteriaceae", "Pseudomonadaceae")) %>%
  ggplot(aes(group = Family)) +
  geom_density(aes(x = MaxGrowthRate, y = ..density.., col = Family), lwd = 2) +
  facet_wrap(CarbonSource~., scale = "free", ncol = 3) +
  theme_bw() +
  theme(legend.position = "top")
