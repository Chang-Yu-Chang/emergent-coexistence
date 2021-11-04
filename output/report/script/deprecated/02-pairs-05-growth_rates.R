#' Compare isolate monoculture growth curve charateristics (e.g., growth rate,
#' maximal OD, and lag time) to outcome of pairwise competition

# Match isolate growth rates to pairwise interactions ----
growth_rate_variables <- colnames(isolates) %>% grep("_", ., value = T) # 3 growth curve variables, 11 carbon sources
pairs_growth_rate <- fread(root$find_file("data/temp/pairs_interaction.csv"))

for (i in 1:length(growth_rate_variables)) {
  # Match isolates growth rate variables
  pairs_growth_rate <- pairs_growth_rate %>%
    invnet::match_isolates_pairs(isolates_df = isolates,
      variable_name = growth_rate_variables[i]) %>%
    as.data.table()

  # Difference in growth rate variables; always the value in isolate 1 subtracted by value in isolate 2
  pairs_growth_rate[,paste0("Pair_", growth_rate_variables[i])] <-
    pairs_growth_rate[,get(paste0(growth_rate_variables[i], 1))] - pairs_growth_rate[,get(paste0(growth_rate_variables[i], 2))]
}



# Tidy df
pairs_growth_rate <- pairs_growth_rate  %>%
  select(Community, Isolate1, Isolate2, InteractionType, From, To, starts_with("Pair")) %>%
  mutate(Community = ordered(Community, communities_name)) %>%
  gather("variable", "PairGrowthRateVariable", -(Community:To)) %>%
  separate(variable, c("temp", "CarbonSource", "GrowthRateVariable")) %>%
  select(-temp) %>%
  arrange(Community, Isolate1, Isolate2, CarbonSource) %>%
  as_tibble()

pairs_growth_rate <- mutate(pairs_growth_rate, PairGrowthRateVariable = abs(PairGrowthRateVariable))



# Plot  ----
"
1. make R function for plotting histogram
2. swith the sign in differnce in pairwise value
"

interaction_type <- c("win", "coexistence", "lose", "bistability", "neutrality", "self", "undefined")
myColor = c("#73C966", "#557BAA", "#DD6E74", "#EECF6D", "#8650C4", "black", "grey80")

## All in one figure
p_pairs_growth_rate_hist <-
  pairs_growth_rate %>%
  filter(!(GrowthRateVariable == "LagTime" & (PairGrowthRateVariable > 500 | PairGrowthRateVariable < -500))) %>%
  ggplot(aes(group = InteractionType)) +
  geom_histogram(aes(x = PairGrowthRateVariable, fill = InteractionType), col = "black") +
  scale_fill_manual(values = myColor) +
  facet_grid(CarbonSource~GrowthRateVariable, scale = "free") +
  theme_bw() +
  labs("")

p_pairs_growth_rate_density <-
  pairs_growth_rate %>%
  filter(!(GrowthRateVariable == "LagTime" & (PairGrowthRateVariable > 500 | PairGrowthRateVariable < -500))) %>%
  ggplot(aes(group = InteractionType)) +
  geom_density(aes(x = PairGrowthRateVariable, y = ..density.., col = InteractionType), lwd = 1) +
  scale_color_manual(values = myColor) +
  facet_grid(CarbonSource~GrowthRateVariable, scale = "free") +
  theme_bw() +
  labs("")



## Plot by carbon source
carbon_source <- unique(pairs_growth_rate$CarbonSource)
p_pairs_growth_rate_hist_list <- rep(list(NA), length(carbon_source))
p_pairs_growth_rate_density_list <- rep(list(NA), length(carbon_source))

for (i in 1:length(carbon_source)) {
  p_pairs_growth_rate_hist_list[[i]] <- pairs_growth_rate %>%
    filter(CarbonSource == carbon_source[i]) %>%
    filter(InteractionType %in% c("exclusion", "coexistence")) %>%
    ggplot() +
    geom_histogram(aes(x = PairGrowthRateVariable, fill = InteractionType), col = "black") +
    scale_fill_manual(values = myColor) +
    facet_grid(CarbonSource~GrowthRateVariable, scale = "free") +
    theme_bw() +
    ggtitle(carbon_source[i]) +
    labs("")
}

for (i in 1:length(carbon_source)) {
  p_pairs_growth_rate_density_list[[i]] <- pairs_growth_rate %>%
    filter(CarbonSource == carbon_source[i]) %>%
    filter(InteractionType %in% c("exclusion", "coexistence")) %>%
    ggplot(aes(group = InteractionType)) +
    geom_density(aes(x = PairGrowthRateVariable, y = ..density.., col = InteractionType), lwd = 1) +
    scale_color_manual(values = myColor) +
    facet_wrap(GrowthRateVariable~., scale = "free") +
    theme_bw() +
    ggtitle(carbon_source[i]) +
    labs("")
}










