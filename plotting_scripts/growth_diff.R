# Plot the difference in growth rates on glucose and acetate
library(tidyverse)
library(cowplot)

isolates_growth <- read_csv(here::here("data/raw/growth_rate/Growthcurver.csv"))

# Dissimilarity using the growth rates on different acids
isolates_growth_w <- isolates_growth %>%
    separate(col = SID, sep = "_", into  = c("ID", "CS"), convert = T) %>%
    select(-SangerID, -Family) %>%
    select(ID, CS, RMid) %>%
    filter(CS %in% c("D-Glucose", "Acetate", "D-Lactate", "Succinate", "Gluconate", "2-Ketogluconate")) %>%
    pivot_wider(names_from = CS, values_from = RMid)
ID <- isolates_growth_w$ID

## Sign determined by the overall sum of the growth rates on acids
isolates_growth_total <- isolates_growth_w %>%
    pivot_longer(-ID) %>%
    group_by(ID) %>%
    summarize(r_total = sum(value))

pairs_temp <- combn(ID, 2) %>% t() %>% as_tibble() %>% setNames(c("ID1", "ID2"))
pairs_r_total <- pairs_temp %>%
    bind_rows(setNames(pairs_temp, c("ID2", "ID1"))) %>%
    left_join(rename_with(isolates_growth_total, ~ paste0(., "1"))) %>%
    left_join(rename_with(isolates_growth_total, ~ paste0(., "2"))) %>%
    mutate(ID1_has_higher_r_total = r_total1 > r_total2, r_total_d = r_total1 - r_total2) %>%
    {.}

## Distance between the two IDs on CSs
isolates_growth_dis <- isolates_growth_w %>%
    # Exclude glucose
    select(-ID, -`D-Glucose`) %>%
    cluster::daisy(metric="euclidean") %>%
    as.matrix() %>% as_tibble() %>%
    mutate(ID1 = ID) %>%
    pivot_longer(cols = -ID1, names_to = "ID2", names_transform = list(ID2 = as.integer), values_to = "r_dissim") %>%
    mutate(ID2 = ID[ID2]) %>%
    # Add signs
    left_join(pairs_r_total) %>%
    filter(!is.na(r_total1)) %>%
    mutate(r_dissim = ifelse(ID1_has_higher_r_total, -r_dissim, r_dissim))

## Scatterplot of r_glu_d against +-|r_dis|, for all isoaltes
isolates_growth_dis %>%
    left_join(rename_with(select(isolates_growth_w, ID, `D-Glucose`), ~ paste0(., "1"))) %>%
    left_join(rename_with(select(isolates_growth_w, ID, `D-Glucose`), ~ paste0(., "2"))) %>%
    mutate(r_glu_d = `D-Glucose1` - `D-Glucose2`) %>%
    ggplot() +
    geom_point(aes(x = r_glu_d, y = r_dissim), size = 2, shape = 1) +
    theme_classic() +
    theme(legend.title = element_blank(), legend.position = "top", legend.direction = "horizontal") +
    guides(shape = "none") +
    labs(x = expression(r[A_glu] - r[B_glu]), y = expression(r[dis_acids]))
