#' Perform rarefaction on T0 samples (the 11 inoclumn at T0)
#'

library(tidyverse)
library(cowplot)
library(broom)
source(here::here("analysis/00-metadata.R"))

communities <- read_csv(paste0(folder_data, "temp/00c-communities.csv"), show_col_types = F) %>%
    mutate(Community = factor(Community, Community))
communities_abundance <- read_csv(paste0(folder_data, "raw/community_ESV/Emergent_Comunity_Data.csv"), show_col_types = F) %>%
    filter(Carbon_Source == "Glucose" | Carbon_Source == "Original") %>%
    mutate(Community = factor(paste0("C", Inoculum, "R", Replicate), paste0("C", rep(1:12, each = 8), "R", rep(1:8, 12))))%>%
    arrange(Community, Family, Transfer, ESV)

communities_abundance_T0 <- communities_abundance %>%
    filter(Transfer == 0) %>%
    filter(Relative_Abundance > 0.0001)

communities_abundance_T0 %>%
    group_by(Inoculum) %>%
    count() %>%
    pull(n) %>%
    range

communities_abundance_T0_wider <- communities_abundance_T0 %>%
    mutate(AboundanceCount = round(Relative_Abundance * 10^4)) %>%
    select(ESV_ID, Inoculum, AboundanceCount) %>%
    pivot_wider(names_from = ESV_ID, values_from = AboundanceCount, values_fill = 0) %>%
    arrange(Inoculum)

# inoculum <- 1
# sampling_size <- 2000
# n_samples <- 100
rarefy_samples <- function (n_samples, sampling_size, inoculum) {
    rarefied_richness <- rep(NA, n_samples)
    for (i in 1:n_samples) {
        rarefied_richness[i] <-
            sample(x = names(communities_abundance_T0_wider[-1]),
                   size = sampling_size,
                   prob = communities_abundance_T0_wider[inoculum,-1],
                   replace = T) %>%
            unique() %>%
            length()
    }
    return(rarefied_richness)
}
communities_rarefaction <- tibble(SamplingSize = c(seq(10, 100, 10), seq(100, 1000, 100), seq(1000, 10000, 1000))) %>%
    slice(rep(1:n(), 11)) %>%
    mutate(Inoculum = rep(1:11, each = n()/11)) %>%
    mutate(RarefiedRichness = NA)

for (i in 1:nrow(communities_rarefaction)) {
    communities_rarefaction$RarefiedRichness[i] <- mean(
        rarefy_samples(n_samples = 100, sampling_size = communities_rarefaction$SamplingSize[i], inoculum = communities_rarefaction$Inoculum[i])
    )
    cat(" ", i)
}

communities_rarefaction %>% write_csv(paste0(folder_data, "temp/18-communities_rarefaction.csv"))
communities_abundance_T0 %>% write_csv(paste0(folder_data, "temp/18-communities_abundance_T0.csv"))










