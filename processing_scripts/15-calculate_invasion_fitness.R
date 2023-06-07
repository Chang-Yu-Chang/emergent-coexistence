#' This script calculates the invasion communities_fitness for ESVs in the 26 glucose communities that have temporal data
#'
#' 1. calculate communities_fitness for all ESVs
#' 2. find the stable ESVs
#' 3. find the transient ESVs with >=5 data points
#' 4. find the traients ESVs with >=3 data points
#' 4. bind the fitness files. Output:
#'  temp/15-eq_freq_stable.csv. 99 stable ESVs
#'  temp/15-eq_freq_transient.csv. 46 transient ESVs with >=5 data points
#'  temp/15-eq_freq_transient2.csv. 110 transient ESVs with >=3 data points

library(tidyverse)
library(cowplot)
library(broom)
source(here::here("processing_scripts/00-metadata.R"))

communities_abundance <- read_csv(paste0(folder_data, "temp/14-communities_abundance.csv"), show_col_types = F) %>%
    mutate(Community = factor(paste0("C", Inoculum, "R", Replicate), paste0("C", rep(1:12, each = 8), "R", rep(1:8, 12))))

# 1. Calculate invasion communities_fitness ----
communities_abundance_temporal <- communities_abundance %>%
    filter(Transfer != 0) %>%
    # Filter for those that has temporal data; either from inocula 2 and 6, or replicate 4
    filter(Inoculum %in% c(2,6) | Replicate == 4) %>%
    select(Community, Transfer, ESV_ID, Relative_Abundance) %>%
    arrange(Community, Transfer, ESV_ID)

# A complete list of ESV per transfer
communities_abundance_temporal_complete <- communities_abundance_temporal %>%
    distinct(Community, ESV_ID) %>%
    slice(rep(1:n(), each = 12)) %>%
    mutate(Transfer = rep(1:12, n()/12))
nrow(distinct(communities_abundance_temporal_complete, Community, ESV_ID)) # 755 unique ESVs
nrow(communities_abundance_temporal_complete) # 755 ESVs * 12 transfers = 9060 rows

# Calculate communities_fitness as ln(x_i/x_{i-1})
communities_fitness <- communities_abundance_temporal %>%
    right_join(communities_abundance_temporal_complete) %>%
    group_by(Community, ESV_ID) %>%
    arrange(Community, ESV_ID, Transfer) %>%
    mutate(Fitness = log(lead(Relative_Abundance) / Relative_Abundance))


# 2. Find stable ESVs ----
# ESVs that have data from the last four transfers T9-12. Three exceptions: C4R4 and C9R4 do not have T11, and C5R4 does not have T9
ESV_stable <- communities_abundance_temporal %>%
    filter(Transfer %in% c(9:12)) %>%
    pivot_wider(id_cols = c(Community, ESV_ID), names_from = Transfer, names_prefix = "T", values_from = Relative_Abundance) %>%
    filter(
        # Species with complete presence at T9-12
        (!is.na(T9) & !is.na(T10) & !is.na(T11) & !is.na(T12) ) |
            # Community C4R4, C9R4 have missing time points at T11. For these communities, add back the species
            (Community %in% c("C4R4", "C9R4") & (!is.na(T9) & !is.na(T10) & !is.na(T12))) |
            # Community C5R4 have missing time points at T9. For these communities, add back the species
            (Community %in% c("C5R4") & (!is.na(T10) & !is.na(T11) & !is.na(T12)))
    ) %>%
    distinct(Community, ESV_ID) %>%
    mutate(CommunityESV = paste0(Community, ESV_ID))
nrow(ESV_stable) # 99 ESVs are considered "stable" ESVs

fitness_stable <- communities_fitness %>%
    drop_na() %>%
    mutate(CommunityESV = paste0(Community, ESV_ID)) %>%
    filter(CommunityESV %in% ESV_stable$CommunityESV)

# Correlation
cor_test_long <- function (long_data) cor.test(long_data$Relative_Abundance, long_data$Fitness, method = "spearman", exact = F, alternative = c("less"))
tb_cor_stable <- fitness_stable %>%
    nest(data = c(-Community, -ESV_ID)) %>%
    mutate(fit = map(data, ~ cor_test_long(.x)),
           tidied = map(fit, tidy)) %>%
    unnest(tidied)
ESV_sig_stable <- tb_cor_stable %>%
    select(Community, ESV_ID, estimate, p.value) %>%
    filter(p.value < 0.05, estimate < 0) %>%
    mutate(CommunityESV = paste0(Community, ESV_ID)) %>%
    mutate(Significance = "p<0.05")

nrow(tb_cor_stable) # 99 stable ESVs
nrow(ESV_sig_stable) # 52 shows significant negative correlation
range(ESV_sig_stable$estimate) # range of correlation [-0.95, -0.53]

# Empirical equilibrium frequency
communities_eq_freq_stable <- communities_abundance_temporal %>%
    mutate(CommunityESV = paste0(Community, ESV_ID)) %>%
    filter(Transfer %in% 9:12) %>%
    filter(CommunityESV %in% ESV_stable$CommunityESV) %>%
    arrange(Community, ESV_ID, CommunityESV, Transfer) %>%  # 387 rows. 99 * 4 - 9 = 387. 9 ESVs from C4R4, C5R4, C9R4 have only three time points
    group_by(Community, ESV_ID, CommunityESV) %>%
    summarize(EmpiricalEqAbundance = mean(Relative_Abundance)) # 99 ESVs

tb_cor_stable %>% filter(Community %in% c("C4R4", "C5R4", "C9R4")) %>% nrow() # 9 ESVs from C4R4, C5R4, C9R4 have only three time points

# Linear model predicted equilibrium frequency
xintercept_stable <- fitness_stable %>%
    nest(data = c(-Community, -ESV_ID, -CommunityESV)) %>%
    mutate(fit = map(data, ~ lm(Fitness ~ Relative_Abundance, data = .x)),
           tidied = map(fit, tidy)) %>%
    unnest(tidied) %>%
    select(Community, ESV_ID, CommunityESV, term, estimate) %>%
    pivot_wider(names_from = term, values_from = estimate) %>%
    # X intercept is the predicted equilibrium when fitness=0, meaning -b/a
    mutate(PredictedEqAbundance = -`(Intercept)` / Relative_Abundance) %>%
    select(Community, ESV_ID, CommunityESV, PredictedEqAbundance, Slope = Relative_Abundance)

table(xintercept_stable$Slope < 0) # 95 ESVs have negative slopes, 4 have positive slope

# 3. find transient ESVs with >=5 data points ----
ESV_transient <- communities_fitness %>%
    drop_na() %>%
    mutate(CommunityESV = paste0(Community, ESV_ID)) %>%
    # Not in the list of stable ESVs
    filter(!(CommunityESV %in% ESV_stable$CommunityESV)) %>%
    group_by(CommunityESV, Community, ESV_ID) %>%
    count() %>%
    # At least five data points
    filter(n>=5) %>%
    # Remove the artifact C10R4 Stenotrophomonas
    filter(CommunityESV != "C10R4Stenotrophomonas")
nrow(ESV_transient) # 46 ESVs that are considered "transient"

fitness_transient <- communities_fitness %>%
    drop_na() %>%
    mutate(CommunityESV = paste0(Community, ESV_ID)) %>%
    filter(CommunityESV %in% ESV_transient$CommunityESV)

# Correlation
tb_cor_transient <- fitness_transient %>%
    nest(data = c(-Community, -ESV_ID)) %>%
    mutate(fit = map(data, ~ cor_test_long(.x)),
           tidied = map(fit, tidy)) %>%
    unnest(tidied)

ESV_sig_transient <- tb_cor_transient %>%
    select(Community, ESV_ID, estimate, p.value) %>%
    filter(p.value < 0.05, estimate < 0) %>%
    mutate(CommunityESV = paste0(Community, ESV_ID)) %>%
    mutate(Significance = "p<0.05")

nrow(tb_cor_transient) # 46 transient ESVs
nrow(ESV_sig_transient) # 10 transient ESVs shows significant negative correlation
range(ESV_sig_transient$estimate) # range of correlation [-1, -0.714]

# Linear model predicted equilibrium frequency
xintercept_transient <- fitness_transient %>%
    nest(data = c(-Community, -ESV_ID, -CommunityESV)) %>%
    mutate(fit = map(data, ~ lm(Fitness ~ Relative_Abundance, data = .x)),
           tidied = map(fit, tidy)) %>%
    unnest(tidied) %>%
    select(Community, ESV_ID, CommunityESV, term, estimate) %>%
    pivot_wider(names_from = term, values_from = estimate) %>%
    # X intercept is the predicted equilibrium when fitness=0, meaning -b/a
    mutate(PredictedEqAbundance = -`(Intercept)` / Relative_Abundance) %>%
    select(Community, ESV_ID, CommunityESV, PredictedEqAbundance, Slope = Relative_Abundance)

# 4. find transient ESVs with >=3 data points ----
ESV_transient2 <- communities_fitness %>%
    drop_na() %>%
    mutate(CommunityESV = paste0(Community, ESV_ID)) %>%
    filter(!(CommunityESV %in% ESV_stable$CommunityESV)) %>%
    group_by(CommunityESV, Community, ESV_ID) %>%
    count() %>%
    # Include the ESVs with >=3 data points
    filter(n>=3) %>%
    # Remove the artifact C10R4 Stenotrophomonas
    filter(CommunityESV != "C10R4Stenotrophomonas")

fitness_transient2 <- communities_fitness %>%
    drop_na() %>%
    mutate(CommunityESV = paste0(Community, ESV_ID)) %>%
    filter(CommunityESV %in% ESV_transient2$CommunityESV)

# Correlation
tb_cor_transient2 <- fitness_transient2 %>%
    nest(data = c(-Community, -ESV_ID)) %>%
    mutate(fit = map(data, ~ cor_test_long(.x)),
           tidied = map(fit, tidy)) %>%
    unnest(tidied)

ESV_sig_transient2 <- tb_cor_transient2 %>%
    select(Community, ESV_ID, estimate, p.value) %>%
    filter(p.value < 0.05, estimate < 0) %>%
    mutate(CommunityESV = paste0(Community, ESV_ID)) %>%
    mutate(Significance = "p<0.05")

nrow(tb_cor_transient2) # 110 transient ESVs
nrow(ESV_sig_transient2) # 21 transient ESVs shows significant negative correlation
range(ESV_sig_transient2$estimate) # range of correlation [-1, -0.714]

# Linear model predicted equilibrium frequency
xintercept_transient2 <- fitness_transient2 %>%
    nest(data = c(-Community, -ESV_ID, -CommunityESV)) %>%
    mutate(fit = map(data, ~ lm(Fitness ~ Relative_Abundance, data = .x)),
           tidied = map(fit, tidy)) %>%
    unnest(tidied) %>%
    select(Community, ESV_ID, CommunityESV, term, estimate) %>%
    pivot_wider(names_from = term, values_from = estimate) %>%
    # X intercept is the predicted equilibrium when fitness=0, meaning -b/a
    mutate(PredictedEqAbundance = -`(Intercept)` / Relative_Abundance) %>%
    select(Community, ESV_ID, CommunityESV, PredictedEqAbundance, Slope = Relative_Abundance) %>%
    mutate(ESVType = "transient")

# 5. Save fitness data ----
eq_freq_stable <- communities_eq_freq_stable %>%
    left_join(select(ESV_sig_stable, CommunityESV, Significance, rho = estimate)) %>%
    left_join(xintercept_stable) %>%
    mutate(ESVType = "stable")
fitness_stable <- fitness_stable %>% mutate(ESVType = "stable")
nrow(eq_freq_stable) # 99 ESVs
write_csv(eq_freq_stable, paste0(folder_data, 'temp/15-eq_freq_stable.csv'))
write_csv(fitness_stable, paste0(folder_data, 'temp/15-fitness_stable.csv'))

eq_freq_transient <- xintercept_transient %>%
    left_join(select(ESV_sig_transient, CommunityESV, Significance, rho = estimate)) %>%
    mutate(ESVType = "transient")
fitness_transient <- fitness_transient %>% mutate(ESVType = "transient")
nrow(eq_freq_transient) # 46 ESV
write_csv(eq_freq_transient, paste0(folder_data, 'temp/15-eq_freq_transient.csv'))
write_csv(fitness_transient, paste0(folder_data, 'temp/15-fitness_transient.csv'))


# For Fig S5
eq_freq_transient2 <- xintercept_transient2 %>%
    left_join(select(ESV_sig_transient, CommunityESV, Significance, rho = estimate)) %>%
    mutate(ESVType = "transient")
fitness_transient2 <- fitness_transient2 %>% mutate(ESVType = "transient")
nrow(eq_freq_transient2) # 110 ESVs

write_csv(eq_freq_transient2, paste0(folder_data, 'temp/15-eq_freq_transient2.csv'))
write_csv(fitness_transient2, paste0(folder_data, 'temp/15-fitness_transient2.csv'))
























