#' Read isolate monoculture growth rates in various carbon sources
library(tidyverse)
library(growthcurver)
library(janitor)
isolates_ID_match <- read_csv(here::here("data/temp/isolates_ID_match.csv"))

# Byproduct measurement on glucose. Data from Sylvie
isolates_byproduct <- read_csv("~/Dropbox/lab/invasion-network/data/raw/growth_rate/By_Products_Glucose.csv") %>%
    select(OD620_16h = OD620, ID = SangerID, Glucose_perc, acetate_mM, succinate_mM, lactate_mM, gluconate_mM, ketogluconate_mM)
isolates_byproduct_time <- read_csv("~/Dropbox/lab/invasion-network/data/raw/growth_rate/Estrela_2021_isolates_ph_OAs.csv") %>%
    select(ID = SangerID, Time = time_hours, OD620, pH, Glucose_perc, acetate_mM, succinate_mM, lactate_mM)

# Growth rate. Growthcurver----
# Growth curve data from Sylvie
isolates_curves1 <- read_csv("~/Dropbox/lab/invasion-network/data/raw/growth_rate/raw_gcurves_all_sylvies.csv") %>%
    select(Date = date, ID = seq, Well = well, CS = csource, Time = t, OD620 = abs) %>%
    filter(ID %in% isolates_ID_match$ID) %>%
    mutate(CS = str_replace(CS, "Dlactate", "lactate")) %>%
    # Use the latest date measured replicates and use one replicate
    group_by(ID, CS, Time) %>%
    filter(Date == max(Date), Well == unique(Well)[1])
isolates_curves_list1 <- isolates_curves1 %>%
    nest(GrowthCurve = c(Time, OD620)) %>%
    rowwise() %>%
    # growthcurver fit
    mutate(GrowthCurveOutcome = SummarizeGrowth(GrowthCurve$Time, GrowthCurve$OD620) %>% list())
isolates_growthcurver1 <- isolates_curves_list1 %>%
    rowwise() %>%
    mutate(r = GrowthCurveOutcome$vals$r, sigma = GrowthCurveOutcome$vals$sigma) %>%
    group_by(ID, CS) %>%
    slice(1) %>%
    select(ID, CS, r) %>%
    pivot_wider(names_from = CS, values_from = r, names_glue = "r_{CS}_curver")

# Growth curve from Jean
isolates_curves2 <- read_csv("~/Dropbox/lab/invasion-network/data/raw/growth_rate/20GC_Data.csv") %>%
    select(ID = SangerID, CS, Time, OD620) %>%
    filter(ID %in% isolates_ID_match$ID) %>%
    mutate(CS = tolower(CS)) %>%
    mutate(CS = str_replace(CS, "d-", "") %>% str_replace("l-", "") %>% str_replace("2-", "")) %>%
    filter(CS %in% c("glucose", "lactate", "acetate", "succinate"))
isolates_curves_list2 <- isolates_curves2 %>%
    nest(GrowthCurve = c(Time, OD620)) %>%
    rowwise() %>%
    mutate(GrowthCurveOutcome = SummarizeGrowth(GrowthCurve$Time, GrowthCurve$OD620) %>% list()) # growthcurver fit
isolates_growthcurver2 <- isolates_curves_list2 %>%
    rowwise() %>%
    mutate(r = GrowthCurveOutcome$vals$r, sigma = GrowthCurveOutcome$vals$sigma) %>%
    select(ID, CS, r) %>%
    pivot_wider(names_from = CS, values_from = r, names_glue = "r_{CS}_curver")



# Growth rates using the end time point 12, 16, 28 hr-----
calculate_r <- function(N0, N1, T0, T1) (log10(N1) - log10(N0)) / (T1 - T0)
isolates_curves_T0 <- isolates_curves1 %>%
    group_by(ID, Well, CS) %>%
    filter(Time == min(Time)) %>%
    mutate(T0 = ifelse(Time == min(Time), 0, Time), N0 = OD620) %>%
    select(-Time, -OD620)
isolates_growth <- isolates_curves1 %>%
    filter(Time %in% c(12, 16, 28)) %>%
    group_by(ID, CS) %>%
    arrange(ID, CS) %>%
    mutate(T1 = Time, N1 = OD620) %>%
    left_join(isolates_curves_T0) %>%
    # Calculate r
    mutate(r = calculate_r(N0, N1, T0, T1)) %>%
    select(Date, ID, Well, CS, Time = T1, r) %>%
    # Remove contamination
    filter(r>0) %>%
    # Average
    group_by(ID, CS, Time) %>%
    summarize(r = mean(r)) %>%
    pivot_wider(names_from = c(Time, CS), values_from = r, names_glue = "r_{CS}_{Time}hr") %>%
    ungroup()

# Jean's growth rate data. Use the fitted Rmid ----
# Growth rate data from Jean
isolates_growth_mid <- read_csv("~/Dropbox/lab/invasion-network/data/raw/growth_rate/Growthcurver.csv")
isolates_growth_w_mid <- isolates_growth_mid %>%
    separate(col = SID, sep = "_", into  = c("ID", "CS"), convert = T) %>%
    select(-SangerID, -Family) %>%
    select(ID, CS, RMid) %>%
    filter(CS %in% c("D-Glucose", "Acetate", "D-Lactate", "Succinate", "Gluconate", "2-Ketogluconate")) %>%
    # names to lower case
    mutate(CS = CS %>% sub("2-", "", .) %>% sub("D-", "", .) %>% tolower()) %>%
    pivot_wider(names_from = CS, values_from = RMid, names_glue = "r_{CS}_midhr")


# Sylvie's growth rate. Use Rmax ----
# Growth rate data from Sylvie
isolates_growth_max <- read_csv("~/Dropbox/lab/invasion-network/data/raw/growth_rate/Estrela_2021_isolates_grmax.csv")
isolates_growth_w_max <- isolates_growth_max %>%
    select(ID = SangerID, CS = cs, gr_max) %>%
    pivot_wider(names_from = CS, values_from = gr_max, names_glue = "r_{CS}_maxhr")



# Isolate preference----
## Match the most secreted CS of the dominant strain to the most preferred CS of the subdominant
find_max <- function (x) {
    if (all(x==0)) return(NaN)
    return(x == max(x))
}
isolates_byproduct_time_peak <- isolates_byproduct_time %>%
    #mutate(Time = paste0(Time, "hr")) %>%
    group_by(ID) %>%
    summarize(acetate_peak = Time[find_max(acetate_mM)],
              lactate_peak = Time[find_max(lactate_mM)],
              succinate_peak = Time[find_max(succinate_mM)]) %>%
    {.}
isolates_byproduct_time_w <- isolates_byproduct_time %>%
    select(ID, Time, ends_with("mM")) %>%
    pivot_wider(names_from = Time, values_from = ends_with("mM")) %>%
    left_join(isolates_byproduct_time_peak)

x = rep(NA, nrow(isolates_byproduct_time_w))
counter = 0
for (i in 1:nrow(isolates_byproduct_time_w)) {
    t_ace_peak <- isolates_byproduct_time_w$acetate_peak[i]
    t_lac_peak <- isolates_byproduct_time_w$lactate_peak[i]
    t_suc_peak <- isolates_byproduct_time_w$succinate_peak[i]
    test_vec <- sort(c(acetate = t_ace_peak, lactate = t_lac_peak, succinate = t_suc_peak), decreasing = F)
    # Skip NA
    if (any(is.na(c(t_ace_peak, t_lac_peak, t_suc_peak)))) next
    # Skip those that does not consume the acids
    if (all(c(t_ace_peak, t_lac_peak, t_suc_peak) == 48)) next
    # Sequential preference
    if (length(unique(test_vec)) == 3) {
        print(test_vec)
        x[i] <- names(test_vec)[1]
        counter = counter + 1
    }
    # Co-utilization on less preferred ones
    if (length(unique(test_vec)) == 2 & duplicated(test_vec)[3]) {
        print(test_vec)
        x[i] <- names(test_vec)[1]
        counter = counter + 1
    }
    # Co-utilization on more preferred ones
    if (length(unique(test_vec)) == 2 & duplicated(test_vec)[2]) {
        print(test_vec)
        temp <- isolates_byproduct_time_w
        t1 <- test_vec[1]
        t2 <- ifelse(t1 == 16, 28, ifelse(t1 == 28, 48, NA))
        # Compare the fold changes
        x_1_t1 <- temp[i, paste0(names(test_vec[1]), "_mM_", t1)]
        x_1_t2 <- temp[i, paste0(names(test_vec[1]), "_mM_", t2)]
        x_2_t1 <- temp[i, paste0(names(test_vec[2]), "_mM_", t1)]
        x_2_t2 <- temp[i, paste0(names(test_vec[2]), "_mM_", t2)]
        if (x_1_t1/x_1_t2 > x_2_t1/x_2_t2) x[i] <- names(test_vec)[1]
        if (x_1_t1/x_1_t2 < x_2_t1/x_2_t2) x[i] <- names(test_vec)[2]
        if (x_1_t2 == 0 | x_2_t2 == 0) x[i] <- names(test_vec)[ifelse(x_1_t1 - x_1_t2 > x_2_t1 - x_2_t2, 1, 2)]
        counter = counter + 1
    }

    # Co-utilization on three CS
    if (length(unique(test_vec)) == 1) {
        print(test_vec)
        temp <- isolates_byproduct_time_w
        t1 <- test_vec[1]
        t2 <- ifelse(t1 == 16, 28, ifelse(t1 == 28, 48, NA))
        x_1_t1 <- temp[i, paste0(names(test_vec[1]), "_mM_", t1)]
        x_1_t2 <- temp[i, paste0(names(test_vec[1]), "_mM_", t2)]
        x_2_t1 <- temp[i, paste0(names(test_vec[2]), "_mM_", t1)]
        x_2_t2 <- temp[i, paste0(names(test_vec[2]), "_mM_", t2)]
        x_3_t1 <- temp[i, paste0(names(test_vec[3]), "_mM_", t1)]
        x_3_t2 <- temp[i, paste0(names(test_vec[3]), "_mM_", t2)]

        if (x_1_t1/x_1_t2 > x_2_t1/x_2_t2 & x_1_t1/x_1_t2 > x_3_t1/x_3_t2) names(test_vec)[1]
        if (x_2_t1/x_2_t2 > x_3_t1/x_3_t2 & x_2_t1/x_2_t2 > x_3_t1/x_3_t2) names(test_vec)[2]
        if (x_3_t1/x_3_t2 > x_2_t1/x_2_t2 & x_3_t1/x_3_t2 > x_1_t1/x_1_t2) names(test_vec)[3]

        counter = counter + 1

    }
}

isolates_preferred <- tibble(ID = isolates_byproduct_time_w$ID, PreferredCS = x)
isolates_secreted <- isolates_byproduct_time_w %>%
    select(ID, ends_with("16")) %>%
    pivot_longer(cols = -ID, names_to = c("CS", "Time"), names_pattern = "(.*)_mM_(.*)") %>%
    group_by(ID) %>%
    filter(value == max(value)) %>%
    ungroup() %>%
    select(ID, SecretedCS = CS)
isolates_preference <- isolates_preferred %>% left_join(isolates_secreted)



## Manual assignment for crappy isolates
### From succinate
isolates_preference$PreferredCS[isolates_preference$ID == 160] <- "lactate"
isolates_preference$PreferredCS[isolates_preference$ID == 169] <- "lactate"
isolates_preference$PreferredCS[isolates_preference$ID == 259] <- "lactate"
isolates_preference$PreferredCS[isolates_preference$ID == 270] <- "lactate"
isolates_preference$PreferredCS[isolates_preference$ID == 276] <- "lactate"
isolates_preference$PreferredCS[isolates_preference$ID == 278] <- "lactate"
isolates_preference$PreferredCS[isolates_preference$ID == 288] <- "lactate"
isolates_preference$PreferredCS[isolates_preference$ID == 289] <- "lactate"
isolates_preference$PreferredCS[isolates_preference$ID == 326] <- "lactate"
isolates_preference$PreferredCS[isolates_preference$ID == 332] <- "lactate"
isolates_preference$PreferredCS[isolates_preference$ID == 334] <- "lactate" # check
isolates_preference$PreferredCS[isolates_preference$ID == 348] <- "lactate"
isolates_preference$PreferredCS[isolates_preference$ID == 364] <- "lactate"
isolates_preference$PreferredCS[isolates_preference$ID == 439] <- "lactate"
isolates_preference$PreferredCS[isolates_preference$ID == 460] <- "lactate"

### From lactate
isolates_preference$PreferredCS[isolates_preference$ID == 338] <- "succinate"

### To both
isolates_preference$PreferredCS[isolates_preference$ID == 280] <- "lactate and succinate"

### From non
isolates_preference$PreferredCS[isolates_preference$ID == 304] <- "acetate"
isolates_preference$PreferredCS[isolates_preference$ID == 350] <- "acetate"

### From acetate
isolates_preference$PreferredCS[isolates_preference$ID == 352] <- "succinate" # check
isolates_preference$PreferredCS[isolates_preference$ID == 358] <- "succinate" # check
isolates_preference$PreferredCS[isolates_preference$ID == 432] <- "acetate" # check
isolates_preference$PreferredCS[isolates_preference$ID == 440] <- "lactate"
isolates_preference$PreferredCS[isolates_preference$ID == 462] <- "succinate"

### none
isolates_preference$PreferredCS[isolates_preference$ID == 274] <- NA
isolates_preference$PreferredCS[isolates_preference$ID == 290] <- NA
isolates_preference$PreferredCS[isolates_preference$ID == 449] <- NA
isolates_preference$PreferredCS[isolates_preference$ID == 453] <- NA

isolates_preference$PreferredCS[isolates_preference$ID == 316] <- "acetate"
isolates_preference$PreferredCS[isolates_preference$ID == 338] <- "lactate"
isolates_preference$PreferredCS[isolates_preference$ID == 440] <- "acetate"
isolates_preference$PreferredCS[isolates_preference$ID == 523] <- "acetate"

# Isolate pH ----
isolates_pH <- isolates_byproduct_time %>%
    select(ID, Time, pH) %>%
    pivot_wider(names_from = Time, values_from = pH, names_glue = "pH_{Time}hr")
# Total amount of acid secretion ----
isolates_byproduct_time_sum <- isolates_byproduct_time %>%
    group_by(ID, Time) %>%
    # Weighted by CS number of carbons
    mutate(ByproductSum = sum(2 * acetate_mM, 4 * succinate_mM, 3* lactate_mM)) %>%
    select(ID, Time, ByproductSum) %>%
    pivot_wider(names_from = Time, values_from = ByproductSum, names_glue = "X_sum_{Time}hr")
# Acid secretion
isolates_acids <- isolates_byproduct_time %>%
    select(ID, Time, ends_with("mM")) %>%
    rename_with(~ paste0("X_", sub("_mM", "", .)), ends_with("_mM")) %>%
    pivot_wider(names_from = Time, values_from = starts_with("X_"), names_glue = "{.value}_{Time}hr")
# OD
isolates_OD <- isolates_curves2 %>%
    mutate(TimeHr = case_when(
        abs((Time-12)) == min(abs((Time-12))) ~ 12,
        abs((Time-16)) == min(abs((Time-16)))  ~ 16,
        abs((Time-28)) == min(abs((Time-28)))  ~ 28,
        abs((Time-48)) == min(abs((Time-48)))  ~ 48
    )) %>%
    filter(!is.na(TimeHr)) %>%
    select(ID, CS, Time = TimeHr, OD = OD620) %>%
    pivot_wider(names_from = c(CS, Time), values_from = OD, names_glue = "OD_{CS}_{Time}hr")
# Leakiness
isolates_leakiness <- isolates_byproduct_time %>%
    select(ID, Time, gluConc = Glucose_perc, X_acetate = acetate_mM, X_succinate = succinate_mM, X_lactate = lactate_mM) %>%
    mutate(X_sum = 2 * X_acetate + 3 * X_lactate + 4 * X_succinate) %>%
    # Convert the glucose from mass fraction to mM per carbon
    mutate(gluConsumption = (0.2 - gluConc) / 180 / 0.1 * 1000 * 6) %>%
    mutate(leakiness = X_sum / gluConsumption) %>%
    filter(Time != 0) %>%
    select(ID, Time, leakiness) %>%
    pivot_wider(names_from = Time, values_from = leakiness, names_glue = "leakiness_{Time}hr") %>%
    # Replace the leakiness of 0 glucose consumption time point by NA
    mutate_all(function(x) ifelse(is.infinite(x), NA, x))


# Combine growth rate, OD, pH, preference, secretion data
isolates_growth_traits <- isolates_ID_match %>%
    left_join(isolates_growthcurver2) %>% # fitted growth curve r from Jean's 20GC_Data.csv
    left_join(isolates_growth_w_mid) %>% # fitted r from Jean's Growthcurver.csv
    left_join(isolates_growth_w_max) %>% # r from Sylvie's Estrela_2021_isolates_grmax
    left_join(isolates_growth) %>%  # r from Sylvie's raw_gcurves_all_sylvies.csv
    left_join(isolates_acids) %>% # from Sylvie's Estrela_2021_isolates_ph_OAs.csv
    left_join(isolates_byproduct_time_sum) %>% # from Sylvie's Estrela_2021_isolates_ph_OAs.csv
    left_join(isolates_OD) %>% # from Jean's 20GC_Data.csv
    left_join(isolates_pH) %>% # from Sylvie's Estrela_2021_isolates_ph_OAs.csv
    left_join(isolates_preference) %>% # from Sylvie's Estrela_2021_isolates_ph_OAs.csv
    left_join(isolates_leakiness) # from Sylvie's Estrela_2021_isolates_ph_OAs.csv



write_csv(isolates_growth_traits, file = here::here("data/temp/isolates_growth_traits.csv"))
write_csv(isolates_curves1, file = here::here("data/output/isolates_curves1.csv"))
write_csv(isolates_curves2, file = here::here("data/output/isolates_curves2.csv"))
