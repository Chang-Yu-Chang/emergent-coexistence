#' Determine the isolate acids preference
library(tidyverse)
library(cowplot)

# Data ----
isolates <- read_csv(here::here("data/output/isolates.csv"))
# Byproduct measurement on glucose. Data from Sylvie
isolates_byproduct_time <- read_csv(here::here("data/raw/growth_rate/Estrela_2021_isolates_ph_OAs.csv")) %>%
    select(ID = SangerID, Time = time_hours, OD620, pH, Glucose_perc, acetate_mM, succinate_mM, lactate_mM)


# Match the most secreted CS of the dominant strain to the most preferred CS of the subdominant
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
counter
x

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

# none
isolates_preference$PreferredCS[isolates_preference$ID == 274] <- NA
isolates_preference$PreferredCS[isolates_preference$ID == 290] <- NA
isolates_preference$PreferredCS[isolates_preference$ID == 449] <- NA
isolates_preference$PreferredCS[isolates_preference$ID == 453] <- NA


isolates_preference$PreferredCS[isolates_preference$ID == 316] <- "acetate"
isolates_preference$PreferredCS[isolates_preference$ID == 338] <- "lactate"
isolates_preference$PreferredCS[isolates_preference$ID == 440] <- "acetate"
isolates_preference$PreferredCS[isolates_preference$ID == 523] <- "acetate"

write_csv(isolates_preference, file = here::here("data/temp/isolates_preference.csv"))


## Acids concentration over time
coeff <- 5
p <- isolates_byproduct_time %>%
    #filter(ID %in% c(160, 165)) %>%
    mutate(ID = factor(ID)) %>%
    # group_by(ID) %>%
    # select(-pH, -OD620, -Glucose_perc) %>%
    # pivot_longer(cols = c(-ID, -Time), names_to = "CS", values_to = "Conc") %>%
    # group_by(ID, CS) %>%
    # Remove a CS if it is not produced in any of the time point sampled
    filter(!all(acetate_mM == 0), !all(lactate_mM == 0), !all(succinate_mM == 0)) %>%
    ggplot() +
    geom_line(aes(x = Time, y = acetate_mM, color = "acetate")) +
    geom_point(aes(x = Time, y = acetate_mM, color = "acetate")) +
    geom_line(aes(x = Time, y = lactate_mM * coeff, color = "lactate")) +
    geom_point(aes(x = Time, y = lactate_mM * coeff, color = "lactate")) +
    geom_line(aes(x = Time, y = succinate_mM * coeff, color = "succinate")) +
    geom_point(aes(x = Time, y = succinate_mM * coeff, color = "succinate")) +
    geom_text(data = isolates_preference, aes(label = PreferredCS), x = 8, y = Inf, size = 2, vjust = 1) +
    scale_y_continuous(name = "acetate (mM)", sec.axis = sec_axis(~./coeff, name = "lactate (mM) or succinate (mM)")) +
    scale_color_manual(values = c("acetate" = "red", "lactate" = "green", "succinate" = "blue")) +
    #geom_point(aes(x = Time, y = Conc, color = CS, group = CS)) +
    facet_wrap(.~ID, scales = "free_y") +
    theme_bw() +
    theme(legend.position = "top") +
    #guides(color = "none") +
    labs(x = "Time", color = "")

ggsave(here::here("plots/Fig_pairs-byproduct_conc.png"), plot = p, width = 20 , height = 10)


# Acetate and pH over time
coeff <- 3
p <- isolates_byproduct_time %>%
    mutate(ID = factor(ID)) %>%
    ggplot() +
    geom_line(aes(x = Time, y = pH*coeff, color = "pH")) +
    geom_point(aes(x = Time, y = pH*coeff, color = "pH")) +
    geom_line(aes(x = Time, y = acetate_mM, color = "acetate")) +
    geom_point(aes(x = Time, y = acetate_mM, color = "acetate")) +
    #geom_text(data = isolates_preference, aes(label = PreferredCS), x = 8, y = Inf, size = 2, vjust = 1) +
    scale_y_continuous(name = "acetate (mM)", sec.axis = sec_axis(~./coeff, name = "pH")) +
    scale_color_manual(values = c("acetate" = "red", "pH" = "grey")) +
    facet_wrap(.~ID, scales = "free_y") +
    theme_bw() +
    theme(legend.position = "top") +
    labs(x = "Time", color = "")

ggsave(here::here("plots/Fig_pairs-pH.png"), plot = p, width = 20 , height = 10)




if (FALSE) {
    which(isolates_byproduct_time_w$ID %in% c(316, 338, 440, 523))
    i = 50
    isolates_byproduct_time_w$ID[6]
    t_ace_peak <- isolates_byproduct_time_w$acetate_peak[i]
    t_lac_peak <- isolates_byproduct_time_w$lactate_peak[i]
    t_suc_peak <- isolates_byproduct_time_w$succinate_peak[i]
    test_vec <- sort(c(acetate = t_ace_peak, lactate = t_lac_peak, succinate = t_suc_peak), decreasing = F)
    isolates_byproduct_time_w[i,] %>%
        pivot_longer(cols = contains("mM"), names_to = c("CS", "Time"), names_pattern = "(.*)_mM_(.*)")

}




