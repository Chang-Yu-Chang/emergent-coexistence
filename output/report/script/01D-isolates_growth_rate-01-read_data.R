#' Read isolate monoculture growth rates in various carbon sources
library(tidyverse)

isolates_ID_match <- read_csv(here::here("data/temp/isolates_ID_match.csv"))
by_product_glucose <- read_csv(here::here("data/raw/growth_rate/By_Products_Glucose.csv")) %>% mutate(ID = SangerID)
growth_curver <- read_csv(here::here("data/raw/growth_rate/Growthcurver.csv")) %>% mutate(ID = SangerID)

isolates_glu <- isolates_ID_match %>% left_join(by_product_glucose, by = "ID")
isolates_growth <- isolates_ID_match %>% left_join(growth_curver, by = "ID")

write_csv(isolates_glu, here::here("data/temp/isolates_glu.csv"))
write_csv(isolates_growth, here::here("data/output/isolates_growth.csv"))


if (FALSE) {
    # Tidy up the variable names ----
    isolates_growth_rate <- d.fit %>%
        select(strain, cs, date, max_gr, lag, max.od) %>%
        setNames(c("ExpID", "CarbonSource", "DateMeasurement", "MaxGrowthRate", "LagTime", "MaxOD")) %>%
        # Subset for my isolates
        filter(ExpID %in% isolates_ID_match$ExpID) %>%
        as_tibble()

    ## Remove the duplicate measurement that were done twice on different days
    isolates_growth_rate <- isolates_growth_rate %>%
        group_by(ExpID, CarbonSource) %>%
        # For isolates with more than two measurement, take smaller maximal growth rates and max OD to prevent using contaminated samples
        summarize(MaxGrowthRate = min(MaxGrowthRate), LagTime = LagTime[1], MaxOD = min(MaxOD))

    ## Join isolate growth rates and isolate ID
    isolates_growth_rate_ID <- left_join(isolates_ID_match, isolates_growth_rate, by = "ExpID") %>% as_tibble()

    # Spread the isolates df so that it only has 68 rows (the number of isolates) ----
    isolates_growth_rate_spread <- isolates_growth_rate %>%
        gather(variable, value, -(ExpID:CarbonSource)) %>%
        unite(temp, CarbonSource, variable) %>%
        spread(temp, value)

    ## Join spread df of growth rates and isolate ID
    isolates_growth_rate_spread_ID <- left_join(isolates_ID_match, isolates_growth_rate_spread, by = "ExpID") %>% as_tibble()

    #
}


