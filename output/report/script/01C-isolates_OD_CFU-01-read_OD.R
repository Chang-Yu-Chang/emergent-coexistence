# Read OD data
library(tidyverse)

# Read raw files and set variable types
OD_B2 <- read_csv("~/Dropbox/lab/emergent-coexistence/data/raw/OD/OD_B2.csv", col_types = cols()) %>%
    mutate(Batch = "B2", Layout = as.character(Layout), Isolate1Freq = as.character(Isolate1Freq), Isolate2Freq = as.character(Isolate2Freq))
OD_C2 <- read_csv("~/Dropbox/lab/emergent-coexistence/data/raw/OD/OD_C2.csv", col_types = cols()) %>% # Data from batch C C11R1 plate is also included
    mutate(Batch = "C2", Isolate1 = as.character(Isolate1), Isolate2 = as.character(Isolate2))
OD_D <- read_csv("~/Dropbox/lab/emergent-coexistence/data/raw/OD/OD_D.csv", col_types = cols()) %>%
    mutate(Batch = "D", Isolate1 = as.character(Isolate1), Isolate2 = as.character(Isolate2))
communities <- read_csv("~/Dropbox/lab/emergent-coexistence/data/output/communities.csv", col_types = cols())

# Remove the contaminated data: B2 community C11R1 isolate1 and C2 C11R2 isolate 13 (streplococcus contamination)
OD_B2 <- OD_B2 %>%
    filter(!(Community == "C11R1" & Isolate1 == 1)) %>%
    filter(!(Community == "C11R1" & Isolate2 == 1))

OD_C2 <- OD_C2 %>%
    filter(!(Community == "C11R2" & Isolate1 == 13)) %>%
    filter(!(Community == "C11R2" & Isolate2 == 13))


# Join OD data from different batches
OD <- OD_B2 %>%
    bind_rows(OD_C2) %>%
    bind_rows(OD_D) %>%
    # Select essential variables for analysis. Remove rows of experimental blanks
    select(Batch, Community, Transfer, Isolate1, Isolate2, Isolate1Freq, Isolate2Freq, MixIsolate, MixPlate, Wavelength, Abs) %>%
    filter(Isolate1 != "blank") %>%
    mutate(Community = factor(Community, communities$Community)) %>%
    arrange(Community)

# Some pairs dont have full wavelength measurement
table(OD$Wavelength)
table(OD$Batch, OD$Community)

# There are 68 isolates that I used in the pairwise competition experiment.
OD %>% distinct(Community, Isolate1)


# Subset single culture that are not mixed
isolates_OD <- OD %>%
    filter(Isolate1 == Isolate2, MixIsolate == F) %>%
    mutate(Isolate = as.numeric(Isolate1)) %>% select(-Isolate1) %>%
    arrange(Community, Transfer, Isolate) %>%
    group_by(Community, Isolate, Transfer, Wavelength) %>%
    summarize(AbsMean = mean(Abs), AbsSd = sd(Abs)) %>% ungroup() %>%
    mutate(Isolate = ordered(Isolate, 1:13))

# Replace C11R2 isolates 2 T8 OD by T7 OD because of contamination at T8
temp <- isolates_OD %>%
    filter(Community == "C11R2", Isolate == 2, Transfer == 7) %>%
    mutate(Transfer = 8)
isolates_OD <- isolates_OD %>%
    filter(!(Community == "C11R2" & Isolate == 2 & Transfer == 8)) %>%
    bind_rows(temp) %>%
    arrange(Community, Isolate, Transfer, Wavelength)

# Set negative OD values to 0
isolates_OD$AbsMean[isolates_OD$AbsMean <=0] <- 0

#
isolates_OD %>% filter(Wavelength == 620, Transfer == 8)
write_csv(OD, "~/Dropbox/lab/emergent-coexistence/data/temp/OD.csv")
write_csv(isolates_OD, "~/Dropbox/lab/emergent-coexistence/data/temp/isolates_OD.csv")




