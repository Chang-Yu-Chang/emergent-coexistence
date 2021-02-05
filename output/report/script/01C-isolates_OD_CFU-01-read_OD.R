# Read OD data from temp files
library(tidyverse)
library(data.table)

# Read raw files ----
OD_B2 <- fread(here::here("data/raw/OD/OD_B2.csv"))
OD_C2 <- fread(here::here("data/raw/OD/OD_C2.csv")) # Data from batch C C11R1 plate is also included
OD_D <- fread(here::here("data/raw/OD/OD_D.csv"))

# Remove the contaminated data: B2 community C11R1 isolate1
OD_B2 <- filter(OD_B2, !(Community == "C11R1" & Isolate1 == 1)) %>% filter(!(Community == "C11R1" & Isolate2 == 1))

# Merge OD data from different batches.
OD <- rbind(cbind(Experiment = "transitivity_B2", OD_B2),
  cbind(Experiment = "transitivity_C2", OD_C2),
  cbind(Experiment = "transitivity_D", OD_D))

# Select essential variables for analysis. Filter out row of experimental blanks.
OD <- OD %>%
  select(Community, Transfer, Isolate1, Isolate2, Isolate1Freq, Isolate2Freq, MixIsolate, MixPlate, Wavelength, Abs, Experiment) %>%
  filter(Isolate1 != "blank") %>%
  mutate(Community = factor(Community, communities_name)) %>% arrange(Community)

# Some pairs dont have full wavelength measurement
# table(OD$Wavelength)
# table(OD$Experiment, OD$Community)

# unique(OD[,c("Community", "Isolate1")]) %>% nrow()
# There are 68 isolates that I used in the pairwise competition experiment.

# Isolate ----
# Subset single culture that are not mixed
isolates_OD <- filter(OD, Isolate1 == Isolate2, MixIsolate == F) %>%
  mutate(Isolate = as.numeric(Isolate1)) %>% select(-Isolate1) %>%
  arrange(Community, Transfer, Isolate) %>%
  group_by(Community, Isolate, Transfer, Wavelength) %>%
  summarize(AbsMean = mean(Abs), AbsSd = sd(Abs)) %>% ungroup() %>%
  mutate(Isolate = ordered(Isolate, 1:13))

# Replace C11R2 isolates 2 T8 OD by T7 OD because of contamination at T8
isolates_OD <- as.data.table(isolates_OD)
temp <-  isolates_OD[Community == "C11R2" & Isolate == 2 & Transfer == 7, ]
isolates_OD[Community == "C11R2" & Isolate == 2 & Transfer == 8, ":="(AbsMean = temp$AbsMean, AbsSd = temp$AbsSd)]


# Remove C11R2 isolate 13 (streplococcus) contamination
isolates_OD <- isolates_OD %>% filter(!(Community == "C11R2" & Isolate == 13)) %>%
  filter(Transfer == 8, Wavelength == 620)

# Set OD values that are smaller than 0 into 0
isolates_OD$AbsMean[isolates_OD$AbsMean <=0] <- 0

#
isolates_OD %>% filter(Wavelength == 620, Transfer == 8)
fwrite(OD, here::here("data/temp/OD.csv"))
fwrite(isolates_OD, file = here::here("data/temp/isolates_OD.csv"))




