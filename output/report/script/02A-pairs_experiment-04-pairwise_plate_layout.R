#' Match the isolate and pairs to the frozen stock saved in 96-well deep-well plates
library(tidyverse)

communities <- read_csv(here::here("data/output/communities.csv"))
# Plate layout in data.frame form ----
## Batch B2, plate 933 ----
plate_B2_933 <- data.frame(stringsAsFactors = F,
  Experiment = rep("Transitivity_B2", 96),
  PlateLayout = rep("933", 96),
  Well = paste0(rep(LETTERS[1:8], each = 12), sprintf("%02d", rep(1:12, 8))),
  MixIsolate = rep(TRUE, 96),
  Community = c(
    rep("C11R1", 9), rep("C8R4", 3),
    rep("C11R1", 9), rep("C8R4", 3),
    rep("C11R1", 9), rep("C8R4", 3),
    rep("C11R1", 9), rep("blank", 3),
    rep("C11R1", 9), rep("C10R2", 3),
    rep("C11R1", 9), rep("C10R2", 3),
    rep("C11R1", 9), rep("C10R2", 3),
    rep("C11R1", 9), rep("blank", 3)
  ),
  # Isolate1 by row
  Isolate1 = c(
    rep(1, 12), rep(2, 12), rep(3, 12),
    rep(4, 9), rep("blank", 3),
    rep(5, 9), rep(1, 3), rep(6, 9), rep(2, 3), rep(7, 9), rep(3, 3),
    rep(8, 9), rep("blank", 3)
  ),
  # Isolate2 by column
  Isolate2 = c(
    rep(c(1:9, 1:3), 3),
    1:9, rep("blank", 3),
    rep(c(1:9, 1:3), 3),
    1:9, rep("blank", 3)
  )
)

plate_B2_933$MixIsolate[plate_B2_933$Well %in% c(
  paste0("D", sprintf("%02d", 10:12)), # Blank
  paste0("H", sprintf("%02d", 10:12)) # Blank
)] <- FALSE


## Batch B2, plate 444 ----
plate_B2_444 <- data.frame(stringsAsFactors = F,
  Experiment = rep("Transitivity_B2", 96),
  PlateLayout = rep("444", 96),
  Well = paste0(rep(LETTERS[1:8], each = 12), sprintf("%02d", rep(1:12, 8))),
  MixIsolate = rep(TRUE, 96),
  Community = c(
    rep("C11R1", 9), rep("C8R4", 3),
    rep("C11R1", 9), rep("C10R2", 3),
    rep("blank", 12),
    rep("C2R6", 4), rep("C2R8", 4), rep("C7R1", 4),
    rep("C2R6", 4), rep("C2R8", 4), rep("C7R1", 4),
    rep("C2R6", 4), rep("C2R8", 4), rep("C7R1", 4),
    rep("C2R6", 4), rep("C2R8", 4), rep("C7R1", 4),
    rep("C2R6", 4), rep("C2R8", 4), rep("C7R1", 4)
  ),
  # Isolate1 by row
  Isolate1 = c(
    rep(9, 9), 1:3,
    1:9, 1:3,
    rep("blank", 12),
    rep(1:4, 3),
    rep(1:4, each = 12)
  ),
  # Isolate2 by column
  Isolate2 = c(
    1:9, 1:3, 1:9, 1:3,
    rep("blank", 12),
    rep(1:4, 15)
  )
)


plate_B2_444$MixIsolate[plate_B2_444$Well %in% c(
  paste0("A", sprintf("%02d", 10:12)),
  paste0("B", sprintf("%02d", 1:12)),
  paste0("C", sprintf("%02d", 1:12)), # Blank
  paste0("D", sprintf("%02d", 1:12))
)] <- FALSE


## Batch C, plate C11R1 ----
plate_C_C11R1 <- data.frame(stringsAsFactors = F,
  Experiment = rep("Transitivity_C", 96),
  PlateLayout = rep("C11R1", 96),
  Well = paste0(rep(LETTERS[1:8], each = 12), sprintf("%02d", rep(1:12, 8))),
  MixIsolate = rep(TRUE, 96),
  Community = c(
    rep("C11R1", 6 * 12 + 9),
    rep("blank", 15)
  ),
  # Isolate1 by row
  Isolate1 = c(
    rep(c(rep(1, 9), rep(5, 3)), 4),
    1:9, rep(6, 3),
    1:9, rep(7, 3),
    1:9,
    rep("blank", 15)
  ),
  # Isolate2 by column
  Isolate2 = c(
    1:9, rep(6, 3),
    1:9, rep(7, 3),
    1:9, rep(6, 3),
    1:9, rep(7, 3),
    rep(1, 9), rep(5, 3),
    rep(1, 9), rep(5, 3),
    1:9,
    rep("blank", 15)
  )
)

plate_C_C11R1$MixIsolate[plate_C_C11R1$Well %in% c(
  paste0("G", sprintf("%02d", 1:12)), # G10-G12 are blank
  paste0("H", sprintf("%02d", 1:12))
)] <- FALSE


## Batch C2, plate 13A ----
plate_C2_13A <- data.frame(stringsAsFactors = F,
  Experiment = rep("Transitivity_C2", 96),
  PlateLayout = rep("13A", 96),
  Well = paste0(rep(LETTERS[1:8], each = 12), sprintf("%02d", rep(1:12, 8))),
  MixIsolate = rep(TRUE, 96),
  Community = c(
    rep("C11R2", 96)
  ),
  # Isolate1 by row
  Isolate1 = c(
    rep(1:8, each = 12)
  ),
  # Isolate2 by column
  Isolate2 = c(
    rep(1:12, 8)
  )
)


## Batch C2, plate 13B ----
plate_C2_13B <- data.frame(stringsAsFactors = F,
  Experiment = rep("Transitivity_C2", 96),
  PlateLayout = rep("13B", 96),
  Well = paste0(rep(LETTERS[1:8], each = 12), sprintf("%02d", rep(1:12, 8))),
  MixIsolate = rep(TRUE, 96),
  Community = c(
    rep("C11R2", 5 * 12 + 2),
    rep("blank", 10),
    rep("C11R2", 2 * 12)
  ),
  # Isolate1 by row
  Isolate1 = c(
    rep(9:13, each = 12),
    13, 13,
    rep("blank", 10),
    1:12,
    1:12
  ),
  # Isolate2 by column
  Isolate2 = c(
    rep(1:12, 5),
    13, 13,
    rep("blank", 10),
    rep(13, 12),
    1:12
  )
)

plate_C2_13B$MixIsolate[plate_C2_13B$Well %in% c(
  paste0("F", sprintf("%02d", 2:12)), # F03-G12 are blank
  paste0("H", sprintf("%02d", 1:12))
)] <- FALSE




## Batch D, plate 75 ----
plate_D_75 <- data.frame(stringsAsFactors = F,
  Experiment = rep("Transitivity_D", 96),
  PlateLayout = rep("75", 96),
  Well = paste0(rep(LETTERS[1:8], each = 12), sprintf("%02d", rep(1:12, 8))),
  MixIsolate = rep(TRUE, 96),
  Community = c(
    rep(c(rep("C1R7", 7), rep("C11R5", 5)), 6),
    rep("C1R7", 7), rep("blank", 5),
    rep("C1R7", 7), rep("blank", 5)
  ),
  # Isolate1 by row
  Isolate1 = c(
    rep(1:5, each = 12),
    rep(6, 7), 1:5,
    rep(7, 7), rep("blank", 5),
    1:7, rep("blank", 5)
  ),
  # Isolate2 by column
  Isolate2 = c(
    rep(c(1:7, 1:5), 6),
    1:7, rep("blank", 5),
    1:7, rep("blank", 5)
  )
)

plate_D_75$MixIsolate[plate_D_75$Well %in% c(
  paste0("F", sprintf("%02d", 8:12)),
  paste0("G", sprintf("%02d", 8:12)),
  paste0("H", sprintf("%02d", 1:12))
)] <- FALSE

## Batch D, plate 5543 ----
plate_D_5543 <- data.frame(stringsAsFactors = F,
  Experiment = rep("Transitivity_D", 96),
  PlateLayout = rep("5543", 96),
  Well = paste0(rep(LETTERS[1:8], each = 12), sprintf("%02d", rep(1:12, 8))),
  MixIsolate = rep(TRUE, 96),
  Community = c(
    rep(c(rep("C1R4", 6), rep("C1R6", 6)), 5),
    rep(c(rep("C1R2", 8), rep("C4R1", 4)), 2),
    rep("blank", 4),
    rep("C1R2", 4), rep("C4R1", 4)
  ),
  # Isolate1 by row
  Isolate1 = c(
    rep(1:5, each = 12),
    rep(1, 4), rep(3, 4), rep(1, 4),
    rep(2, 4), rep(4, 4), rep(2, 4),
    rep("blank", 4),
    1:4, rep(3, 4)
  ),
  # Isolate2 by column
  Isolate2 = c(
    1:5, 1, 1:5, 1,
    1:5, 2, 1:5, 2,
    1:5, 3, 1:5, 3,
    1:5, 4, 1:5, 4,
    1:5, 5, 1:5, 5,
    rep(1:4, 2), 1:3, 1,
    rep(1:4, 2), 1:3, 2,
    rep("blank", 4),
    1:4, 1:3, 3
  )
)

plate_D_5543$MixIsolate[plate_D_5543$Well %in% c(
  paste0(rep(LETTERS[1:5], each = 2), sprintf("%02d", c(6, 12))),
  paste0("H", sprintf("%02d", 1:8)),
  paste0(LETTERS[6:8], sprintf("%02d", 12))
)] <- FALSE









# Make merged df ----
plate_B2_933_P1 <- mutate(plate_B2_933, Plate = "P1", Isolate1Freq = 50, Isolate2Freq = 50)
plate_B2_933_P2 <- mutate(plate_B2_933, Plate = "P2", Isolate1Freq = 95, Isolate2Freq = 5)
plate_B2_444_P1 <- mutate(plate_B2_444, Plate = "P1", Isolate1Freq = 50, Isolate2Freq = 50)
plate_B2_444_P2 <- mutate(plate_B2_444, Plate = "P2",Isolate1Freq = 95, Isolate2Freq = 5)
plate_C_C11R1_P1 <- mutate(plate_C_C11R1, Isolate1Freq = 50, Isolate2Freq = 50)
plate_C_C11R1_P1$Isolate1Freq[plate_C_C11R1_P1$Well %in% paste0(rep(LETTERS[3:6], each = 12), sprintf("%02d", rep(1:12, 4)))] <- 95
plate_C_C11R1_P1$Isolate2Freq[plate_C_C11R1_P1$Well %in% paste0(rep(LETTERS[3:6], each = 12), sprintf("%02d", rep(1:12, 4)))] <- 5
plate_C_C11R1_P1 <- mutate(plate_C_C11R1_P1, Plate = "P1")
plate_C2_13A_P1 <- mutate(plate_C2_13A, Plate = "P1", Isolate1Freq = 50, Isolate2Freq = 50)
plate_C2_13A_P2 <- mutate(plate_C2_13A, Plate = "P2", Isolate1Freq = 95, Isolate2Freq = 5)
plate_C2_13B_P1 <- mutate(plate_C2_13B, Plate = "P1",Isolate1Freq = 50, Isolate2Freq = 50)
plate_C2_13B_P2 <- mutate(plate_C2_13B, Plate = "P2",Isolate1Freq = 95, Isolate2Freq = 5)
plate_D_75_P1 <- mutate(plate_D_75, Plate = "P1", Isolate1Freq = 50, Isolate2Freq = 50)
plate_D_75_P2 <- mutate(plate_D_75, Plate = "P2", Isolate1Freq = 95, Isolate2Freq = 5)
plate_D_5543_P1 <- mutate(plate_D_5543, Plate = "P1", Isolate1Freq = 50, Isolate2Freq = 50)
plate_D_5543_P2 <- mutate(plate_D_5543, Plate = "P2", Isolate1Freq = 95, Isolate2Freq = 5)

plates <- rbind(plate_B2_933_P1, plate_B2_933_P2, plate_B2_444_P1, plate_B2_444_P2, plate_C_C11R1_P1, plate_C2_13A_P1, plate_C2_13A_P2,
  plate_C2_13B_P1, plate_C2_13B_P2, plate_D_75_P1, plate_D_75_P2, plate_D_5543_P1, plate_D_5543_P2)


# Match ambiguous pairs to well position on plates ----
pairs_ambiguous <- read_csv(here::here("data/temp/pairs_ambiguous.csv"))

## Switch the isolate1 and isolate2 since the P1 is 50:50 and
## rows and columns in P2 P3 are for 95 and 5 respectively
pairs_ambiguous[pairs_ambiguous$Isolate1Freq == 5, c("Isolate1", "Isolate2", "Isolate1Freq", "Isolate2Freq")] <- pairs_ambiguous[pairs_ambiguous$Isolate1Freq == 5, c("Isolate2", "Isolate1", "Isolate2Freq", "Isolate1Freq")]

## Mutate the isolate variables
pairs_ambiguous <- pairs_ambiguous %>% mutate(Isolate1 = factor(Isolate1), Isolate2 = factor(Isolate2))

## Join plates
pairs_ambiguous_on_DW96 <- pairs_ambiguous %>%
  left_join(plates, by = c("Community", "Isolate1", "Isolate2", "Experiment", "Isolate1Freq", "Isolate2Freq")) %>%
  mutate(Plate = ifelse((Isolate1Freq != 50 & PlateLayout != "C11R1") , "P2", "P1"),
    Community = ordered(Community, communities$Community)) %>%
  distinct(Community, Isolate1, Isolate2, Isolate1Freq, Isolate2Freq, .keep_all = T) %>%
  select(-MixIsolate)




# Across-community and random assembly ----
# Read data
isolates <- read_csv(here::here("data/output/isolates.csv"))
community_names <- isolates %>% filter(str_detect(Community, "Ass")) %>% pull(Community) %>% unique
#community_names_sizes <- rep(8, 4)
myColor <- c(AcrAss1 = "#ED6A5A", AcrAss2 = "#53A2BE", RanAss1 = "#FFD23F", RanAss2 = "#2CA58D", blank = "#BFBFBF", EP = "purple")
plate_names <- c("AD_P1", "BD_P1", "C_P1", "AD_P2", "BD_P2", "C_P2")
plate_names_DNA <- c("AB_P1", "CD_P1", "AD_P2", "BD_P2", "C_P2")

# Plate layout ----
well_names <- paste0(rep(LETTERS[1:8], 12), sprintf("%02d", rep(1:12, each = 8)))
random_community_names <- c("AcrAss1", "AcrAss2", "RanAss1", "RanAss2")

layout_AD_P1 <- data.frame(
    PlateLayout = "AD",
    MixPlate = "P1",
    Well = well_names,
    Community = c(rep("AcrAss1", 64), rep("RanAss2", 32)),
    Isolate1 = rep(1:8, 12),
    Isolate2 = c(rep(1:8, each = 8), rep(1:4, each = 8)),
    MixIsolate = rep(T, 96),
    Isolate1Freq = rep(0.5, 96),
    Isolate2Freq = rep(0.5, 96)
)

layout_AD_P2 <- layout_AD_P1 %>%
    mutate(Isolate1Freq = rep(0.95, 96), Isolate2Freq = rep(0.05, 96),
           MixPlate = "P2")

layout_AD_P3 <- layout_AD_P2 %>% mutate(MixPlate = "P3")

layout_BD_P1 <- data.frame(
    PlateLayout = "BD",
    MixPlate = "P1",
    Well = well_names,
    Community = c(rep("AcrAss2", 64), rep("RanAss2", 32)),
    Isolate1 = rep(1:8, 12),
    Isolate2 = c(rep(1:8, each = 8), rep(5:8, each = 8)),
    MixIsolate = rep(T, 96),
    Isolate1Freq = rep(0.5, 96),
    Isolate2Freq = rep(0.5, 96)
)

layout_BD_P2 <- layout_BD_P1 %>%
    mutate(Isolate1Freq = rep(0.95, 96), Isolate2Freq = rep(0.05, 96),
           MixPlate = "P2")

layout_BD_P3 <- layout_BD_P2 %>% mutate(MixPlate = "P3")

layout_C_P1 <- data.frame(
    PlateLayout = "C",
    MixPlate = "P1",
    Well = well_names,
    Community = c(rep("RanAss1", 64), rep("blank", 16), rep("EP", 4), rep("blank", 4), rep("EP", 4), rep("blank", 4)),
    Isolate1 = c(rep(1:8, 8), rep(NA, 16), rep(1, 4), rep(NA, 4), rep(1, 3), 2, rep(NA, 4)),
    Isolate2 = c(rep(1:8, each = 8), rep(NA, 16), rep(2, 3), 1, rep(NA, 4), rep(2, 4), rep(NA, 4)),
    MixIsolate = c(rep(T, 64), rep(F, 32)),
    Isolate1Freq = c(rep(0.5, 64), rep(NA, 16), 0.5, 0.95, 0.05, rep(NA, 5), 0.5, 0.95, 0.05, rep(NA, 5)),
    Isolate2Freq = c(rep(0.5, 64), rep(NA, 16), 0.5, 0.05, 0.95, rep(NA, 5), 0.5, 0.95, 0.05, rep(NA, 5))
)


layout_C_P2 <- data.frame(
    PlateLayout = "C",
    MixPlate = "P2",
    Well = well_names,
    Community = c(rep("RanAss1", 64), rep(random_community_names[1:4], each = 8)),
    Isolate1 = c(rep(1:8, 8), rep(1:8, 4)),
    Isolate2 = c(rep(1:8, each = 8), rep(1:8, 4)),
    MixIsolate = c(rep(T, 64), rep(F, 32)),
    Isolate1Freq = c(rep(0.95, 64), rep(NA, 32)),
    Isolate2Freq = c(rep(0.05, 64), rep(NA, 32))
)

layout_C_P3 <- layout_C_P2 %>% mutate(MixPlate = "P3")



# Merge plate layouts ----
plates_random <- bind_rows(layout_AD_P1, layout_AD_P2, layout_AD_P3,
                           layout_BD_P1, layout_BD_P2, layout_BD_P3,
                           layout_C_P1, layout_C_P2, layout_C_P3) %>%
    as_tibble()


# Plot the plate layout
plate_layout_name <- plates_random %>%
    arrange(MixPlate, PlateLayout) %>%
    unite("PlateMixLayout", PlateLayout, MixPlate) %>%
    pull(PlateMixLayout) %>% unique

p_random_network_plates_list <- rep(list(NA), length(plate_layout_name))
names(p_random_network_plates_list) <- plate_layout_name

for (i in 1:length(plate_layout_name)) {
    # Unite the plate and well names
    temp <- plates_random %>%
        unite("PlateMixLayout", PlateLayout, MixPlate) %>%
        filter(PlateMixLayout == plate_layout_name[i]) %>%
        mutate(FillLabel = Community, Isolate1Freq = Isolate1Freq * 100, Isolate2Freq = Isolate2Freq * 100) %>%
        unite("temp1", Isolate1, Isolate2, sep = "_") %>%
        unite("temp2", Isolate1Freq, Isolate2Freq, sep = ":") %>%
        unite("TextLabel", temp1, temp2, sep = "\n") %>%
        select(Community, Well, FillLabel, TextLabel)

    # Monoculture
    temp$TextLabel[grepl("NA", temp$TextLabel)] <- substr(temp$TextLabel, 1, 1)[grepl("NA", temp$TextLabel)]

    # Plot plate
    # p_random_network_plates_list[[i]] <- temp %>%
    #     draw_plate_from_df(fill_legend = F, annotation = T) +
    #     scale_fill_manual(values = myColor)
}



write_csv(plates_random, here::here("data/output/plates_random.csv"))
write_csv(plates, here::here("data/output/plates.csv"))











