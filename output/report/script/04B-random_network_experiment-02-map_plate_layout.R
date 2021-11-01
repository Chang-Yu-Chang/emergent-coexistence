#' Map the plate layout of random network

# Read data ----
isolates_assembly <- read_csv(here::here("data/temp/isolates_assembly.csv"))
community_names <- unique(isolates_assembly$AssemblyCommunity)
community_names_sizes <- rep(8, 4)
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
random_assembly_plates <- bind_rows(layout_AD_P1, layout_AD_P2, layout_AD_P3,
  layout_BD_P1, layout_BD_P2, layout_BD_P3,
  layout_C_P1, layout_C_P2, layout_C_P3) %>%
  as_tibble()


# Plot the plate layout
plate_layout_name <- random_assembly_plates %>%
  arrange(MixPlate, PlateLayout) %>%
  unite("PlateMixLayout", PlateLayout, MixPlate) %>%
  pull(PlateMixLayout) %>% unique

p_random_network_plates_list <- rep(list(NA), length(plate_layout_name))
names(p_random_network_plates_list) <- plate_layout_name

for (i in 1:length(plate_layout_name)) {
  # Unite the plate and well names
  temp <- random_assembly_plates %>%
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
  p_random_network_plates_list[[i]] <- temp %>%
    draw_plate_from_df(fill_legend = F, annotation = T) +
    scale_fill_manual(values = myColor)
}





# Old code -----

if (FALSE) {


  # Plate layout in DNA extraction ----
  temp_vector <- NULL; temp_vector2 <- NULL; temp_vector3 <- NULL; temp_vector4 <- NULL
  for (i in 1:7) {
    # Upper triangle
    temp_vector <- c(temp_vector, paste0(1:i, "_", rep(i+1, i)))
    temp_vector2 <- c(temp_vector2, paste0(LETTERS[1:i], sprintf("%02d", rep(i+1, i))))
    # Lower triangle
    temp_vector3 <- c(temp_vector3, paste0((8-i):1, "_", rep((9-i), 8-i)))
    temp_vector4 <- c(temp_vector4, paste0(LETTERS[(i+1):8], sprintf("%02d", rep(i, (8-i)))))
  }

  # AB P1 ----
  DNA_AB_P1_AcrAss1 <- tibble(
    FillLabel = "AcrAss1",
    TextLabel = paste0(temp_vector, "\n50:50"),
    Well = temp_vector2
  )

  DNA_AB_P1_AcrAss2 <- tibble(
    FillLabel = "AcrAss2",
    TextLabel = paste0(temp_vector3, "\n50:50"),
    Well = temp_vector4
  )

  DNA_AB_P1_EP <- tibble(
    FillLabel = "EP",
    TextLabel = c(paste0(rep("E", 6), "_", rep("P", 6), "\n", rep(c("50:50", "95:5", "5:95"), each = 2)), "E", "P"),
    Well = paste0(rep(LETTERS[1:4], each = 2), sprintf("%02d", rep(11:12, 4)))
  )

  p_DNA_AB_P1 <- bind_rows(DNA_AB_P1_AcrAss1, DNA_AB_P1_AcrAss2, DNA_AB_P1_EP) %>%
    draw_plate_from_df(fill_legend = F, annotation = T) +
    scale_fill_manual(values = myColor)

  # CD P1 ----
  DNA_CD_P1_RanAss1 <- tibble(
    FillLabel = "RanAss1",
    TextLabel = paste0(temp_vector, "\n50:50"),
    Well = temp_vector2
  )

  DNA_CD_P1_RanAss2 <- tibble(
    FillLabel = "RanAss2",
    TextLabel = paste0(temp_vector3, "\n50:50"),
    Well = temp_vector4
  )

  p_DNA_CD_P1 <- bind_rows(DNA_CD_P1_RanAss1, DNA_CD_P1_RanAss2) %>%
    draw_plate_from_df(fill_legend = F, annotation = T) +
    scale_fill_manual(values = myColor)

  # AD P2
  p_DNA_AD_P2 <- p_AD_P2
  # BD P2
  p_DNA_BD_P2 <- p_BD_P2
  # C P2
  p_DNA_C_P2 <- p_C_P2


  # Merge plate
  plates_DNA_list <- list(p_DNA_AB_P1, p_DNA_CD_P1, p_DNA_AD_P2, p_DNA_BD_P2, p_DNA_C_P2)
  names(plates_DNA_list) <- paste0("DNA_", c("AB_P1", "CD_P1", "AD_P2", "BD_P2", "C_P2"))
  p_plates_DNA <- cowplot::plot_grid(plotlist = plates_DNA_list, ncol = 3, labels = c(plate_names_DNA[1:2], "", plate_names_DNA[3:5]))


  # Save plate figure ----
  pdf(here::here("output/protocol/tab_fig/protocol_20191205_random_network_T2-plate_layout.pdf"), width = 18, height = 10); p_plates; dev.off()
  pdf(here::here("output/protocol/tab_fig/protocol_20191205_random_network_T2-DNA_plate_layout.pdf"), width = 18, height = 10); p_plates_DNA; dev.off()
  save(plates_list, plates_DNA_list, p_plates, p_plates_DNA, file = here::here("output/protocol/tab_fig/protocol_20191205_random_network_T2-plate_layout.Rdata"))

}
