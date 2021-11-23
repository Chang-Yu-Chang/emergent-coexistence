#' Read OD data from random network experiment

# Read data
random_assembly_plates <- fread(root$find_file("data/temp/random_assembly_plates.csv"))


# Read OD data ----
## R function for reading raw OD data
OD_file_read <- function (
  file_name_pattern = "20191205",
  folder_directory = "~/ownCloud/Chang-Yu/random_network/"
) {
  # Read raw spec output
  temp_df <- list.files(folder_directory, pattern = file_name_pattern, full.names = T) %>%
    lapply(fread) %>%
    rbindlist() %>%
    select(Plate, Well, Wavelength, Abs)

  # Subset
  ## Blank
  OD_blank <- filter(temp_df, grepl("blank", Plate)) %>%
    mutate(AbsBlank = Abs) %>%
    select(Well, Wavelength, AbsBlank) %>%
    as_tibble()

  ## Sample
  OD_sample <- filter(temp_df, !grepl("blank", Plate)) %>%  as_tibble()

  # Substract the OD of sample by blank's OD
  OD_sample <- OD_sample %>%
    left_join(OD_blank, by = c("Well", "Wavelength")) %>%
    mutate(Abs = Abs - AbsBlank) %>%
    select(-AbsBlank) %>%
    separate(Plate, into = c("Transfer", "PlateLayout", "MixPlate"), sep = "_")

  return(OD_sample)
}


## Read OD raw data
OD_T1 <- OD_file_read(file_name_pattern = "20191205", folder_directory = "~/ownCloud/Chang-Yu/random_network/")
OD_T2 <- OD_file_read(file_name_pattern = "20191207", folder_directory = "~/ownCloud/Chang-Yu/random_network/")
OD_T3 <- OD_file_read(file_name_pattern = "20191209", folder_directory = "~/ownCloud/Chang-Yu/random_network/")
OD <- bind_rows(OD_T1, OD_T2, OD_T3)


# Join OD and plate layout ----
OD <- OD %>%
  left_join(random_assembly_plates, by = c("PlateLayout", "MixPlate", "Well")) %>%
  mutate(Isolate1Freq = factor(Isolate1Freq, c(0.05, 0.5, 0.95)),
    Isolate2Freq = factor(Isolate2Freq, c(0.05, 0.5, 0.95))) %>%
  as_tibble()

temp_index <- which(OD$Isolate2 > OD$Isolate1)
OD[temp_index, c("Isolate1", "Isolate2")] <- OD[temp_index, c("Isolate2", "Isolate1")]
OD[temp_index, c("Isolate1Freq", "Isolate2Freq")] <- OD[temp_index, c("Isolate2Freq", "Isolate1Freq")]




