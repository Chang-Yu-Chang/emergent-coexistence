library(tidyverse)
library(here)

source(here::here("processing_scripts/00-metadata.R"))

b <- batch_names[1]
list_folders[1]
paste0(folder_pipeline, "images/", b, "-", list_folders[1], "/")

cf <- function (x) length(list.files(x))

B2_ori <- paste0(folder_pipeline, "images/", b, "-", "00-original", "/")
cf(B2_ori)


read_csv(here("image_scripts/mapping_files/00-list_image_mapping-B2.csv")) # 138
read_csv(here("image_scripts/mapping_files/00-list_image_mapping-C.csv")) # 24
read_csv(here("image_scripts/mapping_files/00-list_image_mapping-C2.csv")) # 135
read_csv(here("image_scripts/mapping_files/00-list_image_mapping-D.csv")) # 180

138+24+135+180 # 477 coculture images
