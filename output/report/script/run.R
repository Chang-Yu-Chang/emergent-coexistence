#' Generate zsh script for running all R scripts saved in `output/report/script/`
library(tidyverse)

list.files(here::here("output/report/script/"), ".R") %>%
    str_subset("^((?!run).)*$") %>%
    paste0("echo  '${Running script}', ", ., "; Rscript ", .) %>%
    cat(file = here::here("output/report/script/run.sh"), sep = "\n")

