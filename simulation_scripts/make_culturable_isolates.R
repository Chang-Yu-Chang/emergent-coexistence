#' Filter the culturable isolates

suppressWarnings(suppressMessages(library(tidyverse)))
suppressWarnings(suppressMessages(library(data.table)))

args = commandArgs(trailingOnly = T)
input_set <- fread(args[[1]])
#input_set <- fread("../data/raw/simulation/mapping_files/input_set_simple_medium.csv")
temp <- args[[1]] %>% strsplit("/") %>% `[[`(1)
treatment <- sub("input_set_", "", temp[length(temp)]) %>% sub(".csv", "", .) # simple_medium
input_set_monoculture <- input_set %>% filter(grepl("monoculture", exp_id))

cat("\nMaking list of culturable isolates")
for (i in 1:nrow(input_set_monoculture)) {
    seed <- input_set_monoculture$seed[i]
    scale <- input_set_monoculture$scale[i]
    output_dir <- paste0(input_set_monoculture$output_dir[i])
    cat("\nexp_id = ", input_set_monoculture$exp_id[i])

    df_monoculture <- fread(paste0(output_dir, input_set_monoculture$exp_id[i], "_composition.txt"))
    df_monoculture_culturable <- df_monoculture %>%
        filter(Transfer == max(Transfer), Type == "consumer") %>%
        filter(Abundance >= 1/scale) %>%
        mutate(Seed = i) %>%
        select(Seed, Type, ID, Abundance)
    cat(",\tnumber of culturable isolates = ", nrow(df_monoculture_culturable))
    # arrange(desc(Abundance))
    fwrite(df_monoculture_culturable, file = paste0(output_dir, sub("-\\d+$", "", input_set_monoculture$exp_id[i]), "_culturable-", seed, ".txt"))
}
