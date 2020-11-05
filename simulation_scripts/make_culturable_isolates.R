#' Filter the culturable isolates

suppressWarnings(suppressMessages(library(tidyverse)))
suppressWarnings(suppressMessages(library(data.table)))

args = commandArgs(trailingOnly = T)
input_independent <- fread(args[[1]])
#input_independent <- fread("../data/raw/simulation/mapping_files/input_independent.csv")
input_independent_monoculture <- input_independent %>% filter(grepl("monoculture", exp_id))

cat("\nMaking list of culturable isolates")
for (i in 1:nrow(input_independent_monoculture)) {
    seed <- input_independent_monoculture$seed[i]
    scale <- input_independent_monoculture$scale[i] 
    output_dir <- input_independent_monoculture$output_dir[i]
    cat("\nexp_id = ", input_independent_monoculture$exp_id[i])
    
    df_monoculture <- fread(paste0(output_dir, "monoculture-", seed, "_composition.txt"))
    df_monoculture_culturable <- df_monoculture %>% 
        filter(Transfer == max(Transfer), Type == "consumer") %>%
        filter(Abundance >= 1/scale) %>%
        mutate(Seed = i) %>% 
        select(Seed, Type, ID, Abundance)
    cat(",\tnumber of culturable isolates = ", nrow(df_monoculture_culturable))
    # arrange(desc(Abundance))
    fwrite(df_monoculture_culturable, file = paste0(output_dir, "monoculture-culturable-", seed, ".txt"))
}
