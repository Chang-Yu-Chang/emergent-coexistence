# Check the download files on cluster before download
suppressWarnings(suppressMessages(library(tidyverse)))
suppressWarnings(suppressMessages(library(data.table)))


test_file_existence <- function (mapping_file_directory, mapping_file_name, simulation_file_directory) {
    stopifnot(!is.na(list.files(mapping_file_directory, mapping_file_name)))
    cat("\n", paste0(rep("-", getOption("width")), collapse = ""))
    cat("\nInput csv: ", mapping_file_name)
    #cat("\nInput csv: ", paste0(mapping_file_directory, mapping_file_name))
    input_csv <- fread(paste0(mapping_file_directory, mapping_file_name))
    x <- paste0(simulation_file_directory, input_csv$exp_id, "_composition.txt")
    if (sum(file.exists(x)) == nrow(input_csv)) {
        cat("\tAll files exist")
    } else {
        cat("\tNumber of rows: ", nrow(input_csv))
        cat("\t", length(input_csv$exp_id[!file.exists(x)]), " files missing")
    }
}

mapping_files <- c("input_set_simple_medium", "input_set_rich_medium",
                   "input_synthetic_simple_medium", "input_synthetic_rich_medium")

# mapping_file_directory <- "/home/cc2553/project/community-selection/wrapper/"
# simulation_file_directory <- "/home/cc2553/project/community-selection/data/"
mapping_file_directory <- "/home/cc2553/project/invasion-network/wrapper/"
simulation_file_directory <- "/home/cc2553/project/invasion-network/data/"
#mapping_file_directory <- "~/Desktop/Lab/invasion-network/data/raw/simulation/mapping_files/"
#simulation_file_directory <- "~/Desktop/Lab/invasion-network/data/raw/simulation/"


for (i in 1:length(mapping_files)) {
    mfd <- paste0(mapping_file_directory, sub("input_", "", mapping_files[i]), "/")
    sfd <- paste0(simulation_file_directory, sub("input_", "", mapping_files[i]), "/")

    test_file_existence(#mapping_file_directory = paste0(mapping_file_directory, sub("input_", "", mapping_files[i]), "/"),
        mapping_file_directory = mfd,
        paste0(mapping_files[i], ".csv"),
        # simulation_file_directory = paste0(simulation_file_directory, sub("input_", "", mapping_files[i]), "/"))
        simulation_file_directory = sfd)
}
cat("\n")





