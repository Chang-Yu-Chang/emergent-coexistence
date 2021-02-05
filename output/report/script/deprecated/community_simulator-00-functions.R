#' The script contains two R functions
#'
#' 1. `make_list_of_algorithm()`: read the lines in the list of algorithms saved in a single .py file
#' and output a csv that contains the author names and  algorithm names
#' 2. `read_bind_csv()`: read and rbind the simulated result
#'



# Read the lines in the python files ----
# output a csv that contains the author names and  algorithm names
make_list_of_algorithm <- function (
  protocol_function = "selection_function",
  file_name = "script/community_simulator-03A-community_selection_student.py" # python file that saves the algorithms
){

  # Read line in python files
  algorithm_names <-
    read_lines(file_name) %>%
    grep("^def\\s|----", ., value = T)

  # Write the dataframe that containts names of student selection function into output txt. I use repexp in R.
  num_func <- sum(!grepl("# ", algorithm_names)) # Number of functions
  author_names <- grep("# ", algorithm_names, value = T) %>% gsub("# | -*$", "", .) # Author names
  temp_index <- c(which(grepl("# ", algorithm_names)), length(algorithm_names)+1)
  temp_number = (temp_index[-1] - temp_index[-length(temp_index)])-1

  # Function index and function name cleanup
  func_index = rep(1:length(author_names), temp_number)
  func_names = grep("def ", algorithm_names, value = T)  %>%
    gsub("def ", "", .) %>%
    gsub("\\(\\w.+", "", .)

  # Make the algorithm into dataframe
  algorithm_list <-
    data.frame(Algorithm = rep(protocol_function, num_func), # Variable saved for different algorithms, e.g. migration, community_function...
               Author = author_names[func_index],
               AlgorithmName = func_names) %>%
    mutate(AlgorithmID = 0:(nrow(.)-1))

  # Write the list into csv
  fwrite(algorithm_list, paste0("script/data/", protocol_function, "_list.csv"), sep = ",")
  return(algorithm_list)
}




# Read and rbind the simulated result csv files ----
read_bind_csv <- function (community_phenotype) {
  # Read csv
  result_list <- list.files("script/data/", paste0("[0-9]+-", community_phenotype))
  temp_result <- rep(list(NA), length(result_list)) # Make empty R list
  for (f in 1:length(result_list)) temp_result[[f]] <- fread(paste0("script/data/", result_list[f])) # Read csv
  result <- rbindlist(temp_result) # Bind the R list

  # Sort algorithms
  list_algorithm <- fread("script/data/list_algorithm.csv")
  result <- result %>%
    mutate(SelectionFunction = ordered(SelectionFunction, level = list_algorithm$AlgorithmName[list_algorithm$Algorithm == "selection_function"]),
           MigrationFunction = ordered(MigrationFuntion, level = list_algorithm$AlgorithmName[list_algorithm$Algorithm == "migration_function"]))


  # Write the result
  fwrite(result, file = paste0("script/data/binded-", community_phenotype, ".csv"))
}



# Read community and resource composition by transfer

read_bind_transfer_csv <- function (community_phenotype) {
  # Csv file list
  temp_file <- list.files("script/data/passages/", pattern = community_phenotype)
  temp_result <- rep(list(NA), length(temp_file)) # Make empty R list

  #
  temp_file_df <- temp_file %>%
    gsub(".csv", "", .) %>%
    strsplit("-") %>% unlist %>%
    matrix(ncol = 9, byrow = T) %>% as.data.frame() %>%
    setNames(c("CommunityPhenotypeID", "CommunityPhenotype",
               "SelectionFunctionID", "SelectionFunction",
               "MigrationFuntionID", "MigrationFunction",
               "Replicate", "Passage", "Type")) # Type indicates species or resources

  for (i in 1:length(temp_file)) {
    if (grepl("N.csv", temp_file[i])) { # Species
      temp_result[[i]] <-
        fread(paste0("script/data/passages/", temp_file[i])) %>%
        gather(key = "Well", value = "Abundance") %>%
        mutate(TypeID = rep(1:96, 210)) %>%
        cbind(temp_file_df[i,])

    } else if (grepl("R.csv", temp_file[i])) {
      temp_result[[i]] <-
        fread(paste0("script/data/passages/", temp_file[i])) %>%
        gather(key = "Well", value = "Abundance") %>%
        mutate(TypeID = rep(1:96, 90)) %>%
        cbind(temp_file_df[i,])

    }
  }

  result <- rbindlist(temp_result)

  # Write the result
  fwrite(result, file = paste0("script/data/binded-passage-", community_phenotype, ".csv"))
}



