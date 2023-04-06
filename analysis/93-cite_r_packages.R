citeRpacks <- function(pkg_list, filename, RStudio = FALSE) {

    #ht to https://stackoverflow.com/questions/2470248/write-lines-of-text-to-a-file-in-r for sink()
    for (i in 1:length(pkg_list)) {
        sink(file = paste(filename, ".bib", sep = ""), append = T)
        writeLines(toBibtex(citation(package = pkg_list[i])))
        sink()
    }

  if(RStudio) {

    c <- RStudio.Version()$citation

    sink(file = paste(filename, ".bib", sep = ""), append = T)
    writeLines(paste("@Manual{,",
                     "\n   title = {", c$title, "},",
                     "\n   author = {{", c$author, "}},",
                     "\n   organization = {", c$organization, "},",
                     "\n   address = {", c$address, "},",
                     "\n   year = {", c$year, "},",
                     "\n   note = {", RStudio.Version()$version, "}",
                     "\n   url = {", c$url, "},",
                     "\n }",
                     sep = ""))
    sink()

  }

}

my_packs <- c("ggplot2", "base", "rRDP", "Biostrings", "EBImage", "randomForest", "caret", "igraph")
citeRpacks(my_packs, "citeR", RStudio = TRUE)


