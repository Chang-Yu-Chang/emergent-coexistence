#' This scripts reads and aligns raw sanger sequences from Genewiz

if (FALSE) {
    install.packages("BiocManager") # Package for managing and installing Bioconductor packages
    BiocManager::install("sangeranalyseR")
}

library(tidyverse)
library(sangeranalyseR)
source(here::here("analysis/00-metadata.R"))

alignment1 <- sangeranalyseR::SangerAlignment(
    ABIF_Directory = paste0(folder_data, "raw/sanger/sanger_seq_16S_communities"),
    REGEX_SuffixForward = "F.ab1$",
    REGEX_SuffixReverse = "R.ab1$"
)

# This generates a merged sanger fa file: sanger_contigs_alignment.fa
sangeranalyseR::writeFasta(
    alignment1,
    outputDir = paste0(folder_data, "raw/sanger/sanger_seq_16S_communities"),
    compress = FALSE,
    compression_level = NA,
    selection = "all"
)

#
consensus1 <- seqinr::read.fasta(paste0(folder_data, "raw/sanger/sanger_seq_16S_communities/Sanger_contigs_alignment.fa"), seqtype = "DNA")
consensus1 <- consensus1 %>% lapply(function(x) {
    contig <- x %>% as.character()
    contig[contig=="a"] <- "A"
    contig[contig=="t"] <- "T"
    contig[contig=="g"] <- "G"
    contig[contig=="c"] <- "C"
    contig <- paste(contig, collapse = "")
    contig <- gsub("-", "", contig)
    return(data.frame(Sequence = contig))
})

df_seq <- bind_rows(consensus1, .id = "temp") %>%
    tidyr::separate(col = "temp", into = c("ID", "t1", "t2", "t3")) %>%
    dplyr::select(ID, Sequence) %>%
    dplyr::mutate(ConsensusLength = nchar(Sequence))

write_csv(df_seq, paste0(folder_data, "temp/11-isolates_16S.csv"))






