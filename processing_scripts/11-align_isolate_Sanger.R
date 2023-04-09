#' This scripts reads and aligns raw sanger sequences from Genewiz

if (FALSE) {
    install.packages("BiocManager") # Package for managing and installing Bioconductor packages
    BiocManager::install("sangeranalyseR")
}

library(tidyverse)
library(sangeranalyseR)
source(here::here("processing_scripts/00-metadata.R"))

alignment1 <- sangeranalyseR::SangerAlignment(
    ABIF_Directory = paste0(folder_data, "raw/sanger/sanger_seq_16S"),
    REGEX_SuffixForward = "F.ab1$",
    REGEX_SuffixReverse = "R.ab1$"
)

# This generates a merged sanger fa file: sanger_contigs_alignment.fa
sangeranalyseR::writeFasta(
    alignment1,
    outputDir = paste0(folder_data, "raw/sanger/sanger_seq_16S"),
    compress = FALSE,
    compression_level = NA,
    selection = "all"
)

consensus1 <- seqinr::read.fasta(paste0(folder_data, "raw/sanger/sanger_seq_16S/Sanger_contigs_alignment.fa"), seqtype = "DNA")
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

seq_C2R6 <- SangerRead(
    printLevel = "SangerRead",
    inputSource = "ABIF",
    readFeature = "Forward Read",
    readFileName = paste0(folder_data, "raw/sanger/2.6.A.5-27F.ab1")
)

iso_c2r6 <- tibble(ID = "2.6.A.5", Sequence = seq_C2R6@primarySeq)

isolates_16S <- bind_rows(df_seq, iso_c2r6)

write_csv(isolates_16S, paste0(folder_data, "temp/11-isolates_16S.csv"))

















