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

# Read the two isolates
# C2R6 isolate 4: CASEU_RN5_Dec2021_30-623993045_ab1/B2_T7_444_P2_28-27F.ab1
# C10R2 isolate 3: CASEU_RN5_Dec2021_30-623993045_ab1/B2_T7_444_P2_90-27F.ab1
if (FALSE) {
    plates <- read_csv(paste0(folder_data, "../old/data_old/raw/Sanger/CASEU_six_plates/plates_six_plates.csv"), show_col_types = F)
    plates %>%
        mutate(WellID = rep(1:96, 6)) %>%
        filter(
            (Community == "C10R2" & Isolate1 == Isolate2 & Isolate1 == 3 & MixIsolate == F)  |
                (Community == "C2R6" & Isolate1 == Isolate2 & Isolate1 == 4 & MixIsolate == F)
        ) %>%
        select(Community, Isolate1, WellID)

}

seq_C2R6 <- SangerRead(
    printLevel = "SangerRead",
    inputSource = "ABIF",
    readFeature = "Forward Read",
    readFileName = paste0(folder_data, "../old/data_old/raw/Sanger/CASEU_six_plates/CASEU_RN5_Dec2021_30-623993045_ab1/B2_T7_444_P2_28-27F.ab1")
)


seq_C10R2 <- SangerRead(
    printLevel = "SangerRead",
    inputSource = "ABIF",
    readFeature = "Forward Read",
    readFileName = paste0(folder_data, "../old/data_old/raw/Sanger/CASEU_six_plates/CASEU_RN5_Dec2021_30-623993045_ab1/B2_T7_444_P2_90-27F.ab1")
)

isolates_16S_two <- tibble(
#    ExpID = c("10.2.C.4", "2.6.A.5"),
    ID = c("10.2.C.4", "2.6.A.5"),
#    Community = c("C2R6", "C10R2"),
#    Isolate = c("4", "3"),
    Sequence = c(as.character(seq_C2R6@primarySeq), as.character(seq_C10R2@primarySeq))
)

#
isolates_16S <- bind_rows(df_seq, isolates_16S_two)

write_csv(isolates_16S, paste0(folder_data, "temp/11-isolates_16S.csv"))

















