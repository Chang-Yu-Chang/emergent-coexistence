# Read and align raw Sanger sequences from Genewiz

#' For `R >= 3.5.0`, use `BiocManager::install()` instead of `biocLite()` to
#' install Bioconductor-maintained packages. This may take a few seconds.

if (FALSE) {
    install.packages("BiocManager") # Package for managing and installing Bioconductor packages
    BiocManager::install("Biostrings") # DNA string format
    BiocManager::install("DECIPHER") # Dependent packages of sangeranalyseR
    BiocManager::install("sangeranalyseR")
    BiocManager::install("sangerseqR") # Align forward and forward Sanger sequences
    BiocManager::install("rRDP") # ribosomal database project
    BiocManager::install("rRDPData") # 16S Reference database for RDP
    BiocManager::install("msa") # Multiple sequence alignment
    BiocManager::install("ggtree") # Make figure of phylogenetic tree
    install.packages("ape") # Analyses of Phylogenetics and Evolution
}


library(DECIPHER)
library(sangeranalyseR)
library(sangerseqR)
library(rRDP)
library(rRDPData)
library(ggtree)
library(tidyverse)
library(sangeranalyseR)

#
alignment1 <- sangeranalyseR::SangerAlignment(ABIF_Directory = "~/Dropbox/lab/emergent-coexistence/data/raw/Sanger/sanger_seq_16S_20180702",
                                              REGEX_SuffixForward = "F.ab1$",
                                              REGEX_SuffixReverse = "R.ab1$")
sangeranalyseR::writeFasta(alignment1,
                           outputDir = "~/Dropbox/lab/emergent-coexistence/data/raw/Sanger/sanger_seq_16S_20180702",
                           compress = FALSE,
                           compression_level = NA,
                           selection = "all")

alignment2 <- sangeranalyseR::SangerAlignment(ABIF_Directory = "~/Dropbox/lab/emergent-coexistence/data/raw/Sanger/sanger_seq_16S_20180726",
                                              REGEX_SuffixForward = "F.ab1$",
                                              REGEX_SuffixReverse = "R.ab1$")

sangeranalyseR::writeFasta(alignment2,
                           outputDir = "~/Dropbox/lab/emergent-coexistence/data/raw/Sanger/sanger_seq_16S_20180726",
                           compress = FALSE,
                           compression_level = NA,
                           selection = "all")

consensus1 <- seqinr::read.fasta("~/Dropbox/lab/emergent-coexistence/data/raw/Sanger/sanger_seq_16S_20180702/Sanger_contigs_alignment.fa", seqtype = "DNA")
consensus2 <- seqinr::read.fasta("~/Dropbox/lab/emergent-coexistence/data/raw/Sanger/sanger_seq_16S_20180726/Sanger_contigs_alignment.fa", seqtype = "DNA")

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
consensus2 <- consensus2 %>% lapply(function(x) {
    contig <- x %>% as.character()
    contig[contig=="a"] <- "A"
    contig[contig=="t"] <- "T"
    contig[contig=="g"] <- "G"
    contig[contig=="c"] <- "C"
    contig <- paste(contig, collapse = "")
    contig <- gsub("-", "", contig)
    return(data.frame(Sequence = contig))
})

df_seq <- data.table::rbindlist(consensus1, idcol = "temp") %>%
    dplyr::bind_rows(data.table::rbindlist(consensus2, idcol = "temp")) %>%
    tidyr::separate(col = "temp", into = c("ID", "t1", "t2", "t3")) %>%
    dplyr::select(ID, Sequence) %>%
    dplyr::mutate(ConsensusLength = nchar(Sequence))

write_csv(df_seq, file = "~/Dropbox/lab/emergent-coexistence/data/temp/isolates_16S.csv")






