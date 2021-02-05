# Read and align raw Sanger sequences from Genewiz

library(sangeranalyseR)

#
alignment1 <- sangeranalyseR::SangerAlignment(parentDirectory = here::here("data/raw/Sanger/sanger_seq_16S_20180702"),
                              suffixForwardRegExp = "F.ab1$",
                              suffixReverseRegExp = "R.ab1$")
sangeranalyseR::writeFasta(alignment1,
           outputDir = here::here("data/raw/Sanger/sanger_seq_16S_20180702"),
           compress = FALSE,
           compression_level = NA,
           selection = "all")

alignment2 <- sangeranalyseR::SangerAlignment(parentDirectory = here::here("data/raw/Sanger/sanger_seq_16S_20180726"),
                              suffixForwardRegExp = "F.ab1$",
                              suffixReverseRegExp = "R.ab1$")

sangeranalyseR::writeFasta(alignment2,
           outputDir = here::here("data/raw/Sanger/sanger_seq_16S_20180726"),
           compress = FALSE,
           compression_level = NA,
           selection = "all")

consensus1 <- seqinr::read.fasta(here::here("data/raw/Sanger/sanger_seq_16S_20180702/Sanger_contigs_alignment.fa"), seqtype = "DNA")
consensus2 <- seqinr::read.fasta(here::here("data/raw/Sanger/sanger_seq_16S_20180726/Sanger_contigs_alignment.fa"), seqtype = "DNA")

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

df_seq <-
  data.table::rbindlist(consensus1, idcol = "temp") %>%
  dplyr::bind_rows(data.table::rbindlist(consensus2, idcol = "temp")) %>%
  tidyr::separate(col = "temp", into = c("ID", "t1", "t2", "t3")) %>%
  dplyr::select(ID, Sequence) %>%
  dplyr::mutate(ConsensusLength = nchar(Sequence))

data.table::fwrite(df_seq, file = here::here("data/temp/", "isolates_16S.csv"))

if (FALSE) {

  # Create the list of abi file names ----
  f <- c(list.files(here::here("data/raw/Sanger/sanger_seq_16S_20180702")),
         list.files(here::here("data/raw/Sanger/sanger_seq_16S_20180726")))

  fwd_list <- grep("F.ab1", f, value = T)
  rev_list <- grep("R.ab1", f, value = T)

  # Merge forward and backard sequences from each isolates. Loop through all isolates ----
  rss <- rep(list(NA), length(rev_list)) # Readsets
  mrs <- rep(list(NA), length(rev_list))  # Merged reads
  length(rev_list)

  for (i in 1:length(rev_list)) {
    iso_ID <- as.numeric(substr(fwd_list[i], 1, 3))

    if (iso_ID <= 369) {
      rss[[i]] <-
        sangeranalyseR::make.readset(
          split(paste0(root$find_file("data/raw/Sanger/sanger_seq_16S_20180702/"), fwd_list[i]), 1),
          split(paste0(root$find_file("data/raw/Sanger/sanger_seq_16S_20180702/"), rev_list[i]), 1))
    } else if (iso_ID >= 370) {
      rss[[i]] <-
        sangeranalyseR::make.readset(
          split(paste0(root$find_file("data/raw/Sanger/sanger_seq_16S_20180726/"), fwd_list[i]), 1),
          split(paste0(root$find_file("data/raw/Sanger/sanger_seq_16S_20180726/"), rev_list[i]), 1))
    }

    # Merged readset
    mrs[[i]] <- merge.reads(rss[[i]]$readset)

    cat(i)
  }


  # Extract only consensus sequences ----
  seqs <- rep(NA, length(mrs))

  for (i in 1:length(mrs)) {
    if (!is.na(mrs[[i]][1])) {
      seqs[i] <- mrs[[i]]$consensus %>% as.character()
    } else next
  }

  # Name ----
  names_seq <- substr(fwd_list, 1, 3)
  names(rss) <- names_seq # Readsets
  names(mrs) <- names_seq # Sequences
  names(seqs) <- names_seq # All sequences

  # Save output file
  save(rss, mrs, seqs, file = root$find_file("data/temp/isolates_sanger_readsets_201807.Rdata"))
}





