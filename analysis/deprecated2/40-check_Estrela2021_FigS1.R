" R script for Figure S1 in Estrela et al (2021)
Functional attractors in microbial community assembly."

#------------------------------------------------------------------------------
#.. Clear workspace
rm(list=ls())

#.. Load packages
library(data.table)
library(ggplot2)
library(dplyr)
library(gridExtra)
library(ggpubr)
library(reshape2)
source(here::here("analysis/00-metadata.R"))
communities <- read_csv(paste0(folder_data, "temp/00c-communities.csv"), col_types = cols())
#------------------------------------------------------------------------------
#.. Assign path to data
#. usern = path to user directory

# path.d<-paste0(usern, "/data/data_main")
# #.. Assign path to plots output
# path.p<- paste0(usern,"/Rplots")
# #.. Assign path to data figures
# path.ds<- paste0(usern,"/data/data_figures")
# #.. Assign path to data processed
# path.dp<- paste0(usern,"/data/data_processed")
#
# #.. Read data
# setwd(path.dp)

# ..Read data table
sequences_abundance <- fread(paste0(folder_data, "raw/test/figS1_sequences_abundance.csv"))
communities_name <- c(paste0("C", 1:12, "Rpool"), paste0("C", rep(1:12, each = 8), "R", rep(1:8, 12)))

# Family colours
colAerof<-"deepskyblue4"
colAlcf<- "orange2"
colComf<- "firebrick"
colEntbf<- "deepskyblue3"
colEntcf<-  "darkolivegreen3"
colMoraxf<- "darkorchid4"
colPseuf<-  "darkorchid2"
colSphingf<-"deeppink3"
colXantf<-"deeppink4"



family_name <- c("Pseudomonadaceae", "Enterobacteriaceae", "Aeromonadaceae",  "Sphingobacteriaceae", "Xanthomonadaceae", "Moraxellaceae", "Alcaligenaceae", "Comamonadaceae", "Other")
family_color <- c(colPseuf, colEntbf, colAerof, colSphingf, colXantf, colMoraxf, colAlcf, colComf, "#989898")
names(family_color) <- family_name

plot_abundance <- function(df, label_x = "Community", label_y = "RelativeAbundance", fill = "CommunityESVID") {
    ggplot(df) +
        geom_bar(aes_string(x = label_x, y = label_y, fill = fill),
                 position = "stack", stat = "identity", col = 1) +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 90))
}

p <- sequences_abundance %>%
    filter(AlignmentType == "local") %>%
    filter(AllowMismatch == Inf) %>%
    filter(BasePairMismatch <= 4) %>%
    mutate(Community = ordered(Community,  communities_name)) %>%
    plot_abundance(label_x = "Community", label_y = " RelativeAbundance", fill = "Family") +
    labs(x = "community", y = "relative abundance") +
    scale_fill_manual(values = family_color) +
    scale_x_discrete(expand=c(0,0)) +
    scale_y_continuous(expand=c(0,0), limits = c(0,1), breaks = seq(0,1, .25)) +
    theme_classic()+
    theme(
        axis.text.x = element_text(size=14, angle=90, vjust=0.5),
        axis.text.y = element_text(size=14),
        axis.title.y =element_text(size=20),
        axis.title.x = element_text(size=20),
        legend.text = element_text(size=12),
        legend.title = element_text(size=12),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.line = element_line(colour = "gray"),
        panel.border = element_rect(colour = "gray", fill=NA, size=0.5))
p
#ggsave(path=path.p, paste0("figS1.pdf"), width= 10, height=4.2, onefile = TRUE)

ds<-sequences_abundance %>%
    filter(AlignmentType == "local") %>%
    filter(AllowMismatch == Inf) %>%
    filter(BasePairMismatch <= 4)
ds[Community == "C4R1"]

dat.figS1<-select(ds, Community, Family,CommunityESVID, RelativeAbundance)
# ## ..Write data table as csv file
# fname<-paste("figS1.csv")
# write.csv(dat.figS1,file.path(path.ds, fname),row.names=FALSE)

# Averaged coverage
sequences_abundance %>%
    filter(AlignmentType == "local") %>%
    filter(AllowMismatch == Inf) %>%
    filter(BasePairMismatch <= 4) %>%
    mutate(Community = ordered(Community,  communities_name)) %>%
    dplyr::group_by(Community) %>%
    #filter(Community %in% communities$Community) %>%
    dplyr::summarize(TotalRelativeAbundance = sum(RelativeAbundance)) %>%
    pull(TotalRelativeAbundance) %>% mean
