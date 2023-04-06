library(tidyverse)
library(cowplot)
library(data.table)
library(readr)
source(here::here("analysis/00-metadata.R"))

#Script 1 for Extract data  from DADA2 output and rarefies to constant read depth.
# Stores into a melted data.frame with standardized columns for subsequent analysis.

rarefy <- function(dat,n =min(colSums(dat))){
    # normalize to sample count
    for(i in 1:ncol(dat)){dat
        old_column = dat[,i,with=FALSE][[1]]
        elements = factor(row.names(dat),levels=row.names(dat))
        sample = table(sample(rep(elements,old_column),n,replace = FALSE))
        dat[,(i) := as.numeric(sample)]
    }
    return(dat)
}
assign_taxonomy_id <- function(dat){
    #Asigns ID to the higest available taxonomic level
    dat$ESV_ID = dat$Genus
    dat$ESV_ID[which(is.na(dat$ESV_ID))] = dat$Family[which(is.na(dat$ESV_ID))]
    dat$ESV_ID[which(is.na(dat$ESV_ID))] = dat$Order[which(is.na(dat$ESV_ID))]
    dat$ESV_ID[which(is.na(dat$ESV_ID))] = dat$Class[which(is.na(dat$ESV_ID))]
    dat$ESV_ID[which(is.na(dat$ESV_ID))] = dat$Phylum[which(is.na(dat$ESV_ID))]
    dat$ESV_ID[which(is.na(dat$ESV_ID))] = dat$Kingdom[which(is.na(dat$ESV_ID))]
    dat$ESV_ID = as.factor(make.unique(dat$ESV_ID))
    #Make NA use ESV_ID to avoid counting them as the same group
    dat[is.na(dat$Genus)]$Genus = as.character(dat[is.na(dat$Genus)]$ESV_ID)
    dat[is.na(dat$Family)]$Family = dat[is.na(dat$Family)]$Genus
    dat[is.na(dat$Order)]$Order = dat[is.na(dat$Order)]$Family
    dat[is.na(dat$Class)]$Class = dat[is.na(dat$Class)]$Order
    dat[is.na(dat$Phylum)]$Phylum = dat[is.na(dat$Phylum)]$Class

    return(dat)
}

# 1. Rarefaction for T0-12, using the minimum number of reads 4397 ----
set.seed(1) # To ensure reproducibility
#Load Data
aux = fread(paste0(folder_data, 'raw/community_ESV/metadata.csv'))
Taxonomy_Data = fread(paste0(folder_data, 'raw/community_ESV/taxonomy.csv'))
Abundance_Data = fread(paste0(folder_data, 'raw/community_ESV/otu_table.csv')) #Data is actually at an ESV level
Abundance_Data = as.matrix(Abundance_Data[,2:ncol(Abundance_Data)])
Abundance_Data = data.table(Abundance_Data)
min(colSums(Abundance_Data))
Abundance_Data = rarefy(Abundance_Data)

Taxonomy_Data = assign_taxonomy_id(Taxonomy_Data) #Assign Taxonomy ID
Abundance_Data = t(t(Abundance_Data)/colSums(Abundance_Data)) #Calculated Relative Abundance of ESVs

#Create Sample_ID from carbon source, community, replicate and trasnfer point
cref = aux$Carbon
cref[is.na(cref)] <- 'Original' #T0 inoculum
comref = parse_number(as.character(aux$Comm))
repref = as.numeric(aux$Rep)
tref = as.numeric(aux$Transfer)
colnames(Abundance_Data)= paste(cref,comref,repref,tref,sep='.')
Abundance_Data = data.table(Abundance_Data)
Abundance_Data$ESV_ID = Taxonomy_Data$ESV_ID

#Now convert matrix into a data.frame using melt
Abundance_Data = melt(Abundance_Data,id= c('ESV_ID'),variable.name ='Sample_ID',value.name = 'Relative_Abundance')
Abundance_Data = Abundance_Data[Relative_Abundance>0,]
#Extract caronb source, community,replicate number and transfer point from sample id

Abundance_Data$Carbon_Source =sapply(Abundance_Data$Sample_ID,function(x) strsplit(as.character(x),'[.]')[[1]][1])
Abundance_Data$Inoculum =sapply(Abundance_Data$Sample_ID,function(x) strsplit(as.character(x),'[.]')[[1]][2])
Abundance_Data$Replicate =sapply(Abundance_Data$Sample_ID,function(x) strsplit(as.character(x),'[.]')[[1]][3])
Abundance_Data$Transfer =sapply(Abundance_Data$Sample_ID,function(x) strsplit(as.character(x),'[.]')[[1]][4])
Abundance_Data$ESV_ID = droplevels(Abundance_Data$ESV_ID)

#Merge with taxonomy
merged_data = merge(Taxonomy_Data,Abundance_Data)

#Save Eveything
fwrite(merged_data, paste0(folder_data, 'temp/41-communities_abundance.csv')) # This file should be identical to raw/Emergent_Comunity_Data.csv

# 2. Rarefaction for T0 only, use a sample of up to 10000 reads -----
#' For each sampling size, repeat for 100 subsamples
communities_rarefaction <- tibble(SamplingSize = c(seq(10, 100, 10), seq(200, 1000, 100), seq(2000, 10000, 1000))) %>%
    mutate(RarefiedSample = NA)

n_subsamples <- 100
rarefy_samples <- function(sampling_size) {
    aux = fread(paste0(folder_data, 'raw/community_ESV/metadata.csv'))
    Taxonomy_Data = fread(paste0(folder_data, 'raw/community_ESV/taxonomy.csv'))
    Abundance_Data = fread(paste0(folder_data, 'raw/community_ESV/otu_table.csv')) #Data is actually at an ESV level

    aux = aux %>% filter(Transfer == 0)
    Abundance_Data = Abundance_Data[,aux$ID, with = F]
    Abundance_Data = as.matrix(Abundance_Data[,1:ncol(Abundance_Data)])
    Abundance_Data = data.table(Abundance_Data)
    colSums(Abundance_Data)
    Abundance_Data = rarefy(Abundance_Data, n = sampling_size)

    Taxonomy_Data = assign_taxonomy_id(Taxonomy_Data) #Assign Taxonomy ID
    Abundance_Data = t(t(Abundance_Data)/colSums(Abundance_Data)) #Calculated Relative Abundance of ESVs

    #Create Sample_ID from carbon source, community, replicate and trasnfer point
    cref = aux$Carbon
    cref[is.na(cref)] <- 'Original' #T0 inoculum
    comref = parse_number(as.character(aux$Comm))
    repref = as.numeric(aux$Rep)
    tref = as.numeric(aux$Transfer)
    colnames(Abundance_Data)= paste(cref,comref,repref,tref,sep='.')
    Abundance_Data = data.table(Abundance_Data)
    Abundance_Data$ESV_ID = Taxonomy_Data$ESV_ID

    #Now convert matrix into a data.frame using melt
    Abundance_Data = melt(Abundance_Data,id= c('ESV_ID'),variable.name ='Sample_ID',value.name = 'Relative_Abundance')
    Abundance_Data = Abundance_Data[Relative_Abundance>0,]
    #Extract caronb source, community,replicate number and transfer point from sample id

    Abundance_Data$Carbon_Source =sapply(Abundance_Data$Sample_ID,function(x) strsplit(as.character(x),'[.]')[[1]][1])
    Abundance_Data$Inoculum =sapply(Abundance_Data$Sample_ID,function(x) strsplit(as.character(x),'[.]')[[1]][2])
    Abundance_Data$Replicate =sapply(Abundance_Data$Sample_ID,function(x) strsplit(as.character(x),'[.]')[[1]][3])
    Abundance_Data$Transfer =sapply(Abundance_Data$Sample_ID,function(x) strsplit(as.character(x),'[.]')[[1]][4])
    Abundance_Data$ESV_ID = droplevels(Abundance_Data$ESV_ID)

    #Merge with taxonomy
    merged_data = merge(Taxonomy_Data,Abundance_Data)
    return(merged_data)
}
count_richness <- function(rarfied_sample) {
    rarfied_sample %>%
        group_by(Inoculum) %>%
        count(name = "Richness")
}


for (i in 1:nrow(communities_rarefaction)) {
    list_subsamples <- rep(list(NA), n_subsamples)
    for (j in 1:n_subsamples) {
        set.seed(j) # To ensure reproducibility for each subsample
        list_subsamples[[j]] <- rarefy_samples(sampling_size = communities_rarefaction$SamplingSize[i]) %>% count_richness() %>% mutate(Subsample = j)
    }

    communities_rarefaction$RarefiedSample[i] <- list(bind_rows(list_subsamples))
    print(i)
}

communities_rarefaction <- communities_rarefaction %>%
    unnest(RarefiedSample) %>%
    group_by(SamplingSize, Inoculum) %>%
    summarize(MeanRichness = mean(Richness), SdRichness = sd(Richness))

ggsave(paste0(folder_data, "temp/41-T0_rarefaction_10000reads.png"), p, width = 4, height = 4)

fwrite(communities_rarefaction, paste0(folder_data, 'temp/41-communities_rarefaction_T0.csv'))

# Save the rarefied data for 10000 reads
set.seed(1)
communities_abundance_T0_10000_reads <- rarefy_samples(10000)
fwrite(communities_abundance_T0_10000_reads, paste0(folder_data, 'temp/41-communities_abundance_T0_10000_reads.csv'))

# T0 samples
aux = fread(paste0(folder_data, 'raw/community_ESV/metadata.csv'))
Taxonomy_Data = fread(paste0(folder_data, 'raw/community_ESV/taxonomy.csv'))
Abundance_Data = fread(paste0(folder_data, 'raw/community_ESV/otu_table.csv')) #Data is actually at an ESV level

aux = aux %>% filter(Transfer == 0)
Abundance_Data = Abundance_Data[,aux$ID, with = F]
Abundance_Data = as.matrix(Abundance_Data[,1:ncol(Abundance_Data)])
Abundance_Data = data.table(Abundance_Data)
range(colSums(Abundance_Data)) # Number of reads
apply(Abundance_Data, 2, function(x) x>0) %>% colSums() # number of unique ESVs













