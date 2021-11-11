rm(list=ls())
library(data.table)
library(ggplot2)
library(ggpubr)
#Fig 1
q = c(0.1,0.9)
q2 = c(0.1,0.9)
seed = seq(1,1)
mdf=expand.grid(q,q2,seed)
colnames(mdf) = c('q','q2','seed')
mdf$exp_id = seq(1:nrow(mdf))
fwrite(mdf,'Experiment_list.csv')
mdf$P
mdf$N
mdf$P_FF =0.0
mdf$P_FR =0.0
mdf$P_RR =0.0
mdf$N_FF =0.0
mdf$N_FR =0.0
mdf$N_RR =0.0
for(i in 1:nrow(mdf)){
  pairs = fread(paste('../Data/Raw/Random_Pairs_',mdf$exp_id[i],'.csv',sep=''))
  pairs = pairs[pairs$Species_1!=pairs$species_2]
  FF = pairs[pairs$Family_1=='F0' &pairs$Family_2 =='F0']
  FR = pairs[pairs$Family_1=='F0' &pairs$Family_2 =='F1']
  RR = pairs[pairs$Family_1=='F1' &pairs$Family_2 =='F1']
  mdf$N[i]= nrow(pairs)
  mdf$P[i] = sum(pairs$Abundance_1 ==0 |pairs$Abundance_2==0.0)/nrow(pairs)
  mdf$N_FF[i]= nrow(FF)
  mdf$P_FF[i] = sum(FF$Abundance_1 ==0 |FF$Abundance_2==0.0)/nrow(FF)
  mdf$N_RR[i]= nrow(RR)
  mdf$P_RR[i] = sum(RR$Abundance_1 ==0 |RR$Abundance_2==0.0)/nrow(RR)
  mdf$N_FR[i]= nrow(FR)
  mdf$P_FR[i] = sum(FR$Abundance_1 ==0.0 |FR$Abundance_2==0.0)/nrow(FR)
}


p1 <-ggplot(mdf) + geom_tile(aes(x=q,y=q2,fill=P)) +ggtitle('All') +
  scale_fill_continuous(limits=c(0,0.5))
p2 <-ggplot(mdf) + geom_tile(aes(x=q,y=q2,fill=P_FF)) +ggtitle('FF') +
  scale_fill_continuous(limits=c(0,0.5))
p3 <-ggplot(mdf) + geom_tile(aes(x=q,y=q2,fill=P_RR)) + ggtitle('RR') +
  scale_fill_continuous(limits=c(0,0.5))
p4 <-ggplot(mdf) + geom_tile(aes(x=q,y=q2,fill=P_FR)) +ggtitle('RF') +
  scale_fill_continuous(limits=c(0,0.5))


ggarrange(p1,p2,p3,p4)
