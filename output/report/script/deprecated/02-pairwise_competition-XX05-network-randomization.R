#' pairwise-competition_05-network-randomization.R
#' Description: this script randomize all networks
#' Input objects: community_name, community_size
#' Output data: data/motif_large_communities_randomization.Rdata and data/motif_small_communities_randomization.Rdata
#' Each of the data has `df_XX_list` that saves the randomizaed networks and `obv_XX_list` that is the observation

b=100 # number of bootstrapping
source("R/network_functions.R") # Call the network functions

# Single large community----
# Randomization
df <- as.data.frame(matrix(NA, nrow=b, ncol=7))
comm = "C11R1"; comm_size <- community_size[community_name==comm]
tour_list <- tournament(pairs = pairs3, comm, comm_size)
net <- network_make(isolates2, pairs3, comm)

for (i in 1:b) {
  temp <- network_randomize(net, comm_size, method="shuffle_interaction") # or switch_sign
  net_rand <- temp$net_rand
  isolate_rank <- temp$isolate_rank
  df[i, ] <- motif_detect(net_rand)
  cat(i)
}


# Find 5th and 95th percentiles
obv <- data.frame(Motif=paste0("V", 1:7), Observation=motif_detect(net)) %>%
  mutate(Fraction=Observation/sum(Observation))
df_percentile <- 
  df %>% mutate(Replicate=1:b) %>%
  gather("Motif", "Count", 1:7) %>%
  group_by(Motif) %>%
  arrange(Motif, desc(Count)) %>%
  slice(c(b*0.05, b*0.95)) %>%
  select(-Replicate) %>%
  mutate(Percentile=c("p5", "p95"),
         Fraction=Count/sum(obv$Observation))


# Many large communities ----

community_large <- c("C1R7", "C11R1", "C11R2")
# Randomization. This step may take a few minutes for 1000 bootstrapping
df_large_list <- rep(list(as.data.frame(matrix(NA, nrow=b, ncol=7))), 3)
obv_large_list <- rep(list(data.frame(Motif=paste0("V", 1:7), Observation=rep(NA, 7))), 3) 
j=1
for (comm in community_large) {
#  comm<-"C1R7"
  # Make network
  comm_size <- community_size[community_name==comm]
  tour_list <- tournament(pairs = pairs3, comm, comm_size)
  net <- network_make(isolates2, pairs3, comm)
  
  # Observation
  obv_large_list[[j]]$Observation <- motif_detect(net)
  obv_large_list[[j]]$Fraction <- obv_large_list[[j]]$Observation/sum(obv_large_list[[j]]$Observation)
  
  # Randomization
  for (i in 1:b) {
    temp <- network_randomize(net, comm_size, method="shuffle_interaction") # or switch_sign
    net_rand <- temp$net_rand
    isolate_rank <- temp$isolate_rank
    df_large_list[[j]][i, ] <- motif_detect(net_rand)
    if(i%%100==0) cat(i); cat(" ")
  }
  j=j+1
  print(comm)
}

if (save_data) save(df_large_list, obv_large_list, b, file="data/motif_large_communities_randomization.Rdata")

# Small communities ----
# Randomization. This step may take a few minutes for 1000 bootstrapping
community_small <- community_name[!community_name%in%community_large]
df_small_list <- rep(list(as.data.frame(matrix(NA, nrow=b, ncol=7))), length(community_small))
obv_small_list <- rep(list(data.frame(Motif=paste0("V", 1:7), Observation=rep(NA, 7))), length(community_small)) 
j=1

for (comm in community_small) {
  # Make network
  comm_size <- community_size[community_name==comm]
  tour_list <- tournament(pairs = pairs3, comm, comm_size)
  net <- network_make(isolates2, pairs3, comm)
  
  # Observation
  obv_small_list[[j]]$Observation <- motif_detect(net)
  obv_small_list[[j]]$Fraction <- obv_small_list[[j]]$Observation/sum(obv_small_list[[j]]$Observation)
  
  # Randomization
  for (i in 1:b) {
    temp <- network_randomize(net, comm_size, method="shuffle_interaction") # or switch_sign
    net_rand <- temp$net_rand
    isolate_rank <- temp$isolate_rank
    df_small_list[[j]][i, ] <- motif_detect(net_rand)
    if(i%%100==0) cat(i)
  }
  j=j+1
  print(comm)
}

if (save_data) save(df_small_list, obv_small_list, b, file="data/motif_small_communities_randomization.Rdata")