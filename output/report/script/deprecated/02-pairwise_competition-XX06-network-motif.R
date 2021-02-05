#' pairwise-competition_06-network-motif.R
#' Description: this script finds the 5th and 95th percentiles and p-value for bootstrapped networks
#' Input objects: data/motif_large_communities_randomization.Rdata and data/motif_small_communities_randomization.Rdata
#' Each of the data has `df_list` that saves the randomizaed networks and `obv_list` that is the observation
#' Output: obv_large, df_percentile_large, df_p_large

source("R/network_functions.R") # Call the network functions

# Many large communities ----
community_large  <- c("C1R7", "C11R1", "C11R2")

# Find 5th and 95th percentiles
load("data/motif_large_communities_randomization.Rdata")
df_percentile_large_list <- rep(list(NA), 3)
for (i in 1:3) {
  df_percentile_large_list[[i]] <- 
    df_large_list[[i]] %>% mutate(Replicate=1:b) %>%
    gather("Motif", "Count", 1:7) %>%
    group_by(Motif) %>%
    arrange(Motif, desc(Count)) %>%
    slice(c(b*0.05, b*0.95)) %>%
    select(-Replicate) %>%
    mutate(Percentile=c("p5", "p95"), Fraction=Count/sum(obv_list[[i]]$Observation))
}

# Merge list
obv_large <- rbindlist(obv_large_list) %>% 
  mutate(Community=rep(factor(c("C1R7", "C11R1", "C11R2"), levels=c("C1R7", "C11R1", "C11R2")), each=7))
df_percentile_large <- rbindlist(df_percentile_large_list) %>% mutate(Community=rep(factor(c("C1R7", "C11R1", "C11R2"), levels=c("C1R7", "C11R1", "C11R2")), each=14))

# Find p
df_boot_large <- rbindlist(df_large_list) %>%
  mutate(Community=rep(factor(c("C1R7", "C11R1", "C11R2"), levels=c("C1R7", "C11R1", "C11R2")), each=b), Boot=rep(1:b, 3)) %>%
  group_by(Community, Boot) %>% gather("Motif", "Count", 1:7) 
prob1 <- rep(list(rep(NA, 7)), 3)
prob2 <- rep(list(rep(NA, 7)), 3)

for (i in 1:3) {
  for (j in 1:7) { 
    temp <- filter(df_boot_large, Community==community_large[i], Motif==paste0("V", j)) %>% pull(Count) %>% sort(decreasing = T)
    prob1[[i]][j] <- which(temp >= obv$Observation[obv$Motif == paste0("V", j) & obv$Community == community_large[i]]) %>% max
    prob2[[i]][j] <- which(temp <= obv$Observation[obv$Motif == paste0("V", j) & obv$Community == community_large[i]]) %>% min
  }
}

df_p_large <- data.frame(Community = rep(community_large, each = 7), Motif = paste0("V", 1:7), 
                   prob1 = (unlist(prob1))/b, prob2 = 1-(unlist(prob2))/b)
df_p_large$prob1[is.infinite(df_p_large$prob1)] <- 0; df_p_large$prob2[is.infinite(df_p_large$prob2)] <- 0
df_p_large$p_value <- round(apply(df_p_large[,c("prob1", "prob2")], 1, min), 4)
df_p_large$max_x <- df_boot %>% ungroup() %>% # y-axis coordinate on historgram
  mutate(Community = factor(Community, community_large)) %>%
  group_by(Community) %>%
  arrange(desc(Count)) %>%
  slice(1) %>% pull(Count) %>% rep(., each = 7)




# Small communities ----
community_small <- community_name[!community_name%in%community_large]

# Merge all 10 small communities into one and count motif distribution
load("data/motif_small_communities_randomization.Rdata")
df_small <- rbindlist(df_small_list) %>%
  mutate(Community=rep(factor(community_small, levels=community_small), each=b),
         Boot=rep(1:b, length(community_small))) %>%
  group_by(Boot) %>%
  gather("Motif", "Count", 1:7) %>%
  group_by(Boot, Motif) %>%
  summarize(Count=sum(Count))

obv_small <- rbindlist(obv_small_list) %>% select(-Fraction) %>%
  mutate(Community=rep(factor(community_small, levels=community_small), each=7)) %>%
  group_by(Motif) %>%
  summarise(Observation=sum(Observation)) %>%
  mutate(Fraction=Observation/sum(Observation))

# Find 5th and 95th percentiles
df_percentile_small <- df_small %>% 
  group_by(Motif) %>%
  arrange(Motif, desc(Count)) %>%
  select(-Boot) %>%
  slice(c(b*0.05, b*0.95)) %>%
  mutate(Percentile=c("p5", "p95"), Fraction=Count/sum(obv_small$Observation))










