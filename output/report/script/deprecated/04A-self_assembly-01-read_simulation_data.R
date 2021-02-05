#' Read the simulated result
#' 1. Self-assembly
#' 2. Monocultures in the speceis pool
#' 3. Pairs of self-assembled communities
#' 4. Pairs of randomly assembled communities

# Read data ----
file_list <- list.files(root$find_file("output/report/data/"))

## Self assembly
df_self_assembly <-
  grep(file_list, pattern = "self_assembly-community", value = T) %>%
  paste0(root$find_file("output/report/data/"), .) %>%
  as.list() %>% lapply(fread) %>% rbindlist() %>%
  mutate(ID = factor(ID)) %>%
  # Community defined by well
  mutate(Community = Well) %>%
  as_tibble()

## Monocultures
df_mono <-
  grep(file_list, pattern = "self_assembly-mono", value = T) %>%
  paste0(root$find_file("output/report/data/"), .) %>%
  as.list() %>% lapply(fread) %>% rbindlist() %>%
  mutate(Community = Well) %>%
  mutate(ID = factor(ID)) %>%
  as_tibble()

## Growers on the minimal media
df_grower <- df_mono %>% filter(Transfer == 10, Type == "consumer")

## Pairs
file_list_pair <- grep(file_list, pattern = "pair", value = T) %>%paste0(root$find_file("output/report/data/"), .)
file_list_pair_name <- gsub(root$find_file("output/report/data/"), "", file_list_pair) %>%
  gsub("pair-", "", .) %>%
  gsub("^\\w+-", "", .) %>%
  gsub("-T\\d+.txt", "", .)

df_pair <- file_list_pair %>%
  as.list() %>% setNames(file_list_pair_name) %>%
  lapply(fread) %>%
  # Name by the community
  rbindlist(idcol = "Community") %>%
  select(Assembly, Community, Well, Transfer, Type, ID, Abundance) %>%
  mutate(ID = factor(ID)) %>%
  as_tibble()















# Plot self-assembled communities ----
p_self_assembly <- rep(list(NA), 4)

## Richness
p_self_assembly[[1]] <- df_self_assembly %>%
  group_by(Transfer, Community, Type) %>%
  filter(Type == "consumer") %>%
  summarize(Richness = n()) %>%
  ggplot(aes(x = Transfer, y = Richness, color = Community)) +
  geom_point() + geom_line() +
  scale_x_continuous(breaks = seq(0, 10, 1)) +
  scale_y_continuous(limits = c(0,100)) +
  guides(color = F) +
  theme_bw()

## Absolute biomass in dynamics
p_self_assembly[[2]] <- df_self_assembly %>%
  filter(Community %in% paste0("W", 0:3)) %>%
  ggplot(aes(x = Transfer, y = Abundance, color = ID, group = ID)) +
  geom_line() +
  scale_x_continuous(breaks = seq(0, 10, 1)) +
  facet_grid(Type ~ Community, scale = "free_y") +
  theme_bw() +
  guides(color = F)

## Relative abundance in dynamics
p_self_assembly[[3]] <- df_self_assembly %>%
  filter(Community %in% paste0("W", 0:3)) %>%
  group_by(Transfer, Community, Type) %>%
  mutate(RelativeAbundance = Abundance / sum(Abundance)) %>%
  ggplot() +
  geom_bar(aes(x = Transfer, y = RelativeAbundance, fill = ID, group = ID),
           color = "black", lwd = 0.1, stat = "identity") +
  scale_x_continuous(breaks = seq(0, 10, 1)) +
  facet_grid(Type ~ Community, scale = "free_y") +
  theme_bw() +
  guides(fill = F, color = F)


## Composition at the final transfer
p_self_assembly[[4]] <- df_self_assembly %>%
  filter(Transfer == 10) %>%
  group_by(Transfer, Community, Type) %>%
  mutate(RelativeAbundace = Abundance / sum(Abundance)) %>%
  ggplot() +
  geom_bar(aes(x = Community, y = RelativeAbundace, fill = ID),
           color = "black", lwd = .3, stat = "identity") +
  facet_grid(Type ~., scale = "free") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90)) +
  guides(fill = F, color = F)



# Plot the monocultures ----
p_mono <- rep(list(NA), 2)

## Richness
p_mono[[1]] <- df_mono %>%
  group_by(Transfer, Community, Type) %>%
  filter(Type == "consumer") %>%
  summarize(Richness = n()) %>%
  ggplot(aes(x = Transfer, y = Richness, color = Community)) +
  geom_point() + geom_line() +
  scale_x_continuous(breaks = seq(0, 10, 1)) +
  scale_y_continuous(limits = c(0,1)) +
  guides(color = F) +
  theme_bw()


## Absolute biomass in dynamics
p_mono[[2]] <- df_mono %>%
  filter(Type == "consumer") %>%
  mutate(Community = ordered(Community, paste0("W", 0:209))) %>%
#  filter(Community %in% paste0("W", 0:3)) %>%
  ggplot(aes(x = Transfer, y = Abundance, color = ID, group = ID)) +
  geom_point() + geom_line() +
  scale_x_continuous(breaks = seq(0, 10, 1)) +
#  facet_wrap(Community~.) +
  theme_bw() +
  guides(color = F)


# Compare monoculture to community
grower_list <- df_mono %>%
  filter(Transfer == 10, Type == "consumer") %>%
  pull(ID) %>% as.character()


df_self_assembly %>%
  filter(Transfer == 10, Type == "consumer") %>%
  group_by(Community) %>%
  summarize(sum(!ID %in% grower_list))











