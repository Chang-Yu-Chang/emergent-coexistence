# Generate the input_csv files

output_dir = "~/Dropbox/lab/invasion-network/simulation/data/raw2/"

# Pairs in species pool
n_treatment = 3*3*5
n_rep = 1
input_poolPairs <- tibble(
    output_dir = output_dir,
    exp_id = 1:(n_treatment*n_rep),
    seed = 1,
    q = 0.9,
    q2 = c(0.1, 0.5, 0.9) %>% rep(5) %>% rep(3) %>% rep(n_rep),
    l1 = c(0, 0.1, 0.5, 0.9, 1) %>% rep(3) %>% rep(each = 3) %>% rep(n_rep),
    l2 = 0.1,
    vamp = c(0.5, 1, 1.5) %>% rep(each = 3) %>% rep(each = 5) %>% rep(n_rep),
    S = 100,
    n_pairs = 100
)
#write_csv(input_poolPairs, file = "input_poolPairs.csv")

# Pairs of random Networks
input_randomNetworks <- tibble(
    output_dir = output_dir,
    exp_id = 1:(n_treatment*n_rep),
    seed = 1,
    q = 0.9,
    q2 = c(0.1, 0.5, 0.9) %>% rep(5) %>% rep(3) %>% rep(n_rep),
    l1 = c(0, 0.1, 0.5, 0.9, 1) %>% rep(3) %>% rep(each = 3) %>% rep(n_rep),
    l2 = 0.1,
    vamp = c(0.5, 1, 1.5) %>% rep(each = 3) %>% rep(each = 5) %>% rep(n_rep),
    S = 100
)

#input_randomNetworks <- input_randomNetworks %>% filter(l1 == 0.5, vamp == 1)
write_csv(input_randomNetworks, file = "input_randomNetworks.csv")

# Self assembly
n_treatment = 3*3*5
input_selfAssembly <- tibble(
    output_dir = output_dir,
    exp_id = 1:(n_treatment*n_rep),
    seed = 1,
    q = 0.9,
    q2 = c(0.1, 0.5, 0.9) %>% rep(5) %>% rep(3) %>% rep(n_rep),
    l1 = c(0, 0.1, 0.5, 0.9, 1) %>% rep(3) %>% rep(each = 3) %>% rep(n_rep),
    l2 = 0.1,
    vamp = c(0.5, 1, 1.5) %>% rep(each = 3) %>% rep(each = 5) %>% rep(n_rep),
    S = 100,
    n_communities = 10
)
distinct(input_selfAssembly)

write_csv(input_selfAssembly, file = "input_selfAssembly.csv")

# Community pairs
input_selfAssembly <- read_csv("input_selfAssembly.csv")
input_communityPairs <- input_selfAssembly
write_csv(input_communityPairs, file = "input_communityPairs.csv")









