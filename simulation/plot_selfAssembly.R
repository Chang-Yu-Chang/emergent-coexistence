library(tidyverse)
source("misc.R")

# Generate input csv
n_treatment = 3*4
n_rep = 20
input_selfAssembly <- tibble(
    output_dir = "~/Dropbox/lab/invasion-network/simulation/data/raw/",
    exp_id = 1:(n_rep*n_treatment),
    seed = rep(1:n_rep, each = n_treatment),
    q = 0.9,
    q2 = c(0.1, 0.5, 0.9) %>% rep(4) %>% rep(n_rep),
    vamp = c(0.5, 1, 1.5, 2) %>% rep(each = 3) %>% rep(n_rep),
    S = 100,
    n_communities = 10
)
#write_csv(input_selfAssembly, file = "input_selfAssembly.csv")

n_treatment = 3*3*5
n_rep = 1
input_selfAssembly <- tibble(
    output_dir = "~/Dropbox/lab/invasion-network/simulation/data/raw2/",
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



# Read data
input_selfAssembly <- read_csv("input_selfAssembly.csv")
target_exp_id <- input_selfAssembly %>%
    filter(seed == 1) %>%
    filter(q2 %in% c(0.1, 0.5, 0.9), vamp %in% c(0.5, 1, 1.5)) %>%
    pull(exp_id)
pattern_exp_id <- paste(target_exp_id, collapse = "|") %>% paste0("selfAssembly_(", . , ")")
#pattern_exp_id <- "selfAssembly_\\d+"

# Subset the target exp_id
#df_init <- read_simulation_data("~/Dropbox/lab/invasion-network/simulation/data/raw", paste0(pattern_exp_id, "_init"), cs_output_format = F) %>% mutate(Time = 0)
df_end <- read_simulation_data("~/Dropbox/lab/invasion-network/simulation/data/raw", paste0(pattern_exp_id, "_end")) %>% mutate(Time = 1)

df <- df_end %>%
    left_join(select(input_selfAssembly, exp_id, seed, vamp, q2, l1))

# Community
"
plot and find the decent q and q2 (fixed) and a range of vamp
"
"
I am mixing the cross-community and within community
"

for (i in 1:3) {
    p <- df %>%
        filter(l1 == c(0.1, 0.5, 0.9)[i]) %>%
        ggplot() +
        geom_col(aes(x = Well, y = Abundance, fill = Family, alpha = Species), color = 1, position = "fill") +
        facet_grid(vamp~q2, labeller = label_both) +
        theme_classic() +
        theme(legend.position = "top") +
        guides(alpha = "none") +
        ggtitle(paste0("l1 = ", c(0.1,0.5,0.9)[i]))
    ggsave(paste0("plots/selfAssembly_seed1_l1_", i, ".png"), p, width = 10, height = 10)
}
p
#ggsave("plots/selfAssembly_seed1.png", p, width = 10, height = 10)

df_interaction %>%
    filter(l1 == 0.5, q2 == 0.1, vamp == 1) %>%
    view

















