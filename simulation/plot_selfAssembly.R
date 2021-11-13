library(tidyverse)
source("misc.R")

# Generate input csv
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

#write_csv(input_selfAssembly, file = "input_selfAssembly.csv")



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
df_end <- read_simulation_data(input_selfAssembly$output_dir[1], paste0(pattern_exp_id, "_end")) %>% mutate(Time = 1)

df <- df_end %>% left_join(select(input_selfAssembly, exp_id, seed, vamp, q2, l1))

# Check richness
df %>%
    group_by(Experiment, exp_id, seed, vamp, q2, l1, Well) %>%
    mutate(l1 = factor(l1), q2 = factor(q2)) %>%
    filter(l1 == 0.5) %>%
    summarize(Richness = n()) %>%
    ggplot() +
    geom_jitter(aes(x = q2, y = Richness, color = l1), shape = 1, height = 0, width = 0.1, size = 3) +
    facet_grid(.~vamp) +
    theme_bw()

# Community
"
plot and find the decent q and q2 (fixed) and a range of vamp
"
"
I am mixing the cross-community and within community
"
l = c(0, 0.1, 0.5, 0.9, 1)
for (i in 1:5) {
    p <- df %>%
        filter(l1 == l[i]) %>%
        ggplot() +
        geom_col(aes(x = Well, y = Abundance, fill = Family, alpha = Species), color = 1, position = "fill") +
        facet_grid(vamp~q2, labeller = label_both) +
        theme_classic() +
        theme(legend.position = "top") +
        guides(alpha = "none") +
        ggtitle(paste0("l1 = ", l[i]))
    ggsave(paste0("plots/selfAssembly_seed1_l1_", i, ".png"), p, width = 10, height = 10)
}
p
#ggsave("plots/selfAssembly_seed1.png", p, width = 10, height = 10)










