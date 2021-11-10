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
write_csv(input_selfAssembly, file = "input_selfAssembly.csv")


# Read data
input_selfAssembly <- read_csv("input_selfAssembly.csv")
df_init <- read_simulation_data("~/Dropbox/lab/invasion-network/simulation/data/raw", "selfAssembly_\\d+_init") %>% mutate(Time = 0)
df_end <- read_simulation_data("~/Dropbox/lab/invasion-network/simulation/data/raw", "selfAssembly_\\d+_end") %>% mutate(Time = 1)
# df_init <- read_simulation_data("~/Dropbox/lab/invasion-network/simulation/data/raw", "selfAssembly_9_init") %>% mutate(Time = 0) %>%
#     arrange(exp_id, Well)
# df_end <- read_simulation_data("~/Dropbox/lab/invasion-network/simulation/data/raw", "selfAssembly_9_end") %>% mutate(Time = 1) %>%
#     arrange(exp_id, Well)

df <- df_end %>%
    left_join(select(input_selfAssembly, exp_id, seed, vamp, q2))

# Community
"
plot and find the decent q and q2 (fixed) and a range of vamp
"
p <- df %>%
    filter(seed == 1) %>%
    ggplot() +
    geom_col(aes(x = Well, y = Abundance, fill = Family, alpha = Species), color = 1, position = "fill") +
    facet_grid(vamp~q2) +
    theme_classic() +
    guides(alpha = "none")
p
ggsave("plots/selfAssembly_seed1.png", p, width = 10, height = 10)

#
df %>%
    group_by(exp_id, vamp, q2, Well) %>%
    summarize(Richness = n()) %>%
    group_by(exp_id, vamp, q2) %>%
    summarize(MeanRichness = mean(Richness)) %>%
    view





if (FALSE) {
read_csv("data/raw/1_end.csv", col_types = cols()) %>%
    select(Species = `...2`, starts_with("W")) %>%
    na_if(0) %>%
    pivot_longer(cols = -Species, names_to = "Well", values_to = "Abundance", values_drop_na = T)

}

