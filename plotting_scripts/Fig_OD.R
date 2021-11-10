library(tidyverse)
library(cowplot)

OD <- read_csv(here::here("data/temp/OD.csv"))
communities <- read_csv(here::here("data/output/communities.csv"))

index <- which(OD$Isolate1 > OD$Isolate2)
temp <- OD$Isolate1[index]
OD$Isolate1[index] <- OD$Isolate2[index]
OD$Isolate2[index] <- temp
temp <- OD$Isolate1Freq[index]
OD$Isolate1Freq[index] <- OD$Isolate2Freq[index]
OD$Isolate2Freq[index] <- temp


p_list = rep(list(NA), nrow(communities))
for (i in 1:nrow(communities)) {
    p_list[[i]] <- OD %>%
        filter(Wavelength == 620, MixIsolate == T) %>%
        filter(Community == communities$Community[i]) %>%
        group_by(Community, Isolate1, Isolate2, Isolate1Freq, Transfer) %>%
        summarize(Abs_mean = mean(Abs), Abs_sd = sd(Abs)) %>%
        ggplot(aes(x = Transfer, y = Abs_mean, color = Isolate1Freq)) +
        geom_point() + geom_line() +
        geom_segment(aes(x = Transfer, xend = Transfer, y = Abs_mean - Abs_sd, yend = Abs_mean + Abs_sd)) +
        facet_grid(Isolate1~Isolate2, scales = "free_y") +
        theme_bw() +
        ggtitle(communities$Community[i])
    ggsave(here::here(paste0("plots/Fig_OD_", communities$Community[i], ".png")), plot = p_list[[i]], width = 12, height = 10)
}


