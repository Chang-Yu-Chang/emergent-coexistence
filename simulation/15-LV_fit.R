# Fitting the LV model using a specific pair

library(tidyverse)
library(cowplot)
library(gauseR) # For fitting LV model
library(broom)
source(here::here("analysis/00-metadata.R"))
source(here::here("simulation/01-generate_input.R"))


# 0. parameters ----
input_parameters <- read_csv(here::here("simulation/01-input_parameters.csv"), col_types = cols())
# input_poolPairs <- read_csv(here::here("simulation/03a-input_poolPairs.csv"), col_types = cols())
# input_withinCommunityPairs <- read_csv(here::here("simulation/03b-input_withinCommunityPairs.csv"), col_types = cols())
input_LVPairs <- read_csv(here::here("simulation/03c-input_LVPairs.csv"), col_types = cols())

read_wide_file <- function(x, type = "N") {
    temp <- read_csv(x, col_types = cols(), name_repair = "unique_quiet") %>%
        pivot_longer(cols = starts_with("W"), names_to = "Well", values_to = "Abundance")
    if ("...1" %in% colnames(temp)) {
        if (type == "N") temp <- temp %>% rename(Family = ...1, Species = ...2)
        if (type == "R") temp <- temp %>% rename(Class = ...1, Resource = ...2)
    }

    return(temp)
}


# Read data ----
pairs_abundance <- list.files(input_LVPairs$output_dir[1], pattern = "LVPairs_W0-1-N") %>%
    map(c("N_T\\d+t\\d+\\.csv", "init", "end"), str_subset, string = .) %>%
    reduce(union) %>%
    paste0(input_LVPairs$output_dir[1], .) %>%
    lapply(function(x) {
        time_point <- str_split_i(x, "-", -1) %>%
            str_replace("N_", "") %>%
            str_replace(".csv", "")
        read_wide_file(x) %>%
            mutate(Time = time_point) %>%
            return()
    }) %>%
    bind_rows


pairs_abundance_focal <- pairs_abundance %>%
    filter(Well == "W9")
pairs_abundance_focal %>%
    filter(Time == "end", Abundance > 0)
pairs_abundance_focal_sp <- pairs_abundance_focal %>%
    filter(Time %in% c("end", "init")) %>%
    filter(Time == "init", Abundance > 0)

pairs_abundance_focal_subset <- pairs_abundance_focal %>%
    filter(Species %in% pairs_abundance_focal_sp$Species) %>%
    filter(!(Time %in% c("end", "init"))) %>%
    mutate(Time = factor(Time, paste0("T", rep(1:10, each = 50), "t", rep(1:50, 10)))) %>%
    #separate(col = Time, sep = "t", into = c("Transfer", "Timepoint")) %>%
    #mutate(Timepoint = paste0("t", Timepoint)) %>%
    select(-Well) %>%
    arrange(Time)



# Clean up data
pairs_LV <- pairs_abundance_focal_subset %>%
    #filter(Time %in% paste0("T", rep(1, 50), "t", 1:50)) %>%
    separate(col = Time, sep = "t", into = c("Transfer", "Timepoint"), remove = F) %>%
    mutate(Transfer = as.numeric(str_replace(Transfer, "T", "")), Timepoint = as.numeric(Timepoint)) %>%
    mutate(TimeID = (Transfer-1) * 50 + Timepoint) %>%
    select(-Family) %>%
    mutate(Time = as.numeric(Time)) %>%
    mutate(SpeciesID = paste0("Species", rep(1:2, n()/2)))

# Interaction coefficient over time
p1 <- pairs_LV %>%
    ggplot(aes(x = TimeID, y = Abundance, color = Species, group = Species)) +
    geom_point() +
    geom_line() +
    #scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, 10)) +
    facet_grid(.~Transfer, scales = "free_x", labeller = labeller(Transfer = label_both)) +
    theme_classic() +
    theme(panel.border = element_rect(color = 1, fill = NA)) +
    labs()


# Time window
windows <- tibble(Window = c(1:46, 51:96), Transfer = rep(1:2, each = 46)) %>%
    group_by(Window) %>%
    group_split() %>%
    lapply(function(x) tibble(Window = x$Window[1], Transfer = x$Transfer[1], TimeID = (x$Window[1]):(x$Window[1]+4))) %>%
    bind_rows() %>%
    mutate(Window = as.numeric(Window))


calculate_coeffs <- function (pairs_focal) {
    # get time-lagged observations for each species
    species1_lagged <- get_lag(x = pairs_focal$Species1, time = pairs_focal$Time)
    species2_lagged <- get_lag(x = pairs_focal$Species2, time = pairs_focal$Time)

    # calculate per-capita growth rates
    species1_dNNdt <- percap_growth(x = species1_lagged$x, laggedx = species1_lagged$laggedx, dt = species1_lagged$dt)
    species2_dNNdt <- percap_growth(x = species2_lagged$x, laggedx = species2_lagged$laggedx, dt = species2_lagged$dt)


    # fit linear models to dNNdt, based on average
    # abundances between current and lagged time steps
    speceis1_mod_dat <- tibble(species1_dNNdt = species1_dNNdt,
                               species1 = species1_lagged$laggedx,
                               species2 = species2_lagged$laggedx)
    mod_speceis1 <- lm(species1_dNNdt ~ species1 + species2, data=speceis1_mod_dat)

    speceis2_mod_dat <- tibble(species2_dNNdt = species2_dNNdt,
                               species1 = species1_lagged$laggedx,
                               species2 = species2_lagged$laggedx)
    mod_speceis2 <- lm(species2_dNNdt ~ species1 + species2, data=speceis2_mod_dat)

    # Plot
    bind_rows(
        tidy(mod_speceis1) %>% mutate(Term = c("r1", "a11", "a12")),
        tidy(mod_speceis2) %>% mutate(Term = c("r2", "a21", "a22"))
    ) %>%
        select(Term, Estimate = estimate) %>%
        return()
}

LVfit <- rep(list(NA), length(unique(windows$Window)))

for (i in 1:length(unique(windows$Window))) {
    windows_i <- filter(windows, Window == unique(windows$Window)[i])
    LVfit[[i]] <- pairs_LV %>%
        pivot_wider(id_cols = -Species, names_from = SpeciesID, values_from = Abundance) %>%
        filter(TimeID %in% windows_i$TimeID) %>%
        calculate_coeffs() %>%
        mutate(Window = unique(windows$Window)[i])
}

LVfit <- LVfit %>%
    bind_rows() %>%
    mutate(Window = as.numeric(Window)) %>%
    left_join(distinct(windows, Window, Transfer))

# Interaction coefficients
p2 <- LVfit %>%
    filter(!(Term %in% c("r1", "r2"))) %>%
    ggplot(aes(x = Window, y = Estimate, color = Term, group = Term))+
    geom_point() +
    geom_line() +
    theme_classic() +
    #scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, 10)) +
    facet_grid(.~Transfer, scales = "free_x", labeller = labeller(Transfer = label_both)) +
    theme(panel.border = element_rect(color = 1, fill = NA)) +
    labs()


p <- plot_grid(p1, p2, nrow = 2, axis = "lr", align = "v")

ggsave(here::here("simulation/plots/15-LV_fit_time.png"), plot = p, width = 8, height = 8)


# vector plot
p <- LVfit %>%
    filter(!(Term %in% c("r1", "r2"))) %>%
    group_by(Window) %>%
    arrange(Window, Term) %>%
    pivot_wider(names_from = Term, values_from = Estimate) %>%
    ggplot(aes(x = a12, y = a21)) +
    #ggplot(aes(x = sign(a12) * log(abs(a12)), y = sign(a21) * log(abs(a21))))+
    geom_point() +
    geom_text(aes(label = Window), hjust = 0, nudge_x = 0.1, size = 3) +
    geom_vline(xintercept = 0, linetype = 2) +
    geom_hline(yintercept = 0, linetype = 2) +
    theme_classic() +
    theme() +
    labs()

ggsave(here::here("simulation/plots/15-LV_fit_phase.png"), plot = p, width = 5, height = 5)

if (FALSE) {

data(gause_1934_book_f22)

logistic_data<-gause_1934_book_f22[gause_1934_book_f22$Treatment=="Pa",]

plot(Volume_Species2~Day, logistic_data)

# calculate per-capita growth rate
lagged_data$dNNdt <- percap_growth(x = lagged_data$x, laggedx = lagged_data$laggedx, dt = lagged_data$dt)

# plot relationship
plot(dNNdt~x, lagged_data,
     xlab="Abundance (N)",
     ylab="Per-capita Growth Rate (dN/Ndt)",
     xlim=c(0, 250), ylim=c(0, 1))
abline(h=0, v=0, lty=3)

# fit model to relationship
mod<-lm(dNNdt~x, lagged_data)
abline(mod, lwd=2, col=2)

# label parameters
arrows(25, 0.6, 1, 0.8, length = 0.1, lwd=1.5)
text(25, 0.6, "y-intercept: r", pos=1)

arrows(200, 0.4, 232, 0.01, length = 0.1, lwd=1.5)
text(200, 0.4, "x-intercept: K", pos=3)

arrows(80, predict(mod, newdata=data.frame(x=80)),
       130, predict(mod, newdata=data.frame(x=80)),
       length = 0.1, lwd=1.5)
arrows(130, predict(mod, newdata=data.frame(x=80)),
       130, predict(mod, newdata=data.frame(x=130)),
       length = 0.1, lwd=1.5)
text(130, predict(mod, newdata=data.frame(x=80)), "slope: aii", pos=3)

plot(Volume_Species2~Day, logistic_data, xlab="Day", ylab="P. aurelia Std. Volume")

timelst<-seq(0, 25, by=0.1) # sequence of time values to plot
prediction<-get_logistic(time = timelst, N0 = 0.6, r = 0.8, K=230)

lines(timelst, prediction, lwd=2)

# load data from competition experiment
data(gause_1934_book_f32)

# keep all data - no separate treatments exist for this experiment
predatorpreydata<-gause_1934_book_f32

# get time-lagged observations for each species
prey_lagged<-get_lag(x = predatorpreydata$Individuals_Prey, time = predatorpreydata$Day)
predator_lagged<-get_lag(x = predatorpreydata$Individuals_Predator, time = predatorpreydata$Day)

# calculate per-capita growth rates
prey_dNNdt<-percap_growth(x = prey_lagged$x, laggedx = prey_lagged$laggedx, dt = prey_lagged$dt)
predator_dNNdt<-percap_growth(x = predator_lagged$x,
                              laggedx = predator_lagged$laggedx, dt = predator_lagged$dt)

# fit linear models to dNNdt, based on average
# abundances between current and lagged time steps
prey_mod_dat<-data.frame(prey_dNNdt=prey_dNNdt, prey=prey_lagged$laggedx,
                         predator=predator_lagged$laggedx)
mod_prey<-lm(prey_dNNdt~prey+predator, data=prey_mod_dat)

predator_mod_dat<-data.frame(predator_dNNdt=predator_dNNdt,
                             predator=predator_lagged$laggedx, prey=prey_lagged$laggedx)
mod_predator<-lm(predator_dNNdt~predator+prey, data=predator_mod_dat)

# model summaries
summary(mod_prey)

# extract parameters
# growth rates
r1 <- unname(coef(mod_prey)["(Intercept)"])
r2 <- unname(coef(mod_predator)["(Intercept)"])

# self-limitation
a11 <- unname(coef(mod_prey)["prey"])
a22 <- unname(coef(mod_predator)["predator"])

# effect of Pa on Pc
a12 <- unname(coef(mod_prey)["predator"])
# effect of Pc on Pa
a21 <- unname(coef(mod_predator)["prey"])

# make parameter vector:
parms <- c(r1, r2, a11, a12, a21, a22)
initialN <- c(4, 0.1)
out <- deSolve::ode(y=initialN, times=seq(1, 17, length=100), func=lv_interaction, parms=parms)
matplot(out[,1], out[,-1], type="l",
        xlab="time", ylab="N", col=c("black","red"), lty=c(1,3), lwd=2, ylim=c(0, 60))
legend("topright", c("Pc", "Dn"), col=c(1,2), lwd=2, lty=c(1,3))


# now, plot in points from data
points(predatorpreydata$Day, predatorpreydata$Individuals_Predator , col=2)
points(predatorpreydata$Day, predatorpreydata$Individuals_Prey, col=1)


#load competition data
data("gause_1934_science_f02_03")

#subset out data from species grown in mixture
mixturedat<-gause_1934_science_f02_03[gause_1934_science_f02_03$Treatment=="Mixture",]

#extract time and species data
time<-mixturedat$Day
species<-data.frame(mixturedat$Volume_Species1, mixturedat$Volume_Species2)
colnames(species)<-c("P_caudatum", "P_aurelia")

#run wrapper
gause_out<-gause_wrapper(time=time, species=species)




















}











