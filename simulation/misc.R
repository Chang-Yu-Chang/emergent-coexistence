# misc code

interaction_type <- c("exclusion", "coexistence", "neutrality", "mutual exclusion", "frequency-dependent\ncoexistence")
interaction_color = c("#DB7469", "#557BAA", "#8650C4", "red", "blue")
names(interaction_color) <- interaction_type

read_simulation_data <- function (folder, pattern, cs_output_format = T) {
    list.files(folder, pattern = pattern) %>%
        lapply(function(x) {
            experiment = str_split(x, "_", simplify = T)[1]
            exp_id = str_split(x, "_", simplify = T)[2] %>% as.numeric()
            paste0(folder, "/", x) %>%
                read_csv(col_types = cols()) %>%
                suppressMessages() %>%
                # Edit this ifelse if the species pool size changes
                {if (cs_output_format) select(., Family = `...1`, Species = `...2`, starts_with("W")) else {
                    mutate(., Family = paste0("F", rep(0:2, 100)), Species = paste0("S", 0:299))
                }} %>%
                na_if(0) %>%
                pivot_longer(cols = c(-Family, -Species), names_to = "Well", values_to = "Abundance", values_drop_na = T) %>%
                mutate(Experiment = experiment, exp_id = exp_id) %>%
                return()
        }) %>%
        bind_rows()
}


