# Miscellanuous functions

switch_pairwise_column <- function (df, bypair = T) {
    if (any(is.factor(df$Isolate1))) df$Isolate1 <- as.numeric(df$Isolate1); df$Isolate2 <- as.numeric(df$Isolate2)
    if ("Isolate1FreqPredicted" %in% colnames(df)) {
        if (bypair == T) {
            temp_index <- df$Isolate1 > df$Isolate2
            df[temp_index, c("Isolate1", "Isolate2", "Isolate1Freq", "Isolate2Freq", "Isolate1FreqPredicted", "Isolate2FreqPredicted")] <-
                df[temp_index, c("Isolate2", "Isolate1", "Isolate2Freq", "Isolate1Freq", "Isolate2FreqPredicted", "Isolate1FreqPredicted")]

            df %>% arrange(Isolate1, Isolate2, Isolate1Freq) %>% return()
        } else if (bypair == F) {
            temp_index <- df$Isolate1Freq == 5
            df[temp_index, c("Isolate1", "Isolate2", "Isolate1Freq", "Isolate2Freq", "Isolate1FreqPredicted", "Isolate2FreqPredicted")] <-
                df[temp_index, c("Isolate2", "Isolate1", "Isolate2Freq", "Isolate1Freq", "Isolate2FreqPredicted", "Isolate1FreqPredicted")]

            df %>% arrange(Isolate1Freq, Isolate1, Isolate2) %>% return()
        }
    } else {

        if (bypair == T) {
            temp_index <- df$Isolate1 > df$Isolate2
            df[temp_index, c("Isolate1", "Isolate2", "Isolate1Freq", "Isolate2Freq")] <-
                df[temp_index, c("Isolate2", "Isolate1", "Isolate2Freq", "Isolate1Freq")]

            df %>% arrange(Isolate1, Isolate2, Isolate1Freq) %>% return()
        } else if (bypair == F) {
            temp_index <- df$Isolate1Freq == 5
            df[temp_index, c("Isolate1", "Isolate2", "Isolate1Freq", "Isolate2Freq")] <-
                df[temp_index, c("Isolate2", "Isolate1", "Isolate2Freq", "Isolate1Freq")]

            df %>% arrange(Isolate1Freq, Isolate1, Isolate2) %>% return()
        }
    }
}

