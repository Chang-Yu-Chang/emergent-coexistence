#' This script manually removes the outliner objects from the feature csv files
#'
#' The idea is that 03-green_gesmentation.R applies a general filter for identifying objects.
#' This filter works fine for most of the cases, but for some images it may let
#' some objects that obviously are not colony pass. This script is then to manually
#' apply customized filter to remove these object before they are passed to the
#' feature selection and clustering analysis
library(tidyverse)

# Bath D ----
# C1R2 ----
if (any(object_feature_combined$image_name == "D_T8_C1R2_3")) {
    temp_index <- which(object_feature_combined$image_name == "D_T8_C1R2_3" &
                            object_feature_combined$b.q05 > 0.224)
    object_feature_combined <- object_feature_combined[-temp_index,]
}

if (any(object_feature_combined$image_name == "D_T8_C1R2_4")) {
    temp_index <- which(object_feature_combined$image_name == "D_T8_C1R2_4" &
                            object_feature_combined$b.q005 < 0.15)
    object_feature_combined <- object_feature_combined[-temp_index,]
}

if (any(object_feature_combined$image_name == "D_T8_C1R2_5-95_2_1")) {
    temp_index <- which(object_feature_combined$image_name == "D_T8_C1R2_5-95_2_1" &
                            object_feature_combined$b.mean < 0.1)
    object_feature_combined <- object_feature_combined[-temp_index,]
}

if (any(object_feature_combined$image_name == "D_T8_C1R2_5-95_2_4")) {
    temp_index <- which(object_feature_combined$image_name == "D_T8_C1R2_5-95_2_4" &
                            object_feature_combined$s.area < 1000)
    object_feature_combined <- object_feature_combined[-temp_index,]
}

if (any(object_feature_combined$image_name == "D_T8_C1R2_5-95_2_1")) {
    temp_index <- which(object_feature_combined$image_name == "D_T8_C1R2_5-95_2_1" &
                            object_feature_combined$b.mean < 0.265)
    object_feature_combined <- object_feature_combined[-temp_index,]
}

if (any(object_feature_combined$image_name == "D_T8_C1R2_5-95_3_1")) {
    temp_index <- which(object_feature_combined$image_name == "D_T8_C1R2_5-95_3_1" &
                            object_feature_combined$b.mad > 0.2)
    object_feature_combined <- object_feature_combined[-temp_index,]
}


if (any(object_feature_combined$image_name == "D_T8_C1R2_50-50_1_2")) {
    temp_index <- which(object_feature_combined$image_name == "D_T8_C1R2_50-50_1_2" &
                            object_feature_combined$s.area < 1000)
    object_feature_combined <- object_feature_combined[-temp_index,]
}

if (any(object_feature_combined$image_name == "D_T8_C1R2_50-50_2_3")) {
    temp_index <- which(object_feature_combined$image_name == "D_T8_C1R2_50-50_2_3" &
                            object_feature_combined$b.q095 > 0.7)
    object_feature_combined <- object_feature_combined[-temp_index,]
}

if (any(object_feature_combined$image_name == "D_T8_C1R2_50-50_3_4")) {
    temp_index <- which(object_feature_combined$image_name == "D_T8_C1R2_50-50_3_4" &
                            object_feature_combined$b.sd > 0.3)
    object_feature_combined <- object_feature_combined[-temp_index,]
}

# C1R4 ----

if (any(object_feature_combined$image_name == "D_T8_C1R4_1")) {
    temp_index <- which(object_feature_combined$image_name == "D_T8_C1R4_1" &
                            object_feature_combined$b.mad > 0.03)
    object_feature_combined <- object_feature_combined[-temp_index,]
}

if (any(object_feature_combined$image_name == "D_T8_C1R4_2")) {
    temp_index <- which(object_feature_combined$image_name == "D_T8_C1R4_2" &
                            object_feature_combined$b.mad > 0.06)
    object_feature_combined <- object_feature_combined[-temp_index,]
}

if (any(object_feature_combined$image_name == "D_T8_C1R4_5-95_1_4")) {
    temp_index <- which(object_feature_combined$image_name == "D_T8_C1R4_5-95_1_4" &
                            object_feature_combined$b.mad < 0.03)
    object_feature_combined <- object_feature_combined[-temp_index,]
}


if (any(object_feature_combined$image_name == "D_T8_C1R4_5-95_1_5")) {
    temp_index <- which(object_feature_combined$image_name == "D_T8_C1R4_5-95_1_5" &
                            object_feature_combined$b.q05 > 0.4)
    object_feature_combined <- object_feature_combined[-temp_index,]
}


if (any(object_feature_combined$image_name == "D_T8_C1R4_5-95_3_1")) {
    temp_index <- which(object_feature_combined$image_name == "D_T8_C1R4_5-95_3_1" &
                            object_feature_combined$s.area < 1000)
    object_feature_combined <- object_feature_combined[-temp_index,]
}

if (any(object_feature_combined$image_name == "D_T8_C1R4_5-95_4_2")) {
    temp_index <- which(object_feature_combined$image_name == "D_T8_C1R4_5-95_4_2" &
                            object_feature_combined$b.mean < 0.15)
    object_feature_combined <- object_feature_combined[-temp_index,]
}

# C1R6 ----

if (any(object_feature_combined$image_name == "D_T8_C1R6_1")) {
    temp_index <- which(object_feature_combined$image_name == "D_T8_C1R6_1" &
                            object_feature_combined$s.area < 1000)
    object_feature_combined <- object_feature_combined[-temp_index,]
}

if (any(object_feature_combined$image_name == "D_T0_C1R6_3")) {
    temp_index <- which(object_feature_combined$image_name == "D_T0_C1R6_3" &
                            object_feature_combined$b.tran.sd < 0.05)
    object_feature_combined <- object_feature_combined[-temp_index,]
}

if (any(object_feature_combined$image_name == "D_T8_C1R6_5")) {
    temp_index <- which(object_feature_combined$image_name == "D_T8_C1R6_5" &
                            object_feature_combined$s.area < 500)
    object_feature_combined <- object_feature_combined[-temp_index,]
}

if (any(object_feature_combined$image_name == "D_T8_C1R6_5-95_1_4")) {
    temp_index <- which(object_feature_combined$image_name == "D_T8_C1R6_5-95_1_4" &
                            object_feature_combined$b.mad > 0.13)
    object_feature_combined <- object_feature_combined[-temp_index,]
}

if (any(object_feature_combined$image_name == "D_T8_C1R6_5-95_4_1")) {
    temp_index <- which(object_feature_combined$image_name == "D_T8_C1R6_5-95_4_1" &
                            object_feature_combined$s.area < 500)
    object_feature_combined <- object_feature_combined[-temp_index,]
}

if (any(object_feature_combined$image_name == "D_T8_C1R6_5-95_5_1")) {
    temp_index <- which(object_feature_combined$image_name == "D_T8_C1R6_5-95_5_1" &
                            object_feature_combined$s.area < 500)
    object_feature_combined <- object_feature_combined[-temp_index,]
}

if (any(object_feature_combined$image_name == "D_T8_C1R6_50-50_1_2")) {
    temp_index <- which(object_feature_combined$image_name == "D_T8_C1R6_50-50_1_2" &
                            object_feature_combined$s.area < 500)
    object_feature_combined <- object_feature_combined[-temp_index,]
}

if (any(object_feature_combined$image_name == "D_T8_C1R6_50-50_1_4")) {
    temp_index <- which(object_feature_combined$image_name == "D_T8_C1R6_50-50_1_4" &
                            object_feature_combined$b.mad > 0.13)
    object_feature_combined <- object_feature_combined[-temp_index,]
}

if (any(object_feature_combined$image_name == "D_T8_C1R6_50-50_2_3")) {
    temp_index <- which(object_feature_combined$image_name == "D_T8_C1R6_50-50_2_3" &
                            object_feature_combined$s.area < 500)
    object_feature_combined <- object_feature_combined[-temp_index,]
}

# C1R7 ----

if (any(object_feature_combined$image_name == "D_T8_C1R7_5-95_7_1")) {
    temp_index <- which(object_feature_combined$image_name == "D_T8_C1R7_5-95_7_1" &
                            object_feature_combined$s.area < 500)
    object_feature_combined <- object_feature_combined[-temp_index,]
}

if (any(object_feature_combined$image_name == "D_T8_C1R7_5-95_3_4")) {
    temp_index <- which(object_feature_combined$image_name == "D_T8_C1R7_5-95_3_4" &
                            object_feature_combined$b.tran.mad > 0.2)
    object_feature_combined <- object_feature_combined[-temp_index,]
}
if (any(object_feature_combined$image_name == "D_T8_C1R7_5-95_4_3")) {
    temp_index <- which(object_feature_combined$image_name == "D_T8_C1R7_5-95_4_3" &
                            object_feature_combined$b.tran.mad > 0.075)
    object_feature_combined <- object_feature_combined[-temp_index,]
}
if (any(object_feature_combined$image_name == "D_T8_C1R7_5-95_4_6")) {
    temp_index <- which(object_feature_combined$image_name == "D_T8_C1R7_5-95_4_6" &
                            object_feature_combined$b.tran.mad > 0.075)
    object_feature_combined <- object_feature_combined[-temp_index,]
}

if (any(object_feature_combined$image_name == "D_T8_C1R7_50-50_1_2")) {
    temp_index <- which(object_feature_combined$image_name == "D_T8_C1R7_50-50_1_2" &
                            object_feature_combined$b.q05 < 0.175 )
    object_feature_combined <- object_feature_combined[-temp_index,]
}

if (any(object_feature_combined$image_name == "D_T8_C1R7_50-50_1_3")) {
    temp_index <- which(object_feature_combined$image_name == "D_T8_C1R7_50-50_1_3" &
                            object_feature_combined$b.q005 < -0.2)
    object_feature_combined <- object_feature_combined[-temp_index,]
}


if (any(object_feature_combined$image_name == "D_T8_C1R7_50-50_3_7")) {
    temp_index <- which(object_feature_combined$image_name == "D_T8_C1R7_50-50_3_7" &
                            object_feature_combined$b.mad > 0.5 )
    object_feature_combined <- object_feature_combined[-temp_index,]
}


# C4R1 ----
if (any(object_feature_combined$image_name == "D_T8_C4R1_5-95_1_2")) {
    temp_index <- which(object_feature_combined$image_name == "D_T8_C4R1_5-95_1_2" &
                            object_feature_combined$b.q005 < -0.15)
    object_feature_combined <- object_feature_combined[-temp_index,]
}



# C11R5 ----

if (any(object_feature_combined$image_name == "D_T8_C11R5_5-95_1_5")) {
    temp_index <- which(object_feature_combined$image_name == "D_T8_C11R5_5-95_1_5" &
                            object_feature_combined$s.area < 1000)
    object_feature_combined <- object_feature_combined[-temp_index,]
}

if (any(object_feature_combined$image_name == "D_T8_C11R5_5-95_3_5")) {
    temp_index <- which(object_feature_combined$image_name == "D_T8_C11R5_5-95_3_5" &
                            object_feature_combined$s.area < 500)
    object_feature_combined <- object_feature_combined[-temp_index,]
}

if (any(object_feature_combined$image_name == "D_T8_C11R5_5-95_5_4")) {
    temp_index <- which(object_feature_combined$image_name == "D_T8_C11R5_5-95_5_4" &
                            object_feature_combined$s.area < 500)
    object_feature_combined <- object_feature_combined[-temp_index,]
}

if (any(object_feature_combined$image_name == "D_T8_C11R5_5-95_4_5")) {
    temp_index <- which(object_feature_combined$image_name == "D_T8_C11R5_5-95_4_5" &
                            object_feature_combined$b.tran.q005 < 0.41)
    object_feature_combined <- object_feature_combined[-temp_index,]
}


if (any(object_feature_combined$image_name == "D_T8_C11R5_50-50_4_5")) {
    temp_index <- which(object_feature_combined$image_name == "D_T8_C11R5_50-50_4_5" &
                            object_feature_combined$b.tran.q005 < 0.4)
    object_feature_combined <- object_feature_combined[-temp_index,]
}



# Batch C2 ----
# C11R2 ----

if (any(object_feature_combined$image_name == "C2_T0_C11R2_2")) {
    temp_index <- which(object_feature_combined$image_name == "C2_T0_C11R2_2" &
                            object_feature_combined$b.q005 < -0.1)
    object_feature_combined <- object_feature_combined[-temp_index,]
}

if (any(object_feature_combined$image_name == "C2_T8_C11R2_4")) {
    temp_index <- which(object_feature_combined$image_name == "C2_T8_C11R2_4" &
                            object_feature_combined$b.mad > 0.2)
    object_feature_combined <- object_feature_combined[-temp_index,]
}

if (any(object_feature_combined$image_name == "C2_T8_C11R2_6")) {
    temp_index <- which(object_feature_combined$image_name == "C2_T8_C11R2_6" &
                            object_feature_combined$b.q005 < 0.05)
    object_feature_combined <- object_feature_combined[-temp_index,]
}

if (any(object_feature_combined$image_name == "C2_T8_C11R2_7")) {
    temp_index <- which(object_feature_combined$image_name == "C2_T8_C11R2_7" &
                            object_feature_combined$b.mad < 0.03)
    object_feature_combined <- object_feature_combined[-temp_index,]
}

if (any(object_feature_combined$image_name == "C2_T8_C11R2_9")) {
    temp_index <- which(object_feature_combined$image_name == "C2_T8_C11R2_9" &
                            object_feature_combined$b.q005 < -0.1)
    object_feature_combined <- object_feature_combined[-temp_index,]
}


if (any(object_feature_combined$image_name == "C2_T8_C11R2_5-95_1_9")) {
    temp_index <- which(object_feature_combined$image_name == "C2_T8_C11R2_5-95_1_9" &
                            object_feature_combined$b.mad > 0.06)
    object_feature_combined <- object_feature_combined[-temp_index,]
}

if (any(object_feature_combined$image_name == "C2_T8_C11R2_5-95_2_5")) {
    temp_index <- which(object_feature_combined$image_name == "C2_T8_C11R2_5-95_2_5" &
                            object_feature_combined$b.mad > 0.3)
    object_feature_combined <- object_feature_combined[-temp_index,]
}

if (any(object_feature_combined$image_name == "C2_T8_C11R2_5-95_2_9")) {
    temp_index <- which(object_feature_combined$image_name == "C2_T8_C11R2_5-95_2_9" &
                            object_feature_combined$b.q005 > 0.4)
    object_feature_combined <- object_feature_combined[-temp_index,]
}

if (any(object_feature_combined$image_name == "C2_T8_C11R2_5-95_3_10")) {
    temp_index <- which(object_feature_combined$image_name == "C2_T8_C11R2_5-95_3_10" &
                            object_feature_combined$b.mad > 0.15)
    object_feature_combined <- object_feature_combined[-temp_index,]
}

if (any(object_feature_combined$image_name == "C2_T8_C11R2_5-95_3_12")) {
    temp_index <- which(object_feature_combined$image_name == "C2_T8_C11R2_5-95_3_12" &
                            object_feature_combined$b.q005 < 0.15)
    object_feature_combined <- object_feature_combined[-temp_index,]
}

if (any(object_feature_combined$image_name == "C2_T8_C11R2_5-95_5_8")) {
    temp_index <- which(object_feature_combined$image_name == "C2_T8_C11R2_5-95_5_8" &
                            object_feature_combined$b.q005 < -0.1)
    object_feature_combined <- object_feature_combined[-temp_index,]
}

if (any(object_feature_combined$image_name == "C2_T8_C11R2_5-95_5_12")) {
    temp_index <- which(object_feature_combined$image_name == "C2_T8_C11R2_5-95_5_12" &
                            object_feature_combined$b.q005 < 0)
    object_feature_combined <- object_feature_combined[-temp_index,]
}

if (any(object_feature_combined$image_name == "C2_T8_C11R2_5-95_5_13")) {
    temp_index <- which(object_feature_combined$image_name == "C2_T8_C11R2_5-95_5_13" &
                            object_feature_combined$b.q005 < 0)
    object_feature_combined <- object_feature_combined[-temp_index,]
}

if (any(object_feature_combined$image_name == "C2_T8_C11R2_5-95_7_1")) {
    temp_index <- which(object_feature_combined$image_name == "C2_T8_C11R2_5-95_7_1" &
                            object_feature_combined$b.q005 < 0)
    object_feature_combined <- object_feature_combined[-temp_index,]
}

if (any(object_feature_combined$image_name == "C2_T8_C11R2_5-95_7_3")) {
    temp_index <- which(object_feature_combined$image_name == "C2_T8_C11R2_5-95_7_3" &
                            object_feature_combined$b.q005 < -.5)
    object_feature_combined <- object_feature_combined[-temp_index,]
}

if (any(object_feature_combined$image_name == "C2_T8_C11R2_5-95_7_8")) {
    temp_index <- which(object_feature_combined$image_name == "C2_T8_C11R2_5-95_7_8" &
                            object_feature_combined$b.mad > 0.075)
    object_feature_combined <- object_feature_combined[-temp_index,]
}

if (any(object_feature_combined$image_name == "C2_T8_C11R2_5-95_7_9")) {
    temp_index <- which(object_feature_combined$image_name == "C2_T8_C11R2_5-95_7_9" &
                            object_feature_combined$b.q05 < 0.24)
    object_feature_combined <- object_feature_combined[-temp_index,]
}

if (any(object_feature_combined$image_name == "C2_T8_C11R2_5-95_8_7")) {
    temp_index <- which(object_feature_combined$image_name == "C2_T8_C11R2_5-95_8_7" &
                            object_feature_combined$b.q005 < -0.1)
    object_feature_combined <- object_feature_combined[-temp_index,]
}

if (any(object_feature_combined$image_name == "C2_T8_C11R2_5-95_8_10")) {
    temp_index <- which(object_feature_combined$image_name == "C2_T8_C11R2_5-95_8_10" &
                            object_feature_combined$b.mad > 0.18)
    object_feature_combined <- object_feature_combined[-temp_index,]
}


if (any(object_feature_combined$image_name == "C2_T8_C11R2_5-95_9_6")) {
    temp_index <- which(object_feature_combined$image_name == "C2_T8_C11R2_5-95_9_6" &
                            object_feature_combined$b.q005 < 0.03)
    object_feature_combined <- object_feature_combined[-temp_index,]
}

if (any(object_feature_combined$image_name == "C2_T8_C11R2_5-95_10_1")) {
    temp_index <- which(object_feature_combined$image_name == "C2_T8_C11R2_5-95_10_1" &
                            object_feature_combined$b.q005 > 0.4)
    object_feature_combined <- object_feature_combined[-temp_index,]
}

if (any(object_feature_combined$image_name == "C2_T8_C11R2_5-95_10_12")) {
    temp_index <- which(object_feature_combined$image_name == "C2_T8_C11R2_5-95_10_12" &
                            object_feature_combined$b.mad > 0.15)
    object_feature_combined <- object_feature_combined[-temp_index,]
}

if (any(object_feature_combined$image_name == "C2_T8_C11R2_5-95_11_4")) {
    temp_index <- which(object_feature_combined$image_name == "C2_T8_C11R2_5-95_11_4" &
                            object_feature_combined$b.tran.q05 < 0.4)
    object_feature_combined <- object_feature_combined[-temp_index,]
}

if (any(object_feature_combined$image_name == "C2_T8_C11R2_5-95_11_10")) {
    temp_index <- which(object_feature_combined$image_name == "C2_T8_C11R2_5-95_11_10" &
                            object_feature_combined$b.mad > 0.06)
    object_feature_combined <- object_feature_combined[-temp_index,]
}


if (any(object_feature_combined$image_name == "C2_T8_C11R2_5-95_12_7")) {
    temp_index <- which(object_feature_combined$image_name == "C2_T8_C11R2_5-95_12_7" &
                            object_feature_combined$b.q005 < 0.15)
    object_feature_combined <- object_feature_combined[-temp_index,]
}

if (any(object_feature_combined$image_name == "C2_T8_C11R2_5-95_12_9")) {
    temp_index <- which(object_feature_combined$image_name == "C2_T8_C11R2_5-95_12_9" &
                            object_feature_combined$b.q05 < 0.24)
    object_feature_combined <- object_feature_combined[-temp_index,]
}

if (any(object_feature_combined$image_name == "C2_T8_C11R2_5-95_13_11")) {
    temp_index <- which(object_feature_combined$image_name == "C2_T8_C11R2_5-95_13_11" &
                            (object_feature_combined$b.q005 < -0.75 | object_feature_combined$b.q005 > 0.4))
    object_feature_combined <- object_feature_combined[-temp_index,]
}


if (any(object_feature_combined$image_name == "C2_T8_C11R2_50-50_1_10")) {
    temp_index <- which(object_feature_combined$image_name == "C2_T8_C11R2_50-50_1_10" &
                            object_feature_combined$b.q05 < 0.1)
    object_feature_combined <- object_feature_combined[-temp_index,]
}


if (any(object_feature_combined$image_name == "C2_T8_C11R2_50-50_3_11")) {
    temp_index <- which(object_feature_combined$image_name == "C2_T8_C11R2_50-50_3_11" &
                            object_feature_combined$b.q05 > 0.4)
    object_feature_combined <- object_feature_combined[-temp_index,]
}

if (any(object_feature_combined$image_name == "C2_T8_C11R2_50-50_3_12")) {
    temp_index <- which(object_feature_combined$image_name == "C2_T8_C11R2_50-50_3_12" &
                            object_feature_combined$b.q005 < 0.1)
    object_feature_combined <- object_feature_combined[-temp_index,]
}

if (any(object_feature_combined$image_name == "C2_T8_C11R2_50-50_5_8")) {
    temp_index <- which(object_feature_combined$image_name == "C2_T8_C11R2_50-50_5_8" &
                            object_feature_combined$b.q005 < -0.1)
    object_feature_combined <- object_feature_combined[-temp_index,]
}

if (any(object_feature_combined$image_name == "C2_T8_C11R2_50-50_7_12")) {
    temp_index <- which(object_feature_combined$image_name == "C2_T8_C11R2_50-50_7_12" &
                            object_feature_combined$b.q005 < 0.15)
    object_feature_combined <- object_feature_combined[-temp_index,]
}

if (any(object_feature_combined$image_name == "C2_T8_C11R2_50-50_8_11")) {
    temp_index <- which(object_feature_combined$image_name == "C2_T8_C11R2_50-50_8_11" &
                            object_feature_combined$b.q005 < 0.15)
    object_feature_combined <- object_feature_combined[-temp_index,]
}

if (any(object_feature_combined$image_name == "C2_T8_C11R2_50-50_11_12")) {
    temp_index <- which(object_feature_combined$image_name == "C2_T8_C11R2_50-50_11_12" &
                            object_feature_combined$b.mad > 0.1)
    object_feature_combined <- object_feature_combined[-temp_index,]
}




