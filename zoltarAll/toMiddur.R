#####################################################
# ~ ZFAD: convert all experiments to middur to give to ZOLTAR ~
#
# for eLife reviews
#
# Francois Kroll 2024
# francois@kroll.be
#####################################################


# packages ----------------------------------------------------------------

library(here)
library(FramebyFrame)


# psen1 -------------------------------------------------------------------

rawToMiddur(ffpath=here('210913_PSEN1/210913_12_RAWs.csv'),
            freezing=3,
            burst=200,
            exportOrNo=TRUE)

rawToMiddur(ffpath=here('210913_PSEN1/210913_13_RAWs.csv'),
            freezing=3,
            burst=200,
            exportOrNo=TRUE)


# psen2 -------------------------------------------------------------------

# ! just noticed that psen2PredictRescue used pixelAdjusted data
# but I think decision was to keep raw data,
# just show with pixelAdjust that it does not affect much the phenotype
# should ideally re-do that figure

# here, re-calculate the middur data
rawToMiddur(ffpath=here('210907_PSEN2/210907_12_RAWs.csv'),
            freezing=3,
            burst=200,
            exportOrNo=TRUE)

rawToMiddur(ffpath=here('210907_PSEN2/210907_13_RAWs.csv'),
            freezing=3,
            burst=200,
            exportOrNo=TRUE)


# appa/appb ---------------------------------------------------------------

rawToMiddur(ffpath=here('220524_APPAB/220524_14_RAWs.csv'),
            freezing=3,
            burst=200,
            exportOrNo=TRUE)

rawToMiddur(ffpath=here('220524_APPAB/220524_15_RAWs.csv'),
            freezing=3,
            burst=200,
            exportOrNo=TRUE)


# apoea/apoeb -------------------------------------------------------------

rawToMiddur(ffpath=here('220313_APOEAB_2/220313_15_RAWs.csv'),
            freezing=3,
            burst=200,
            exportOrNo=TRUE)

rawToMiddur(ffpath=here('220516_APOEAB_3/220516_14_RAWs.csv'),
            freezing=3,
            burst=200,
            exportOrNo=TRUE)

rawToMiddur(ffpath=here('220516_APOEAB_3/220516_15_RAWs.csv'),
            freezing=3,
            burst=200,
            exportOrNo=TRUE)


# cd2ap -------------------------------------------------------------------

rawToMiddur(ffpath=here('220725_CD2AP/220725_16_RAWs.csv'),
            freezing=3,
            burst=200,
            exportOrNo=TRUE)

rawToMiddur(ffpath=here('220725_CD2AP/220725_17_RAWs.csv'),
            freezing=3,
            burst=200,
            exportOrNo=TRUE)


# clu ---------------------------------------------------------------------

rawToMiddur(ffpath=here('220601_CLU/220601_16_RAWs.csv'),
            freezing=3,
            burst=200,
            exportOrNo=TRUE)

rawToMiddur(ffpath=here('220601_CLU/220601_17_RAWs.csv'),
            freezing=3,
            burst=200,
            exportOrNo=TRUE)


# sorl1 -------------------------------------------------------------------

rawToMiddur(ffpath=here('220531_SORL1/220531_14_RAWs.csv'),
            freezing=3,
            burst=200,
            exportOrNo=TRUE)

rawToMiddur(ffpath=here('220531_SORL1/220531_15_RAWs.csv'),
            freezing=3,
            burst=200,
            exportOrNo=TRUE)