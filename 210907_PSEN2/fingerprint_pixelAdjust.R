### psen2: fingerprint before/after pixelAdjust


# packages & functions ----------------------------------------------------

# load the Frame-by-Frame package
# library(devtools)
# install_github('francoiskroll/FramebyFrame')

library(FramebyFrame)
library(here)

library(dplyr)
library(tidyr)
library(tibble)

library(openxlsx)


##### PSEN2 #####
### before pixelAdjust ###

# calculate clutch fingerprints
fgs <- calculateFingerprint(paDir=here('210907_PSEN2', 'bhvparams_beforeAdjust'),
                            controlGrp='scr',
                            grporder=NA,
                            skipNight0=TRUE,
                            avgDayNight=TRUE)


# draw fingerprint --------------------------------------------------------

# delete parameters
# night_sunsetStartle
# night_activeboutMin
# because they cannot be used after pixelAdjust and we want to make it easy to compare fingerprints
fgsall <- fgs
fgs <- fgs %>%
  filter(! parameter %in% c('activitySunsetStartle', 'activeboutMin'))

ggFingerprint(fgp=fgs,
              lmePath=here('LMEreports/psen2_f0_LME_beforeAdjust.csv'),
              onlyFish=NA,
              metric='mean',
              grporder=NA,
              controlGrp='scr',
              removeControl=TRUE,
              onlyExp=NA,
              removeParam=NA,
              colours=c('#b1ce92', '#92b471'), # light = box1 = clutch2 // dark = box2 = clutch1
              connectOrNo=TRUE,
              legendOrNo=FALSE,
              ynameOrNo=FALSE,
              ytextOrNo=TRUE,
              xtextOrNo=FALSE,
              xParamNum=TRUE,
              nightBgOrNo=TRUE,
              ymin=-2,
              ymax=2,
              dotSize=0.1,
              lineSize=0.3,
              asteriskSize=4,
              exportOrNo=TRUE,
              exportPath=here('fingerprints/psen2f0_beforeAdjust.pdf'),
              width=100,
              height=28)

####
# for the record, complete PSEN2 fingerprint
ggFingerprint(fgp=fgsall,
              lmePath=here('LMEreports/psen2_f0_LME_beforeAdjust.csv'),
              onlyFish=NA,
              metric='mean',
              grporder=NA,
              controlGrp='scr',
              removeControl=TRUE,
              onlyExp=NA,
              removeParam=NA,
              colours=c('#b1ce92', '#92b471'), # light = box1 = clutch2 // dark = box2 = clutch1
              connectOrNo=TRUE,
              legendOrNo=FALSE,
              ynameOrNo=FALSE,
              ytextOrNo=TRUE,
              xtextOrNo=TRUE,
              xParamNum=FALSE,
              nightBgOrNo=TRUE,
              ymin=-2,
              ymax=2,
              dotSize=0.1,
              lineSize=0.3,
              asteriskSize=4,
              exportOrNo=TRUE,
              exportPath=here('fingerprints/psen2f0_beforeAdjust_ALLPARAMETERS.pdf'),
              width=100,
              height=47)

####


### after pixelAdjust ###
# calculate clutch fingerprints
fgs <- calculateFingerprint(paDir=here('210907_PSEN2', 'bhvparams'),
                            controlGrp='scr',
                            grporder=NA,
                            skipNight0=TRUE,
                            avgDayNight=TRUE)


# draw fingerprint --------------------------------------------------------

# delete parameters
# night_sunsetStartle
# night_activeboutMin
# because they cannot be used after pixelAdjust and we want to make it easy to compare fingerprints
fgs <- fgs %>%
  filter(! parameter %in% c('activitySunsetStartle', 'activeboutMin'))

ggFingerprint(fgp=fgs,
              lmePath=here('LMEreports/psen2_f0_LME.csv'),
              onlyFish=NA,
              metric='mean',
              grporder=NA,
              controlGrp='scr',
              removeControl=TRUE,
              onlyExp=NA,
              removeParam=NA,
              colours=c('#b1ce92', '#92b471'), # light = box1 = clutch2 // dark = box2 = clutch1
              connectOrNo=TRUE,
              legendOrNo=FALSE,
              ynameOrNo=FALSE,
              ytextOrNo=TRUE,
              xtextOrNo=TRUE,
              xParamNum=FALSE,
              nightBgOrNo=TRUE,
              ymin=-2,
              ymax=2,
              dotSize=0.1,
              lineSize=0.3,
              asteriskSize=4,
              exportOrNo=TRUE,
              exportPath=here('fingerprints/psen2f0_afterAdjust.pdf'),
              width=100,
              height=47)