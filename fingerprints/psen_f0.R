#####################################################
# ~ ZFAD: prepare psen1/2 f0 fingerprint ~
#
#
# Francois Kroll 2023
# francois@kroll.be
#####################################################

# last ran 25/04/2023 after new definition of activitySunsetStartle


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



# PSEN1 -------------------------------------------------------------------

# calculate clutch fingerprints
fgs <- calculateFingerprint(paDir=here('210913_PSEN1', 'bhvparams'),
                            controlGrp='scr',
                            grporder=NA,
                            skipNight0=TRUE,
                            avgDayNight=TRUE)


# some edits on fgs
# add gene & type column
fgs <- fgs %>%
  add_column(gene='psen1', .before='parameter')

# write fgs, we will use it in euclideanClutch.R
write.csv(fgs, here('fingerprints', 'psen1Fgp.csv'), row.names=FALSE)


# draw fingerprint --------------------------------------------------------

fgs <- read.csv(here('fingerprints', 'psen1Fgp.csv'))

fgs <- read.csv('~/Dropbox/ZFAD/fingerprints/psen1Fgp.csv')

### psen1
fgs %>%
  ggFingerprint(fgp=.,
                lmePath="~/Dropbox/ZFAD/LMEreports/psen1_f0_LME.csv",
                onlyFish=NA,
                metric='mean',
                grporder=NA,
                controlGrp='scr',
                removeControl=TRUE,
                onlyExp=NA,
                removeParam=NA,
                colours=c('#fcb505', '#fcce92'),
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
                exportPath='~/Dropbox/ZFAD/fingerprints/psen1f0.pdf',
                width=100,
                height=47)


# fingerprint similarity heatmap ------------------------------------------

ggPairwiseHeat(fgp=fgs,
               simScore='cosine',
               grporder='psen1',
               removeControl=TRUE,
               controlGrp='scr',
               minCol=NA,
               maxCol=NA,
               onlyHalf='upper',
               scoreSize=2,
               legendOrNo=FALSE,
               labelsOrNo=FALSE,
               width=36,
               height=36,
               exportPath='~/Dropbox/ZFAD/fingerprints/cos_psen1.pdf')





# PSEN2 -------------------------------------------------------------------

# I compared fingerprints before vs after pixelAdjust
# overall behaves as expected and I think effort is worth it
# some notes:
# pixelAdjust makes every parameter same or *less* significant, except:

# active bout minimum (day & night):
# This is an obvious artefact.
# Not exactly sure what happens, maybe it brings the min px of all SCR bouts below 1, so gets rounded to 0, so min px becomes a quite high px, like 5.

# active bout duration (night):
# It makes SCR bouts shorter because start/end of bout gets below 1, which gets rounded to 0, hence makes bouts a bit shorter.
# This may not be an artefact, as it was already significant and it does not appear as an artefact during the day. Will not worry about it.

# active bout mean (night):
# psen2 becomes more negative. I cannot explain why but not an artefact in day, so will leave it. Seemed already (slightly) decreased before adjust.

## here, final plot: will delete active bout minimum (clear artefact) & activity sunset startle (obviously meaningless now)

# below, using parameter tables *after pixelAdjust*

# calculate clutch fingerprints
fgs <- calculateFingerprint(paDir=here('210907_PSEN2', 'bhvparams'),
                            controlGrp='scr',
                            grporder=NA,
                            skipNight0=TRUE,
                            avgDayNight=TRUE)


# some edits on fgs
# add gene & type column
fgs <- fgs %>%
  add_column(gene='psen2', .before='parameter')

# write fgs, we will use it in euclideanClutch.R
write.csv(fgs, here('fingerprints', 'psen2Fgp.csv'), row.names=FALSE)


# draw fingerprint --------------------------------------------------------

fgs <- read.csv(here('fingerprints', 'psen2Fgp.csv'))

fgs %>%
  ggFingerprint(fgp=.,
                lmePath="~/Dropbox/ZFAD/LMEreports/psen2_f0_LME.csv",
                onlyFish=NA,
                metric='mean',
                grporder=NA,
                controlGrp='scr',
                removeControl=TRUE,
                onlyExp=NA,
                removeParam=c('activeboutMin', 'activitySunsetStartle'),
                colours=c('#E94735', '#F18D82'),
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
                exportPath='~/Dropbox/ZFAD/fingerprints/psen2f0.pdf',
                width=100,
                height=68)


# fingerprint similarity heatmap ------------------------------------------

ggPairwiseHeat(fgp=fgs,
               simScore='cosine',
               grporder='psen2',
               removeControl=TRUE,
               controlGrp='scr',
               minCol=NA,
               maxCol=NA,
               onlyHalf='upper',
               scoreSize=2,
               legendOrNo=FALSE,
               labelsOrNo=FALSE,
               width=36,
               height=36,
               exportPath='~/Dropbox/ZFAD/fingerprints/cos_psen2.pdf')
