#####################################################
# ~ ZFAD: prepare appab f0 fingerprint ~
#
#
# Francois Kroll 2023
# francois@kroll.be
#####################################################

# last ran 25/04/2023 to update after new definition of activitySunsetStartle

# Note, main experiment is 220524_APPAB which 2 gRNAs for APPA / 2 gRNAs for APPB
# but I had also done APPA alone (210927_APPA), see thesis


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


# calculate fingerprints --------------------------------------------------

# calculate clutch fingerprints
fgs <- calculateFingerprint(paDir=here('220524_APPAB', 'bhvparams'),
                            controlGrp=c('scr', 'scr4x'),
                            grporder=NA,
                            skipNight0=TRUE,
                            avgDayNight=TRUE)


# some edits on fgs
# add gene & type column
fgs <- fgs %>%
  add_column(gene='appab', .before='parameter')

# write fgs, we will use it in euclideanClutch.R
write.csv(fgs, here('fingerprints', 'appabFgp.csv'), row.names=FALSE)


# draw fingerprint --------------------------------------------------------

fgs <- read.csv(here('fingerprints', 'appabFgp.csv'))

fgs <- read.csv('~/Dropbox/ZFAD/fingerprints/appabFgp.csv')

### appab
fgs %>%
  ggFingerprint(fgp=.,
                lmePath="~/Dropbox/ZFAD/LMEreports/appab_f0_LME.csv",
                onlyFish=NA,
                metric='mean',
                grporder=NA,
                controlGrp='scr4x',
                removeControl=TRUE,
                onlyExp=NA,
                removeParam=NA,
                colours=c('#E94735', '#F18D82'),
                connectOrNo=TRUE,
                legendOrNo=FALSE,
                ynameOrNo=FALSE,
                ytextOrNo=TRUE,
                xtextOrNo=TRUE,
                xParamNum=FALSE,
                nightBgOrNo=TRUE,
                ymin=-2,
                ymax=1,
                dotSize=0.1,
                lineSize=0.3,
                asteriskSize=4,
                exportOrNo=TRUE,
                exportPath='~/Dropbox/ZFAD/fingerprints/appabf0.pdf',
                width=100,
                height=47)


# fingerprint similarity heatmap ------------------------------------------

ggPairwiseHeat(fgp=fgs,
               simScore='cosine',
               grporder='appab',
               removeControl=TRUE,
               controlGrp='scr4x',
               minCol=NA,
               maxCol=NA,
               onlyHalf='upper',
               scoreSize=2,
               legendOrNo=FALSE,
               labelsOrNo=FALSE,
               width=36,
               height=36,
               exportPath='~/Dropbox/ZFAD/fingerprints/cos_appab.pdf')



# appa --------------------------------------------------------------------
# not sure if will use it

# calculate clutch fingerprints
fgs <- calculateFingerprint(paDir=here('210927_APPA', 'bhvparams'),
                            controlGrp='scr',
                            grporder=NA,
                            skipNight0=TRUE,
                            avgDayNight=TRUE)


# some edits on fgs
# add gene & type column
fgs <- fgs %>%
  add_column(gene='appa', .before='parameter')

# write fgs, we will use it in euclideanClutch.R
write.csv(fgs, here('fingerprints', 'appaFgp.csv'), row.names=FALSE)

fgs <- read.csv(here('fingerprints', 'appaFgp.csv'))
fgs <- read.csv('~/Dropbox/ZFAD/fingerprints/appaFgp.csv')

### appab
fgs %>%
  ggFingerprint(fgp=.,
                lmePath="~/Dropbox/ZFAD/LMEreports/appa_f0_LME.csv",
                onlyFish=NA,
                metric='mean',
                grporder=NA,
                controlGrp='scr',
                removeControl=TRUE,
                onlyExp=NA,
                removeParam=NA,
                colours=c('#E94735', '#F18D82'),
                connectOrNo=TRUE,
                legendOrNo=FALSE,
                ynameOrNo=FALSE,
                ytextOrNo=TRUE,
                xtextOrNo=TRUE,
                xParamNum=FALSE,
                nightBgOrNo=TRUE,
                ymin=-2,
                ymax=1,
                dotSize=0.1,
                lineSize=0.3,
                asteriskSize=4,
                exportOrNo=TRUE,
                exportPath='~/Dropbox/ZFAD/fingerprints/appaf0.pdf',
                width=100,
                height=47)

ggPairwiseHeat(fgp=fgs,
               simScore='cosine',
               grporder='appa',
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
               exportPath='~/Dropbox/ZFAD/fingerprints/cos_appa.pdf')
