#####################################################
# ~ ZFAD: SORL1 fingerprints F0 vs stable ~
#
# fingerprint & correlation heatmap
# 
# Francois Kroll 2023
# francois@kroll.be
#####################################################

# last ran 25/04/2023 to update after new definition of activitySunsetStartle


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

fgs <- calculateFingerprint(paDir=c(here('220531_SORL1', 'bhvparams'),
                                    here('220316_sorl1Stable', 'bhvparams')),
                            controlGrp=c('wt', 'scr'),
                            grporder=NA,
                            skipNight0=TRUE,
                            avgDayNight=TRUE)


# some corrections on fgs
# need to do manually as stable experiments are now called mergeExp_1 and mergeExp_2
fgs <- fgs %>%
  # add gene column
  add_column(gene='sorl1', .before='parameter') %>%
  # add type (f0 or stable) column
  add_column(type=sapply(.$grp, function(g) {
    if (g %in% c('wt', 'het', 'hom')) return('stable')
    else return('f0')
  }), .before='parameter')

# create composite columns gene_type, e.g. sorl1_f0 or clu_stable
fgs <- fgs %>%
  add_column(gene_type=paste(fgs$gene, fgs$type, sep='_'), .after='type')

# write fgs, we will use it in euclideanClutch.R
write.csv(fgs, here('fingerprints', 'sorl1Fgp.csv'), row.names=FALSE)


# draw fingerprint --------------------------------------------------------

fgs <- read.csv(here('fingerprints', 'sorl1Fgp.csv'))

# fgs <- read.csv('~/Dropbox/ZFAD/fingerprints/cluFgp.csv')

### f0 only
fgs %>%
  filter(type=='f0') %>%
  ggFingerprint(fgp=.,
                lmePath="~/Dropbox/ZFAD/LMEreports/sorl1_f0_LME.csv",
                onlyFish=NA,
                metric='mean',
                grporder='sorl1',
                controlGrp='scr',
                removeControl=FALSE,
                onlyExp=NA,
                removeParam=NA,
                colours=c('#F18D82', '#E94735'), # light = box14 = clutch2 n = 27 // dark = box15 = clutch1 n = 47
                connectOrNo=TRUE,
                legendOrNo=FALSE,
                ynameOrNo=FALSE,
                ytextOrNo=TRUE,
                xtextOrNo=TRUE,
                xParamNum=FALSE,
                nightBgOrNo=TRUE,
                ymin=-2,
                ymax=2.3,
                dotSize=0.1,
                lineSize=0.3,
                asteriskSize=5,
                exportOrNo=TRUE,
                exportPath='~/Dropbox/ZFAD/fingerprints/sorl1f0.pdf',
                width=100,
                height=47)
# 27/09/2023 making ymin & ymax same for all late-onset AD so comparable on figure

# version with labels
fgs %>%
  filter(type=='stable') %>%
  ggFingerprint(fgp=.,
                lmePath="~/Dropbox/ZFAD/LMEreports/sorl1_stable_LME.csv",
                onlyFish=NA,
                metric='mean',
                grporder='hom',
                controlGrp='wt',
                removeControl=FALSE,
                onlyExp=NA,
                removeParam=NA,
                colours=c('#316cbb', '#6999d8'),
                connectOrNo=TRUE,
                legendOrNo=FALSE,
                ynameOrNo=FALSE,
                ytextOrNo=TRUE,
                xtextOrNo=TRUE,
                xParamNum=FALSE,
                nightBgOrNo=TRUE,
                ymin=-2,
                ymax=2.3,
                dotSize=0.1,
                lineSize=0.3,
                asteriskSize=5,
                exportOrNo=TRUE,
                exportPath='~/Dropbox/ZFAD/fingerprints/sorl1Stable.pdf',
                width=100,
                height=47)


# fingerprint similarity heatmap ------------------------------------------

ggPairwiseHeat(fgp=fgs,
               simScore='cosine',
               grporder=c('sorl1', 'hom'),
               removeControl=TRUE,
               controlGrp=c('scr', 'wt'),
               minCol=NA,
               maxCol=NA,
               onlyHalf='upper',
               scoreSize=2,
               legendOrNo=FALSE,
               labelsOrNo=FALSE,
               width=33,
               height=33,
               exportPath='~/Dropbox/ZFAD/fingerprints/cos_sorl1.pdf')
