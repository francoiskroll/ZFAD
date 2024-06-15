#####################################################
# ~ ZFAD: prepare cd2ap f0/stable fingerprint ~
#
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


### about APOEAB experiments

## f0

# 220308_APOEAB_1 had many issues (lights turned ON in the room, Zebrabox jumped open, incomplete because drive full)
# so not included in ZFAD folder
# anecdotally,
# box14: less active during the day/more sleep throughout
# box15: more active/less sleep night (strong effect)

# 220313_APOEAB_2: I was not controlling the pump well at the time, gave sharp peaks on activity trace
# so traces not pretty but should be OK to include
# box15 is good Ns
# box14 is N=8 apoeab, probably not worth including

# 220313_APOEAB_3: good experiment

## stable

# 220531_apoeabStable_1
# same clutch Box10/11; mostly good experiment
# Box10: N=3 Awt_Bwt / N=2 Ahom_Bhom; probably pointless to include
# Box11: N=5 Awt_Bwt / N=5 Ahom_Bhom

# 220706_apoeabStable_2
# same clutch Box14/15; good experiment
# Box14: N=7 Awt_Bwt / N=6 Ahom_Bhom
# Box15: N=9 Awt_Bwt / N=5 Ahom_Bhom

# 220706_apoeabStable_3
# same clutch Box16/17; good experiment
# Box16: N=8 Awt_Bwt / N=7 Ahom_Bhom
# Box17: N=0 Awt_Bwt / N=1 Ahom_Bhom, pointless to include


# calculate fingerprints --------------------------------------------------

# calculate clutch fingerprints
fgs <- calculateFingerprint(paDir=c(here('220313_APOEAB_2', 'bhvparams'),
                                    here('220516_APOEAB_3', 'bhvparams'),
                                    here('220531_apoeabStable_1', 'bhvparams'),
                                    here('220706_apoeabStable_2', 'bhvparams'), 
                                    here('220706_apoeabStable_3', 'bhvparams')),
                            controlGrp=c('Awt_Bwt', 'scr4x'),
                            mergeExp1=c('220706_14', '220706_15'), # ***
                            grporder=NA,
                            skipNight0=TRUE,
                            avgDayNight=TRUE)
# *** pooling Z-scores to controls within each box, not pooling raw datapoints
# 220706_14 & 220706_15 tracked the same clutch

# exclude some experiments
# (see comments above for more context)
fgs <- fgs %>%
  # 220531_10 is too low N of Awt_Bwt / Ahom_Bhom
  filter(date_box!='220531_10') %>%
  # 220706_17 is N= 0 Awt_Bwt, so fingerprint is already NA
  filter(date_box!='220706_17')


# some edits on fgs
# add gene & type column
fgs <- fgs %>%
  add_column(gene='apoeab', .before='parameter') %>%
  add_column(type=sapply(.$date_box, function(db) {
    if (db=='mergeExp_1') return('stable')
    if (db=='220531_11') return('stable')
    if (db=='220706_16') return('stable')
    else return('f0')
  }), .before='parameter')

# create composite columns gene_type, e.g. sorl1_f0 or clu_stable
fgs <- fgs %>%
  add_column(gene_type=paste(fgs$gene, fgs$type, sep='_'), .after='type')

# write fgs, we will use it in euclideanClutch.R
write.csv(fgs, here('fingerprints', 'apoeabFgp.csv'), row.names=FALSE)


# draw fingerprint --------------------------------------------------------

fgs <- read.csv(here('fingerprints', 'apoeabFgp.csv'))

fgs <- read.csv('~/Dropbox/ZFAD/fingerprints/apoeabFgp.csv')

### f0 only
fgs %>%
  filter(type=='f0') %>%
  ggFingerprint(fgp=.,
                lmePath="~/Dropbox/ZFAD/LMEreports/apoeab_f0_LME.csv",
                onlyFish=NA,
                metric='mean',
                grporder=NA,
                controlGrp=c('scr4x'),
                removeControl=TRUE,
                onlyExp=NA,
                removeParam=NA,
                colours=c('#E94735', '#F07F72', '#F6B7B0'),
                connectOrNo=TRUE,
                legendOrNo=FALSE,
                ynameOrNo=FALSE,
                ytextOrNo=TRUE,
                xtextOrNo=FALSE,
                xParamNum=TRUE,
                nightBgOrNo=TRUE,
                ymin=-2.2,
                ymax=2.3,
                dotSize=0.1,
                lineSize=0.3,
                asteriskSize=4,
                exportOrNo=TRUE,
                exportPath='~/Dropbox/ZFAD/fingerprints/apoeabf0.pdf',
                width=100,
                height=28)
# 27/09/2023 making ymin & ymax same for all late-onset AD so comparable on figure


### stable only
# version with labels
fgs %>%
  filter(type=='stable') %>%
  ggFingerprint(fgp=.,
                lmePath="~/Dropbox/ZFAD/LMEreports/apoeab_stable_LME.csv",
                onlyFish=NA,
                metric='mean',
                grporder='Ahom_Bhom',
                controlGrp='Awt_Bwt',
                removeControl=TRUE,
                onlyExp=NA,
                removeParam=NA,
                colours=c('#316cbb', '#5c8fd4', '#92b4e2'),
                connectOrNo=TRUE,
                legendOrNo=FALSE,
                ynameOrNo=FALSE,
                ytextOrNo=TRUE,
                xtextOrNo=TRUE,
                xParamNum=FALSE,
                nightBgOrNo=TRUE,
                ymin=-3,
                ymax=2.5,
                dotSize=0.1,
                lineSize=0.3,
                asteriskSize=4,
                exportOrNo=TRUE,
                exportPath='~/Dropbox/ZFAD/fingerprints/apoeabStable.pdf',
                width=100,
                height=47)


# fingerprint similarity heatmap ------------------------------------------

ggPairwiseHeat(fgp=fgs,
               simScore='cosine',
               grporder=c('apoeab', 'Ahom_Bhom'),
               removeControl=TRUE,
               controlGrp=c('scr4x', 'Awt_Bwt'),
               minCol=NA,
               maxCol=NA,
               onlyHalf='upper',
               scoreSize=2,
               legendOrNo=FALSE,
               labelsOrNo=FALSE,
               width=36,
               height=36,
               exportPath='~/Dropbox/ZFAD/fingerprints/cos_apoeab.pdf')