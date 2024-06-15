#####################################################
# ~ ZFAD: prepare clu f0/stable fingerprint ~
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



# calculate fingerprints --------------------------------------------------

fgs <- calculateFingerprint(paDir=c(here('220601_CLU', 'bhvparams'),
                                    here('220906_cluStable_1', 'bhvparams'),
                                    here('220906_cluStable_2', 'bhvparams')),
                            controlGrp=c('wt', 'scr'),
                            mergeExp1=c('220906_14', '220906_15'), # ***
                            mergeExp2=c('220906_16', '220906_17'), # ***
                            grporder=NA,
                            skipNight0=TRUE,
                            avgDayNight=TRUE)
# *** ! mergeExp does not pool raw datapoints, it pools Z-scores to controls within each box
# each set (mergeExp1 / mergeExp2) is two experiments where same clutch was tracked in two parallel boxes


# some corrections on fgs
# need to do manually as stable experiments are now called mergeExp_1 and mergeExp_2
fgs <- fgs %>%
  # add gene column
  add_column(gene='clu', .before='parameter') %>%
  # add type (f0 or stable) column
  add_column(type=sapply(.$date_box, function(db) {
    if (db=='mergeExp_1') return('stable')
    if (db=='mergeExp_2') return('stable')
    else return('f0')
  }), .before='parameter')


# # edit grp column
# # undo factors
# fgs$grp <- as.character(fgs$grp)
# # f0 KOs should always be called ko, not gene name
# fgs[which(fgs$type=='f0' & fgs$grp!='scr'), 'grp'] <- 'ko'
# # consequently, need to re-create date_box_grp and date_box_period_grp
# fgs$date_box_grp <- paste(fgs$date_box, fgs$grp, sep='_')
# fgs$date_box_period_grp <- paste(fgs$date_box, fgs$period, fgs$grp, sep='_')


# create composite columns gene_type, e.g. sorl1_f0 or clu_stable
fgs <- fgs %>%
  add_column(gene_type=paste(fgs$gene, fgs$type, sep='_'), .after='type')

# write fgs, we will use it in euclideanClutch.R
write.csv(fgs, here('fingerprints', 'cluFgp.csv'), row.names=FALSE)


# draw fingerprint --------------------------------------------------------

fgs <- read.csv(here('fingerprints', 'cluFgp.csv'))

# fgs <- read.csv('~/Dropbox/ZFAD/fingerprints/cluFgp.csv')

### f0 only
fgs %>%
  filter(type=='f0') %>%
  ggFingerprint(fgp=.,
                lmePath="~/Dropbox/ZFAD/LMEreports/clu_f0_LME.csv",
                onlyFish=NA,
                metric='mean',
                grporder='clu',
                controlGrp='scr',
                removeControl=FALSE,
                onlyExp=NA,
                removeParam=NA,
                colours=c('#E94735', '#F18D82', '#316cbb', '#6999d8'),
                connectOrNo=TRUE,
                legendOrNo=FALSE,
                ynameOrNo=FALSE,
                ytextOrNo=TRUE,
                xtextOrNo=FALSE,
                xParamNum=TRUE,
                nightBgOrNo=TRUE,
                ymin=-2,
                ymax=2.3,
                dotSize=0.1,
                lineSize=0.3,
                asteriskSize=5,
                exportOrNo=TRUE,
                exportPath='~/Dropbox/ZFAD/fingerprints/cluf0.pdf',
                width=100,
                height=28)
# 27/09/2023 making ymin & ymax same for all late-onset AD so comparable on figure

### stable
fgs %>%
  filter(type=='stable') %>%
  ggFingerprint(fgp=.,
                lmePath="~/Dropbox/ZFAD/LMEreports/clu_stable_LME.csv",
                onlyFish=NA,
                metric='mean',
                grporder=c('hom'),
                controlGrp=c('wt', 'scr'),
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
                ymin=-1.5,
                ymax=1.5,
                dotSize=0.1,
                lineSize=0.3,
                asteriskSize=5,
                exportOrNo=TRUE,
                exportPath='~/Dropbox/ZFAD/fingerprints/cluStablev2.pdf',
                width=100,
                height=47)


# fingerprint similarity heatmap ------------------------------------------

ggPairwiseHeat(fgp=fgs,
               simScore='cosine',
               grporder=c('clu', 'hom'),
               removeControl=TRUE,
               controlGrp=c('wt', 'scr'),
               minCol=NA,
               maxCol=NA,
               onlyHalf='upper',
               scoreSize=2,
               legendOrNo=FALSE,
               labelsOrNo=FALSE,
               width=33,
               height=33,
               exportPath='~/Dropbox/ZFAD/fingerprints/cos_clu.pdf')
