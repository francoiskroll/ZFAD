###################################################
# ~~~ ZFAD: predictive pharmacology for PSEN2 ~~~
#
# using original PSEN2 experiments 210907_PSEN2
#
# Francois Kroll 2022
# francois@kroll.be
###################################################

# will replicate in code what predPharma Shiny app does
# we are interested in drug rescue, but let's also calculate enrichments

## ran again 25/05/2023 because updated pixelAdjust


# packages/functions ------------------------------------------------------

# library(devtools)
# install_github('francoiskroll/FramebyFrame')
library(FramebyFrame)

library(here)
library(dplyr)
library(tidyr)
library(tibble)


# sorry for hardcoded paths, kept predPharma Shiny app project separate
source(here('~/Dropbox/predPharma/drawEnrich_v5.R'))
source(here('~/Dropbox/predPharma/ggEnrich.R'))
source(here('~/Dropbox/predPharma/gglegacyFingerprint.R'))
source(here('~/Dropbox/predPharma/paramsFromMid.R'))

dir.create(here('psen2predictRescue', 'pubplots'))


# convert PSEN2 RAWs data to middur ---------------------------------------
# using adjusted-pixel data

# commented out as time-consuming
# only need to be run once

# mid12 <- rawToMiddur(ffpath=here('210907_PSEN2', '210907_12_RAWsadjusted.csv'),
#                      zebpath=here('210907_PSEN2', '210907_12_13_PSEN2.xls'),
#                      freezing=3,
#                      burst=200,
#                      exportOrNo=TRUE)
# 
# mid13 <- rawToMiddur(ffpath=here('210907_PSEN2', '210907_13_RAWsadjusted.csv'),
#                      zebpath=here('210907_PSEN2', '210907_12_13_PSEN2.xls'),
#                      freezing=3,
#                      burst=200,
#                      exportOrNo=TRUE)

## re-import middur data so can skip this step next time ##
mid12 <- read.csv(here('210907_PSEN2', '210907_12_middur.csv'))
mid13 <- read.csv(here('210907_PSEN2', '210907_13_middur.csv'))


# prepare legacy fingerprint ----------------------------------------------

lfp12 <- legacyFingerprintMid(mid=mid12,
                              genopath=here('210907_PSEN2', '210907_12genotype.txt'),
                              treGrp='psen2',
                              conGrp='scr',
                              nights=c('night1', 'night2'),
                              days=c('day1', 'day2'))

lfp13 <- legacyFingerprintMid(mid=mid13,
                              genopath=here('210907_PSEN2', '210907_13genotype.txt'),
                              treGrp='psen2',
                              conGrp='scr',
                              nights=c('night1', 'night2'),
                              days=c('day1', 'day2'))

# prepare mean fingerprint of the two F0 experiments
# put both fingerprints in a list
lfpL <- list(lfp12, lfp13)
# join using plyr::join_all
lfp <- plyr::join_all(lfpL, by=c('uparam', 'parameter', 'win'), type='left')
# (avoid loading plyr because it replaces functions from dplyr)
colnames(lfp) <- c('uparam', 'parameter', 'win', 'exp1', 'exp2')

lfp <- lfp %>%
  mutate(zsco = rowMeans(select(., starts_with('exp')))) # mean of columns that start with 'exp'


# rank drugs --------------------------------------------------------------

vdbr <- rankDrugDb(legacyFgp=lfp,
                   dbPath='~/Dropbox/predPharma/drugDb.csv', # hard-coded path, sorry!
                   metric='cosine')
write.csv(vdbr, file=here('psen2predictRescue', 'f0mean_rankedDrugs.csv'), row.names=FALSE)


# prepare fingerprints ----------------------------------------------------
# of drugs predicted to rescue PSEN2 phenotype

# about onlyGrp below:
# NA would plot exp1, exp2, and mean fingerprint; zsco only plots mean fingerprint

ggDrugFgp(drugDb='~/Dropbox/predPharma/drugDb.csv',
          dnames=c('Tinidazole', 'TINIDAZOLE'),
          legacyFgp=lfp, 
          onlyGrp=c('zsco', 'Tinidazole', 'TINIDAZOLE'), 
          colours=c('#78ac63', '#c5dbbc', '#697a87'),
          legendOrNo=FALSE,
          ynameOrNo=FALSE,
          ytextOrNo=TRUE,
          xtextOrNo=FALSE,
          paramNumOrNo=FALSE,
          nightBgOrNo=TRUE,
          ymin=-2,
          ymax=3.68,
          exportOrNo=TRUE,
          exportPath=here('psen2predictRescue', 'pubplots', 'f0mean_tinidazole.pdf'),
          width=90,
          height=33)

ggDrugFgp(drugDb='~/Dropbox/predPharma/drugDb.csv',
          dnames=c('PYRAZINAMIDE', 'Pyrazinamide', 'Pyrazinecarboxamide'),
          legacyFgp=lfp, 
          onlyGrp=c('zsco', 'PYRAZINAMIDE', 'Pyrazinamide', 'Pyrazinecarboxamide'), # NA would plot exp1, exp2, and mean fingerprint; zsco only plots mean fingerprint
          colours=c('#db5072', '#ea96ab', '#f2c1cd', '#697a87'),
          legendOrNo=FALSE,
          ynameOrNo=TRUE,
          ytextOrNo=TRUE,
          xtextOrNo=FALSE,
          nightBgOrNo=TRUE,
          ymin=-1.5,
          ymax=1.5,
          exportOrNo=TRUE,
          exportPath=here('psen2predictRescue', 'pubplots', 'f0mean_pyrazinamide.pdf'),
          width=90,
          height=40)

# note, there is another entry for "oxcarbazepine" but difficult to plot on top
ggDrugFgp(drugDb='~/Dropbox/predPharma/drugDb.csv',
          dnames=c('OXCARBAZEPINE', 'oxcarbazepine'),
          legacyFgp=lfp, 
          onlyGrp=c('zsco', 'OXCARBAZEPINE', 'oxcarbazepine'), # NA would plot exp1, exp2, and mean fingerprint; zsco only plots mean fingerprint
          colours=c('#bbd0f6', '#437ce6', '#697a87'),
          legendOrNo=FALSE,
          ynameOrNo=TRUE,
          ytextOrNo=TRUE,
          xtextOrNo=FALSE,
          nightBgOrNo=TRUE,
          ymin=-2,
          ymax=2.1,
          exportOrNo=TRUE,
          exportPath=here('psen2predictRescue', 'pubplots', 'f0mean_oxcarbazepine.pdf'),
          width=120,
          height=50)

## fenoprofen
# exists:
# "FENOPROFEN" cos=-0.63
# "Fenoprofen calcium salt dihydrate" cos=-0.74
# can plot both, should not be too much

ggDrugFgp(drugDb='~/Dropbox/predPharma/drugDb.csv',
          dnames=c('FENOPROFEN', 'Fenoprofen calcium salt dihydrate'),
          legacyFgp=lfp, 
          onlyGrp=c('mean fgp', 'FENOPROFEN', 'Fenoprofen calcium salt dihydrate'),
          colours=c('#417dcd', '#adc7e9', '#697a87'),
          legendOrNo=FALSE,
          ynameOrNo=FALSE,
          ytextOrNo=TRUE,
          xtextOrNo=TRUE,
          paramNumOrNo=TRUE,
          nightBgOrNo=TRUE,
          ymin=-1.7,
          ymax=2.2,
          exportOrNo=TRUE,
          exportPath=here('psen2predictRescue', 'pubplots', 'f0mean_fenoprofen.pdf'),
          width=90,
          height=35)

## betamethasone
# exists:
# "Betamethasone" cos=-0.74
# "Betamethasone" cos=-0.08; extreme fingerprints, not sure what to do with it
# "Betamethasone Valerate" cos=-0.20; I think same idea as propionate below, added a small fatty acid, this one with 5C, propionate is 3C
# "*(D)Betamethasone" cos=-0.53; not sure if this is the normal isomer or not
# "BETAMETHASONE 17,21-DIPROPIONATE" cos=-0.56
# all seem valid to include, but it is too many
# I think drop valerate/dipropionate. Can mention in Results or Legend.
# I think drop Betamethasone.1 too, I think two drugs behaving +/- same way is really the maximum to grasp

# full version for reference:
ggDrugFgp(drugDb='~/Dropbox/predPharma/drugDb.csv',
          dnames=c('Betamethasone', 'Betamethasone Valerate', '*(D)Betamethasone', 'BETAMETHASONE 17,21-DIPROPIONATE'),
          legacyFgp=lfp, 
          onlyGrp=c('zsco', 'Betamethasone', 'Betamethasone.1', 'Betamethasone Valerate', '*(D)Betamethasone', 'BETAMETHASONE 17,21-DIPROPIONATE'),
          colours=NA,
          legendOrNo=TRUE,
          ynameOrNo=TRUE,
          ytextOrNo=TRUE,
          xtextOrNo=FALSE,
          nightBgOrNo=TRUE,
          ymin=-2,
          ymax=2.1,
          exportOrNo=FALSE,
          exportPath=NA,
          width=NA,
          height=NA)

ggDrugFgp(drugDb='~/Dropbox/predPharma/drugDb.csv',
          dnames=c('Betamethasone', '*(D)Betamethasone'),
          legacyFgp=lfp, 
          onlyGrp=c('zsco', 'Betamethasone', '*(D)Betamethasone'), # ! second Betamethasone gets called Betamethasone.1, here we do not mention it so gets excluded
          colours=c('#fee29c', '#fcb505', '#697a87'),
          legendOrNo=FALSE,
          ynameOrNo=FALSE,
          ytextOrNo=TRUE,
          xtextOrNo=TRUE,
          paramNumOrNo=TRUE,
          nightBgOrNo=TRUE,
          ymin=-1.7,
          ymax=2.23,
          exportOrNo=TRUE,
          exportPath=here('psen2predictRescue', 'pubplots', 'f0mean_betamethasone.pdf'),
          width=90,
          height=35)

ggDrugFgp(drugDb='~/Dropbox/predPharma/drugDb.csv',
          dnames=c('Betamethasone', 'BETAMETHASONE 17,21-DIPROPIONATE'),
          legacyFgp=lfp, 
          onlyGrp=c('mean fgp', 'Betamethasone', 'BETAMETHASONE 17,21-DIPROPIONATE'), # ! second Betamethasone gets called Betamethasone.1, here we do not mention it so gets excluded
          colours=c('#fcb505', '#fee29c', '#697a87'),
          legendOrNo=FALSE,
          ynameOrNo=FALSE,
          ytextOrNo=TRUE,
          xtextOrNo=TRUE,
          paramNumOrNo=TRUE,
          nightBgOrNo=TRUE,
          ymin=-1.7,
          ymax=2.23,
          exportOrNo=TRUE,
          exportPath=here('psen2predictRescue', 'pubplots', 'f0mean_betamethasone_v2.pdf'),
          width=90,
          height=35)


ggDrugFgp(drugDb='~/Dropbox/predPharma/drugDb.csv',
          dnames=c('Diflunisal', 'DIFLUNISAL'),
          legacyFgp=lfp, 
          onlyGrp=c('Diflunisal', 'DIFLUNISAL'),
          colours=c('#fee29c', '#fcb505', '#697a87'),
          legendOrNo=FALSE,
          ynameOrNo=FALSE,
          ytextOrNo=TRUE,
          xtextOrNo=TRUE,
          paramNumOrNo=TRUE,
          nightBgOrNo=TRUE,
          ymin=-1.7,
          ymax=2.23,
          exportOrNo=TRUE,
          exportPath=here('230315_sorl1Diflunisal', 'plots', 'diflunisal.pdf'),
          width=90,
          height=35)
