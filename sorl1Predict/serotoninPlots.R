#####################################################
# ~ ZFAD: plots for serotonin transporter prediction of SORL1 ~
#
#
# Francois Kroll 2023
# francois@kroll.be
#####################################################


# packages/functions ------------------------------------------------------

library(FramebyFrame)
library(here)

source('~/Dropbox/predPharma/drawEnrich_v5.R')
source('~/Dropbox/predPharma/ggEnrich.R')
source('~/Dropbox/predPharma/paramsFromMid.R')
source('~/Dropbox/predPharma/gglegacyFingerprint.R')

dir.create(here('sorl1Predict', 'plots'))


# calculate legacy fingerprint --------------------------------------------

# first need middur files
# probably created previously, but can re-create them here as it does not take too long (few minutes)
# mid14 <- rawToMiddur(ffpath=here('220531_SORL1', '220531_14_RAWs.csv'),
#                      zebpath=here('220531_SORL1', '220531_14_15_sorl1.xls'),
#                      freezing=3,
#                      burst=200,
#                      exportOrNo=TRUE)
# 
# mid15 <- rawToMiddur(ffpath=here('220531_SORL1', '220531_15_RAWs.csv'),
#                      zebpath=here('220531_SORL1', '220531_14_15_sorl1.xls'),
#                      freezing=3,
#                      burst=200,
#                      exportOrNo=TRUE)

## re-import middur data so can skip this step next time ##
mid14 <- read.csv(here('220531_SORL1', '220531_14_middur.csv'))
mid15 <- read.csv(here('220531_SORL1', '220531_15_middur.csv'))


# prepare legacy fingerprint ----------------------------------------------

lfp14 <- legacyFingerprintMid(mid=mid14,
                              genopath=here('220531_SORL1', '220531_14genotype.txt'),
                              treGrp='sorl1',
                              conGrp='scr',
                              nights=c('night1', 'night2'),
                              days=c('day1', 'day2'))

lfp15 <- legacyFingerprintMid(mid=mid15,
                              genopath=here('220531_SORL1', '220531_15genotype.txt'),
                              treGrp='sorl1',
                              conGrp='scr',
                              nights=c('night1', 'night2'),
                              days=c('day1', 'day2'))

# prepare mean fingerprint of the two F0 experiments
# put both fingerprints in a list
lfpL <- list(lfp14, lfp15)
# join using plyr::join_all
lfp <- plyr::join_all(lfpL, by=c('uparam', 'parameter', 'win'), type='left')
# (avoid loading plyr because it replaces functions from dplyr)

lfp <- lfp %>%
  mutate(zsco = rowMeans(select(., matches('^\\d{6}_\\d{2}')))) # mean of columns that have name YYMMDD_BX
# regex is start / 6 digits (YYMMDD) / _ / 2 digits (BX)
# column called zsco is now mean of Z-scores



# rank drugs --------------------------------------------------------------

vdbr <- rankDrugDb(legacyFgp=lfp,
                   dbPath='~/Dropbox/predPharma/drugDb.csv', # hard-coded path, sorry!
                   metric='cosine')
write.csv(vdbr, file=here('sorl1Predict', 'f0mean_rankedDrugs.csv'), row.names=FALSE)




# barcode TTD SLC6A4 ------------------------------------------------------

ggBarcode(vdbr=vdbr,
          namesPath='~/Dropbox/predPharma/compounds.csv',
          annotationPath='~/Dropbox/predPharma/TTDtargets.csv',
          annotation='TTDtargets',
          testAnnotation='T27812',
          minScore=NA,
          barwidth1=0.5,
          barwidth2=25,
          exportPath=here('sorl1Predict', 'plots', 'bc_TTDslc6a4.pdf'),
          width=100,
          height=25)



# barcode Depression ------------------------------------------------------

ggBarcode(vdbr=vdbr,
          namesPath='~/Dropbox/predPharma/compounds.csv',
          annotationPath='~/Dropbox/predPharma/TTDindications.csv',
          annotation='indications',
          testAnnotation='Depression',
          minScore=NA,
          barwidth1=0.5,
          barwidth2=25,
          exportPath=here('sorl1Predict', 'plots', 'bc_TTDdepression.pdf'),
          width=100,
          height=25)



# draws plot --------------------------------------------------------------

# SLC6A4
ggDraws(vdbr=vdbr,
        namesPath='~/Dropbox/predPharma/compounds.csv',
        annotationPath='~/Dropbox/predPharma/TTDtargets.csv',
        annotation='TTDtargets',
        testAnnotation='T27812',
        ndraws=100000,
        minScore=NA,
        whichRank='rankeq',
        exportPath=here('sorl1Predict', 'plots', 'draws_TTDslc6a4.pdf'),
        width=55,
        height=42)

# for KEGG serotonergic synapse
ggDraws(vdbr=vdbr,
        namesPath='~/Dropbox/predPharma/compounds.csv',
        annotationPath='~/Dropbox/predPharma/TTDkegg.csv',
        annotation='KEGG',
        testAnnotation='hsa04726',
        ndraws=100000,
        minScore=NA,
        whichRank='rankeq',
        exportPath=here('sorl1Predict', 'plots', 'draws_serotonergicsynapse.pdf'),
        width=55,
        height=42)

# dopamine D4 receptor
ggDraws(vdbr=vdbr,
        namesPath='~/Dropbox/predPharma/compounds.csv',
        annotationPath='~/Dropbox/predPharma/TTDtargets.csv',
        annotation='TTDtargets',
        testAnnotation='T24983',
        ndraws=100000,
        minScore=NA,
        whichRank='rankeq',
        exportPath=here('sorl1Predict', 'plots', 'draws_TTDdopamineD4receptor.pdf'),
        width=55,
        height=42)

# dopamine D2 receptor
ggDraws(vdbr=vdbr,
        namesPath='~/Dropbox/predPharma/compounds.csv',
        annotationPath='~/Dropbox/predPharma/TTDtargets.csv',
        annotation='TTDtargets',
        testAnnotation='T67162',
        ndraws=100000,
        minScore=NA,
        whichRank='rankeq',
        exportPath=here('sorl1Predict', 'plots', 'draws_TTDdopamineD2receptor.pdf'),
        width=55,
        height=42)


# common settings ---------------------------------------------------------

namesPath <- here('annotateDrugDb', 'compounds.csv')
annotationDir <- here('annotateDrugDb/')
whichRank <- 'rankeq'
minScore <- 900
minNex <- 3
ndraws <- 100000
alphaThr <- 0.05
maxPval <- 0.2
statsDir <- here('predPharma', 'sorl1All/')

dir.create(statsDir, showWarnings=FALSE)



# some plots --------------------------------------------------------------



vdbr <- rankDrugDb(legacyFgp=lfp,
                   dbPath=here('annotateDrugDb', 'drugDb.csv'), 
                   metric='cosine')

ggBarcode(vdbr=vdbr,
          namesPath=namesPath,
          annotationPath=here('annotateDrugDb', 'zebrafishSTITCH.csv'),
          annotation='zebrafishSTITCH',
          testAnnotation='ENSDARP00000109199',
          minScore=minScore,
          exportPath=here('predPharma', 'sorl1All', 'f0mean_hsp701.pdf'),
          width=200,
          height=100)

ggBarcode(vdbr=vdbr,
          namesPath=namesPath,
          annotationPath=here('annotateDrugDb', 'zebrafishSTITCH.csv'),
          annotation='zebrafishSTITCH',
          testAnnotation='ENSDARP00000052988',
          minScore=minScore,
          exportPath=here('predPharma', 'sorl1All', 'f0mean_ache.pdf'),
          width=200,
          height=100)

ggBarcode(vdbr=vdbr,
          namesPath=namesPath,
          annotationPath=here('annotateDrugDb', 'zebrafishSTITCH.csv'),
          annotation='zebrafishSTITCH',
          testAnnotation='ENSDARP00000070269',
          minScore=minScore,
          exportPath=here('predPharma', 'sorl1All', 'f0mean_f2.pdf'),
          width=200,
          height=100)

ggBarcode(vdbr=vdbr,
          namesPath=namesPath,
          annotationPath=here('annotateDrugDb', 'zebrafishSTITCH.csv'),
          annotation='zebrafishSTITCH',
          testAnnotation='ENSDARP00000085242',
          minScore=minScore,
          exportPath=here('predPharma', 'sorl1All', 'f0mean_kcnh6a.pdf'),
          width=200,
          height=100)

ggBarcode(vdbr=vdbr,
          namesPath=namesPath,
          annotationPath=here('annotateDrugDb', 'zebrafishSTITCH.csv'),
          annotation='zebrafishSTITCH',
          testAnnotation='ENSDARP00000085307',
          minScore=minScore,
          exportPath=here('predPharma', 'sorl1All', 'f0mean_kcnh7.pdf'),
          width=200,
          height=100)

ggBarcode(vdbr=vdbr,
          namesPath=namesPath,
          annotationPath=here('annotateDrugDb', 'zebrafishSTITCH.csv'),
          annotation='zebrafishSTITCH',
          testAnnotation='ENSDARP00000009558',
          minScore=minScore,
          exportPath=here('predPharma', 'sorl1All', 'f0mean_kcnh2a.pdf'),
          width=200,
          height=100)

ggBarcode(vdbr=vdbr,
          namesPath=namesPath,
          annotationPath=here('annotateDrugDb', 'zebrafishSTITCH.csv'),
          annotation='zebrafishSTITCH',
          testAnnotation='ENSDARP00000081336',
          minScore=minScore,
          exportPath=here('predPharma', 'sorl1All', 'f0mean_slc6a4a.pdf'),
          width=200,
          height=100)

### stable mean ###
lfp <- legacyFingerprintMEAN(matPaths=c(here('220316_sorl1Stable', 'legacyMiddur', '220316_14', '220316_14.mat'),
                                        here('220316_sorl1Stable', 'legacyMiddur', '220316_15', '220316_15.mat')),
                             conGrp='wt',
                             treGrp='hom',
                             days=c(2,3),
                             nights=c(2,3))
vdbr <- rankDrugDb(legacyFgp=lfp,
                   dbPath=here('annotateDrugDb', 'drugDb.csv'), 
                   metric='cosine')
ggBarcode(vdbr=vdbr,
          namesPath=namesPath,
          annotationPath=here('annotateDrugDb', 'zebrafishSTITCH.csv'),
          annotation='zebrafishSTITCH',
          testAnnotation='ENSDARP00000109199',
          minScore=minScore,
          exportPath=here('predPharma', 'sorl1All', 'stamean_hsp701.pdf'),
          width=200,
          height=100)
ggBarcode(vdbr=vdbr,
          namesPath=namesPath,
          annotationPath=here('annotateDrugDb', 'zebrafishSTITCH.csv'),
          annotation='zebrafishSTITCH',
          testAnnotation='ENSDARP00000052988',
          minScore=minScore,
          exportPath=here('predPharma', 'sorl1All', 'stamean_ache.pdf'),
          width=200,
          height=100)
ggBarcode(vdbr=vdbr,
          namesPath=namesPath,
          annotationPath=here('annotateDrugDb', 'zebrafishSTITCH.csv'),
          annotation='zebrafishSTITCH',
          testAnnotation='ENSDARP00000070269',
          minScore=minScore,
          exportPath=here('predPharma', 'sorl1All', 'stamean_f2.pdf'),
          width=200,
          height=100)