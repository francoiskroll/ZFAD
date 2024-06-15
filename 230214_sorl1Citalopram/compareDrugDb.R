#####################################################
# ~ ZFAD: does citalopram treatment correlate with data in drug database? ~
#
#
# Francois Kroll 2024
# francois@kroll.be
#####################################################


# packages ----------------------------------------------------------------

library(FramebyFrame)

library(here)

# sorry for hardcoded paths, kept predPharma Shiny app project separate
source(here('~/Dropbox/predPharma/drawEnrich_v5.R'))
source(here('~/Dropbox/predPharma/ggEnrich.R'))
source(here('~/Dropbox/predPharma/gglegacyFingerprint.R'))
source(here('~/Dropbox/predPharma/paramsFromMid.R'))


# does citalopram treatment correlate with citalopram in drugdb? ----------

# convert to middur dataset
# mid14 <- rawToMiddur(ffpath=here('230214_sorl1Citalopram/230214_14_RAWs.csv'),
#                      freezing=3,
#                      burst=200,
#                      exportOrNo=TRUE)
# wrote 230214_sorl1Citalopram/230214_14_middur.csv

# mid15 <- rawToMiddur(ffpath=here('230214_sorl1Citalopram/230214_15_RAWs.csv'),
#                      freezing=3,
#                      burst=200,
#                      exportOrNo=TRUE)
# wrote 230214_sorl1Citalopram/230214_15_middur.csv

mid14 <- read.csv(here('230214_sorl1Citalopram/230214_14_middur.csv'))
mid15 <- read.csv(here('230214_sorl1Citalopram/230214_15_middur.csv'))

# calculate legacy fingerprint for scr + 1 uM
lfp14_1um <- legacyFingerprintMid(mid=mid14,
                                  genopath=here('230214_sorl1Citalopram/230214_14genotype.txt'),
                                  treGrp='scr_1',
                                  conGrp='scr_0',
                                  nights=c('night1', 'night2'),
                                  days=c('day1', 'day2'))

lfp15_1um <- legacyFingerprintMid(mid=mid15,
                                  genopath=here('230214_sorl1Citalopram/230214_15_1031genotype.txt'),
                                  treGrp='scr_1',
                                  conGrp='scr_0',
                                  nights=c('night1', 'night2'),
                                  days=c('day1', 'day2'))

# calculate legacy fingerprint for scr + 10 uM
lfp14_10um <- legacyFingerprintMid(mid=mid14,
                                   genopath=here('230214_sorl1Citalopram/230214_14genotype.txt'),
                                   treGrp='scr_10',
                                   conGrp='scr_0',
                                   nights=c('night1', 'night2'),
                                   days=c('day1', 'day2'))

lfp15_10um <- legacyFingerprintMid(mid=mid15,
                                   genopath=here('230214_sorl1Citalopram/230214_15_1031genotype.txt'),
                                   treGrp='scr_10',
                                   conGrp='scr_0',
                                   nights=c('night1', 'night2'),
                                   days=c('day1', 'day2'))

# for each dose, prepare mean fingerprint of the two clutches

## 1 uM
# put both 1 uM fingerprints in a list
lfpL <- list(lfp14_1um, lfp15_1um)
# join using plyr::join_all
lfp1um <- plyr::join_all(lfpL, by=c('uparam', 'parameter', 'win'), type='left')
# (avoid loading plyr because it replaces functions from dplyr)
colnames(lfp1um) <- c('uparam', 'parameter', 'win', 'exp1', 'exp2')

lfp1um <- lfp1um %>%
  mutate(zsco = rowMeans(select(., starts_with('exp')))) # mean of columns that start with 'exp'

## 10 uM
# put both 10 uM fingerprints in a list
lfpL <- list(lfp14_10um, lfp15_10um)
# join using plyr::join_all
lfp10um <- plyr::join_all(lfpL, by=c('uparam', 'parameter', 'win'), type='left')
# (avoid loading plyr because it replaces functions from dplyr)
colnames(lfp10um) <- c('uparam', 'parameter', 'win', 'exp1', 'exp2')

lfp10um <- lfp10um %>%
  mutate(zsco = rowMeans(select(., starts_with('exp')))) # mean of columns that start with 'exp'

## now put both doses in one
lfp <- data.frame(lfp1um$uparam, lfp1um$parameter, lfp1um$win, lfp1um$zsco, lfp10um$zsco)
colnames(lfp) <- c('uparam', 'parameter', 'win', 'zsco1um', 'zsco10um')

ggDrugFgp(drugDb='~/Dropbox/predPharma/drugDb.csv',
          dnames=c('citalopram', 'Citalopram Hydrobromide'), # there is also "ESCITALOPRAM OXALATE", but gives a quite different fingerprint
          legacyFgp=lfp, 
          onlyGrp=c('zsco1um', 'zsco10um', 'citalopram', 'Citalopram Hydrobromide', 'Citalopram Hydrobromide.1'), 
          colours=c('#5a6974', '#84949f', '#b4bdc4', '#78ac63', '#bbd4ae'),
          legendOrNo=FALSE,
          ynameOrNo=FALSE,
          ytextOrNo=TRUE,
          xtextOrNo=TRUE,
          paramNumOrNo=TRUE,
          nightBgOrNo=TRUE,
          ymin=-2.3,
          ymax=7.2,
          exportOrNo=TRUE,
          exportPath=here('230214_sorl1Citalopram', 'plots', 'fgp_citalopramVsDb.pdf'),
          width=103.5,
          height=33)


# eLife reviews: which other SSRI to try? ---------------------------------

# Jason by email recommended fluoxetine, fluvoxamine, paroxetine, sertraline
# will only try one

### FLUOXETINE ###
# in drugdb:
# fluoxetine
# fluoxetine

# FLUOXETINE
# FLUOXETINE

# "(R)"-fluoxetine
# (S)-fluoxetine
# S-(+)-Fluoxetine Hydrochloride
# *(D)R-(-)-Fluoxetine Hydrochloride

# Fluoxetine hydrochloride
# Fluoxetine Hydrochloride

# total 10

# will compare each time with citalopram 10 uM
ggDrugFgp(drugDb='~/Dropbox/predPharma/drugDb.csv',
          dnames=c('"(R)"-fluoxetine', 'fluoxetine', '(S)-fluoxetine', 'S-(+)-Fluoxetine Hydrochloride',
                   '*(D)R-(-)-Fluoxetine Hydrochloride', 'Fluoxetine hydrochloride', 'Fluoxetine Hydrochloride', 'FLUOXETINE'),
          legacyFgp=lfp, 
          onlyGrp=c('zsco10um',
                    'fluoxetine', 'fluoxetine.1',
                    'FLUOXETINE', 'FLUOXETINE.1',
                    '"(R)"-fluoxetine', '(S)-fluoxetine', 'S-(+)-Fluoxetine Hydrochloride', '*(D)R-(-)-Fluoxetine Hydrochloride',
                    'Fluoxetine hydrochloride', 'Fluoxetine Hydrochloride'), 
          colours=c(rep('#b3bcc3', 10), '#EE7163'),
          legendOrNo=TRUE,
          ynameOrNo=FALSE,
          ytextOrNo=TRUE,
          xtextOrNo=TRUE,
          paramNumOrNo=TRUE,
          nightBgOrNo=TRUE,
          ymin=-2.3,
          ymax=7.2,
          exportOrNo=TRUE,
          exportPath=here('230214_sorl1Citalopram', 'plots', 'fgp_fluoxetine.pdf'),
          width=103.5,
          height=33)


### FLUVOXAMINE ###
# in drugdb:
# fluvoxamine

# Fluvoxamine Maleate
# Fluvoxamine Maleate

# Fluvoxamine maleate
ggDrugFgp(drugDb='~/Dropbox/predPharma/drugDb.csv',
          dnames=c('fluvoxamine', 'Fluvoxamine Maleate', 'Fluvoxamine maleate'),
          legacyFgp=lfp, 
          onlyGrp=c('zsco10um',
                    'fluvoxamine', 'Fluvoxamine Maleate', 'Fluvoxamine Maleate.1', 'Fluvoxamine maleate'), 
          colours=NA,
          legendOrNo=TRUE,
          ynameOrNo=FALSE,
          ytextOrNo=TRUE,
          xtextOrNo=TRUE,
          paramNumOrNo=TRUE,
          nightBgOrNo=TRUE,
          ymin=-2.3,
          ymax=7.2,
          exportOrNo=TRUE,
          exportPath=here('230214_sorl1Citalopram', 'plots', 'fgp_fluvoxamine.pdf'),
          width=103.5,
          height=33)


### PAROXETINE ###
# PAROXETINE HYDROCHLORIDE
# paroxetine
# Paroxetine Hydrochloride Hemihydrate (Mw = 374.83)
# Paroxetine Hydrochloride
ggDrugFgp(drugDb='~/Dropbox/predPharma/drugDb.csv',
          dnames=c('PAROXETINE HYDROCHLORIDE', 'paroxetine', 'Paroxetine Hydrochloride Hemihydrate (Mw = 374.83)', 'Paroxetine Hydrochloride'),
          legacyFgp=lfp, 
          onlyGrp=c('zsco10um',
                    'PAROXETINE HYDROCHLORIDE', 'paroxetine', 'Paroxetine Hydrochloride Hemihydrate (Mw = 374.83)', 'Paroxetine Hydrochloride'), 
          colours=c(rep('#b3bcc3', 4), '#EE7163'),
          legendOrNo=TRUE,
          ynameOrNo=FALSE,
          ytextOrNo=TRUE,
          xtextOrNo=TRUE,
          paramNumOrNo=TRUE,
          nightBgOrNo=TRUE,
          ymin=-2.3,
          ymax=7.2,
          exportOrNo=TRUE,
          exportPath=here('230214_sorl1Citalopram', 'plots', 'fgp_paroxetine.pdf'),
          width=103.5,
          height=33)


### SERTRALINE ###
# Sertraline hydrochloride
# Sertraline
# SERTRALINE HYDROCHLORIDE
ggDrugFgp(drugDb='~/Dropbox/predPharma/drugDb.csv',
          dnames=c('Sertraline hydrochloride', 'Sertraline', 'SERTRALINE HYDROCHLORIDE'),
          legacyFgp=lfp, 
          onlyGrp=c('zsco10um',
                    'Sertraline hydrochloride', 'Sertraline', 'SERTRALINE HYDROCHLORIDE'), 
          colours=c(rep('#b3bcc3', 3), '#EE7163'),
          legendOrNo=TRUE,
          ynameOrNo=FALSE,
          ytextOrNo=TRUE,
          xtextOrNo=TRUE,
          paramNumOrNo=TRUE,
          nightBgOrNo=TRUE,
          ymin=-2.3,
          ymax=7.2,
          exportOrNo=TRUE,
          exportPath=here('230214_sorl1Citalopram', 'plots', 'fgp_sertraline.pdf'),
          width=103.5,
          height=33)

# most reproducible is clearly fluvoxamine