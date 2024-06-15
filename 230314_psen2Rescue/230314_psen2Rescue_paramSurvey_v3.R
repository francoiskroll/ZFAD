#####################################################
# ~ ZFAD: psen2Rescue, plot to survey rescue or "side effect" parameter by parameter ~
#
# Francois Kroll 2023
# francois@kroll.be
#####################################################

# idea is to have one plot for each drug
# and plot Z-scores parameter by parameter
# essentially just another version of the fingerprint (exactly the same datapoints, in fact)
# but for each parameter, having next to each other: Z-score of PSEN2 + DMSO // Z-score of PSEN2 + drug
# 1/ if Z-score decreases when adding drug, this is a *rescue*
# 2/ if Z-score stays around the same when adding drug, this is a *symptom not addressed*
# 3/ if Z-score was low and then increases when adding drug, this is a *side effect*
# could colour categories


# packages ----------------------------------------------------------------

library(here)

# library(devtools)
# install_github('francoiskroll/FramebyFrame')
library(FramebyFrame)

library(dplyr)
library(tidyr)
library(tibble)
library(data.table)
library(ggplot2)
library(ggbeeswarm)

# in v2 of this script, experimented with a few different solutions
# I think best result was to pool boxes and use LME statistics from both boxes together
# here, going more directly to this

source(here('230314_psen2Rescue', 'paramSurvey.R'))
source(here('230314_psen2Rescue', 'paramSurvey_v2.R'))


# calculate clutch fingerprints -------------------------------------------

# reference (origin of the multidimensional space) is SCR + DMSO
fgcl <- calculateFingerprint(paDir=here('230314_psen2Rescue', 'bhvparams'),
                             controlGrp='SCR.DMSO',
                             singleFish=FALSE,
                             grporder=c('SCR.DMSO', 'PSEN2.DMSO', 'PSEN2.Tinidazole', 'PSEN2.Betamethasone', 'PSEN2.Fenoprofen'),
                             skipNight0=TRUE)


### before pixelAdjust
fgclBA <- calculateFingerprint(paDir=here('230314_psen2Rescue', 'bhvparams_beforeAdjust'),
                               controlGrp='SCR.DMSO',
                               singleFish=FALSE,
                               grporder=c('SCR.DMSO', 'PSEN2.DMSO', 'PSEN2.Tinidazole', 'PSEN2.Betamethasone', 'PSEN2.Fenoprofen'),
                               skipNight0=TRUE)



# generate paramSurvey plots ----------------------------------------------

# see script paramSurvey for the function


### before pixelAdjust
ggParamSurvey(fgcl=fgclBA,
              lmePath=here('230314_psen2Rescue', 'psen2Rescue_LME_beforeAdjust.csv'),
              conGrp='PSEN2.DMSO',
              treGrp='PSEN2.Tinidazole',
              ytextOrNo=TRUE,
              exportPath=here('230314_psen2Rescue', 'pubplots', 'paramSurvey_tinidazole_v2beforeAdj.pdf'),
              width=60,
              height=105)

ggParamSurvey(fgcl=fgclBA,
              lmePath=here('230314_psen2Rescue', 'psen2Rescue_LME_beforeAdjust.csv'),
              conGrp='PSEN2.DMSO',
              treGrp='PSEN2.Fenoprofen',
              ytextOrNo=FALSE,
              exportPath=here('230314_psen2Rescue', 'pubplots', 'paramSurvey_fenoprofen_v2beforeAdj.pdf'),
              width=45.8,
              height=105)

ggParamSurvey(fgcl=fgclBA,
              lmePath=here('230314_psen2Rescue', 'psen2Rescue_LME_beforeAdjust.csv'),
              conGrp='PSEN2.DMSO',
              treGrp='PSEN2.Betamethasone',
              ytextOrNo=FALSE,
              exportPath=here('230314_psen2Rescue', 'pubplots', 'paramSurvey_betamethasone_v2beforeAdj.pdf'),
              width=45.8,
              height=105)

# temporary for Montpellier poster
ggParamSurvey(fgcl=fgclBA,
              lmePath=here('230314_psen2Rescue', 'psen2Rescue_LME_beforeAdjust.csv'),
              conGrp='PSEN2.DMSO',
              treGrp='PSEN2.Betamethasone',
              ytextOrNo=FALSE,
              exportPath=here('230314_psen2Rescue', 'pubplots', 'paramSurvey_betamethasone_v2beforeAdjTMP.pdf'),
              width=45.8,
              height=130)
