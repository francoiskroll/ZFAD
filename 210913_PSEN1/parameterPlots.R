#####################################################
# ~ ZFAD: make simple parameter plots of exp 210913_PSEN1 for main figure of ZFAD paper ~
#
#
# Francois Kroll 2023
# francois@kroll.be
#####################################################


# packages ----------------------------------------------------------------

# library(devtools)
# install_github('francoiskroll/FramebyFrame')

library(FramebyFrame)
library(here)

source(here('utilities', 'percEffectSize.R'))

### DAY ###

# common settings
width <- 34
height <- 50

dir.create(here('210913_PSEN1', 'pubplots'))



# activityTotalPx ---------------------------------------------------------

ggParameter(pa=c(here('210913_PSEN1', 'bhvparams', 'activityPercentageTimeActive_210913_12.csv'),
                 here('210913_PSEN1', 'bhvparams', 'activityPercentageTimeActive_210913_13.csv')),
            grporder=c('psen1', 'scr'),
            onlyDayorNight='night',
            skipNight0=TRUE,
            colours=c('#fcb505', '#697a87'),
            fainterExp=TRUE,
            faintMax=0.35,
            ymin=NA,
            ymax=NA,
            dotSize=0.3,
            connectMean=TRUE,
            legendOrNo=FALSE,
            xtextOrNo=FALSE,
            ynameOrNo=FALSE,
            yunitOrNo=TRUE,
            splitBy='experiment',
            titleOrNo=FALSE,
            statsOrNo=FALSE,
            width=width,
            height=height,
            exportPath=here('210913_PSEN1', 'pubplots', 'actTimeActive_night.pdf'))

percEffectSize(pa=c(here('210913_PSEN1', 'bhvparams', 'activityPercentageTimeActive_210913_12.csv'),
                    here('210913_PSEN1', 'bhvparams', 'activityPercentageTimeActive_210913_13.csv')),
               lmePath=here('LMEreports', 'psen1_f0_LME.csv'),
               grporder=c('psen1', 'scr'),
               dayornight='night')


# activeboutNum -----------------------------------------------------------

ggParameter(pa=c(here('210913_PSEN1', 'bhvparams', 'activeboutNum_210913_12.csv'),
                 here('210913_PSEN1', 'bhvparams', 'activeboutNum_210913_13.csv')),
            grporder=c('psen1', 'scr'),
            onlyDayorNight='night',
            skipNight0=TRUE,
            colours=c('#fcb505', '#697a87'),
            fainterExp=TRUE,
            faintMax=0.35,
            ymin=NA,
            ymax=NA,
            dotSize=0.3,
            connectMean=TRUE,
            legendOrNo=FALSE,
            xtextOrNo=FALSE,
            ynameOrNo=FALSE,
            yunitOrNo=FALSE,
            splitBy='experiment',
            titleOrNo=FALSE,
            statsOrNo=FALSE,
            width=width,
            height=height,
            exportPath=here('210913_PSEN1', 'pubplots', 'actboutNum_night.pdf'))

percEffectSize(pa=c(here('210913_PSEN1', 'bhvparams', 'activeboutNum_210913_12.csv'),
                    here('210913_PSEN1', 'bhvparams', 'activeboutNum_210913_13.csv')),
               lmePath=here('LMEreports', 'psen1_f0_LME.csv'),
               grporder=c('psen1', 'scr'),
               dayornight='night')



# sunset startle ----------------------------------------------------------

ggParameter(pa=c(here('210913_PSEN1', 'bhvparams', 'activitySunsetStartle_210913_12.csv'),
                 here('210913_PSEN1', 'bhvparams', 'activitySunsetStartle_210913_13.csv')),
            grporder=c('psen1', 'scr'),
            onlyDayorNight='night',
            skipNight0=TRUE,
            colours=c('#fcb505', '#697a87'),
            fainterExp=TRUE,
            faintMax=0.35,
            ymin=NA,
            ymax=NA,
            dotSize=0.3,
            connectMean=TRUE,
            legendOrNo=FALSE,
            xtextOrNo=FALSE,
            ynameOrNo=FALSE,
            yunitOrNo=TRUE,
            splitBy='experiment',
            titleOrNo=FALSE,
            statsOrNo=FALSE,
            width=width,
            height=height,
            exportPath=here('210913_PSEN1', 'pubplots', 'actStartle_night.pdf'))

percEffectSize(pa=c(here('210913_PSEN1', 'bhvparams', 'activitySunsetStartle_210913_12.csv'),
                    here('210913_PSEN1', 'bhvparams', 'activitySunsetStartle_210913_13.csv')),
               lmePath=here('LMEreports', 'psen1_f0_LME.csv'),
               grporder=c('psen1', 'scr'),
               dayornight='night')



# total sleep -------------------------------------------------------------

ggParameter(pa=c(here('210913_PSEN1', 'bhvparams', 'sleepHours_210913_12.csv'),
                 here('210913_PSEN1', 'bhvparams', 'sleepHours_210913_13.csv')),
            grporder=c('psen1', 'scr'),
            onlyDayorNight='night',
            skipNight0=TRUE,
            colours=c('#fcb505', '#697a87'),
            fainterExp=TRUE,
            faintMax=0.35,
            ymin=NA,
            ymax=NA,
            dotSize=0.3,
            connectMean=TRUE,
            legendOrNo=FALSE,
            xtextOrNo=FALSE,
            ynameOrNo=FALSE,
            yunitOrNo=TRUE,
            splitBy='experiment',
            titleOrNo=FALSE,
            statsOrNo=FALSE,
            width=width,
            height=height,
            exportPath=here('210913_PSEN1', 'pubplots', 'sleepTot_night.pdf'))

percEffectSize(pa=c(here('210913_PSEN1', 'bhvparams', 'sleepHours_210913_12.csv'),
                    here('210913_PSEN1', 'bhvparams', 'sleepHours_210913_13.csv')),
               lmePath=here('LMEreports', 'psen1_f0_LME.csv'),
               grporder=c('psen1', 'scr'),
               dayornight='night')

