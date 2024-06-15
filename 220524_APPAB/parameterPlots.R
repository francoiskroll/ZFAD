#####################################################
# ~ ZFAD: make simple parameter plots of exp 220601_CLU for main figure of ZFAD paper ~
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

source(here('utilities/percEffectSize.R'))


### DAY ###

# common settings
width <- 38
height <- 50

dir.create(here('220524_APPAB', 'pubplots'))


ggParameter(pa=c(here('220524_APPAB', 'bhvparams', 'activityPercentageTimeActive_220524_14.csv'),
                 here('220524_APPAB', 'bhvparams', 'activityPercentageTimeActive_220524_15.csv')),
            grporder=c('appab', 'scr4x'),
            onlyDayorNight='day',
            skipNight0=TRUE,
            colours=c('#E94735', '#697a87'),
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
            exportPath=here('220524_APPAB', 'pubplots', 'actTimeActive_day.pdf'))

ggParameter(pa=c(here('220524_APPAB', 'bhvparams', 'activeboutLength_220524_14.csv'),
                 here('220524_APPAB', 'bhvparams', 'activeboutLength_220524_15.csv')),
            grporder=c('appab', 'scr4x'),
            onlyDayorNight='day',
            skipNight0=TRUE,
            colours=c('#E94735', '#697a87'),
            fainterExp=TRUE,
            faintMax=0.35,
            ymin=0.115,
            ymax=0.227,
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
            exportPath=here('220524_APPAB', 'pubplots', 'actboutLength_day.pdf'))


ggParameter(pa=c(here('220524_APPAB', 'bhvparams', 'activeboutLength_220524_14.csv'),
                 here('220524_APPAB', 'bhvparams', 'activeboutLength_220524_15.csv')),
            grporder=c('appab', 'scr4x'),
            onlyDayorNight='night',
            skipNight0=TRUE,
            colours=c('#E94735', '#697a87'),
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
            exportPath=here('220524_APPAB', 'pubplots', 'actboutLength_night.pdf'))

ggParameter(pa=c(here('220524_APPAB', 'bhvparams', 'sleepHours_220524_14.csv'),
                 here('220524_APPAB', 'bhvparams', 'sleepHours_220524_15.csv')),
            grporder=c('appab', 'scr4x'),
            onlyDayorNight='night',
            skipNight0=TRUE,
            colours=c('#E94735', '#697a87'),
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
            exportPath=here('220524_APPAB', 'pubplots', 'sleepTot_night.pdf'))


### effect size
percEffectSize(pa=c(here('220524_APPAB', 'bhvparams', 'activityPercentageTimeActive_220524_14.csv'),
                    here('220524_APPAB', 'bhvparams', 'activityPercentageTimeActive_220524_15.csv')),
               lmePath=here('LMEreports/appab_f0_LME.csv'),
               grporder=c('appab', 'scr4x'),
               dayornight='day')

percEffectSize(pa=c(here('220524_APPAB', 'bhvparams', 'activeboutLength_220524_14.csv'),
                    here('220524_APPAB', 'bhvparams', 'activeboutLength_220524_15.csv')),
               lmePath=here('LMEreports/appab_f0_LME.csv'),
               grporder=c('appab', 'scr4x'),
               dayornight='day')

percEffectSize(pa=c(here('220524_APPAB', 'bhvparams', 'activeboutLength_220524_14.csv'),
                    here('220524_APPAB', 'bhvparams', 'activeboutLength_220524_15.csv')),
               lmePath=here('LMEreports/appab_f0_LME.csv'),
               grporder=c('appab', 'scr4x'),
               dayornight='night')

percEffectSize(pa=c(here('220524_APPAB', 'bhvparams', 'sleepHours_220524_14.csv'),
                    here('220524_APPAB', 'bhvparams', 'sleepHours_220524_15.csv')),
               lmePath=here('LMEreports/appab_f0_LME.csv'),
               grporder=c('appab', 'scr4x'),
               dayornight='night')
