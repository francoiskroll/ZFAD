#####################################################
# ~ ZFAD: make simple parameter plots of exp 220725_CD2AP for main figure of ZFAD paper ~
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


### DAY ###

# common settings
width <- 34
height <- 50

dir.create(here('220725_CD2AP', 'pubplots'))


ggParameter(pa=c(here('220725_CD2AP', 'bhvparams', 'activityPercentageTimeActive_220725_16.csv'),
                 here('220725_CD2AP', 'bhvparams', 'activityPercentageTimeActive_220725_17.csv')),
            grporder=c('cd2ap', 'scr'),
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
            exportPath=here('220725_CD2AP', 'pubplots', 'actTimeActive_night.pdf'))
### effect size
percEffectSize(pa=c(here('220725_CD2AP', 'bhvparams', 'activityPercentageTimeActive_220725_16.csv'),
                    here('220725_CD2AP', 'bhvparams', 'activityPercentageTimeActive_220725_17.csv')),
               lmePath=here('LMEreports/cd2ap_f0_LME.csv'),
               grporder=c('cd2ap', 'scr'),
               dayornight='night')


ggParameter(pa=c(here('220725_CD2AP', 'bhvparams', 'activityTotalPx_220725_16.csv'),
                 here('220725_CD2AP', 'bhvparams', 'activityTotalPx_220725_17.csv')),
            grporder=c('cd2ap', 'scr'),
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
            exportPath=here('220725_CD2AP', 'pubplots', 'actTotPx_night.pdf'))
### effect size
percEffectSize(pa=c(here('220725_CD2AP', 'bhvparams', 'activityTotalPx_220725_16.csv'),
                    here('220725_CD2AP', 'bhvparams', 'activityTotalPx_220725_17.csv')),
               lmePath=here('LMEreports/cd2ap_f0_LME.csv'),
               grporder=c('cd2ap', 'scr'),
               dayornight='night')

ggParameter(pa=c(here('220725_CD2AP', 'bhvparams', 'activeboutNum_220725_16.csv'),
                 here('220725_CD2AP', 'bhvparams', 'activeboutNum_220725_17.csv')),
            grporder=c('cd2ap', 'scr'),
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
            exportPath=here('220725_CD2AP', 'pubplots', 'actboutNum_night.pdf'))
### effect size
percEffectSize(pa=c(here('220725_CD2AP', 'bhvparams', 'activeboutNum_220725_16.csv'),
                    here('220725_CD2AP', 'bhvparams', 'activeboutNum_220725_17.csv')),
               lmePath=here('LMEreports/cd2ap_f0_LME.csv'),
               grporder=c('cd2ap', 'scr'),
               dayornight='night')

ggParameter(pa=c(here('220725_CD2AP', 'bhvparams', 'sleepHours_220725_16.csv'),
                 here('220725_CD2AP', 'bhvparams', 'sleepHours_220725_17.csv')),
            grporder=c('cd2ap', 'scr'),
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
            exportPath=here('220725_CD2AP', 'pubplots', 'sleeptotal_night.pdf'))
### effect size
percEffectSize(pa=c(here('220725_CD2AP', 'bhvparams', 'sleepHours_220725_16.csv'),
                    here('220725_CD2AP', 'bhvparams', 'sleepHours_220725_17.csv')),
               lmePath=here('LMEreports/cd2ap_f0_LME.csv'),
               grporder=c('cd2ap', 'scr'),
               dayornight='night')
