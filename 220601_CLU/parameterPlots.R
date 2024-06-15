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


### DAY ###

# common settings
width <- 34
height <- 50

dir.create(here('220601_CLU', 'pubplots'))


ggParameter(pa=c(here('220601_CLU', 'bhvparams', 'activityPercentageTimeActive_220601_16.csv'),
                 here('220601_CLU', 'bhvparams', 'activityPercentageTimeActive_220601_17.csv')),
            grporder=c('clu', 'scr'),
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
            exportPath=here('220601_CLU', 'pubplots', 'actTimeActive_day.pdf'))

ggParameter(pa=c(here('220601_CLU', 'bhvparams', 'activeboutMean_220601_16.csv'),
                 here('220601_CLU', 'bhvparams', 'activeboutMean_220601_17.csv')),
            grporder=c('clu', 'scr'),
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
            exportPath=here('220601_CLU', 'pubplots', 'actboutMean_day.pdf'))


ggParameter(pa=c(here('220601_CLU', 'bhvparams', 'activityPercentageTimeActive_220601_16.csv'),
                 here('220601_CLU', 'bhvparams', 'activityPercentageTimeActive_220601_17.csv')),
            grporder=c('clu', 'scr'),
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
            exportPath=here('220601_CLU', 'pubplots', 'actTimeActive_night.pdf'))
### effect size
percEffectSize(pa=c(here('220601_CLU', 'bhvparams', 'activityPercentageTimeActive_220601_16.csv'),
                    here('220601_CLU', 'bhvparams', 'activityPercentageTimeActive_220601_17.csv')),
               lmePath=here('LMEreports/clu_f0_LME.csv'),
               grporder=c('clu', 'scr'),
               dayornight='night')


ggParameter(pa=c(here('220601_CLU', 'bhvparams', 'sleepHours_220601_16.csv'),
                 here('220601_CLU', 'bhvparams', 'sleepHours_220601_17.csv')),
            grporder=c('clu', 'scr'),
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
            exportPath=here('220601_CLU', 'pubplots', 'sleeptotal_night.pdf'))
### effect size
percEffectSize(pa=c(here('220601_CLU', 'bhvparams', 'sleepHours_220601_16.csv'),
                    here('220601_CLU', 'bhvparams', 'sleepHours_220601_17.csv')),
               lmePath=here('LMEreports/clu_f0_LME.csv'),
               grporder=c('clu', 'scr'),
               dayornight='night')


### effect size
percEffectSize(pa=c(here('220601_CLU', 'bhvparams', 'activeboutNum_220601_16.csv'),
                    here('220601_CLU', 'bhvparams', 'activeboutNum_220601_17.csv')),
               lmePath=here('LMEreports/clu_f0_LME.csv'),
               grporder=c('clu', 'scr'),
               dayornight='night')
