#####################################################
# ~ ZFAD: make simple parameter plots of exp APOEAB for main figure of ZFAD paper ~
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
width <- 38 # widen a bit compare to other scripts because 3 clutches (also, 4 plots and not 5 so have a little bit more space)
height <- 50

dir.create(here('220516_APOEAB_3', 'pubplots'))

# pasting here notes from apoeab_f0Stable.R

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

### effect size
percEffectSize(pa=c(here('220313_APOEAB_2', 'bhvparams', 'activityPercentageTimeActive_220313_15.csv'),
                    here('220516_APOEAB_3', 'bhvparams', 'activityPercentageTimeActive_220516_14.csv'),
                    here('220516_APOEAB_3', 'bhvparams', 'activityPercentageTimeActive_220516_15.csv')),
               lmePath=here('LMEreports/apoeab_f0_LME.csv'),
               grporder=c('apoeab', 'scr4x'),
               dayornight='night')


ggParameter(pa=c(here('220313_APOEAB_2', 'bhvparams', 'activeboutMean_220313_15.csv'),
                 here('220516_APOEAB_3', 'bhvparams', 'activeboutMean_220516_14.csv'),
                 here('220516_APOEAB_3', 'bhvparams', 'activeboutMean_220516_15.csv')),
            grporder=c('apoeab', 'scr4x'),
            onlyDayorNight='day',
            skipNight0=TRUE,
            colours=c('#e94735', '#697a87'),
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
            exportPath=here('220516_APOEAB_3', 'pubplots', 'actboutMean_day.pdf'))
### effect size
percEffectSize(pa=c(here('220313_APOEAB_2', 'bhvparams', 'activeboutMean_220313_15.csv'),
                    here('220516_APOEAB_3', 'bhvparams', 'activeboutMean_220516_14.csv'),
                    here('220516_APOEAB_3', 'bhvparams', 'activeboutMean_220516_15.csv')),
               lmePath=here('LMEreports/apoeab_f0_LME.csv'),
               grporder=c('apoeab', 'scr4x'),
               dayornight='day')


ggParameter(pa=c(here('220313_APOEAB_2', 'bhvparams', 'activeboutMax_220313_15.csv'),
                 here('220516_APOEAB_3', 'bhvparams', 'activeboutMax_220516_14.csv'),
                 here('220516_APOEAB_3', 'bhvparams', 'activeboutMax_220516_15.csv')),
            grporder=c('apoeab', 'scr4x'),
            onlyDayorNight='day',
            skipNight0=TRUE,
            colours=c('#e94735', '#697a87'),
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
            exportPath=here('220516_APOEAB_3', 'pubplots', 'actboutMax_day.pdf'))
### effect size
percEffectSize(pa=c(here('220313_APOEAB_2', 'bhvparams', 'activeboutMax_220313_15.csv'),
                    here('220516_APOEAB_3', 'bhvparams', 'activeboutMax_220516_14.csv'),
                    here('220516_APOEAB_3', 'bhvparams', 'activeboutMax_220516_15.csv')),
               lmePath=here('LMEreports/apoeab_f0_LME.csv'),
               grporder=c('apoeab', 'scr4x'),
               dayornight='day')

ggParameter(pa=c(here('220313_APOEAB_2', 'bhvparams', 'activeboutNum_220313_15.csv'),
                 here('220516_APOEAB_3', 'bhvparams', 'activeboutNum_220516_14.csv'),
                 here('220516_APOEAB_3', 'bhvparams', 'activeboutNum_220516_15.csv')),
            grporder=c('apoeab', 'scr4x'),
            onlyDayorNight='night',
            skipNight0=TRUE,
            colours=c('#e94735', '#697a87'),
            fainterExp=TRUE,
            faintMax=0.35,
            ymin=0,
            ymax=20000,
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
            exportPath=here('220516_APOEAB_3', 'pubplots', 'actboutNum_night.pdf'))
# ! ymax excludes some dots !
### effect size
percEffectSize(pa=c(here('220313_APOEAB_2', 'bhvparams', 'activeboutNum_220313_15.csv'),
                    here('220516_APOEAB_3', 'bhvparams', 'activeboutNum_220516_14.csv'),
                    here('220516_APOEAB_3', 'bhvparams', 'activeboutNum_220516_15.csv')),
               lmePath=here('LMEreports/apoeab_f0_LME.csv'),
               grporder=c('apoeab', 'scr4x'),
               dayornight='night')

ggParameter(pa=c(here('220313_APOEAB_2', 'bhvparams', 'sleepHours_220313_15.csv'),
                 here('220516_APOEAB_3', 'bhvparams', 'sleepHours_220516_14.csv'),
                 here('220516_APOEAB_3', 'bhvparams', 'sleepHours_220516_15.csv')),
            grporder=c('apoeab', 'scr4x'),
            onlyDayorNight='night',
            skipNight0=TRUE,
            colours=c('#e94735', '#697a87'),
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
            exportPath=here('220516_APOEAB_3', 'pubplots', 'sleephrs_night.pdf'))
### effect size
percEffectSize(pa=c(here('220313_APOEAB_2', 'bhvparams', 'sleepHours_220313_15.csv'),
                    here('220516_APOEAB_3', 'bhvparams', 'sleepHours_220516_14.csv'),
                    here('220516_APOEAB_3', 'bhvparams', 'sleepHours_220516_15.csv')),
               lmePath=here('LMEreports/apoeab_f0_LME.csv'),
               grporder=c('apoeab', 'scr4x'),
               dayornight='night')
