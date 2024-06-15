#####################################################
# ~ ZFAD: make simple parameter plots of exp 220531_SORL1 for main figure of ZFAD paper ~
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

# activityPercentageTimeActive --------------------------------------------

ggParameter(pa=c(here('220531_SORL1', 'bhvparams', 'activityPercentageTimeActive_220531_15.csv'),
                 here('220531_SORL1', 'bhvparams', 'activityPercentageTimeActive_220531_14.csv')),
            grporder=c('sorl1', 'scr'),
            onlyDayorNight='day',
            colours=c('#417dcd', '#697a87'),
            fainterExp=TRUE,
            faintMax=0.35,
            ymin=0,
            ymax=25,
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
            exportPath=here('220531_SORL1', 'pubplots', 'actPerc_day.pdf'))
### effect size
percEffectSize(pa=c(here('220531_SORL1', 'bhvparams', 'activityPercentageTimeActive_220531_15.csv'),
                    here('220531_SORL1', 'bhvparams', 'activityPercentageTimeActive_220531_14.csv')),
               lmePath=here('LMEreports/sorl1_f0_LME.csv'),
               grporder=c('sorl1', 'scr'),
               dayornight='day')


# activeboutLength --------------------------------------------------------

ggParameter(pa=c(here('220531_SORL1', 'bhvparams', 'activeboutLength_220531_15.csv'),
                 here('220531_SORL1', 'bhvparams', 'activeboutLength_220531_14.csv')),
            grporder=c('sorl1', 'scr'),
            onlyDayorNight='day',
            colours=c('#417dcd', '#697a87'),
            fainterExp=TRUE,
            faintMax=0.35,
            ymin=0.10,
            ymax=0.25,
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
            exportPath=here('220531_SORL1', 'pubplots', 'actboutLength_day.pdf'))



# activeboutMean ----------------------------------------------------------

ggParameter(pa=c(here('220531_SORL1', 'bhvparams', 'activeboutMean_220531_15.csv'),
                 here('220531_SORL1', 'bhvparams', 'activeboutMean_220531_14.csv')),
            grporder=c('sorl1', 'scr'),
            onlyDayorNight='day',
            colours=c('#417dcd', '#697a87'),
            fainterExp=TRUE,
            faintMax=0.35,
            ymin=10,
            ymax=23,
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
            exportPath=here('220531_SORL1', 'pubplots', 'actboutMean_day.pdf'))



# activeboutNum -----------------------------------------------------------

ggParameter(pa=c(here('220531_SORL1', 'bhvparams', 'activeboutNum_220531_15.csv'),
                 here('220531_SORL1', 'bhvparams', 'activeboutNum_220531_14.csv')),
            grporder=c('sorl1', 'scr'),
            onlyDayorNight='day',
            colours=c('#417dcd', '#697a87'),
            fainterExp=TRUE,
            faintMax=0.35,
            ymin=5000,
            ymax=65000,
            dotSize=0.3,
            connectMean=TRUE,
            legendOrNo=FALSE,
            xtextOrNo=FALSE,
            ynameOrNo=FALSE,
            yunitOrNo=FALSE, # can guess the unit in this guess, saves some space
            splitBy='experiment',
            titleOrNo=FALSE,
            statsOrNo=FALSE,
            width=width+2, # need a bit larger because the Y axis labels take a lot of space
            height=height,
            exportPath=here('220531_SORL1', 'pubplots', 'actboutNum_day.pdf'))
### effect size
percEffectSize(pa=c(here('220531_SORL1', 'bhvparams', 'activeboutNum_220531_15.csv'),
                    here('220531_SORL1', 'bhvparams', 'activeboutNum_220531_14.csv')),
               lmePath=here('LMEreports/sorl1_f0_LME.csv'),
               grporder=c('sorl1', 'scr'),
               dayornight='day')


# sleepHours --------------------------------------------------------------

ggParameter(pa=c(here('220531_SORL1', 'bhvparams', 'sleepHours_220531_15.csv'),
                 here('220531_SORL1', 'bhvparams', 'sleepHours_220531_14.csv')),
            grporder=c('sorl1', 'scr'),
            onlyDayorNight='day',
            colours=c('#417dcd', '#697a87'),
            fainterExp=TRUE,
            faintMax=0.35,
            ymin=0,
            ymax=7,
            dotSize=0.3,
            connectMean=TRUE,
            legendOrNo=FALSE,
            xtextOrNo=FALSE,
            ynameOrNo=FALSE,
            yunitOrNo=TRUE,
            splitBy='experiment',
            titleOrNo=FALSE,
            statsOrNo=FALSE,
            width=width, # need a bit larger because the Y axis labels take a lot of space
            height=height,
            exportPath=here('220531_SORL1', 'pubplots', 'sleepHours_day.pdf'))



### NIGHT ###


# activityPercentageTimeActive --------------------------------------------

ggParameter(pa=c(here('220531_SORL1', 'bhvparams', 'activityPercentageTimeActive_220531_15.csv'),
                 here('220531_SORL1', 'bhvparams', 'activityPercentageTimeActive_220531_14.csv')),
            grporder=c('sorl1', 'scr'),
            onlyDayorNight='night',
            skipNight0=TRUE,
            colours=c('#417dcd', '#697a87'),
            fainterExp=TRUE,
            faintMax=0.35,
            ymin=0,
            ymax=2,
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
            exportPath=here('220531_SORL1', 'pubplots', 'actPerc_night.pdf'))
### effect size
percEffectSize(pa=c(here('220531_SORL1', 'bhvparams', 'activityPercentageTimeActive_220531_15.csv'),
                    here('220531_SORL1', 'bhvparams', 'activityPercentageTimeActive_220531_14.csv')),
               lmePath=here('LMEreports/sorl1_f0_LME.csv'),
               grporder=c('sorl1', 'scr'),
               dayornight='night')


# activeboutDuration ------------------------------------------------------

ggParameter(pa=c(here('220531_SORL1', 'bhvparams', 'activeboutLength_220531_15.csv'),
                 here('220531_SORL1', 'bhvparams', 'activeboutLength_220531_14.csv')),
            grporder=c('sorl1', 'scr'),
            onlyDayorNight='night',
            skipNight0=TRUE,
            colours=c('#417dcd', '#697a87'),
            fainterExp=TRUE,
            faintMax=0.35,
            ymin=0.12,
            ymax=0.25,
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
            exportPath=here('220531_SORL1', 'pubplots', 'actboutDur_night.pdf'))
### effect size
percEffectSize(pa=c(here('220531_SORL1', 'bhvparams', 'activeboutLength_220531_15.csv'),
                    here('220531_SORL1', 'bhvparams', 'activeboutLength_220531_14.csv')),
               lmePath=here('LMEreports/sorl1_f0_LME.csv'),
               grporder=c('sorl1', 'scr'),
               dayornight='night')

# sleepHours --------------------------------------------------------------

ggParameter(pa=c(here('220531_SORL1', 'bhvparams', 'sleepHours_220531_15.csv'),
                 here('220531_SORL1', 'bhvparams', 'sleepHours_220531_14.csv')),
            grporder=c('sorl1', 'scr'),
            onlyDayorNight='night',
            skipNight0=TRUE,
            colours=c('#417dcd', '#697a87'),
            fainterExp=TRUE,
            faintMax=0.35,
            ymin=4,
            ymax=9.9,
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
            exportPath=here('220531_SORL1', 'pubplots', 'sleepHours_night.pdf'))
### effect size
percEffectSize(pa=c(here('220531_SORL1', 'bhvparams', 'sleepHours_220531_15.csv'),
                    here('220531_SORL1', 'bhvparams', 'sleepHours_220531_14.csv')),
               lmePath=here('LMEreports/sorl1_f0_LME.csv'),
               grporder=c('sorl1', 'scr'),
               dayornight='night')



# sleepNumNaps ------------------------------------------------------------

ggParameter(pa=c(here('220531_SORL1', 'bhvparams', 'sleepNumNaps_220531_15.csv'),
                 here('220531_SORL1', 'bhvparams', 'sleepNumNaps_220531_14.csv')),
            grporder=c('sorl1', 'scr'),
            onlyDayorNight='night',
            skipNight0=TRUE,
            colours=c('#417dcd', '#697a87'),
            fainterExp=TRUE,
            faintMax=0.35,
            ymin=40,
            ymax=230,
            dotSize=0.3,
            connectMean=TRUE,
            legendOrNo=FALSE,
            xtextOrNo=FALSE,
            ynameOrNo=FALSE,
            yunitOrNo=FALSE, # can guess units here
            splitBy='experiment',
            titleOrNo=FALSE,
            statsOrNo=FALSE,
            width=width-1, # minus because we removed units from Y axis
            height=height,
            exportPath=here('220531_SORL1', 'pubplots', 'sleepNumNaps_night.pdf'))
### effect size
percEffectSize(pa=c(here('220531_SORL1', 'bhvparams', 'sleepNumNaps_220531_15.csv'),
                    here('220531_SORL1', 'bhvparams', 'sleepNumNaps_220531_14.csv')),
               lmePath=here('LMEreports/sorl1_f0_LME.csv'),
               grporder=c('sorl1', 'scr'),
               dayornight='night')


# sleepNapDuration --------------------------------------------------------

ggParameter(pa=c(here('220531_SORL1', 'bhvparams', 'sleepNapDuration_220531_15.csv'),
                 here('220531_SORL1', 'bhvparams', 'sleepNapDuration_220531_14.csv')),
            grporder=c('sorl1', 'scr'),
            onlyDayorNight='night',
            skipNight0=TRUE,
            colours=c('#417dcd', '#697a87'),
            fainterExp=TRUE,
            faintMax=0.35,
            ymin=1.95,
            ymax=6.7,
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
            exportPath=here('220531_SORL1', 'pubplots', 'sleepNapDur_night.pdf'))


### effect size
percEffectSize(pa=c(here('220531_SORL1', 'bhvparams', 'activeboutNum_220531_15.csv'),
                    here('220531_SORL1', 'bhvparams', 'activeboutNum_220531_14.csv')),
               lmePath=here('LMEreports/sorl1_f0_LME.csv'),
               grporder=c('sorl1', 'scr'),
               dayornight='night')
