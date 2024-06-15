#####################################################
# ~ ZFAD: make simple parameter plots of exp 210907_PSEN2 for main figure of ZFAD paper ~
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

# ! plotting below in order BOX13, BOX12
# clutch1 is BOX13
# clutch2 is BOX12



# activityTotalPx ---------------------------------------------------------

ggParameter(pa=c(here('210907_PSEN2', 'bhvparams_beforeAdjust', 'activityTotalPx_210907_13.csv'),
                 here('210907_PSEN2', 'bhvparams_beforeAdjust', 'activityTotalPx_210907_12.csv')),
            grporder=c('psen2', 'scr'),
            onlyDayorNight='day',
            colours=c('#78ac63', '#697a87'),
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
            exportPath=here('210907_PSEN2', 'pubplots', 'actPx_day_noadjust.pdf'))

### effect size
percEffectSize(pa=c(here('210907_PSEN2', 'bhvparams_beforeAdjust', 'activityTotalPx_210907_13.csv'),
                    here('210907_PSEN2', 'bhvparams_beforeAdjust', 'activityTotalPx_210907_12.csv')),
               lmePath=here('LMEreports/psen2_f0_LME_beforeAdjust.csv'),
               grporder=c('psen2', 'scr'),
               dayornight='day')


# activeboutMean ----------------------------------------------------------

ggParameter(pa=c(here('210907_PSEN2', 'bhvparams_beforeAdjust', 'activeboutMean_210907_13.csv'),
                 here('210907_PSEN2', 'bhvparams_beforeAdjust', 'activeboutMean_210907_12.csv')),
            grporder=c('psen2', 'scr'),
            onlyDayorNight='day',
            colours=c('#78ac63', '#697a87'),
            fainterExp=TRUE,
            faintMax=0.35,
            ymin=5,
            ymax=20,
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
            exportPath=here('210907_PSEN2', 'pubplots', 'actboutMean_day_noadjust.pdf'))
### effect size
percEffectSize(pa=c(here('210907_PSEN2', 'bhvparams_beforeAdjust', 'activeboutMean_210907_13.csv'),
                    here('210907_PSEN2', 'bhvparams_beforeAdjust', 'activeboutMean_210907_12.csv')),
               lmePath=here('LMEreports/psen2_f0_LME_beforeAdjust.csv'),
               grporder=c('psen2', 'scr'),
               dayornight='day')



# activeboutNum -----------------------------------------------------------

ggParameter(pa=c(here('210907_PSEN2', 'bhvparams_beforeAdjust', 'activeboutNum_210907_13.csv'),
                 here('210907_PSEN2', 'bhvparams_beforeAdjust', 'activeboutNum_210907_12.csv')),
            grporder=c('psen2', 'scr'),
            onlyDayorNight='day',
            colours=c('#78ac63', '#697a87'),
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
            width=width+2, # need a bit larger because the Y axis labels take a lot of space
            height=height,
            exportPath=here('210907_PSEN2', 'pubplots', 'actboutNum_day_noadjust.pdf'))
### effect size
percEffectSize(pa=c(here('210907_PSEN2', 'bhvparams_beforeAdjust', 'activeboutNum_210907_13.csv'),
                    here('210907_PSEN2', 'bhvparams_beforeAdjust', 'activeboutNum_210907_12.csv')),
               lmePath=here('LMEreports/psen2_f0_LME_beforeAdjust.csv'),
               grporder=c('psen2', 'scr'),
               dayornight='day')



# sleepHours --------------------------------------------------------------

ggParameter(pa=c(here('210907_PSEN2', 'bhvparams_beforeAdjust', 'sleepHours_210907_13.csv'),
                 here('210907_PSEN2', 'bhvparams_beforeAdjust', 'sleepHours_210907_12.csv')),
            grporder=c('psen2', 'scr'),
            onlyDayorNight='day',
            colours=c('#78ac63', '#697a87'),
            fainterExp=TRUE,
            faintMax=0.35,
            ymin=0,
            ymax=8.2,
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
            exportPath=here('210907_PSEN2', 'pubplots', 'sleepHours_day_noadjust.pdf'))
### effect size
percEffectSize(pa=c(here('210907_PSEN2', 'bhvparams_beforeAdjust', 'sleepHours_210907_13.csv'),
                    here('210907_PSEN2', 'bhvparams_beforeAdjust', 'sleepHours_210907_12.csv')),
               lmePath=here('LMEreports/psen2_f0_LME_beforeAdjust.csv'),
               grporder=c('psen2', 'scr'),
               dayornight='day')



# sleepNumNaps ------------------------------------------------------------

ggParameter(pa=c(here('210907_PSEN2', 'bhvparams_beforeAdjust', 'sleepNumNaps_210907_13.csv'),
                 here('210907_PSEN2', 'bhvparams_beforeAdjust', 'sleepNumNaps_210907_12.csv')),
            grporder=c('psen2', 'scr'),
            onlyDayorNight='day',
            colours=c('#78ac63', '#697a87'),
            fainterExp=TRUE,
            faintMax=0.35,
            ymin=0,
            ymax=220,
            dotSize=0.3,
            connectMean=TRUE,
            legendOrNo=FALSE,
            xtextOrNo=FALSE,
            ynameOrNo=FALSE,
            yunitOrNo=FALSE, # can guess unit when counts
            splitBy='experiment',
            titleOrNo=FALSE,
            statsOrNo=FALSE,
            width=width,
            height=height,
            exportPath=here('210907_PSEN2', 'pubplots', 'sleepNum_day_noadjust.pdf'))
### effect size
percEffectSize(pa=c(here('210907_PSEN2', 'bhvparams_beforeAdjust', 'sleepNumNaps_210907_13.csv'),
                    here('210907_PSEN2', 'bhvparams_beforeAdjust', 'sleepNumNaps_210907_12.csv')),
               lmePath=here('LMEreports/psen2_f0_LME_beforeAdjust.csv'),
               grporder=c('psen2', 'scr'),
               dayornight='day')
