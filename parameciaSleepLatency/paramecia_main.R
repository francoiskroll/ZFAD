#####################################################
# ~ ZFAD: parameciaSleepLatency ~
#
# trying to find a good illustration of sleep latency survival plots, maybe paramecia data from PhD thesis would work well
#
# Francois Kroll 2023
# francois@kroll.be
#####################################################



# packages ----------------------------------------------------------------

# library(devtools)
# install_github('francoiskroll/FramebyFrame')

library(FramebyFrame)
library(here)
library(ggplot2)



# vpSorter ----------------------------------------------------------------

# vpSorter(ffDir=here('parameciaSleepLatency', '210325_10_11_WhatmanDark_rawoutput/'),
#          zebpath=here('parameciaSleepLatency', '210325_10_11_WhatmanDark.xls'),
#          boxGen=2,
#          twoBoxMode=TRUE,
#          boxnum=1,
#          zt0='09:00:00',
#          dayduration=14)



# sleep latency -----------------------------------------------------------

# by day/night
# night0 is full
multiBehaviourParameter(parameters='sleepLatency',
                        ffpath=here('parameciaSleepLatency', '210325_10_RAWs.csv'),
                        genopath=here('parameciaSleepLatency', '210325_10genotype.txt'),
                        skipNight0=FALSE)


ggParameter(pa=here('parameciaSleepLatency', 'bhvparams', 'sleepLatency_210325_10.csv'),
            grporder=NA,
            skipNight0=FALSE,
            colours=NA,
            ymin=0,
            ymax=25,
            dotSize=1.2,
            legendOrNo=TRUE,
            xtextOrNo=TRUE,
            ynameOrNo=TRUE,
            yunitOrNo=FALSE,
            titleOrNo=FALSE,
            statsOrNo=FALSE,
            width=120,
            height=60,
            exportPath=here('parameciaSleepLatency', 'sleepLatency.pdf'))

ggSleepLatencySurvival(pa=here('parameciaSleepLatency', 'bhvparams', 'sleepLatency_210325_10.csv'),
                       grporder=NA,
                       skipNight0=FALSE,
                       colours=NA,
                       legendOrNo=FALSE,
                       xtextOrNo=TRUE,
                       xnameOrNo=TRUE,
                       ynameOrNo=TRUE,
                       nightBgOrNo=TRUE,
                       xmaxh=0.5,
                       detailsOrNo=TRUE, 
                       exportOrNo=TRUE,
                       exportDir=here('parameciaSleepLatency'),
                       width=55,
                       height=45,
                       dayduration=14)


ggSleepLatencySurvival(pa='~/Dropbox/ZFAD/parameciaSleepLatency/bhvparams/sleepLatency_210325_10.csv',
                       grporder=c('wt_Whatman_Dark', 'wt_Whatman_parameciaWell_Dark'),
                       skipNight0=FALSE,
                       colours=c('#697a87', '#EE7163'),
                       legendOrNo=FALSE,
                       xtextOrNo=TRUE,
                       xnameOrNo=FALSE,
                       ynameOrNo=TRUE,
                       nightBgOrNo=TRUE,
                       xmaxh=0.5,
                       detailsOrNo=TRUE, 
                       exportOrNo=TRUE,
                       exportDir='~/Dropbox/ZFAD/parameciaSleepLatency/',
                       width=80,
                       height=45,
                       dayduration=14)
