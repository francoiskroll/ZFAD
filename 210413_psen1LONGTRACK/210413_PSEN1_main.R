#####################################################
# ~ ZFAD: main script for analysis/plotting experiment 210413_PSEN1 ~
#
#
# Francois Kroll 2024
# francois@kroll.be
#####################################################


# packages ----------------------------------------------------------------

library(here)

# library(devtools)
# install_github('francoiskroll/FramebyFrame')
library(FramebyFrame)

library(data.table)


# add time columns --------------------------------------------------------

# 210413_13_RAWs.csv was created with code before FramebyFrame package
# add time columns to it
# (only need to do it once)

# addTimetoRAWs(ffpath=here('210413_psen1LONGTRACK', '210413_13_RAWs.csv'),
#               zebpath=here('210413_psen1LONGTRACK', '210413_12_13_PSEN1.xls'),
#               dayduration=14)

# writes 210413_13_RAWsv2.csv



# plot trace --------------------------------------------------------------

### v4 -- colours to go with sleep parameters for new version of FramebyFrame figure with zooms
ggActivityTraceByGroup(ffpath=here('data', '210413_psen1LONGTRACK', '210413_13_RAWsv2.csv'),
                       genopath=here('210413_psen1LONGTRACK', '210413_13genotype.txt'),
                       smoothOrNo=TRUE,
                       smooth_nsecs=30*60,
                       binOrNo=TRUE,
                       bin_nsecs=10*60,
                       grporder='scr',
                       tracecols='#417dcd',
                       ribboncols='#a3bbdb',
                       linethick=0.4,
                       xname='',
                       yname='activity (sum of Δ px/10 min)',
                       xtextOrNo=TRUE,
                       ytextOrNo=TRUE,
                       xmajorOrNo=FALSE,
                       ymin=0,
                       ymax=30000,
                       xstart=20,
                       xstop=210,
                       trimstart=5,
                       trimstop=217,
                       nightBgOrNo=FALSE,
                       sunlinesOrNo=FALSE,
                       legendOrNo=FALSE,
                       exportOrNo=TRUE,
                       exportPath=here('210413_psen1LONGTRACK', 'plots', 'tracebygroup_v3.pdf'),
                       width=112,
                       height=50)

ff <- fread(here('data', '210413_psen1LONGTRACK/210413_13_RAWsv2.csv'))
# tracking started 2021-04-13 18:18:09
# tacking finished 2021-04-22 10:48:24
# which is 208.5 hr
# first 9 AM is 14.7 hr
# second 9 AM is 38.7 hr
# etc.