#####################################################
# ~ ZFAD: traces 210907_PSEN2 ~
#
#
# Francois Kroll 2023
# francois@kroll.be
#####################################################


# scr = injected with 3x non-targeting gRNAs
# psen2 = psen2 F0 knockouts


# packages & functions ----------------------------------------------------

# load the Frame-by-Frame package

# library(devtools)
# install_github('francoiskroll/FramebyFrame')

library(FramebyFrame)
library(here)

# create a new folder to store the plots we make
dir.create(here('210907_PSEN2', 'pubplots'), showWarnings=FALSE)


# ACTIVITY TRACES ---------------------------------------------------------

# box12
ggActivityTraceByGroup(ffpath=here('data', '210907_PSEN2', '210907_12_RAWs.csv'),
                       genopath=here('210907_PSEN2', '210907_12genotype.txt'),
                       smoothOrNo=TRUE,
                       smooth_nsecs=30*60,
                       binOrNo=TRUE,
                       bin_nsecs=10*60,
                       grporder=c('psen2', 'scr'),
                       tracecols=c('#78ac63', '#697a87'),
                       ribboncols=c('#bbd4ae', '#b3bcc3'),
                       linethick=0.4,
                       xname='',
                       yname='activity (sum of Δ px/10 min)',
                       xtextOrNo=FALSE,
                       ytextOrNo=TRUE,
                       ymin=0,
                       ymax=40000,
                       xstart=24,
                       xstop=69.8,
                       trimstart=24,
                       trimstop=72,
                       nightBgOrNo=TRUE,
                       legendOrNo=FALSE,
                       exportOrNo=TRUE,
                       exportPath=here('210907_PSEN2', 'pubplots', '210913_12_tracebygroup_noadjust.pdf'),
                       width=68,
                       height=50)

# box13
ggActivityTraceByGroup(ffpath=here('data', '210907_PSEN2', '210907_13_RAWs.csv'),
                       genopath=here('210907_PSEN2', '210907_13genotype.txt'),
                       smoothOrNo=TRUE,
                       smooth_nsecs=30*60,
                       binOrNo=TRUE,
                       bin_nsecs=10*60,
                       grporder=c('psen2', 'scr'),
                       tracecols=c('#78ac63', '#697a87'),
                       ribboncols=c('#bbd4ae', '#b3bcc3'),
                       linethick=0.4,
                       xname='',
                       yname='activity (sum of Δ px/10 min)',
                       xtextOrNo=FALSE,
                       ytextOrNo=TRUE,
                       ymin=0,
                       ymax=40000,
                       xstart=24,
                       xstop=69.8,
                       trimstart=24,
                       trimstop=72,
                       nightBgOrNo=TRUE,
                       legendOrNo=FALSE,
                       exportOrNo=TRUE,
                       exportPath=here('210907_PSEN2', 'pubplots', '210913_13_tracebygroup_noadjust.pdf'),
                       width=68,
                       height=50)



# SLEEP TRACES ------------------------------------------------------------

# box12
ggSleepTraceByGroup(ffpath=here('data', '210907_PSEN2', '210907_12_RAWs.csv'),
                    genopath=here('210907_PSEN2', '210907_12genotype.txt'),
                    epo_min=10,
                    smoothOrNo=TRUE, 
                    smooth_npoints=5,
                    grporder=c('psen2', 'scr'),
                    tracecols=c('#78ac63', '#697a87'),
                    ribboncols=c('#bbd4ae', '#b3bcc3'),
                    linethick=0.4,
                    xname='',
                    yname='sleep (min/10 min)',
                    xtextOrNo=FALSE,
                    ytextOrNo=TRUE,
                    ymin=0,
                    ymax=10,
                    xstart=24,
                    xstop=69.8,
                    trimstart=24,
                    trimstop=72,
                    nightBgOrNo=TRUE,
                    legendOrNo=FALSE,
                    exportOrNo=TRUE,
                    exportPath=here('210907_PSEN2', 'pubplots', '210913_12_sleepbygroup_noadjust.pdf'),
                    width=68,
                    height=50)

# box13
ggSleepTraceByGroup(ffpath=here('data', '210907_PSEN2', '210907_13_RAWs.csv'),
                    genopath=here('210907_PSEN2', '210907_13genotype.txt'),
                    epo_min=10,
                    smoothOrNo=TRUE, 
                    smooth_npoints=5,
                    grporder=c('psen2', 'scr'),
                    tracecols=c('#78ac63', '#697a87'),
                    ribboncols=c('#bbd4ae', '#b3bcc3'),
                    linethick=0.4,
                    xname='',
                    yname='sleep (min/10 min)',
                    xtextOrNo=FALSE,
                    ytextOrNo=TRUE,
                    ymin=0,
                    ymax=10,
                    xstart=24,
                    xstop=69.8,
                    trimstart=24,
                    trimstop=72,
                    nightBgOrNo=TRUE,
                    legendOrNo=FALSE,
                    exportOrNo=TRUE,
                    exportPath=here('210907_PSEN2', 'pubplots', '210913_13_sleepbygroup_noadjust.pdf'),
                    width=68,
                    height=50)