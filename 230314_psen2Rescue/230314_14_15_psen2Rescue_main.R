#####################################################
# ~ ZFAD: 230314_14_15_psen2Rescue ~
#
# main script
#
# Francois Kroll 2023
# francois@kroll.be
#####################################################


# packages & functions ----------------------------------------------------

library(here)

# install/load FramebyFrame package

# library(devtools)
# install_github('francoiskroll/FramebyFrame')
library(FramebyFrame)

dir.create(here('230314_psen2Rescue', 'plots'), showWarnings=FALSE)
dir.create(here('230314_psen2Rescue', 'pubplots'), showWarnings=FALSE)


# vpSorter ----------------------------------------------------------------

# below ran on Rihel lab common computer
vpSorter(ffDir=here('230314_psen2Rescue', '230314_14_15_psen2Rescue_rawoutput'),
         zebpath=here('230314_psen2Rescue', '230314_14_15_psen2Rescue.xls'),
         boxGen=2,
         twoBoxMode=TRUE,
         boxnum=1,
         zt0='09:00:00',
         dayduration=14,
         exportXlsOrNo=FALSE)

vpSorter(ffDir=here('230314_psen2Rescue', '230314_14_15_psen2Rescue_rawoutput'),
         zebpath=here('230314_psen2Rescue', '230314_14_15_psen2Rescue.xls'),
         boxGen=2,
         twoBoxMode=TRUE,
         boxnum=2,
         zt0='09:00:00',
         dayduration=14,
         exportXlsOrNo=FALSE)


# genotype file -----------------------------------------------------------

genotypeGenerator(plateMap=here('230314_psen2Rescue', '230314_14_genotypeMap.xlsx'))
genotypeGenerator(plateMap=here('230314_psen2Rescue', '230314_15_genotypeMap.xlsx'))




# adjust pixel to counteract fainter pigmentation of psen2 F0 -------------
# approach first developed for experiment 210907_PSEN2

# we decrease the delta pixels of SCR larvae by a bit

adjustPixel(ffpath='~/Dropbox/ZFAD/230314_psen2Rescue/230314_14_RAWs.csv',
            genopath='~/Dropbox/ZFAD/230314_psen2Rescue/230314_14genotype.txt',
            grpL='SCR.DMSO',
            grpS='PSEN2.DMSO',
            scale=NA,
            round='down',
            dayduration=14)
# scaling ratio = 0.97; p=0.4153

adjustPixel(ffpath='~/Dropbox/ZFAD/230314_psen2Rescue/230314_15_RAWs.csv',
            genopath='~/Dropbox/ZFAD/230314_psen2Rescue/230314_15genotype.txt',
            grpL='SCR.DMSO',
            grpS='PSEN2.DMSO',
            scale=NA,
            round='down',
            dayduration=14)
# scaling ratio = 0.95; p=0.1157

# checking with ggActivityTraceByGroup below; looks good


# frame rate --------------------------------------------------------------

ggFramerate(ffpath=here('230314_psen2Rescue', '230314_14_RAWs.csv'),
            subsample=TRUE,
            subsample_by=1000,
            xstart=0,
            xstop=0,
            ymin=0,
            ymax=50,
            sunlines=FALSE,
            dayduration=14,
            xname='hours since first 9 AM',
            yname='frames-per-second',
            exportOrNo=TRUE,
            width=75,
            height=55,
            exportPath=here('230314_psen2Rescue', 'plots', '230314_14_framerate.pdf'))

ggFramerate(ffpath=here('230314_psen2Rescue', '230314_15_RAWs.csv'),
            subsample=TRUE,
            subsample_by=1000,
            xstart=0,
            xstop=0,
            ymin=0,
            ymax=50,
            sunlines=FALSE,
            dayduration=14,
            xname='hours since first 9 AM',
            yname='frames-per-second',
            exportOrNo=TRUE,
            width=75,
            height=55,
            exportPath=here('230314_psen2Rescue', 'plots', '230314_15_framerate.pdf'))



# activity grid -----------------------------------------------------------

### box14 ###

# after adjust
ggActivityTraceGrid(ffpath=here('230314_psen2Rescue', '230314_14_RAWs.csv'),
                    genopath=here('230314_psen2Rescue', '230314_14genotype.txt'),
                    smoothOrNo=TRUE,
                    smooth_nsecs=30*60,
                    binOrNo=TRUE,
                    bin_nsecs=10*60,
                    tracecols=NA,
                    linethick=0.4,
                    ymin=0,
                    ymax=60000,
                    xstart=0,
                    xstop=0,
                    trimstart=0,
                    trimstop=0,
                    nightBgOrNo=TRUE,
                    ncol=12,
                    nrow=8,
                    exportOrNo=TRUE,
                    exportPath=here('230314_psen2Rescue', 'plots', '230314_14_activitygrid.pdf'),
                    width=255,
                    height=171)

# before adjust
ggActivityTraceGrid(ffpath=here('230314_psen2Rescue', '230314_14_RAWsadjusted.csv'),
                    genopath=here('230314_psen2Rescue', '230314_14genotype.txt'),
                    smoothOrNo=TRUE,
                    smooth_nsecs=30*60,
                    binOrNo=TRUE,
                    bin_nsecs=10*60,
                    tracecols=NA,
                    linethick=0.4,
                    ymin=0,
                    ymax=60000,
                    xstart=0,
                    xstop=0,
                    trimstart=0,
                    trimstop=0,
                    nightBgOrNo=TRUE,
                    ncol=12,
                    nrow=8,
                    exportOrNo=TRUE,
                    exportPath=here('230314_psen2Rescue', 'plots', '230314_14_activitygrid_adjusted.pdf'),
                    width=255,
                    height=171)

### box15 ###

# after adjust
ggActivityTraceGrid(ffpath=here('230314_psen2Rescue', '230314_15_RAWs.csv'),
                    genopath=here('230314_psen2Rescue', '230314_15genotype.txt'),
                    smoothOrNo=TRUE,
                    smooth_nsecs=30*60,
                    binOrNo=TRUE,
                    bin_nsecs=10*60,
                    tracecols=NA,
                    linethick=0.4,
                    ymin=0,
                    ymax=60000,
                    xstart=0,
                    xstop=0,
                    trimstart=0,
                    trimstop=0,
                    nightBgOrNo=TRUE,
                    ncol=12,
                    nrow=8,
                    exportOrNo=TRUE,
                    exportPath=here('230314_psen2Rescue', 'plots', '230314_15_activitygrid.pdf'),
                    width=255,
                    height=171)

# before adjust
ggActivityTraceGrid(ffpath=here('230314_psen2Rescue', '230314_15_RAWsadjusted.csv'),
                    genopath=here('230314_psen2Rescue', '230314_15genotype.txt'),
                    smoothOrNo=TRUE,
                    smooth_nsecs=30*60,
                    binOrNo=TRUE,
                    bin_nsecs=10*60,
                    tracecols=NA,
                    linethick=0.4,
                    ymin=0,
                    ymax=60000,
                    xstart=0,
                    xstop=0,
                    trimstart=0,
                    trimstop=0,
                    nightBgOrNo=TRUE,
                    ncol=12,
                    nrow=8,
                    exportOrNo=TRUE,
                    exportPath=here('230314_psen2Rescue', 'plots', '230314_15_activitygrid_adjusted.pdf'),
                    width=255,
                    height=171)


# BOX14: activity trace by group ------------------------------------------

### DMSO only, i.e. original phenotype
ggActivityTraceByGroup(ffpath=here('230314_psen2Rescue', '230314_14_RAWs.csv'),
                       genopath=here('230314_psen2Rescue', '230314_14genotype.txt'),
                       smoothOrNo=TRUE,
                       smooth_nsecs=30*60,
                       binOrNo=TRUE,
                       bin_nsecs=10*60,
                       grporder=c('SCR.DMSO', 'PSEN2.DMSO'),
                       tracecols=c('#697a87', '#cb2a20'),
                       ribboncols=c('#b3bcc3', '#e4928f'),
                       linethick=0.4,
                       xname='hours since ZT0',
                       yname='activity (sum of Δ px/10 min)',
                       xtextOrNo=FALSE,
                       ytextOrNo=TRUE,
                       ymin=0,
                       ymax=50000,
                       xstart=0,
                       xstop=0,
                       trimstart=0,
                       trimstop=0,
                       nightBgOrNo=TRUE,
                       legendOrNo=TRUE,
                       exportOrNo=TRUE,
                       exportPath=here('230314_psen2Rescue', 'plots', '230314_14_tracebygroup_DMSO.pdf'),
                       width=75,
                       height=55)

ggActivityTraceByGroup(ffpath='~/Dropbox/ZFAD/230314_psen2Rescue/230314_14_RAWsadjusted.csv',
                       genopath='~/Dropbox/ZFAD/230314_psen2Rescue/230314_14genotype.txt',
                       smoothOrNo=TRUE,
                       smooth_nsecs=30*60,
                       binOrNo=TRUE,
                       bin_nsecs=10*60,
                       grporder=c('SCR.DMSO', 'PSEN2.DMSO'),
                       tracecols=c('#697a87', '#cb2a20'),
                       ribboncols=c('#b3bcc3', '#e4928f'),
                       linethick=0.4,
                       xname='hours since ZT0',
                       yname='activity (sum of Δ px/10 min)',
                       xtextOrNo=FALSE,
                       ytextOrNo=TRUE,
                       ymin=0,
                       ymax=50000,
                       xstart=0,
                       xstop=0,
                       trimstart=0,
                       trimstop=0,
                       nightBgOrNo=TRUE,
                       legendOrNo=TRUE,
                       exportOrNo=TRUE,
                       exportPath='~/Dropbox/ZFAD/230314_psen2Rescue/plots/230314_14_DMSO_traceAfterAdjust.pdf',
                       width=75,
                       height=55)

### Tinidazole treatment
ggActivityTraceByGroup(ffpath=here('230314_psen2Rescue', '230314_14_RAWs.csv'),
                       genopath=here('230314_psen2Rescue', '230314_14genotype.txt'),
                       smoothOrNo=TRUE,
                       smooth_nsecs=30*60,
                       binOrNo=TRUE,
                       bin_nsecs=10*60,
                       grporder=c('SCR.DMSO', 'PSEN2.DMSO', 'PSEN2.Tinidazole'),
                       tracecols=c('#697a87', '#cb2a20', '#fcb505'),
                       ribboncols=c('#b3bcc3', '#e4928f', '#fcd98a'),
                       linethick=0.4,
                       xname='hours since ZT0',
                       yname='activity (sum of Δ px/10 min)',
                       xtextOrNo=FALSE,
                       ytextOrNo=TRUE,
                       ymin=0,
                       ymax=50000,
                       xstart=0,
                       xstop=0,
                       trimstart=0,
                       trimstop=0,
                       nightBgOrNo=TRUE,
                       legendOrNo=TRUE,
                       exportOrNo=TRUE,
                       exportPath=here('230314_psen2Rescue', 'plots', '230314_14_tracebygroup_Tinidazole.pdf'),
                       width=75,
                       height=55)
# works during the day?

### Betamethasone treatment

# POLISHED PLOT
# promising

ggActivityTraceByGroup(ffpath=here('230314_psen2Rescue', '230314_14_RAWs.csv'),
                       genopath=here('230314_psen2Rescue', '230314_14genotype.txt'),
                       smoothOrNo=TRUE,
                       smooth_nsecs=30*60,
                       binOrNo=TRUE,
                       bin_nsecs=10*60,
                       grporder=c('SCR.DMSO', 'PSEN2.DMSO', 'PSEN2.Betamethasone'),
                       tracecols=c('#78ac63', '#697a87', '#fcb505'),
                       ribboncols=c('#bbd4ae', '#b3bcc3', '#fcd98a'),
                       # tracecols=c('#697a87', '#78ac63', '#fcb505'),
                       # ribboncols=c('#b3bcc3', '#bbd4ae', '#fcd98a'),
                       linethick=0.4,
                       xname='',
                       yname='',
                       xtextOrNo=FALSE,
                       ytextOrNo=TRUE,
                       ymin=0,
                       ymax=60000,
                       xstart=24,
                       xstop=69.8,
                       trimstart=24,
                       trimstop=72,
                       nightBgOrNo=TRUE,
                       legendOrNo=FALSE,
                       exportOrNo=TRUE,
                       exportPath=here('230314_psen2Rescue', 'pubplots', 'box14_acttrace_betamethasone_beforeAdj.pdf'),
                       width=68,
                       height=50)

### Fenoprofen treatment
ggActivityTraceByGroup(ffpath=here('230314_psen2Rescue', '230314_14_RAWs.csv'),
                       genopath=here('230314_psen2Rescue', '230314_14genotype.txt'),
                       smoothOrNo=TRUE,
                       smooth_nsecs=30*60,
                       binOrNo=TRUE,
                       bin_nsecs=10*60,
                       grporder=c('SCR.DMSO', 'PSEN2.DMSO', 'PSEN2.Fenoprofen'),
                       tracecols=c('#697a87', '#cb2a20', '#fcb505'),
                       ribboncols=c('#b3bcc3', '#e4928f', '#fcd98a'),
                       linethick=0.4,
                       xname='hours since ZT0',
                       yname='activity (sum of Δ px/10 min)',
                       xtextOrNo=FALSE,
                       ytextOrNo=TRUE,
                       ymin=0,
                       ymax=50000,
                       xstart=0,
                       xstop=0,
                       trimstart=0,
                       trimstop=0,
                       nightBgOrNo=TRUE,
                       legendOrNo=TRUE,
                       exportOrNo=TRUE,
                       exportPath=here('230314_psen2Rescue', 'plots', '230314_14_tracebygroup_Fenoprofen.pdf'),
                       width=75,
                       height=55)
# wrong direction!



# BOX15: activity trace by group ------------------------------------------

### DMSO only, i.e. original phenotype
ggActivityTraceByGroup(ffpath=here('230314_psen2Rescue', '230314_15_RAWs.csv'),
                       genopath=here('230314_psen2Rescue', '230314_15genotype.txt'),
                       smoothOrNo=TRUE,
                       smooth_nsecs=30*60,
                       binOrNo=TRUE,
                       bin_nsecs=10*60,
                       grporder=c('SCR.DMSO', 'PSEN2.DMSO'),
                       tracecols=c('#697a87', '#cb2a20'),
                       ribboncols=c('#b3bcc3', '#e4928f'),
                       linethick=0.4,
                       xname='hours since ZT0',
                       yname='activity (sum of Δ px/10 min)',
                       xtextOrNo=FALSE,
                       ytextOrNo=TRUE,
                       ymin=0,
                       ymax=50000,
                       xstart=0,
                       xstop=0,
                       trimstart=0,
                       trimstop=0,
                       nightBgOrNo=TRUE,
                       legendOrNo=TRUE,
                       exportOrNo=TRUE,
                       exportPath=here('230314_psen2Rescue', 'plots', '230314_15_tracebygroup_DMSO.pdf'),
                       width=90,
                       height=55)

ggActivityTraceByGroup(ffpath='~/Dropbox/ZFAD/230314_psen2Rescue/230314_15_RAWsadjusted.csv',
                       genopath='~/Dropbox/ZFAD/230314_psen2Rescue/230314_15genotype.txt',
                       smoothOrNo=TRUE,
                       smooth_nsecs=30*60,
                       binOrNo=TRUE,
                       bin_nsecs=10*60,
                       grporder=c('SCR.DMSO', 'PSEN2.DMSO'),
                       tracecols=c('#697a87', '#cb2a20'),
                       ribboncols=c('#b3bcc3', '#e4928f'),
                       linethick=0.4,
                       xname='hours since ZT0',
                       yname='activity (sum of Δ px/10 min)',
                       xtextOrNo=FALSE,
                       ytextOrNo=TRUE,
                       ymin=0,
                       ymax=50000,
                       xstart=0,
                       xstop=0,
                       trimstart=0,
                       trimstop=0,
                       nightBgOrNo=TRUE,
                       legendOrNo=TRUE,
                       exportOrNo=TRUE,
                       exportPath='~/Dropbox/ZFAD/230314_psen2Rescue/plots/230314_15_DMSO_traceAfterAdjust.pdf',
                       width=90,
                       height=55)

### Tinidazole treatment
ggActivityTraceByGroup(ffpath=here('230314_psen2Rescue', '230314_15_RAWs.csv'),
                       genopath=here('230314_psen2Rescue', '230314_15genotype.txt'),
                       smoothOrNo=TRUE,
                       smooth_nsecs=30*60,
                       binOrNo=TRUE,
                       bin_nsecs=10*60,
                       grporder=c('SCR.DMSO', 'PSEN2.DMSO', 'PSEN2.Tinidazole'),
                       tracecols=c('#697a87', '#cb2a20', '#fcb505'),
                       ribboncols=c('#b3bcc3', '#e4928f', '#fcd98a'),
                       linethick=0.4,
                       xname='hours since ZT0',
                       yname='activity (sum of Δ px/10 min)',
                       xtextOrNo=FALSE,
                       ytextOrNo=TRUE,
                       ymin=0,
                       ymax=50000,
                       xstart=0,
                       xstop=0,
                       trimstart=0,
                       trimstop=0,
                       nightBgOrNo=TRUE,
                       legendOrNo=TRUE,
                       exportOrNo=TRUE,
                       exportPath=here('230314_psen2Rescue', 'plots', '230314_15_tracebygroup_Tinidazole.pdf'),
                       width=75,
                       height=55)
# works during the night? (looked more like day in BOX14)

### Betamethasone treatment
ggActivityTraceByGroup(ffpath=here('230314_psen2Rescue', '230314_15_RAWs.csv'),
                       genopath=here('230314_psen2Rescue', '230314_15genotype.txt'),
                       smoothOrNo=TRUE,
                       smooth_nsecs=30*60,
                       binOrNo=TRUE,
                       bin_nsecs=10*60,
                       grporder=c('SCR.DMSO', 'PSEN2.DMSO', 'PSEN2.Betamethasone'),
                       tracecols=c('#697a87', '#cb2a20', '#fcb505'),
                       ribboncols=c('#b3bcc3', '#e4928f', '#fcd98a'),
                       linethick=0.4,
                       xname='hours since ZT0',
                       yname='activity (sum of Δ px/10 min)',
                       xtextOrNo=FALSE,
                       ytextOrNo=TRUE,
                       ymin=0,
                       ymax=60000,
                       xstart=0,
                       xstop=0,
                       trimstart=0,
                       trimstop=0,
                       nightBgOrNo=TRUE,
                       legendOrNo=FALSE,
                       exportOrNo=TRUE,
                       exportPath=here('230314_psen2Rescue', 'plots', '230314_15_tracebygroup_Betamethasone.pdf'),
                       width=105,
                       height=60)
# really good!

# POLISHED VERSION

### before adjust
# ! delete from folder acttc... and clear Env
ggActivityTraceByGroup(ffpath=here('230314_psen2Rescue', '230314_15_RAWs.csv'),
                       genopath=here('230314_psen2Rescue', '230314_15genotype.txt'),
                       smoothOrNo=TRUE,
                       smooth_nsecs=30*60,
                       binOrNo=TRUE,
                       bin_nsecs=10*60,
                       grporder=c('SCR.DMSO', 'PSEN2.DMSO', 'PSEN2.Betamethasone'),
                       tracecols=c('#78ac63', '#697a87', '#fcb505'),
                       ribboncols=c('#bbd4ae', '#b3bcc3', '#fcd98a'),
                       # tracecols=c('#697a87', '#78ac63', '#fcb505'),
                       # ribboncols=c('#b3bcc3', '#bbd4ae', '#fcd98a'),
                       linethick=0.4,
                       xname='',
                       yname='',
                       xtextOrNo=FALSE,
                       ytextOrNo=TRUE,
                       ymin=0,
                       ymax=60000,
                       xstart=24,
                       xstop=69.8,
                       trimstart=24,
                       trimstop=72,
                       nightBgOrNo=TRUE,
                       legendOrNo=FALSE,
                       exportOrNo=TRUE,
                       exportPath=here('230314_psen2Rescue', 'pubplots', 'box15_acttrace_betamethasone_beforeAdj.pdf'),
                       width=68,
                       height=50)

### Fenoprofen treatment
ggActivityTraceByGroup(ffpath=here('230314_psen2Rescue', '230314_15_RAWs.csv'),
                       genopath=here('230314_psen2Rescue', '230314_15genotype.txt'),
                       smoothOrNo=TRUE,
                       smooth_nsecs=30*60,
                       binOrNo=TRUE,
                       bin_nsecs=10*60,
                       grporder=c('SCR.DMSO', 'PSEN2.DMSO', 'PSEN2.Fenoprofen'),
                       tracecols=c('#697a87', '#cb2a20', '#fcb505'),
                       ribboncols=c('#b3bcc3', '#e4928f', '#fcd98a'),
                       linethick=0.4,
                       xname='hours since ZT0',
                       yname='activity (sum of Δ px/10 min)',
                       xtextOrNo=FALSE,
                       ytextOrNo=TRUE,
                       ymin=0,
                       ymax=50000,
                       xstart=0,
                       xstop=0,
                       trimstart=0,
                       trimstop=0,
                       nightBgOrNo=TRUE,
                       legendOrNo=TRUE,
                       exportOrNo=TRUE,
                       exportPath=here('230314_psen2Rescue', 'plots', '230314_15_tracebygroup_Fenoprofen.pdf'),
                       width=75,
                       height=55)
# seems to go halfway there?


# BOX14: sleep trace by group ---------------------------------------------

ggSleepTraceByGroup(ffpath=here('230314_psen2Rescue', '230314_14_RAWs.csv'),
                    genopath=here('230314_psen2Rescue', '230314_14genotype.txt'),
                    epo_min=10,
                    smoothOrNo=TRUE,
                    smooth_npoints=5,
                    grporder=c('SCR.DMSO', 'PSEN2.DMSO'),
                    tracecols=c('#697a87', '#cb2a20'),
                    ribboncols=c('#b3bcc3', '#e4928f'),
                    linethick=0.4,
                    xname='',
                    yname='',
                    xtextOrNo=FALSE,
                    ytextOrNo=TRUE,
                    ymin=0,
                    ymax=10,
                    xstart=0,
                    xstop=0,
                    trimstart=0,
                    trimstop=0,
                    nightBgOrNo=TRUE,
                    legendOrNo=FALSE,
                    exportOrNo=TRUE,
                    exportPath=here('230314_psen2Rescue', 'plots', '230314_14_sleepbygroup_DMSO.pdf'),
                    width=75,
                    height=55)

ggSleepTraceByGroup(ffpath=here('230314_psen2Rescue', '230314_14_RAWs.csv'),
                    genopath=here('230314_psen2Rescue', '230314_14genotype.txt'),
                    epo_min=10,
                    smoothOrNo=TRUE,
                    smooth_npoints=5,
                    grporder=c('SCR.DMSO', 'PSEN2.DMSO', 'PSEN2.Tinidazole'),
                    tracecols=c('#697a87', '#cb2a20', '#fcb505'),
                    ribboncols=c('#b3bcc3', '#e4928f', '#fcd98a'),
                    linethick=0.4,
                    xname='',
                    yname='',
                    xtextOrNo=FALSE,
                    ytextOrNo=TRUE,
                    ymin=0,
                    ymax=10,
                    xstart=0,
                    xstop=0,
                    trimstart=0,
                    trimstop=0,
                    nightBgOrNo=TRUE,
                    legendOrNo=FALSE,
                    exportOrNo=TRUE,
                    exportPath=here('230314_psen2Rescue', 'plots', '230314_14_sleepbygroup_Tinidazole.pdf'),
                    width=75,
                    height=55)

ggSleepTraceByGroup(ffpath=here('230314_psen2Rescue', '230314_14_RAWs.csv'),
                    genopath=here('230314_psen2Rescue', '230314_14genotype.txt'),
                    epo_min=10,
                    smoothOrNo=TRUE,
                    smooth_npoints=5,
                    grporder=c('SCR.DMSO', 'PSEN2.DMSO', 'PSEN2.Betamethasone'),
                    tracecols=c('#78ac63', '#697a87', '#fcb505'),
                    ribboncols=c('#bbd4ae', '#b3bcc3', '#fcd98a'),
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
                    exportPath=here('230314_psen2Rescue', 'pubplots', 'box14_sleepTrace_betamethasone_beforeAdj.pdf'),
                    width=68,
                    height=50)


# BOX15: sleep trace by group ---------------------------------------------

ggSleepTraceByGroup(ffpath=here('230314_psen2Rescue', '230314_15_RAWs.csv'),
                    genopath=here('230314_psen2Rescue', '230314_15genotype.txt'),
                    epo_min=10,
                    smoothOrNo=TRUE,
                    smooth_npoints=5,
                    grporder=c('SCR.DMSO', 'PSEN2.DMSO'),
                    tracecols=c('#697a87', '#cb2a20'),
                    ribboncols=c('#b3bcc3', '#e4928f'),
                    linethick=0.4,
                    xname='',
                    yname='',
                    xtextOrNo=FALSE,
                    ytextOrNo=TRUE,
                    ymin=0,
                    ymax=10,
                    xstart=0,
                    xstop=0,
                    trimstart=0,
                    trimstop=0,
                    nightBgOrNo=TRUE,
                    legendOrNo=FALSE,
                    exportOrNo=TRUE,
                    exportPath=here('230314_psen2Rescue', 'plots', '230314_15_sleepbygroup_DMSO.pdf'),
                    width=75,
                    height=55)

ggSleepTraceByGroup(ffpath=here('230314_psen2Rescue', '230314_15_RAWs.csv'),
                    genopath=here('230314_psen2Rescue', '230314_15genotype.txt'),
                    epo_min=10,
                    smoothOrNo=TRUE,
                    smooth_npoints=5,
                    grporder=c('SCR.DMSO', 'PSEN2.DMSO', 'PSEN2.Tinidazole'),
                    tracecols=c('#697a87', '#cb2a20', '#fcb505'),
                    ribboncols=c('#b3bcc3', '#e4928f', '#fcd98a'),
                    linethick=0.4,
                    xname='',
                    yname='',
                    xtextOrNo=FALSE,
                    ytextOrNo=TRUE,
                    ymin=0,
                    ymax=10,
                    xstart=0,
                    xstop=0,
                    trimstart=0,
                    trimstop=0,
                    nightBgOrNo=TRUE,
                    legendOrNo=FALSE,
                    exportOrNo=TRUE,
                    exportPath=here('230314_psen2Rescue', 'plots', '230314_15_sleepbygroup_Tinidazole.pdf'),
                    width=75,
                    height=55)

## polished version ##
# no pixelAdjust
ggSleepTraceByGroup(ffpath=here('230314_psen2Rescue', '230314_15_RAWs.csv'),
                    genopath=here('230314_psen2Rescue', '230314_15genotype.txt'),
                    epo_min=10,
                    smoothOrNo=TRUE,
                    smooth_npoints=5,
                    grporder=c('SCR.DMSO', 'PSEN2.DMSO', 'PSEN2.Betamethasone'),
                    tracecols=c('#78ac63', '#697a87', '#fcb505'),
                    ribboncols=c('#bbd4ae', '#b3bcc3', '#fcd98a'),
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
                    exportPath=here('230314_psen2Rescue', 'pubplots', 'box15_sleepTrace_betamethasone_beforeAdj.pdf'),
                    width=68,
                    height=50)

ggSleepTraceByGroup(ffpath=here('230314_psen2Rescue', '230314_15_RAWs.csv'),
                    genopath=here('230314_psen2Rescue', '230314_15genotype.txt'),
                    epo_min=10,
                    smoothOrNo=TRUE,
                    smooth_npoints=5,
                    grporder=c('SCR.DMSO', 'PSEN2.DMSO', 'PSEN2.Fenoprofen'),
                    tracecols=c('#697a87', '#cb2a20', '#fcb505'),
                    ribboncols=c('#b3bcc3', '#e4928f', '#fcd98a'),
                    linethick=0.4,
                    xname='',
                    yname='',
                    xtextOrNo=FALSE,
                    ytextOrNo=TRUE,
                    ymin=0,
                    ymax=10,
                    xstart=0,
                    xstop=0,
                    trimstart=0,
                    trimstop=0,
                    nightBgOrNo=TRUE,
                    legendOrNo=FALSE,
                    exportOrNo=TRUE,
                    exportPath=here('230314_psen2Rescue', 'plots', '230314_15_sleepbygroup_Fenoprofen.pdf'),
                    width=75,
                    height=55)


# behavioural parameters --------------------------------------------------

multiBehaviourParameter(parameters='all',
                        ffpath=c(here('230314_psen2Rescue', '230314_14_RAWs.csv'),
                                 here('230314_psen2Rescue', '230314_15_RAWs.csv')),
                        genopath=c(here('230314_psen2Rescue', '230314_14genotype.txt'),
                                   here('230314_psen2Rescue', '230314_15genotype.txt')),
                        dayduration=14)

# parameter grid ----------------------------------------------------------

### BOX14
ggParameterGrid(paDir=here('230314_psen2Rescue', 'bhvparams'),
                grporder=c('SCR.DMSO', 'PSEN2.DMSO', 'PSEN2.Tinidazole', 'PSEN2.Betamethasone', 'PSEN2.Fenoprofen'),
                skipNight0=TRUE,
                colours=c('#697a87', '#cb2a20', '#78ac63', '#fcb505', '#417dcd'),
                legendOrNo=FALSE,
                ynameOrNo=FALSE,
                yunitOrNo=TRUE,
                xtextOrNo=TRUE,
                titleOrNo=TRUE,
                nightBgOrNo=TRUE,
                statsOrNo=TRUE,
                ncol=5,
                nrow=4,
                width=800,
                height=400,
                exportPath=here('230314_psen2Rescue', 'plots', '230314_14_grid.pdf'))



# fingerprints ------------------------------------------------------------

# rescue should appear as 'flattening' the fingerprint

calculateFingerprint(paDir=here('230314_psen2Rescue', 'bhvparams'),
                     controlGrp=c('SCR.DMSO'),
                     grporder=c('SCR.DMSO', 'PSEN2.DMSO', 'PSEN2.Tinidazole', 'PSEN2.Betamethasone', 'PSEN2.Fenoprofen'),
                     skipNight0=TRUE)

### BOX14
ggFingerprint(fgp=here('230314_psen2Rescue', 'fingerprint.csv'),
              grporder=c('SCR.DMSO', 'PSEN2.DMSO', 'PSEN2.Tinidazole', 'PSEN2.Betamethasone', 'PSEN2.Fenoprofen'),
              controlGrp='SCR.DMSO',
              removeControl=FALSE,
              onlyExp='230314_14',
              colours=c('#697a87', '#cb2a20', '#78ac63', '#fcb505', '#417dcd'),
              legendOrNo=TRUE,
              ynameOrNo=TRUE,
              ytextOrNo=TRUE,
              xtextOrNo=FALSE,
              xParamNum=TRUE,
              nightBgOrNo=TRUE,
              ymin=-4,
              ymax=4,
              exportOrNo=TRUE,
              exportPath=here('230314_psen2Rescue', 'plots', 'fgp14.pdf'),
              width=200,
              height=100)

### BOX15
ggFingerprint(fgp=here('230314_psen2Rescue', 'fingerprint.csv'),
              grporder=c('SCR.DMSO', 'PSEN2.DMSO', 'PSEN2.Tinidazole', 'PSEN2.Betamethasone', 'PSEN2.Fenoprofen'),
              controlGrp='SCR.DMSO',
              removeControl=FALSE,
              onlyExp='230314_15',
              colours=c('#697a87', '#cb2a20', '#78ac63', '#fcb505', '#417dcd'),
              legendOrNo=TRUE,
              ynameOrNo=TRUE,
              ytextOrNo=TRUE,
              xtextOrNo=TRUE,
              xParamNum=FALSE,
              nightBgOrNo=TRUE,
              ymin=-4,
              ymax=4,
              exportOrNo=TRUE,
              exportPath=here('230314_psen2Rescue', 'plots', 'fgp15.pdf'),
              width=200,
              height=100)
# can see quite clearly see which parameters are rescued and which are not

# could imagine counting *specificity* as "of those parameters which were *not* affected, how many are affected by the drug"
# and efficiency or something similar as how many affected parameters are rescued
# similar as ON/OFF target score



# statistics --------------------------------------------------------------

LMEreport(paDir=here('230314_psen2Rescue', 'bhvparams'),
          grporder=c('SCR.DMSO', 'PSEN2.DMSO', 'PSEN2.Tinidazole', 'PSEN2.Betamethasone', 'PSEN2.Fenoprofen'),
          skipNight0=TRUE,
          exportPath=here('230314_psen2Rescue', 'LMEreport.csv'))

LMEreport(paDir=here('230314_psen2Rescue', 'bhvparams14'),
          grporder=c('SCR.DMSO', 'PSEN2.DMSO', 'PSEN2.Tinidazole', 'PSEN2.Betamethasone', 'PSEN2.Fenoprofen'),
          skipNight0=TRUE,
          exportPath=here('230314_psen2Rescue', 'LMEreport14.csv'))
