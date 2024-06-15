#####################################################
# ~ ZFAD: main script for experiment 210913_PSEN1~
#
#
# Francois Kroll 2023
# francois@kroll.be
#####################################################


# packages & functions ----------------------------------------------------

library(FramebyFrame)


# quality check -----------------------------------------------------------

# box12
ggFramerate(ffpath=here('210913_12_RAWs.csv'),
            zebpath=here('210913_12_13_PSEN1.xls'),
            zebDeprecatedFormat=FALSE,
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
            exportPath=here('plots/210913_12_framerate.pdf'))

# box13
ggFramerate(ffpath=here('210913_13_RAWs.csv'),
            zebpath=here('210913_12_13_PSEN1.xls'),
            zebDeprecatedFormat=FALSE,
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
            exportPath=here('plots/210913_13_framerate.pdf'))


# ACTIVITY GRIDS ----------------------------------------------------------

ggActivityTraceGrid(ffpath=here('210913_PSEN1', '210913_12_RAWs.csv'),
                    genopath=here('210913_PSEN1', '210913_12genotype.txt'),
                    dayduration=14,
                    smoothOrNo=TRUE,
                    smooth_nsecs=30*60,
                    binOrNo=TRUE,
                    bin_nsecs=10*60,
                    onlyWell=NA,
                    tracecols=NA,
                    linethick=0.4,
                    ymin=0,
                    ymax=60000,
                    xstart=0,
                    xstop=0,
                    trimstart=0,
                    trimstop=0,
                    xmajorOrNo=FALSE,
                    ymajorOrNo=TRUE,
                    markTimes=NA,
                    nightBgOrNo=TRUE,
                    sunlinesOrNo=FALSE,
                    ncol=12,
                    nrow=8,
                    exportOrNo=TRUE,
                    exportPath=here('210913_PSEN1', 'plots', '210913_12_activitygrid.pdf'),
                    width=255,
                    height=171)

ggActivityTraceGrid(ffpath=here('210913_PSEN1', '210913_13_RAWs.csv'),
                    genopath=here('210913_PSEN1', '210913_13genotype.txt'),
                    dayduration=14,
                    smoothOrNo=TRUE,
                    smooth_nsecs=30*60,
                    binOrNo=TRUE,
                    bin_nsecs=10*60,
                    onlyWell=NA,
                    tracecols=NA,
                    linethick=0.4,
                    ymin=0,
                    ymax=60000,
                    xstart=0,
                    xstop=0,
                    trimstart=0,
                    trimstop=0,
                    xmajorOrNo=FALSE,
                    ymajorOrNo=TRUE,
                    markTimes=NA,
                    nightBgOrNo=TRUE,
                    sunlinesOrNo=FALSE,
                    ncol=12,
                    nrow=8,
                    exportOrNo=TRUE,
                    exportPath=here('210913_PSEN1', 'plots', '210913_13_activitygrid.pdf'),
                    width=255,
                    height=171)


# activity traces ---------------------------------------------------------

# box12
ggActivityTraceByGroup(ffpath=here('210913_PSEN1', '210913_12_RAWs.csv'),
                       genopath=here('210913_PSEN1', '210913_12genotype.txt'),
                       dayduration=14,
                       smoothOrNo=TRUE,
                       smooth_nsecs=30*60,
                       binOrNo=TRUE,
                       bin_nsecs=10*60,
                       grporder=c('psen1', 'scr'),
                       tracecols=c('#fcb505', '#697a87'),
                       ribboncols=c('#fcd98a', '#b3bcc3'),
                       linethick=0.4,
                       xname='',
                       yname='',
                       xtextOrNo=FALSE,
                       ytextOrNo=TRUE,
                       ymin=0,
                       ymax=40000,
                       xstart=24,
                       xstop=72,
                       trimstart=24,
                       trimstop=72,
                       xmajorOrNo=FALSE,
                       ymajorOrNo=TRUE,
                       sunlinesOrNo=FALSE,
                       legendOrNo=FALSE,
                       exportOrNo=TRUE,
                       exportPath=here('210913_PSEN1', 'plots', '210913_12_tracebygroup.pdf'),
                       width=75,
                       height=55)

# box13
ggActivityTraceByGroup(ffpath=here('210913_PSEN1', '210913_13_RAWs.csv'),
                       genopath=here('210913_PSEN1', '210913_13genotype.txt'),
                       dayduration=14,
                       smoothOrNo=TRUE,
                       smooth_nsecs=30*60,
                       binOrNo=TRUE,
                       bin_nsecs=10*60,
                       grporder=c('psen1', 'scr'),
                       tracecols=c('#fcb505', '#697a87'),
                       ribboncols=c('#fcd98a', '#b3bcc3'),
                       linethick=0.4,
                       xname='',
                       yname='',
                       xtextOrNo=FALSE,
                       ytextOrNo=TRUE,
                       ymin=0,
                       ymax=40000,
                       xstart=24,
                       xstop=72,
                       trimstart=24,
                       trimstop=72,
                       xmajorOrNo=FALSE,
                       ymajorOrNo=TRUE,
                       sunlinesOrNo=FALSE,
                       legendOrNo=FALSE,
                       exportOrNo=TRUE,
                       exportPath=here('210913_PSEN1', 'plots', '210913_13_tracebygroup.pdf'),
                       width=75,
                       height=55)


# sleep traces ------------------------------------------------------------

ggSleepTraceByGroup(ffpath=here('210913_13_RAWs.csv'),
                    genopath=here('210913_13genotype.txt'),
                    zebpath=here('210913_12_13_PSEN1.xls'),
                    zebDeprecatedFormat=FALSE,
                    zthr_min=1,
                    epo_min=10,
                    dayduration=14,
                    smoothOrNo=TRUE, 
                    smooth_npoints=5,
                    ribbon='sem',
                    grporder=c('psen1', 'scr'),
                    tracecols=c('#fcb505', '#697a87'),
                    ribboncols=c('#fcd98a', '#b3bcc3'),
                    linethick=0.4,
                    xname='',
                    yname='',
                    xtextOrNo=FALSE,
                    ytextOrNo=TRUE,
                    ymin=0,
                    ymax=8,
                    xstart=24,
                    xstop=72,
                    trimstart=24,
                    trimstop=72,
                    xmajorOrNo=FALSE,
                    ymajorOrNo=TRUE,
                    sunlinesOrNo=FALSE,
                    legendOrNo=FALSE,
                    exportOrNo=TRUE,
                    exportPath=here('plots/210913_13_sleepbygroup.pdf'),
                    width=75,
                    height=55)





# CALCULATE ALL PARAMETERS ------------------------------------------------

multiBehaviourParameter(parameter='all',
                        ffpath=c(here('210913_PSEN1', '210913_12_RAWs.csv'), here('210913_PSEN1', '210913_13_RAWs.csv')),
                        genopath=c(here('210913_PSEN1', '210913_12genotype.txt'), here('210913_PSEN1', '210913_13genotype.txt')),
                        zebpath=here('210913_PSEN1', '210913_12_13_PSEN1.xls'),
                        zebDeprecatedFormat=FALSE,
                        zthr_min=1,
                        dayduration=14)



# PARAMETER GRID ----------------------------------------------------------

# for LMEreport
ggParameterGrid(paDir=here('210913_PSEN1', 'bhvparams/'),
                statsReport=TRUE,
                grporder=c('scr', 'psen1'),
                skipNight0=TRUE,
                poolDayNight=FALSE,
                onlyExp=NA,
                onlyDayorNight=NA,
                onlyWin=NA,
                colours=c('#fcb505', '#697a87'),
                ymin=NA,
                ymax=NA,
                legendOrNo=FALSE,
                ynameOrNo=FALSE,
                yunitOrNo=TRUE,
                xtextOrNo=FALSE,
                titleOrNo=TRUE,
                nightBgOrNo=TRUE,
                statsOrNo=TRUE,
                ncol=5,
                nrow=4,
                width=500,
                height=230,
                exportPath=here('210913_PSEN1', 'plots/210913_grid.pdf'))

# for plot
ggParameterGrid(paDir=here('210913_PSEN1', 'bhvparams/'),
                statsReport=FALSE, # note here
                grporder=c('psen1', 'scr'),
                skipNight0=TRUE,
                poolDayNight=FALSE,
                onlyExp=NA,
                onlyDayorNight=NA,
                onlyWin=NA,
                colours=c('#fcb505', '#697a87'),
                ymin=NA,
                ymax=NA,
                legendOrNo=FALSE,
                ynameOrNo=FALSE,
                yunitOrNo=TRUE,
                xtextOrNo=FALSE,
                titleOrNo=TRUE,
                nightBgOrNo=TRUE,
                statsOrNo=TRUE,
                ncol=5,
                nrow=4,
                width=500,
                height=230,
                exportPath=here('210913_PSEN1', 'plots/210913_grid.pdf'))



# SLEEP LATENCY PLOTS -----------------------------------------------------

dir.create(here('plots/sleeplatency/'))

ggSleepLatencyGrid(pa=c(here('bhvparams/sleepLatency_210913_12.csv'), here('bhvparams/sleepLatency_210913_13.csv')),
                   grporder=c('psen1', 'scr'),
                   skipNight0=TRUE,
                   colours=c('#fcb505', '#697a87'),
                   legendOrNo=FALSE,
                   xtextOrNo=TRUE,
                   xnameOrNo=TRUE,
                   ynameOrNo=FALSE,
                   nightBgOrNo=TRUE,
                   xmaxDay=14,
                   xmaxNight=3,
                   detailsOrNo=FALSE,
                   exportDir=here('plots/sleeplatency/'),
                   width=150,
                   height=100,
                   dayduration=14)