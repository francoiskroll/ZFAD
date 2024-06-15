### main script experiment 220516_APOEAB ###

# scr4x = injected with 4x non-targeting gRNAs
# apoeab = apoea & apoeb double F0 knockouts


# packages & functions ----------------------------------------------------

# load the Frame-by-Frame package

library(FramebyFrame)
library(here)

# create a new folder to store the plots we make
dir.create(here('plots'))

# QUALITY CHECKS ----------------------------------------------------------

### Framerate ###
# box14
ggFramerate(ffpath=here('220516_14_RAWs.csv'),
            zebpath=here('220516_14_15_apoeab.xls'),
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
            exportPath=here('plots/220516_14_framerate.pdf'))

# box15
ggFramerate(ffpath=here('220516_15_RAWs.csv'),
            zebpath=here('220516_14_15_apoeab.xls'),
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
            exportPath=here('plots/220516_15_framerate.pdf'))


# ACTIVITY GRIDS ----------------------------------------------------------

ggActivityTraceGrid(ffpath=here('220516_APOEAB_3', '220516_14_RAWs.csv'),
                    genopath=here('220516_APOEAB_3', '220516_14genotype.txt'),
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
                    exportPath=here('220516_APOEAB_3', 'plots', '220516_14_activitygrid.pdf'),
                    width=255,
                    height=171)

ggActivityTraceGrid(ffpath=here('220516_APOEAB_3', '220516_15_RAWs.csv'),
                    genopath=here('220516_APOEAB_3', '220516_15genotype.txt'),
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
                    exportPath=here('220516_APOEAB_3', 'plots', '220516_15_activitygrid.pdf'),
                    width=255,
                    height=171)



# ACTIVITY TRACES ---------------------------------------------------------

# box12
ggActivityTraceByGroup(ffpath=here('220516_14_RAWs.csv'),
                       genopath=here('220516_14genotype.txt'),
                       zebpath=here('220516_14_15_apoeab.xls'),
                       zebDeprecatedFormat=FALSE,
                       dayduration=14,
                       smoothOrNo=TRUE,
                       smooth_nsecs=30*60,
                       binOrNo=TRUE,
                       bin_nsecs=10*60,
                       ribbon='sem',
                       grporder=c('apoeab', 'scr4x'),
                       tracecols=c('#db5072', '#697a87'),
                       ribboncols=c('#eca6b8', '#b3bcc3'),
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
                       exportPath=here('plots/220516_14_tracebygroup.pdf'),
                       width=75,
                       height=55)

# box15
ggActivityTraceByGroup(ffpath=here('220516_15_RAWs.csv'),
                       genopath=here('220516_15genotype.txt'),
                       zebpath=here('220516_14_15_apoeab.xls'),
                       zebDeprecatedFormat=FALSE,
                       dayduration=14,
                       smoothOrNo=TRUE,
                       smooth_nsecs=30*60,
                       binOrNo=TRUE,
                       bin_nsecs=10*60,
                       ribbon='sem',
                       grporder=c('apoeab', 'scr4x'),
                       tracecols=c('#db5072', '#697a87'),
                       ribboncols=c('#eca6b8', '#b3bcc3'),
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
                       exportPath=here('plots/220516_15_tracebygroup.pdf'),
                       width=75,
                       height=55)



# SLEEP TRACES ------------------------------------------------------------

# box14
ggSleepTraceByGroup(ffpath=here('220516_14_RAWs.csv'),
                    genopath=here('220516_14genotype.txt'),
                    zebpath=here('220516_14_15_apoeab.xls'),
                    zebDeprecatedFormat=FALSE,
                    zthr_min=1,
                    epo_min=10,
                    dayduration=14,
                    smoothOrNo=TRUE, 
                    smooth_npoints=5,
                    ribbon='sem',
                    grporder=c('apoeab', 'scr4x'),
                    tracecols=c('#db5072', '#697a87'),
                    ribboncols=c('#eca6b8', '#b3bcc3'),
                    linethick=0.4,
                    xname='',
                    yname='',
                    xtextOrNo=FALSE,
                    ytextOrNo=TRUE,
                    ymin=0,
                    ymax=10,
                    xstart=24,
                    xstop=72,
                    trimstart=24,
                    trimstop=72,
                    xmajorOrNo=FALSE,
                    ymajorOrNo=TRUE,
                    sunlinesOrNo=FALSE,
                    legendOrNo=FALSE,
                    exportOrNo=TRUE,
                    exportPath=here('plots/220516_14_sleepbygroup.pdf'),
                    width=75,
                    height=55)

# box15
ggSleepTraceByGroup(ffpath=here('220516_15_RAWs.csv'),
                    genopath=here('220516_15genotype.txt'),
                    zebpath=here('220516_14_15_apoeab.xls'),
                    zebDeprecatedFormat=FALSE,
                    zthr_min=1,
                    epo_min=10,
                    dayduration=14,
                    smoothOrNo=TRUE, 
                    smooth_npoints=5,
                    ribbon='sem',
                    grporder=c('apoeab', 'scr4x'),
                    tracecols=c('#db5072', '#697a87'),
                    ribboncols=c('#eca6b8', '#b3bcc3'),
                    linethick=0.4,
                    xname='',
                    yname='',
                    xtextOrNo=FALSE,
                    ytextOrNo=TRUE,
                    ymin=0,
                    ymax=10,
                    xstart=24,
                    xstop=72,
                    trimstart=24,
                    trimstop=72,
                    xmajorOrNo=FALSE,
                    ymajorOrNo=TRUE,
                    sunlinesOrNo=FALSE,
                    legendOrNo=FALSE,
                    exportOrNo=TRUE,
                    exportPath=here('plots/220516_15_sleepbygroup.pdf'),
                    width=75,
                    height=55)


# CALCULATE ALL PARAMETERS ------------------------------------------------


multiBehaviourParameter(parameter='all',
                        ffpath=c(here('220516_APOEAB_3', '220516_14_RAWs.csv'), here('220516_APOEAB_3', '220516_15_RAWs.csv')),
                        genopath=c(here('220516_APOEAB_3', '220516_14genotype.txt'), here('220516_APOEAB_3', '220516_15genotype.txt')),
                        zebpath=here('220516_APOEAB_3', '220516_14_15_apoeab.xls'),
                        zebDeprecatedFormat=FALSE,
                        zthr_min=1,
                        dayduration=14)


# PARAMETER GRID ----------------------------------------------------------

ggParameterGrid(paDir=here('bhvparams/'),
                statsReport=TRUE,
                grporder=c('apoeab', 'scr4x'),
                skipNight0=TRUE,
                poolDayNight=FALSE,
                onlyExp=NA,
                onlyDayorNight=NA,
                onlyWin=NA,
                colours=c('#db5072', '#697a87'),
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
                exportPath=here('plots/220516_grid.pdf'))

# together with APOEAB_2 box2
# i.e. 220313_15 + 220516_14 + 220516_15

# for LMEreport
ggParameterGrid(paDir=c(here('220313_APOEAB_2', 'bhvparams'), here('220516_APOEAB_3', 'bhvparams')),
                statsReport=TRUE,
                grporder=c('scr4x', 'apoeab'),
                skipNight0=TRUE,
                poolDayNight=FALSE,
                onlyExp=NA,
                onlyDayorNight=NA,
                onlyWin=NA,
                colours=c('#db5072', '#697a87'),
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
                exportPath=here('220516_APOEAB_3', 'plots', '220516_grid_exp23.pdf'))

# for plot
ggParameterGrid(paDir=c(here('220313_APOEAB_2', 'bhvparams'), here('220516_APOEAB_3', 'bhvparams')),
                statsReport=FALSE,
                grporder=c('apoeab', 'scr4x'),
                skipNight0=TRUE,
                poolDayNight=FALSE,
                onlyExp=NA,
                onlyDayorNight=NA,
                onlyWin=NA,
                colours=c('#db5072', '#697a87'),
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
                exportPath=here('220516_APOEAB_3', 'plots', '220516_grid_exp23.pdf'))


# SLEEP LATENCY PLOTS -----------------------------------------------------


ggSleepLatencyGrid(pa=c(here('bhvparams/sleepLatency_220516_14.csv'), here('bhvparams/sleepLatency_220516_15.csv')),
                   grporder=c('apoeab', 'scr4x'),
                   skipNight0=TRUE,
                   colours=c('#db5072', '#697a87'),
                   legendOrNo=FALSE,
                   xtextOrNo=TRUE,
                   xnameOrNo=TRUE,
                   ynameOrNo=FALSE,
                   nightBgOrNo=TRUE,
                   xmaxDay=14,
                   xmaxNight=3,
                   detailsOrNo=FALSE,
                   exportDir=here('plots/'),
                   width=150,
                   height=100,
                   dayduration=14)