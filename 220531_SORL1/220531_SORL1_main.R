#####################################################
# ~ ZFAD: main script experiment 220531_SORL1 ~
#
# ... ... ...
# Francois Kroll 2023
# francois@kroll.be
#####################################################

# scr = injected with 3x non-targeting gRNAs
# sorl1 = sorl1 F0 knockouts


# packages & functions ----------------------------------------------------

# load the Frame-by-Frame package
library(FramebyFrame)
library(here)

# create a new folder to store the plots we make
dir.create(here('220531_SORL1', 'plots'), showWarnings=FALSE)
dir.create(here('220531_SORL1', 'pubplots'), showWarnings=FALSE)


# QUALITY CHECKS ----------------------------------------------------------

### Framerate ###
# box14
ggFramerate(ffpath=here('220531_14_RAWs.csv'),
            zebpath=here('220531_14_15_sorl1.xls'),
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
            exportPath=here('plots/220531_14_framerate.pdf'))

# box15
ggFramerate(ffpath=here('220531_15_RAWs.csv'),
            zebpath=here('220531_14_15_sorl1.xls'),
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
            exportPath=here('plots/220531_15_framerate.pdf'))


# ACTIVITY GRIDS ----------------------------------------------------------

ggActivityTraceGrid(ffpath=here('220531_14_RAWs.csv'),
                    genopath=here('220531_14genotype.txt'),
                    zebpath=here('220531_14_15_sorl1.xls'),
                    zebDeprecatedFormat=FALSE,
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
                    exportPath=here('plots/220531_14_activitygrid.pdf'),
                    width=255,
                    height=171)

ggActivityTraceGrid(ffpath=here('220531_15_RAWs.csv'),
                    genopath=here('220531_15genotype.txt'),
                    zebpath=here('220531_14_15_sorl1.xls'),
                    zebDeprecatedFormat=FALSE,
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
                    exportPath=here('plots/220531_15_activitygrid.pdf'),
                    width=255,
                    height=171)



# ACTIVITY TRACES ---------------------------------------------------------

# box14
ggActivityTraceByGroup(ffpath=here('220531_SORL1', '220531_14_RAWs.csv'),
                       genopath=here('220531_SORL1', '220531_14genotype.txt'),
                       smoothOrNo=TRUE,
                       smooth_nsecs=30*60,
                       binOrNo=TRUE,
                       bin_nsecs=10*60,
                       grporder=c('sorl1', 'scr'),
                       tracecols=c('#417dcd', '#697a87'),
                       ribboncols=c('#a3bbdb', '#b3bcc3'),
                       linethick=0.4,
                       xname='',
                       yname='activity (sum of Δ px/10 min)',
                       xtextOrNo=FALSE,
                       ytextOrNo=TRUE,
                       ymin=0,
                       ymax=50000,
                       xstart=24,
                       xstop=69.8,
                       trimstart=24,
                       trimstop=72,
                       nightBgOrNo=TRUE,
                       legendOrNo=FALSE,
                       exportOrNo=TRUE,
                       exportPath=here('220531_SORL1', 'pubplots', '220531_14_tracebygroup.pdf'),
                       width=68,
                       height=50)

# box15
ggActivityTraceByGroup(ffpath=here('220531_SORL1', '220531_15_RAWs.csv'),
                       genopath=here('220531_SORL1', '220531_15genotype.txt'),
                       smoothOrNo=TRUE,
                       smooth_nsecs=30*60,
                       binOrNo=TRUE,
                       bin_nsecs=10*60,
                       grporder=c('sorl1', 'scr'),
                       tracecols=c('#417dcd', '#697a87'),
                       ribboncols=c('#a3bbdb', '#b3bcc3'),
                       linethick=0.4,
                       xname='',
                       yname='activity (sum of Δ px/10 min)',
                       xtextOrNo=FALSE,
                       ytextOrNo=TRUE,
                       ymin=0,
                       ymax=50000,
                       xstart=24,
                       xstop=69.8,
                       trimstart=24,
                       trimstop=72,
                       nightBgOrNo=TRUE,
                       legendOrNo=FALSE,
                       exportOrNo=TRUE,
                       exportPath=here('220531_SORL1', 'pubplots', '220531_15_tracebygroup.pdf'),
                       width=68,
                       height=50)



# SLEEP TRACES ------------------------------------------------------------

# box14
ggSleepTraceByGroup(ffpath=here('220531_SORL1', '220531_14_RAWs.csv'),
                    genopath=here('220531_SORL1', '220531_14genotype.txt'),
                    epo_min=10,
                    smoothOrNo=TRUE, 
                    smooth_npoints=5,
                    ribbon='sem',
                    grporder=c('sorl1', 'scr'),
                    tracecols=c('#417dcd', '#697a87'),
                    ribboncols=c('#a3bbdb', '#b3bcc3'),
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
                    exportPath=here('220531_SORL1', 'pubplots', '220531_14_sleepbygroup.pdf'),
                    width=68,
                    height=50)

# box15
ggSleepTraceByGroup(ffpath=here('220531_SORL1', '220531_15_RAWs.csv'),
                    genopath=here('220531_SORL1', '220531_15genotype.txt'),
                    epo_min=10,
                    smoothOrNo=TRUE, 
                    smooth_npoints=5,
                    ribbon='sem',
                    grporder=c('sorl1', 'scr'),
                    tracecols=c('#417dcd', '#697a87'),
                    ribboncols=c('#a3bbdb', '#b3bcc3'),
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
                    exportPath=here('220531_SORL1', 'pubplots', '220531_15_sleepbygroup.pdf'),
                    width=68,
                    height=50)



# CALCULATE ALL PARAMETERS ------------------------------------------------


multiBehaviourParameter(parameter='all',
                        ffpath=c(here('220531_SORL1', '220531_14_RAWs.csv'), here('220531_SORL1', '220531_15_RAWs.csv')),
                        genopath=c(here('220531_SORL1', '220531_14genotype.txt'), here('220531_SORL1', '220531_15genotype.txt')),
                        zebpath=here('220531_SORL1', '220531_14_15_sorl1.xls'),
                        zebDeprecatedFormat=FALSE,
                        zthr_min=1,
                        dayduration=14)


# PARAMETER GRID ----------------------------------------------------------

# for LME report
ggParameterGrid(paDir=here('220531_SORL1', 'bhvparams'),
                statsReport=TRUE,
                grporder=c('scr', 'sorl1'),
                skipNight0=TRUE,
                poolDayNight=FALSE,
                onlyExp=NA,
                onlyDayorNight=NA,
                onlyWin=NA,
                colours=c('#417dcd', '#697a87'),
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
                exportPath=here('220531_SORL1', 'plots', '220531_grid.pdf'))

# for plot
ggParameterGrid(paDir=here('220531_SORL1', 'bhvparams'),
                statsReport=FALSE, # note here
                grporder=c('sorl1', 'scr'),
                skipNight0=TRUE,
                poolDayNight=FALSE,
                onlyExp=NA,
                onlyDayorNight=NA,
                onlyWin=NA,
                colours=c('#417dcd', '#697a87'),
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
                exportPath=here('220531_SORL1', 'plots', '220531_grid.pdf'))

# for ppt
ggParameter(pa=c(here('220531_SORL1', 'bhvparams', 'activityPercentageTimeActive_220531_14.csv'),
                 here('220531_SORL1', 'bhvparams', 'activityPercentageTimeActive_220531_15.csv')),
            grporder=c('sorl1', 'scr'),
            skipNight0=TRUE,
            poolDayNight=FALSE,
            onlyDayorNight='day',
            colours=c('#417dcd', '#697a87'),
            ymin=NA,
            ymax=NA,
            legendOrNo=FALSE,
            xtextOrNo=FALSE,
            ynameOrNo=FALSE,
            yunitOrNo=TRUE,
            titleOrNo=TRUE,
            blankTitle=TRUE,
            nightBgOrNo=TRUE,
            statsOrNo=TRUE,
            silent=TRUE,
            width=60,
            height=75,
            exportPath=here('220531_SORL1', 'plots', 'percentageTimeActive_day.pdf'))

ggParameter(pa=c(here('220531_SORL1', 'bhvparams', 'sleepNumNaps_220531_14.csv'),
                 here('220531_SORL1', 'bhvparams', 'sleepNumNaps_220531_15.csv')),
            grporder=c('sorl1', 'scr'),
            skipNight0=TRUE,
            poolDayNight=FALSE,
            onlyDayorNight='night',
            colours=c('#417dcd', '#697a87'),
            ymin=NA,
            ymax=NA,
            legendOrNo=FALSE,
            xtextOrNo=FALSE,
            ynameOrNo=FALSE,
            yunitOrNo=TRUE,
            titleOrNo=TRUE,
            blankTitle=TRUE,
            nightBgOrNo=TRUE,
            statsOrNo=TRUE,
            silent=TRUE,
            width=60,
            height=75,
            exportPath=here('220531_SORL1', 'plots', 'sleepNumNaps_night.pdf'))


# SLEEP LATENCY PLOTS -----------------------------------------------------


ggSleepLatencyGrid(pa=c(here('bhvparams/sleepLatency_220531_14.csv'), here('bhvparams/sleepLatency_220531_15.csv')),
                   grporder=c('sorl1', 'scr'),
                   skipNight0=TRUE,
                   colours=c('#417dcd', '#697a87'),
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