### main script experiment 220524_APPAB ###

# scr = injected with 4x non-targeting gRNAs
# appab = appab F0 knockouts


# packages & functions ----------------------------------------------------

# load the Frame-by-Frame package

library(FramebyFrame)

# create a new folder to store the plots we make
dir.create(here('plots'))


# QUALITY CHECKS ----------------------------------------------------------

### Framerate ###
# box14
ggFramerate(ffpath=here('220524_14_RAWs.csv'),
            zebpath=here('220524_14_15_appab.xls'),
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
            exportPath=here('plots/220524_14_framerate.pdf'))

# box15
ggFramerate(ffpath=here('220524_15_RAWs.csv'),
            zebpath=here('220524_14_15_appab.xls'),
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
            exportPath=here('plots/220524_15_framerate.pdf'))


# ACTIVITY GRIDS ----------------------------------------------------------

ggActivityTraceGrid(ffpath=here('220524_APPAB', '220524_14_RAWs.csv'),
                    genopath=here('220524_APPAB', '220524_14genotype.txt'),
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
                    exportPath=here('220524_APPAB', 'plots', '220524_14_activitygrid.pdf'),
                    width=255,
                    height=171)

ggActivityTraceGrid(ffpath=here('220524_APPAB', '220524_15_RAWs.csv'),
                    genopath=here('220524_APPAB', '220524_15genotype.txt'),
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
                    exportPath=here('220524_APPAB', 'plots', '220524_15_activitygrid.pdf'),
                    width=255,
                    height=171)



# ACTIVITY TRACES ---------------------------------------------------------

# box14
ggActivityTraceByGroup(ffpath=here('220524_14_RAWs.csv'),
                       genopath=here('220524_14genotype.txt'),
                       zebpath=here('220524_14_15_appab.xls'),
                       zebDeprecatedFormat=FALSE,
                       dayduration=14,
                       smoothOrNo=TRUE,
                       smooth_nsecs=30*60,
                       binOrNo=TRUE,
                       bin_nsecs=10*60,
                       ribbon='sem',
                       grporder=c('appab', 'scr4x'),
                       tracecols=c('#cb2a20', '#697a87'),
                       ribboncols=c('#e4928f', '#b3bcc3'),
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
                       exportPath=here('plots/220524_14_tracebygroup.pdf'),
                       width=75,
                       height=55)

# box15
ggActivityTraceByGroup(ffpath=here('220524_15_RAWs.csv'),
                       genopath=here('220524_15genotype.txt'),
                       zebpath=here('220524_14_15_appab.xls'),
                       zebDeprecatedFormat=FALSE,
                       dayduration=14,
                       smoothOrNo=TRUE,
                       smooth_nsecs=30*60,
                       binOrNo=TRUE,
                       bin_nsecs=10*60,
                       ribbon='sem',
                       grporder=c('appab', 'scr4x'),
                       tracecols=c('#cb2a20', '#697a87'),
                       ribboncols=c('#e4928f', '#b3bcc3'),
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
                       exportPath=here('plots/220524_15_tracebygroup.pdf'),
                       width=75,
                       height=55)



# SLEEP TRACES ------------------------------------------------------------

# box14
ggSleepTraceByGroup(ffpath=here('220524_14_RAWs.csv'),
                    genopath=here('220524_14genotype.txt'),
                    zebpath=here('220524_14_15_appab.xls'),
                    zebDeprecatedFormat=FALSE,
                    zthr_min=1,
                    epo_min=10,
                    dayduration=14,
                    smoothOrNo=TRUE, 
                    smooth_npoints=5,
                    ribbon='sem',
                    grporder=c('appab', 'scr4x'),
                    tracecols=c('#cb2a20', '#697a87'),
                    ribboncols=c('#e4928f', '#b3bcc3'),
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
                    exportPath=here('plots/220524_14_sleepbygroup.pdf'),
                    width=75,
                    height=55)

# box15
ggSleepTraceByGroup(ffpath=here('220524_15_RAWs.csv'),
                    genopath=here('220524_15genotype.txt'),
                    zebpath=here('220524_14_15_appab.xls'),
                    zebDeprecatedFormat=FALSE,
                    zthr_min=1,
                    epo_min=10,
                    dayduration=14,
                    smoothOrNo=TRUE, 
                    smooth_npoints=5,
                    ribbon='sem',
                    grporder=c('appab', 'scr4x'),
                    tracecols=c('#cb2a20', '#697a87'),
                    ribboncols=c('#e4928f', '#b3bcc3'),
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
                    exportPath=here('plots/220524_15_sleepbygroup.pdf'),
                    width=75,
                    height=55)




# CALCULATE ALL PARAMETERS ------------------------------------------------


multiBehaviourParameter(parameter='all',
                        ffpath=c(here('220524_APPAB', '220524_14_RAWs.csv'), here('220524_APPAB', '220524_15_RAWs.csv')),
                        genopath=c(here('220524_APPAB', '220524_14genotype.txt'), here('220524_APPAB', '220524_15genotype.txt')),
                        zebpath=here('220524_APPAB', '220524_14_15_appab.xls'),
                        zebDeprecatedFormat=FALSE,
                        zthr_min=1,
                        dayduration=14)


# PARAMETER GRID ----------------------------------------------------------

# for LME report
ggParameterGrid(paDir=here('220524_APPAB', 'bhvparams/'),
                statsReport=TRUE,
                grporder=c('scr4x', 'appab'),
                skipNight0=TRUE,
                poolDayNight=FALSE,
                onlyExp=NA,
                onlyDayorNight=NA,
                onlyWin=NA,
                colours=c('#cb2a20', '#697a87'),
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
                exportPath=here('220524_APPAB', 'plots', '220524_grid.pdf'))

# for plot
ggParameterGrid(paDir=here('220524_APPAB', 'bhvparams'),
                statsReport=FALSE, # note here
                grporder=c('appab', 'scr4x'),
                skipNight0=TRUE,
                poolDayNight=FALSE,
                onlyExp=NA,
                onlyDayorNight=NA,
                onlyWin=NA,
                colours=c('#cb2a20', '#697a87'),
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
                exportPath=here('220524_APPAB', 'plots', '220524_grid.pdf'))


# SLEEP LATENCY PLOTS -----------------------------------------------------

dir.create(here('plots/sleeplatency/'))

# day
ggSleepLatencySurvival(pa=c(here('bhvparams/sleepLatency_220524_14.csv'), here('bhvparams/sleepLatency_220524_15.csv')),
                       grporder=c('appab', 'scr4x'),
                       skipNight0=FALSE,
                       onlyDayorNight='day',
                       onlyWin=NA,
                       colours=c('#cb2a20', '#697a87'),
                       legendOrNo=FALSE,
                       xtextOrNo=TRUE,
                       xnameOrNo=TRUE,
                       ynameOrNo=TRUE,
                       nightBgOrNo=TRUE,
                       xmaxh=14,
                       detailsOrNo=FALSE,
                       exportDir=here('plots/sleeplatency/'),
                       width=80,
                       height=50,
                       dayduration=14)
# night
ggSleepLatencySurvival(pa=c(here('bhvparams/sleepLatency_220524_14.csv'), here('bhvparams/sleepLatency_220524_15.csv')),
                       grporder=c('appab', 'scr4x'),
                       skipNight0=FALSE,
                       onlyDayorNight='night',
                       onlyWin=NA,
                       colours=c('#cb2a20', '#697a87'),
                       legendOrNo=FALSE,
                       xtextOrNo=TRUE,
                       xnameOrNo=TRUE,
                       ynameOrNo=TRUE,
                       nightBgOrNo=TRUE,
                       xmaxh=3,
                       detailsOrNo=FALSE,
                       exportDir=here('plots/sleeplatency/'),
                       width=80,
                       height=50,
                       dayduration=14)