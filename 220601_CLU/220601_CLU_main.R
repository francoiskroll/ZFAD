### main script experiment 220601_CLU ###

# scr = injected with 3x non-targeting gRNAs
# clu = clu F0 knockouts


# packages & functions ----------------------------------------------------

# load the Frame-by-Frame package

# source('D:\\FramebyFrame\\package\\FramebyFrame.R')
source('~/Dropbox/phd/framebyframe/FramebyFrame.R')

expdir <- here('220601_CLU/')

# create a new folder to store the plots we make
dir.create(here(expdir, 'plots'))


# QUALITY CHECKS ----------------------------------------------------------

### Framerate ###
# box16
ggFramerate(ffpath=here('220601_16_RAWs.csv'),
            zebpath=here('220601_16_17_clu.xls'),
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
            exportPath=here('plots/220601_16_framerate.pdf'))

# box17
ggFramerate(ffpath=here('220601_17_RAWs.csv'),
            zebpath=here('220601_16_17_clu.xls'),
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
            exportPath=here('plots/220601_17_framerate.pdf'))


# ACTIVITY GRIDS ----------------------------------------------------------

ggActivityTraceGrid(ffpath=here('220601_16_RAWs.csv'),
                    genopath=here('220601_16genotype.txt'),
                    zebpath=here('220601_16_17_clu.xls'),
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
                    exportPath=here('plots/220601_16_activitygrid.pdf'),
                    width=255,
                    height=171)

ggActivityTraceGrid(ffpath=here('220601_17_RAWs.csv'),
                    genopath=here('220601_17genotype.txt'),
                    zebpath=here('220601_16_17_clu.xls'),
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
                    exportPath=here('plots/220601_17_activitygrid.pdf'),
                    width=255,
                    height=171)


# ACTIVITY TRACES ---------------------------------------------------------

# box16
ggActivityTraceByGroup(ffpath=here('220601_16_RAWs.csv'),
                       genopath=here('220601_16genotype.txt'),
                       zebpath=here('220601_16_17_clu.xls'),
                       zebDeprecatedFormat=FALSE,
                       dayduration=14,
                       smoothOrNo=TRUE,
                       smooth_nsecs=30*60,
                       binOrNo=TRUE,
                       bin_nsecs=10*60,
                       ribbon='sem',
                       grporder=c('scr', 'clu'),
                       tracecols=c('#982150', '#697a87'),
                       ribboncols=c('#cb8fa8', '#b3bcc3'),
                       linethick=0.4,
                       xname='',
                       yname='',
                       xtextOrNo=FALSE,
                       ytextOrNo=TRUE,
                       ymin=0,
                       ymax=50000,
                       xstart=24,
                       xstop=72,
                       trimstart=24,
                       trimstop=72,
                       xmajorOrNo=FALSE,
                       ymajorOrNo=TRUE,
                       sunlinesOrNo=FALSE,
                       legendOrNo=FALSE,
                       exportOrNo=TRUE,
                       exportPath=here('plots/220601_16_tracebygroup.pdf'),
                       width=75,
                       height=55)

# box17
ggActivityTraceByGroup(ffpath=here('220601_17_RAWs.csv'),
                       genopath=here('220601_17genotype.txt'),
                       zebpath=here('220601_16_17_clu.xls'),
                       zebDeprecatedFormat=FALSE,
                       dayduration=14,
                       smoothOrNo=TRUE,
                       smooth_nsecs=30*60,
                       binOrNo=TRUE,
                       bin_nsecs=10*60,
                       ribbon='sem',
                       grporder=c('scr', 'clu'),
                       tracecols=c('#982150', '#697a87'),
                       ribboncols=c('#cb8fa8', '#b3bcc3'),
                       linethick=0.4,
                       xname='',
                       yname='',
                       xtextOrNo=FALSE,
                       ytextOrNo=TRUE,
                       ymin=0,
                       ymax=50000,
                       xstart=24,
                       xstop=72,
                       trimstart=24,
                       trimstop=72,
                       xmajorOrNo=FALSE,
                       ymajorOrNo=TRUE,
                       sunlinesOrNo=FALSE,
                       legendOrNo=FALSE,
                       exportOrNo=TRUE,
                       exportPath=here('plots/220601_17_tracebygroup.pdf'),
                       width=75,
                       height=55)



# SLEEP TRACES ------------------------------------------------------------

# box16
ggSleepTraceByGroup(ffpath=here('220601_16_RAWs.csv'),
                    genopath=here('220601_16genotype.txt'),
                    zebpath=here('220601_16_17_clu.xls'),
                    zebDeprecatedFormat=FALSE,
                    zthr_min=1,
                    epo_min=10,
                    dayduration=14,
                    smoothOrNo=TRUE, 
                    smooth_npoints=5,
                    ribbon='sem',
                    grporder=c('scr', 'clu'),
                    tracecols=c('#982150', '#697a87'),
                    ribboncols=c('#cb8fa8', '#b3bcc3'),
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
                    exportPath=here('plots/220601_16_sleepbygroup.pdf'),
                    width=75,
                    height=55)

# box17
ggSleepTraceByGroup(ffpath=here('220601_17_RAWs.csv'),
                    genopath=here('220601_17genotype.txt'),
                    zebpath=here('220601_16_17_clu.xls'),
                    zebDeprecatedFormat=FALSE,
                    zthr_min=1,
                    epo_min=10,
                    dayduration=14,
                    smoothOrNo=TRUE, 
                    smooth_npoints=5,
                    ribbon='sem',
                    grporder=c('scr', 'clu'),
                    tracecols=c('#982150', '#697a87'),
                    ribboncols=c('#cb8fa8', '#b3bcc3'),
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
                    exportPath=here('plots/220601_17_sleepbygroup.pdf'),
                    width=75,
                    height=55)




# CALCULATE ALL PARAMETERS ------------------------------------------------


multiBehaviourParameter(parameter='all',
                        ffpath=c(here('220601_CLU', '220601_16_RAWs.csv'), here('220601_CLU', '220601_17_RAWs.csv')),
                        genopath=c(here('220601_CLU', '220601_16genotype.txt'), here('220601_CLU', '220601_17genotype.txt')),
                        zebpath=here('220601_CLU', '220601_16_17_clu.xls'),
                        zebDeprecatedFormat=FALSE,
                        zthr_min=1,
                        dayduration=14)


# PARAMETER GRID ----------------------------------------------------------

# for LME report
ggParameterGrid(paDir=here('220601_CLU', 'bhvparams'),
                statsReport=TRUE,
                grporder=c('scr', 'clu'),
                skipNight0=TRUE,
                poolDayNight=FALSE,
                onlyExp=NA,
                onlyDayorNight=NA,
                onlyWin=NA,
                colours=c('#982150', '#697a87'),
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
                exportPath=here('220601_CLU', 'plots', '220601_grid.pdf'))

# for plot
ggParameterGrid(paDir=here('220601_CLU', 'bhvparams'),
                statsReport=FALSE, # note here
                grporder=c('clu', 'scr'),
                skipNight0=TRUE,
                poolDayNight=FALSE,
                onlyExp=NA,
                onlyDayorNight=NA,
                onlyWin=NA,
                colours=c('#982150', '#697a87'),
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
                exportPath=here('220601_CLU', 'plots', '220601_grid.pdf'))


# SLEEP LATENCY PLOTS -----------------------------------------------------


ggSleepLatencyGrid(pa=c(here('bhvparams/sleepLatency_220601_16.csv'), here('bhvparams/sleepLatency_220601_17.csv')),
                   grporder=c('clu', 'scr'),
                   skipNight0=TRUE,
                   colours=c('#982150', '#697a87'),
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