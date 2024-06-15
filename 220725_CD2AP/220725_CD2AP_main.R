### main script experiment 220725_16_17_CD2AP ###

# scr = injected with 3x non-targeting gRNAs
# cd2ap = cd2ap F0 knockouts

# experiment was interrupted for 17 min because of a power cut
# here, data after appending part1/part2

# remember the only use of the Zebralab XLS file is to get the experiment start date/time
# so here we should give Zebralab XLS file of part1


# packages & functions ----------------------------------------------------

# load the Frame-by-Frame package

# source('D:\\FramebyFrame\\package\\FramebyFrame.R')
source('~/Dropbox/phd/framebyframe/FramebyFrame.R')

expdir <- here('220725_CD2AP/')

# create a new folder to store the plots we make
dir.create(here('plots'))


# QUALITY CHECKS ----------------------------------------------------------

### Framerate ###
# box16
ggFramerate(ffpath=here('220725_16_RAWsapp.csv'),
            zebpath=here('220725_16_17_CD2AP_part1.xls'),
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
            exportPath=here('plots/220725_16_framerate.pdf'))

# box17
ggFramerate(ffpath=here('220725_17_RAWsapp.csv'),
            zebpath=here('220725_16_17_CD2AP_part1.xls'),
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
            exportPath=here('plots/220725_17_framerate.pdf'))

# we do see the interruption
# note, we would expect one frame at 1 frame per 17 min i.e. 0.0009 frames-per-second
# but it is expected that we do not see it as subsampling


# ACTIVITY GRIDS ----------------------------------------------------------

ggActivityTraceGrid(ffpath=here('220725_CD2AP', '220725_16_RAWs.csv'),
                    genopath=here('220725_CD2AP', '220725_16genotype.txt'),
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
                    exportPath=here('220725_CD2AP', 'plots', '220725_16_activitygrid.pdf'),
                    width=255,
                    height=171)

ggActivityTraceGrid(ffpath=here('220725_CD2AP', '220725_17_RAWs.csv'),
                    genopath=here('220725_CD2AP', '220725_17genotype.txt'),
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
                    exportPath=here('220725_CD2AP', 'plots', '220725_17_activitygrid.pdf'),
                    width=255,
                    height=171)


# ACTIVITY TRACES ---------------------------------------------------------

# box16
ggActivityTraceByGroup(ffpath=here('220725_16_RAWsapp.csv'),
                       genopath=here('220725_16genotype.txt'),
                       zebpath=here('220725_16_17_CD2AP_part1.xls'),
                       zebDeprecatedFormat=FALSE,
                       dayduration=14,
                       smoothOrNo=TRUE,
                       smooth_nsecs=30*60,
                       binOrNo=TRUE,
                       bin_nsecs=10*60,
                       ribbon='sem',
                       grporder=c('cd2ap', 'scr'),
                       tracecols=c('#2aa844', '#697a87'),
                       ribboncols=c('#96d3a1', '#b3bcc3'),
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
                       exportPath=here('plots/220725_16_tracebygroup.pdf'),
                       width=75,
                       height=55)

# box17
ggActivityTraceByGroup(ffpath=here('220725_17_RAWsapp.csv'),
                       genopath=here('220725_17genotype.txt'),
                       zebpath=here('220725_16_17_CD2AP_part1.xls'),
                       zebDeprecatedFormat=FALSE,
                       dayduration=14,
                       smoothOrNo=TRUE,
                       smooth_nsecs=30*60,
                       binOrNo=TRUE,
                       bin_nsecs=10*60,
                       ribbon='sem',
                       grporder=c('cd2ap', 'scr'),
                       tracecols=c('#2aa844', '#697a87'),
                       ribboncols=c('#96d3a1', '#b3bcc3'),
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
                       exportPath=here('plots/220725_17_tracebygroup.pdf'),
                       width=75,
                       height=55)



# SLEEP TRACES ------------------------------------------------------------

# box16
ggSleepTraceByGroup(ffpath=here('220725_16_RAWsapp.csv'),
                    genopath=here('220725_16genotype.txt'),
                    zebpath=here('220725_16_17_CD2AP_part1.xls'),
                    zebDeprecatedFormat=FALSE,
                    zthr_min=1,
                    epo_min=10,
                    dayduration=14,
                    smoothOrNo=TRUE, 
                    smooth_npoints=5,
                    ribbon='sem',
                    grporder=c('cd2ap', 'scr'),
                    tracecols=c('#2aa844', '#697a87'),
                    ribboncols=c('#96d3a1', '#b3bcc3'),
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
                    exportPath=here('plots/220725_16_sleepbygroup.pdf'),
                    width=75,
                    height=55)

# box17
ggSleepTraceByGroup(ffpath=here('220725_17_RAWsapp.csv'),
                    genopath=here('220725_17genotype.txt'),
                    zebpath=here('220725_16_17_CD2AP_part1.xls'),
                    zebDeprecatedFormat=FALSE,
                    zthr_min=1,
                    epo_min=10,
                    dayduration=14,
                    smoothOrNo=TRUE, 
                    smooth_npoints=5,
                    ribbon='sem',
                    grporder=c('cd2ap', 'scr'),
                    tracecols=c('#2aa844', '#697a87'),
                    ribboncols=c('#96d3a1', '#b3bcc3'),
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
                    exportPath=here('plots/220725_17_sleepbygroup.pdf'),
                    width=75,
                    height=55)



# CALCULATE ALL PARAMETERS ------------------------------------------------


multiBehaviourParameter(parameter='all',
                        ffpath=c(here('220725_CD2AP', '220725_16_RAWsapp.csv'), here('220725_CD2AP', '220725_17_RAWsapp.csv')),
                        genopath=c(here('220725_CD2AP', '220725_16genotype.txt'), here('220725_CD2AP', '220725_17genotype.txt')),
                        zebpath=here('220725_CD2AP', '220725_16_17_CD2AP_part1.xls'),
                        zebDeprecatedFormat=FALSE,
                        zthr_min=1,
                        dayduration=14)

# PARAMETER GRID ----------------------------------------------------------

 # for LME report
ggParameterGrid(paDir=here('220725_CD2AP', 'bhvparams'),
                statsReport=TRUE,
                grporder=c('scr', 'cd2ap'),
                skipNight0=TRUE,
                poolDayNight=FALSE,
                onlyExp=NA,
                onlyDayorNight=NA,
                onlyWin=NA,
                colours=c('#2aa844', '#697a87'),
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
                exportPath=here('220725_CD2AP', 'plots', '220725_grid.pdf'))

# for plot
ggParameterGrid(paDir=here('220725_CD2AP', 'bhvparams'),
                statsReport=FALSE, # note here
                grporder=c('cd2ap', 'scr'),
                skipNight0=TRUE,
                poolDayNight=FALSE,
                onlyExp=NA,
                onlyDayorNight=NA,
                onlyWin=NA,
                colours=c('#2aa844', '#697a87'),
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
                exportPath=here('220725_CD2AP', 'plots', '220725_grid.pdf'))


# SLEEP LATENCY PLOTS -----------------------------------------------------


ggSleepLatencyGrid(pa=c(here('bhvparams/sleepLatency_220725_16.csv'), here('bhvparams/sleepLatency_220725_17.csv')),
                   grporder=c('cd2ap', 'scr'),
                   skipNight0=TRUE,
                   colours=c('#2aa844', '#697a87'),
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
