### main script experiment 220313_APOEAB_2 ###

# scr4x = injected with 4x non-targeting gRNAs
# apoeab = apoeab F0 knockouts


# packages & functions ----------------------------------------------------

# load the Frame-by-Frame package

library(FramebyFrame)
library(here)

# create a new folder to store the plots we make
dir.create(here('plots'), showWarnings=FALSE)


# QUALITY CHECKS ----------------------------------------------------------

### Framerate ###
# box15
ggFramerate(ffpath=here('220313_15_RAWs.csv'),
            zebpath=here('220313_14_15_apoeab.xls'),
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
            exportPath=here('plots/220313_14_framerate.pdf'))


# ACTIVITY GRIDS ----------------------------------------------------------

ggActivityTraceGrid(ffpath=here('220313_APOEAB_2', '220313_15_RAWs.csv'),
                    genopath=here('220313_APOEAB_2', '220313_15genotype.txt'),
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
                    exportPath=here('220313_APOEAB_2', 'plots', '220313_15_activitygrid.pdf'),
                    width=255,
                    height=171)



# ACTIVITY TRACES ---------------------------------------------------------

# box15
ggActivityTraceByGroup(ffpath=here('220313_15_RAWs.csv'),
                       genopath=here('220313_15genotype.txt'),
                       zebpath=here('220313_14_15_apoeab.xls'),
                       zebDeprecatedFormat=FALSE,
                       dayduration=14,
                       smoothOrNo=TRUE,
                       smooth_nsecs=30*60,
                       binOrNo=TRUE,
                       bin_nsecs=10*60,
                       ribbon='sem',
                       grporder=c('apoeab', 'scr4x'),
                       tracecols=c('#7e7ba4', '#697a87'),
                       ribboncols=c('#bebdd1', '#b3bcc3'),
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
                       exportPath=here('plots/220313_15_tracebygroup.pdf'),
                       width=75,
                       height=55)



# SLEEP TRACES ------------------------------------------------------------

# box15
ggSleepTraceByGroup(ffpath=here('220313_15_RAWs.csv'),
                    genopath=here('220313_15genotype.txt'),
                    zebpath=here('220313_14_15_apoeab.xls'),
                    zebDeprecatedFormat=FALSE,
                    zthr_min=1,
                    epo_min=10,
                    dayduration=14,
                    smoothOrNo=TRUE, 
                    smooth_npoints=5,
                    ribbon='sem',
                    grporder=c('apoeab', 'scr4x'),
                    tracecols=c('#7e7ba4', '#697a87'),
                    ribboncols=c('#bebdd1', '#b3bcc3'),
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
                    exportPath=here('plots/220313_15_sleepbygroup.pdf'),
                    width=75,
                    height=55)



# CALCULATE ALL PARAMETERS ------------------------------------------------


multiBehaviourParameter(parameter='all',
                        ffpath=here('220313_APOEAB_2', '220313_15_RAWs.csv'),
                        genopath=here('220313_APOEAB_2', '220313_15genotype.txt'),
                        zebpath=here('220313_APOEAB_2', '220313_14_15_apoeab.xls'),
                        zebDeprecatedFormat=FALSE,
                        zthr_min=1,
                        dayduration=14)

# for grid, see 220516_APOEAB_3

# SLEEP LATENCY PLOTS -----------------------------------------------------

ggSleepLatencyGrid(pa=here('bhvparams/sleepLatency_220313_15.csv'),
                   grporder=c('apoeab', 'scr4x'),
                   skipNight0=TRUE,
                   colours=c('#7e7ba4', '#697a87'),
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
                   height=50,
                   dayduration=14)
