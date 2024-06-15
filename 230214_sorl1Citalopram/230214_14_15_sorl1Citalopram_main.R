#####################################################
# ~ ZFAD: main script for experiment 230214_sorl1Citalopram ~
#
# Francois Kroll 2023
# francois@kroll.be
#####################################################

# packages ----------------------------------------------------------------

library(here)

# library(devtools)
# install_github('francoiskroll/FramebyFrame')

library(FramebyFrame)

dir.create(here('230214_sorl1Citalopram', 'plots'), showWarnings=FALSE)
dir.create(here('230214_sorl1Citalopram', 'pubplots'), showWarnings=FALSE)

source(here('utilities', 'percEffectSize.R'))


# vpSorter ----------------------------------------------------------------
# ran on common Rihel lab computer


# genotype file -----------------------------------------------------------

genotypeGenerator(plateMap=here('230214_sorl1Citalopram', '230214_14_genotypeMap.xlsx'))
genotypeGenerator(plateMap=here('230214_sorl1Citalopram', '230214_15_genotypeMap1102.xlsx'))
genotypeGenerator(plateMap=here('230214_sorl1Citalopram', '230214_15_genotypeMap1031.xlsx'))


# frame rate --------------------------------------------------------------

ggFramerate(ffpath=here('230214_sorl1Citalopram', '230214_14_RAWs.csv'),
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
            exportPath=here('230214_sorl1Citalopram', 'plots', '230214_14_framerate.pdf'))

ggFramerate(ffpath=here('230214_sorl1Citalopram', '230214_15_RAWs.csv'),
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
            exportPath=here('230214_sorl1Citalopram', 'plots', '230214_15_framerate.pdf'))
# there is a framerate issue for box15; from ~ 25 to 48



# activity grid -----------------------------------------------------------

ggActivityTraceGrid(ffpath=here('230214_sorl1Citalopram', '230214_14_RAWs.csv'),
                    genopath=here('230214_sorl1Citalopram', '230214_14genotype.txt'),
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
                    exportPath=here('230214_sorl1Citalopram', 'plots', '230214_14_grid.pdf'),
                    width=255,
                    height=171)

ggActivityTraceGrid(ffpath=here('230214_sorl1Citalopram', '230214_15_RAWs.csv'),
                    genopath=here('230214_sorl1Citalopram', '230214_15_1102genotype.txt'),
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
                    exportPath=here('230214_sorl1Citalopram', 'plots', '230214_15_1102actgrid.pdf'),
                    width=255,
                    height=171)

ggActivityTraceGrid(ffpath=here('230214_sorl1Citalopram', '230214_15_RAWs.csv'),
                    genopath=here('230214_sorl1Citalopram', '230214_15_1031genotype.txt'),
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
                    exportPath=here('230214_sorl1Citalopram', 'plots', '230214_15_1031actgrid.pdf'),
                    width=255,
                    height=171)



# box14: activity trace by group ------------------------------------------

ggActivityTraceByGroup(ffpath=here('230214_sorl1Citalopram', '230214_14_RAWs.csv'),
                       genopath=here('230214_sorl1Citalopram', '230214_14genotype.txt'),
                       smoothOrNo=TRUE,
                       smooth_nsecs=30*60,
                       binOrNo=TRUE,
                       bin_nsecs=10*60,
                       grporder=c('sorl1_0', 'scr_0'),
                       tracecols=c('#EE7163', '#697a87'),
                       ribboncols=c('#EE7163', '#697a87'),
                       linethick=0.4,
                       xname='hours since ZT0',
                       yname='activity (sum of Δ px/10 min)',
                       xtextOrNo=FALSE,
                       ytextOrNo=TRUE,
                       ymin=0,
                       ymax=50000,
                       xstart=24,
                       xstop=72,
                       trimstart=24,
                       trimstop=72,
                       nightBgOrNo=TRUE,
                       legendOrNo=FALSE,
                       exportOrNo=TRUE,
                       exportPath=here('230214_sorl1Citalopram', 'plots', '230214_14_tracebygroup_sorl1scr.pdf'),
                       width=75,
                       height=55)

ggActivityTraceByGroup(ffpath=here('230214_sorl1Citalopram', '230214_14_RAWs.csv'),
                       genopath=here('230214_sorl1Citalopram', '230214_14genotype.txt'),
                       smoothOrNo=TRUE,
                       smooth_nsecs=30*60,
                       binOrNo=TRUE,
                       bin_nsecs=10*60,
                       grporder=c('scr_0', 'scr_1', 'scr_10'),
                       tracecols=c('#b4bdc4', '#8d9ca7', '#697a87'),
                       ribboncols=c('#b4bdc4', '#8d9ca7', '#697a87'),
                       linethick=0.4,
                       xname='hours since ZT0',
                       yname='activity (sum of Δ px/10 min)',
                       xtextOrNo=FALSE,
                       ytextOrNo=TRUE,
                       ymin=0,
                       ymax=50000,
                       xstart=24,
                       xstop=72,
                       trimstart=24,
                       trimstop=72,
                       nightBgOrNo=TRUE,
                       legendOrNo=FALSE,
                       exportOrNo=TRUE,
                       exportPath=here('230214_sorl1Citalopram', 'plots', '230214_14_tracebygroup_scr.pdf'),
                       width=75,
                       height=55)

ggActivityTraceByGroup(ffpath=here('230214_sorl1Citalopram', '230214_14_RAWs.csv'),
                       genopath=here('230214_sorl1Citalopram', '230214_14genotype.txt'),
                       smoothOrNo=TRUE,
                       smooth_nsecs=30*60,
                       binOrNo=TRUE,
                       bin_nsecs=10*60,
                       grporder=c('sorl1_0', 'sorl1_1', 'sorl1_10'),
                       tracecols=c('#b4bdc4', '#8d9ca7', '#697a87'),
                       ribboncols=c('#b4bdc4', '#8d9ca7', '#697a87'),
                       linethick=0.4,
                       xname='hours since ZT0',
                       yname='activity (sum of Δ px/10 min)',
                       xtextOrNo=FALSE,
                       ytextOrNo=TRUE,
                       ymin=0,
                       ymax=50000,
                       xstart=24,
                       xstop=72,
                       trimstart=24,
                       trimstop=72,
                       nightBgOrNo=TRUE,
                       legendOrNo=FALSE,
                       exportOrNo=TRUE,
                       exportPath=here('230214_sorl1Citalopram', 'plots', '230214_14_tracebygroup_sorl1.pdf'),
                       width=75,
                       height=55)



# sleep trace by group ----------------------------------------------------

ggSleepTraceByGroup(ffpath=here('230214_sorl1Citalopram', '230214_14_RAWs.csv'),
                    genopath=here('230214_sorl1Citalopram', '230214_14genotype.txt'),
                    epo_min=10,
                    smoothOrNo=TRUE,
                    smooth_npoints=5,
                    grporder=c('sorl1_0', 'scr_0'),
                    tracecols=c('#EE7163', '#697a87'),
                    ribboncols=c('#EE7163', '#697a87'),
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
                    nightBgOrNo=TRUE,
                    legendOrNo=FALSE,
                    exportOrNo=TRUE,
                    exportPath=here('230214_sorl1Citalopram', 'plots', '230214_14_sleepbygroup_sorl1scr.pdf'),
                    width=75,
                    height=55)

ggSleepTraceByGroup(ffpath=here('230214_sorl1Citalopram', '230214_14_RAWs.csv'),
                    genopath=here('230214_sorl1Citalopram', '230214_14genotype.txt'),
                    epo_min=10,
                    smoothOrNo=TRUE,
                    smooth_npoints=5,
                    grporder=c('scr_0', 'scr_1', 'scr_10'),
                    tracecols=c('#b4bdc4', '#8d9ca7', '#697a87'),
                    ribboncols=c('#b4bdc4', '#8d9ca7', '#697a87'),
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
                    nightBgOrNo=TRUE,
                    legendOrNo=FALSE,
                    exportOrNo=TRUE,
                    exportPath=here('230214_sorl1Citalopram', 'plots', '230214_14_sleepbygroup_scr.pdf'),
                    width=75,
                    height=55)

ggSleepTraceByGroup(ffpath=here('230214_sorl1Citalopram', '230214_14_RAWs.csv'),
                    genopath=here('230214_sorl1Citalopram', '230214_14genotype.txt'),
                    epo_min=10,
                    smoothOrNo=TRUE,
                    smooth_npoints=5,
                    grporder=c('sorl1_0', 'sorl1_1', 'sorl1_10'),
                    tracecols=c('#b4bdc4', '#8d9ca7', '#697a87'),
                    ribboncols=c('#b4bdc4', '#8d9ca7', '#697a87'),
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
                    nightBgOrNo=TRUE,
                    legendOrNo=FALSE,
                    exportOrNo=TRUE,
                    exportPath=here('230214_sorl1Citalopram', 'plots', '230214_14_sleepbygroup_sorl1.pdf'),
                    width=75,
                    height=55)


# behavioural parameters --------------------------------------------------

multiBehaviourParameter(parameters='all',
                        ffpath=here('230214_sorl1Citalopram', '230214_14_RAWs.csv'),
                        genopath=here('230214_sorl1Citalopram', '230214_14genotype.txt'),
                        dayduration=14)

multiBehaviourParameter(parameters='all',
                        ffpath=here('230214_sorl1Citalopram', '230214_15_RAWs.csv'),
                        genopath=here('230214_sorl1Citalopram', '230214_15_1102genotype.txt'),
                        dayduration=14)

multiBehaviourParameter(parameters='all',
                        ffpath=here('230214_sorl1Citalopram', '230214_15_RAWs.csv'),
                        genopath=here('230214_sorl1Citalopram', '230214_15_1031genotype.txt'),
                        dayduration=14)

# parameter grid ----------------------------------------------------------

### box14
ggParameterGrid(paDir=here('230214_sorl1Citalopram', 'bhvparams14'),
                grporder=c('scr_0', 'scr_1', 'scr_10', 'sorl1_0', 'sorl1_1', 'sorl1_10'),
                skipNight0=TRUE,
                colours=c('#b4bdc4', '#8d9ca7', '#697a87', '#92b4e2', '#6999d8', '#417dcd'),
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
                exportPath=here('230214_sorl1Citalopram', 'plots', '230214_14_grid.pdf'))

ggParameterGrid(paDir=here('230214_sorl1Citalopram', 'bhvparams15_1031'),
                grporder=c('scr_0', 'scr_1', 'scr_10', 'sorl1_0', 'sorl1_1', 'sorl1_10'),
                skipNight0=TRUE,
                colours=c('#b4bdc4', '#8d9ca7', '#697a87', '#92b4e2', '#6999d8', '#417dcd'),
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
                exportPath=here('230214_sorl1Citalopram', 'plots', '230214_15_1031grid.pdf'))


# calculate fingerprints --------------------------------------------------

calculateFingerprint(paDir=here('230214_sorl1Citalopram', 'bhvparams'),
                     controlGrp=c('scr_0'),
                     grporder=c('scr_0', 'scr_1', 'scr_10'),
                     skipNight0=TRUE)
# changed manually filename fingerprint.csv to fingerprint_scr.csv

calculateFingerprint(paDir=here('230214_sorl1Citalopram', 'bhvparams'),
                     controlGrp=c('sorl1_0'),
                     grporder=c('sorl1_0', 'sorl1_1', 'sorl1_10'),
                     skipNight0=TRUE)
# changed manually filename fingerprint.csv to fingerprint_sorl1.csv


# box14 fingerprint -------------------------------------------------------

ggFingerprint(fgp=c(here('230214_sorl1Citalopram', 'fingerprint_scr.csv'),
                    here('230214_sorl1Citalopram', 'fingerprint_sorl1.csv')),
              grporder=c('scr_1', 'sorl1_1'),
              controlGrp=c('scr_0', 'sorl1_0'),
              removeControl=TRUE,
              onlyExp='230214_14',
              colours=c('#b3bcc3', '#F6B8B1'),
              legendOrNo=FALSE,
              ynameOrNo=FALSE,
              ytextOrNo=TRUE,
              xtextOrNo=FALSE,
              xParamNum=TRUE,
              nightBgOrNo=TRUE,
              ymin=-8,
              ymax=8,
              dotSize=0.1,
              lineSize=0.3,
              exportOrNo=TRUE,
              exportPath=here('230214_sorl1Citalopram', 'pubplots', 'fgp_14dose1.pdf'),
              width=90,
              height=28)

ggFingerprint(fgp=c(here('230214_sorl1Citalopram', 'fingerprint_scr.csv'),
                    here('230214_sorl1Citalopram', 'fingerprint_sorl1.csv')),
              grporder=c('scr_10', 'sorl1_10'),
              controlGrp=c('scr_0', 'sorl1_0'),
              removeControl=TRUE,
              onlyExp='230214_14',
              colours=c('#697a87', '#EE7163'),
              legendOrNo=FALSE,
              ynameOrNo=FALSE,
              ytextOrNo=TRUE,
              xtextOrNo=TRUE,
              xParamNum=FALSE,
              nightBgOrNo=TRUE,
              ymin=-8,
              ymax=8,
              dotSize=0.1,
              lineSize=0.3,
              exportOrNo=TRUE,
              exportPath=here('230214_sorl1Citalopram', 'pubplots', 'fgp_14dose10.pdf'),
              width=90,
              height=47)


# box15/clutch1102 fingerprint --------------------------------------------

ggFingerprint(fgp=c(here('230214_sorl1Citalopram', 'fingerprint_1102scr.csv'),
                    here('230214_sorl1Citalopram', 'fingerprint_1102sorl1.csv')),
              grporder=c('scr_1', 'sorl1_1'),
              controlGrp=c('scr_0', 'sorl1_0'),
              removeControl=TRUE,
              colours=c('#b3bcc3', '#a3bbdb'),
              legendOrNo=TRUE,
              ynameOrNo=TRUE,
              ytextOrNo=TRUE,
              xtextOrNo=TRUE,
              xParamNum=FALSE,
              nightBgOrNo=TRUE,
              ymin=-8,
              ymax=8,
              exportOrNo=TRUE,
              exportPath=here('230214_sorl1Citalopram', 'plots', 'fgp_1102dose1.pdf'),
              width=200,
              height=100)

ggFingerprint(fgp=c(here('230214_sorl1Citalopram', 'fingerprint_1102scr.csv'),
                    here('230214_sorl1Citalopram', 'fingerprint_1102sorl1.csv')),
              grporder=c('scr_10', 'sorl1_10'),
              controlGrp=c('scr_0', 'sorl1_0'),
              removeControl=TRUE,
              colours=c('#697a87', '#417dcd'),
              legendOrNo=TRUE,
              ynameOrNo=TRUE,
              ytextOrNo=TRUE,
              xtextOrNo=TRUE,
              xParamNum=FALSE,
              nightBgOrNo=TRUE,
              ymin=-8,
              ymax=8,
              exportOrNo=TRUE,
              exportPath=here('230214_sorl1Citalopram', 'plots', 'fgp_1102dose2.pdf'),
              width=200,
              height=100)


# box15/clutch1031 fingerprint --------------------------------------------

ggFingerprint(fgp=c(here('230214_sorl1Citalopram', 'fingerprint_scr.csv'),
                    here('230214_sorl1Citalopram', 'fingerprint_sorl1.csv')),
              grporder=c('scr_1', 'sorl1_1'),
              controlGrp=c('scr_0', 'sorl1_0'),
              removeControl=TRUE,
              onlyExp='230214_15',
              colours=c('#b3bcc3', '#F6B8B1'),
              legendOrNo=FALSE,
              ynameOrNo=FALSE,
              ytextOrNo=TRUE,
              xtextOrNo=FALSE,
              xParamNum=TRUE,
              nightBgOrNo=TRUE,
              ymin=-8,
              ymax=8,
              dotSize=0.1,
              lineSize=0.3,
              exportOrNo=TRUE,
              exportPath=here('230214_sorl1Citalopram', 'pubplots', 'fgp_1031dose1.pdf'),
              width=90,
              height=28)

# previously blue shades:
# colours=c('#b3bcc3', '#a3bbdb'),

ggFingerprint(fgp=c(here('230214_sorl1Citalopram', 'fingerprint_scr.csv'),
                    here('230214_sorl1Citalopram', 'fingerprint_sorl1.csv')),
              grporder=c('scr_10', 'sorl1_10'),
              controlGrp=c('scr_0', 'sorl1_0'),
              removeControl=TRUE,
              onlyExp='230214_15',
              colours=c('#697a87', '#EE7163'),
              legendOrNo=FALSE,
              ynameOrNo=FALSE,
              ytextOrNo=TRUE,
              xtextOrNo=TRUE,
              xParamNum=FALSE,
              nightBgOrNo=TRUE,
              ymin=-8,
              ymax=8,
              dotSize=0.1,
              lineSize=0.3,
              exportOrNo=TRUE,
              exportPath=here('230214_sorl1Citalopram', 'pubplots', 'fgp_1031dose10.pdf'),
              width=90,
              height=47)

# previously blue shades:
# colours=c('#697a87', '#417dcd'),


# mean fingerprints -------------------------------------------------------

calculateFingerprint(paDir=here('230214_sorl1Citalopram', 'bhvparams'),
                     controlGrp=c('scr_0'),
                     grporder=c('scr_0', 'scr_1', 'scr_10'),
                     mergeExp1=c('230214_14', '230214_15'),
                     skipNight0=TRUE)
# changed manually filename fingerprint.csv to fingerprintMean_scr.csv

calculateFingerprint(paDir=here('230214_sorl1Citalopram', 'bhvparams'),
                     controlGrp=c('sorl1_0'),
                     grporder=c('sorl1_0', 'sorl1_1', 'sorl1_10'),
                     mergeExp1=c('230214_14', '230214_15'),
                     skipNight0=TRUE)
# changed manually filename fingerprint.csv to fingerprintMean_sorl1.csv

ggFingerprint(fgp=c(here('230214_sorl1Citalopram', 'fingerprintMean_scr.csv'),
                    here('230214_sorl1Citalopram', 'fingerprintMean_sorl1.csv')),
              grporder=c('scr_1', 'sorl1_1'),
              controlGrp=c('scr_0', 'sorl1_0'),
              removeControl=TRUE,
              colours=c('#b3bcc3', '#F6B8B1'),
              legendOrNo=FALSE,
              ynameOrNo=FALSE,
              ytextOrNo=TRUE,
              xtextOrNo=FALSE,
              xParamNum=TRUE,
              nightBgOrNo=TRUE,
              ymin=-8,
              ymax=8,
              dotSize=0.1,
              lineSize=0.3,
              exportOrNo=TRUE,
              exportPath=here('230214_sorl1Citalopram', 'pubplots', 'fgpMean_dose1.pdf'),
              width=90,
              height=28)

ggFingerprint(fgp=c(here('230214_sorl1Citalopram', 'fingerprintMean_scr.csv'),
                    here('230214_sorl1Citalopram', 'fingerprintMean_sorl1.csv')),
              grporder=c('scr_10', 'sorl1_10'),
              controlGrp=c('scr_0', 'sorl1_0'),
              removeControl=TRUE,
              colours=c('#697a87', '#EE7163'),
              legendOrNo=FALSE,
              ynameOrNo=FALSE,
              ytextOrNo=TRUE,
              xtextOrNo=TRUE,
              xParamNum=FALSE,
              nightBgOrNo=TRUE,
              ymin=-8,
              ymax=8,
              dotSize=0.1,
              lineSize=0.3,
              exportOrNo=TRUE,
              exportPath=here('230214_sorl1Citalopram', 'pubplots', 'fgpMean_dose10.pdf'),
              width=90,
              height=47)



# effect sizes ------------------------------------------------------------

### 1 µM ###
# day, scr
percEffectSize(pa=c(here('230214_sorl1Citalopram', 'bhvparams', 'sleepHours_230214_14.csv'),
                    here('230214_sorl1Citalopram', 'bhvparams', 'sleepHours_230214_15.csv')),
               lmePath=here('230214_sorl1Citalopram/LMEscr0.csv'), # i.e. baseline is SCR + DMSO
               grporder=c('scr_1', 'scr_0'),
               dayornight='day')

# day, sorl1
percEffectSize(pa=c(here('230214_sorl1Citalopram', 'bhvparams', 'sleepHours_230214_14.csv'),
                    here('230214_sorl1Citalopram', 'bhvparams', 'sleepHours_230214_15.csv')),
               lmePath=here('230214_sorl1Citalopram/LMEsorl10.csv'), # i.e. baseline is SCR + DMSO
               grporder=c('sorl1_1', 'scr_0'),
               dayornight='day')

# night, scr
percEffectSize(pa=c(here('230214_sorl1Citalopram', 'bhvparams', 'sleepHours_230214_14.csv'),
                    here('230214_sorl1Citalopram', 'bhvparams', 'sleepHours_230214_15.csv')),
               lmePath=here('230214_sorl1Citalopram/LMEscr0.csv'), # i.e. baseline is SCR + DMSO
               grporder=c('scr_1', 'scr_0'),
               dayornight='night')

# night, sorl1
percEffectSize(pa=c(here('230214_sorl1Citalopram', 'bhvparams', 'sleepHours_230214_14.csv'),
                    here('230214_sorl1Citalopram', 'bhvparams', 'sleepHours_230214_15.csv')),
               lmePath=here('230214_sorl1Citalopram/LMEsorl10.csv'), # i.e. baseline is SCR + DMSO
               grporder=c('sorl1_1', 'scr_0'),
               dayornight='night')
# sleepHours	day	scr_0	scr_1	1.145928332	0.565301424 >> 1.1
# sleepHours	night	scr_0	scr_1	0.572542242	0.427001755 >> 0.6

# sleepHours	day	sorl1_0	sorl1_1	1.241892791	0.724146026 >> 1.2
# sleepHours	night	sorl1_0	sorl1_1	1.589954745	0.546813967 >> 1.6

# day both genotypes avg: 1.19 hr
# night both genotypes avg: 1.08 hr


### 10 µM ###
# sleepHours	day	scr_0	scr_10	3.396733325	0.565301424
# sleepHours	night	scr_0	scr_10	2.061143341	0.427001755

# sleepHours	day	sorl1_0	sorl1_10	5.122826975	0.723931541
# sleepHours	night	sorl1_0	sorl1_10	3.861330079	0.546813967


### differential effect ###
# e.g. sleep during the day 
# increase SCR + 10
percEffectSize(pa=c(here('230214_sorl1Citalopram', 'bhvparams', 'sleepHours_230214_14.csv'),
                    here('230214_sorl1Citalopram', 'bhvparams', 'sleepHours_230214_15.csv')),
               lmePath=here('230214_sorl1Citalopram/LMEscr0.csv'), # i.e. baseline is SCR + DMSO
               grporder=c('scr_10', 'scr_0'),
               dayornight='day')
# 116.295 so 2.16x

# increase SORL1 + 10
percEffectSize(pa=c(here('230214_sorl1Citalopram', 'bhvparams', 'sleepHours_230214_14.csv'),
                    here('230214_sorl1Citalopram', 'bhvparams', 'sleepHours_230214_15.csv')),
               lmePath=here('230214_sorl1Citalopram/LMEsorl10.csv'), # i.e. baseline is SORL1 + DMSO
               grporder=c('sorl1_10', 'sorl1_0'),
               dayornight='day')
# 153.8799 so 2.54x



# eLife reviews: check replication of baseline phenotype ------------------

calculateFingerprint(paDir=here('230214_sorl1Citalopram', 'bhvparams'),
                     controlGrp='scr_0',
                     grporder=c('scr_0', 'sorl1_0'),
                     skipNight0=TRUE)
# changed manually filename fingerprint.csv to fingerprint_H2O.csv, i.e. only comparing H2O-treated sorl1 vs scr

ggFingerprint(fgp=here('230214_sorl1Citalopram', 'fingerprint_H2O.csv'),
              grporder='sorl1_0',
              controlGrp='scr_0',
              removeControl=TRUE,
              #onlyExp='230214_14',
              #colours=c('#b3bcc3', '#F6B8B1'),
              legendOrNo=TRUE,
              ynameOrNo=FALSE,
              ytextOrNo=TRUE,
              xtextOrNo=TRUE,
              xParamNum=FALSE,
              nightBgOrNo=TRUE,
              ymin=-2,
              ymax=2,
              dotSize=0.1,
              lineSize=0.3,
              exportOrNo=TRUE,
              exportPath=here('230214_sorl1Citalopram', 'plots', 'fgp_H2O.pdf'),
              width=90,
              height=28)



# eLife reviews: sorl1 vs scr + citalopram --------------------------------

fgp <- calculateFingerprint(paDir=here('230214_sorl1Citalopram', 'bhvparams'),
                            mergeExp1=c('230214_14', '230214_15'),
                            controlGrp='scr_0',
                            grporder=c('scr_0', 'scr_1', 'scr_10', 'sorl1_0', 'sorl1_1', 'sorl1_10'),
                            skipNight0=TRUE)
# changed manually filename fingerprint.csv to fingerprint_SCR0.csv, i.e. only comparing H2O-treated sorl1 vs scr

ggFingerprint(fgp=here('230214_sorl1Citalopram', 'fingerprint_SCR0.csv'),
              grporder=c('scr_0', 'scr_1', 'sorl1_0'),
              controlGrp='scr_0',
              removeControl=TRUE,
              #onlyExp='230214_14',
              #colours=c('#b3bcc3', '#F6B8B1'),
              legendOrNo=TRUE,
              ynameOrNo=FALSE,
              ytextOrNo=TRUE,
              xtextOrNo=TRUE,
              xParamNum=FALSE,
              nightBgOrNo=TRUE,
              ymin=-2,
              ymax=2,
              dotSize=0.1,
              lineSize=0.3,
              exportOrNo=TRUE,
              exportPath=here('230214_sorl1Citalopram', 'plots', 'fgp_SCR0.pdf'),
              width=90,
              height=28)

ggPairwiseHeat(fgp=here('230214_sorl1Citalopram', 'fingerprint_SCR0.csv'),
               metric='mean',
               simScore='cosine',
               grporder=c('scr_0', 'scr_1', 'scr_10', 'sorl1_0'),
               removeControl=TRUE,
               controlGrp='scr_0',
               medianMid=NA,
               minCol=NA,
               maxCol=NA,
               onlyHalf='upper',
               scoreSize=5,
               legendOrNo=TRUE,
               labelsOrNo=TRUE,
               width=36,
               height=36, 
               exportPath=here('230214_sorl1Citalopram', 'plots', 'heatSCR0.pdf'))


### can also bring baseline SORL1 tracking
fgp <- calculateFingerprint(paDir=c(here('220531_SORL1/bhvparams/'),
                                    here('230214_sorl1Citalopram', 'bhvparams')),
                            mergeExp1=c('220531_14', '220531_15'),
                            mergeExp2=c('230214_14', '230214_15'),
                            controlGrp=c('scr', 'scr_0'),
                            grporder=c('scr', 'sorl1', 'scr_0', 'scr_1', 'scr_10', 'sorl1_0', 'sorl1_1', 'sorl1_10'),
                            skipNight0=TRUE)

ggFingerprint(fgp=fgp,
              grporder=c('scr', 'sorl1', 'scr_1', 'scr_10'),
              controlGrp='scr',
              removeControl=TRUE,
              #onlyExp='230214_14',
              #colours=c('#b3bcc3', '#F6B8B1'),
              legendOrNo=TRUE,
              ynameOrNo=FALSE,
              ytextOrNo=TRUE,
              xtextOrNo=TRUE,
              xParamNum=FALSE,
              nightBgOrNo=TRUE,
              ymin=-2,
              ymax=2,
              dotSize=0.1,
              lineSize=0.3,
              exportOrNo=FALSE,
              exportPath=NA,
              width=90,
              height=28)

ggPairwiseHeat(fgp=fgp,
               metric='mean',
               simScore='cosine',
               grporder=c('sorl1', 'scr_1', 'scr_10', 'sorl1_0'),
               removeControl=TRUE,
               controlGrp=c('scr', 'scr_0'),
               medianMid=NA,
               minCol=NA,
               maxCol=NA,
               onlyHalf='upper',
               scoreSize=5,
               legendOrNo=TRUE,
               labelsOrNo=TRUE,
               width=36,
               height=36, 
               exportPath=here('230214_sorl1Citalopram', 'plots', 'heatBASELINE.pdf'))



### more polished fingerprint baseline experiment vs. citalopram experiment
# will merge citalopram experiment clutches because low Ns
# but will keep the two baseline clutches separated
fgp <- calculateFingerprint(paDir=c(here('220531_SORL1/bhvparams/'),
                                    here('230214_sorl1Citalopram', 'bhvparams')),
                            mergeExp2=c('230214_14', '230214_15'),
                            controlGrp=c('scr', 'scr_0'),
                            grporder=c('scr', 'sorl1', 'scr_0', 'scr_1', 'scr_10', 'sorl1_0', 'sorl1_1', 'sorl1_10'),
                            skipNight0=TRUE)

# need to do it in two plots otherwise difficult to see
### PLOT1: compare all sorl1 KO experiments
ggFingerprint(fgp=fgp,
              grporder=c('scr', 'sorl1', 'sorl1_0'),
              controlGrp=c('scr', 'scr_0'),
              removeControl=TRUE,
              colours=c('#F18D82', '#E94735', '#5a6974'), # light = box14 = clutch2 // dark = box15 = clutch1
              legendOrNo=FALSE,
              ynameOrNo=FALSE,
              ytextOrNo=TRUE,
              xtextOrNo=TRUE,
              xParamNum=FALSE,
              nightBgOrNo=TRUE,
              ymin=-1.5,
              ymax=2.2,
              dotSize=0.1,
              lineSize=0.3,
              exportOrNo=TRUE,
              exportPath=here('230214_sorl1Citalopram/pubplots/fgp_SCR0.pdf'),
              width=100,
              height=47)

ggPairwiseHeat(fgp=fgp,
               metric='mean',
               simScore='cosine',
               grporder=c('scr', 'sorl1', 'sorl1_0'),
               removeControl=TRUE,
               controlGrp=c('scr', 'scr_0'),
               medianMid=NA,
               minCol=NA,
               maxCol=NA,
               onlyHalf='upper',
               scoreSize=5,
               legendOrNo=TRUE,
               labelsOrNo=TRUE,
               width=36,
               height=36, 
               exportPath=here('230214_sorl1Citalopram', 'plots', 'heatSCR0.pdf'))

### PLOT2: compare sorl1 KO vs citalopram treatment
# will only keep it within citalopram experiment
ggFingerprint(fgp=fgp,
              grporder=c('scr_0', 'sorl1_0', 'scr_1', 'scr_10'),
              controlGrp='scr_0',
              removeControl=TRUE,
              colours=c('#5a6974', '#96d3a1', '#2aa844'),
              legendOrNo=FALSE,
              ynameOrNo=FALSE,
              ytextOrNo=TRUE,
              xtextOrNo=TRUE,
              xParamNum=FALSE,
              nightBgOrNo=TRUE,
              ymin=-2,
              ymax=4.2,
              dotSize=0.1,
              lineSize=0.3,
              exportOrNo=TRUE,
              exportPath=here('230214_sorl1Citalopram/pubplots/fgp_Sorl1vsCita.pdf'),
              width=100,
              height=47)

ggPairwiseHeat(fgp=fgp,
               metric='mean',
               simScore='cosine',
               grporder=c('scr_0', 'sorl1_0', 'scr_1', 'scr_10'),
               removeControl=TRUE,
               controlGrp='scr_0',
               medianMid=NA,
               minCol=NA,
               maxCol=NA,
               onlyHalf='upper',
               scoreSize=5,
               legendOrNo=TRUE,
               labelsOrNo=TRUE,
               width=36,
               height=36, 
               exportPath=here('230214_sorl1Citalopram', 'plots', 'heat_Sorl1vsCita.pdf'))

### PLOT3: effect of citalopram on SORL1, compared to SCR + H2O
ggFingerprint(fgp=fgp,
              grporder=c('scr_0', 'sorl1_0', 'sorl1_1', 'sorl1_10'),
              controlGrp='scr_0',
              removeControl=TRUE,
              colours=c('#5a6974', '#fcd98a', '#fcb505'),
              legendOrNo=FALSE,
              ynameOrNo=FALSE,
              ytextOrNo=TRUE,
              xtextOrNo=TRUE,
              xParamNum=FALSE,
              nightBgOrNo=TRUE,
              ymin=-2.5,
              ymax=7,
              dotSize=0.1,
              lineSize=0.3,
              exportOrNo=TRUE,
              exportPath=here('230214_sorl1Citalopram/pubplots/fgp_CitalopramonSORL1.pdf'),
              width=102.5,
              height=47)

ggPairwiseHeat(fgp=fgp,
               metric='mean',
               simScore='cosine',
               grporder=c('scr_0', 'sorl1_0', 'scr_1', 'scr_10'),
               removeControl=TRUE,
               controlGrp='scr_0',
               medianMid=NA,
               minCol=NA,
               maxCol=NA,
               onlyHalf='upper',
               scoreSize=5,
               legendOrNo=TRUE,
               labelsOrNo=TRUE,
               width=36,
               height=36, 
               exportPath=here('230214_sorl1Citalopram', 'plots', 'heat_Sorl1vsCita.pdf'))
