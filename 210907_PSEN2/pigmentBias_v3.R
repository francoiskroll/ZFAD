#####################################################
# ~ ZFAD: adjust delta-pixel data to control for fainter pigmentation of PSEN2 larvae ~
#
#
# Francois Kroll 2023
# francois@kroll.be
#####################################################

# we want to know if the Zebrabox camera perceives the psen2 as significantly fainter
# idea here: look at the maximum number of pixels triggered during the startle response when lights off
# this should be close to the total number of pixels of each larva

# then scale the data accordingly

# v3:
# now that FramebyFrame package exists, can simply use activitySunsetStartle parameter
# and we are making adjustPixel a function in FramebyFrame

# note, v2 version was taking 4 seconds after sunset
# current definition in FramebyFrame package is 10 seconds


# packages ----------------------------------------------------------------

# library(devtools)
# install_github('francoiskroll/FramebyFrame')

library(FramebyFrame)

library(here)
library(dplyr)
library(tidyr)
library(tibble)

library(ggplot2)



# is there a difference in activitySunsetStartle? -------------------------
# if there is none, there might be no support for adjusting the pixel data

# Note, issued version v0.9.0 with better definition of parameter activitySunsetStartle
# motivated because I had trouble replicating pigmentBias_v2 when looking only *after* the transition
# clearly the precise frame of the light transition is not perfectly accurate

genotypeGenerator(plateMap='~/Dropbox/ZFAD/210907_PSEN2/210907_12_genotypeMap.xlsx')
genotypeGenerator(plateMap='~/Dropbox/ZFAD/210907_PSEN2/210907_13_genotypeMap.xlsx')

behaviourParameter(parameter='activitySunsetStartle',
                   ffpath=c('~/Dropbox/ZFAD/210907_PSEN2/210907_12_RAWs.csv',
                            '~/Dropbox/ZFAD/210907_PSEN2/210907_13_RAWs.csv'),
                   genopath=c('~/Dropbox/ZFAD/210907_PSEN2/210907_12genotype.txt',
                              '~/Dropbox/ZFAD/210907_PSEN2/210907_13genotype.txt'),
                   skipNight0=FALSE)

# change names to
# activitySunsetStartle_210907_12BEFORE.csv
# activitySunsetStartle_210907_13BEFORE.csv

LMEdaynight(pa=c('~/Dropbox/ZFAD/210907_PSEN2/bhvparams/activitySunsetStartle_210907_12.csv',
                 '~/Dropbox/ZFAD/210907_PSEN2/bhvparams/activitySunsetStartle_210907_13.csv'),
            grporder=c('scr', 'psen2'),
            skipNight0=TRUE,
            silent=FALSE,
            detailsOrNo=FALSE)

ggParameter(pa=c(here('210907_PSEN2', 'bhvparams', 'activitySunsetStartle_210907_13.csv'),
                 here('210907_PSEN2', 'bhvparams', 'activitySunsetStartle_210907_12.csv')),
            grporder=c('psen2', 'scr'),
            onlyDayorNight='night',
            colours=c('#78ac63', '#697a87'),
            fainterExp=TRUE,
            faintMax=0.35,
            ymin=NA,
            ymax=NA,
            dotSize=0.3,
            connectMean=TRUE,
            legendOrNo=FALSE,
            xtextOrNo=FALSE,
            ynameOrNo=FALSE,
            yunitOrNo=TRUE,
            splitBy='experiment',
            titleOrNo=FALSE,
            statsOrNo=FALSE,
            width=70,
            height=50,
            exportPath=here('210907_PSEN2', 'plots', 'sunStartle_night.pdf'))
### yes, as expected, sunsetStartle significantly lower in PSEN2 knockouts


# plotting a small timecourse for illustration ----------------------------

### use function ggZoomFrames originally written for bhvparamFigs

# wdf = typically frame-by-frame dataframe of one window
# expecting column 'exsecs'
# and expecting some columns to be f1, f2, etc.

ggZoomFrames <- function(wdf,
                         fid,
                         from,
                         to,
                         colour='#fcb505',
                         ymax=50,
                         labelActive=FALSE,
                         exportPath,
                         width=80,
                         height=50) {
  
  wdf <- as.data.frame(wdf)
  
  frames <- wdf[from:to , c( which(colnames(wdf)=='exsecs') , which(colnames(wdf)==fid))]
  colnames(frames) <- c('exsecs', 'px')
  
  # add a column to label the active frames
  frames$active <- frames$px>0
  
  # have breaks at multiple of fps
  fps <- round(averageFramerate(frames$exsecs))
  brks <- frames$exsecs[seq(1, nrow(frames), fps)]
  
  if(labelActive) {
    ggz <- ggplot(frames, aes(x=exsecs, y=px)) +
      geom_line(colour=colour, linewidth=0.8) +
      geom_point(aes(fill=active), shape=21, colour=colour, size=1, stroke=0.4) +
      scale_fill_manual(values=c('white', colour))
  } else {
    ggz <- ggplot(frames, aes(x=exsecs, y=px)) +
      geom_line(colour=colour, linewidth=0.8) +
      geom_point(shape=21, fill='white', colour=colour, size=1, stroke=0.4)
  }
  
  # plotting function
  ggz <- ggz +
    #geom_text(aes(label=px)) +
    theme_minimal() +
    theme(
      panel.grid.minor.x=element_blank(),
      panel.grid.minor.y=element_blank(),
      axis.title.x=element_blank(),
      axis.title.y=element_text(size=9),
      axis.text.x=element_text(size=7),
      axis.text.y=element_text(size=7),
      legend.position='none'
    ) +
    coord_cartesian(ylim=c(0, ymax)) +
    scale_x_continuous(breaks=brks, labels=-1:(length(brks)-2)) +
    xlab('seconds') +
    ylab('\u0394 px')
  
  # export
  ggsave(exportPath, ggz, width=width, height=height, units='mm', device=cairo_pdf)
  
  # return
  return(ggz)
  
}


### import some data
ff <- data.table::fread(here('210907_PSEN2', '210907_12_RAWs.csv'))

dn <- splitFramesbyDayNight(tc=ff,
                            path=here('210907_PSEN2', '210907_12_13_PSEN2.xls'))

fps <- averageFramerate(dn$night2$exsecs)


### keep 4 seconds around sunset

# will only keep sunset2, for example

# box1
# sunset 2
# getting light transition frames from 210907_12_lights.csv


### plot example startle responses
ggZoomFrames(wdf=ff,
             fid='f61', # psen2 larva
             from=(5058178-round(fps)),
             to=(5058178+3*round(fps)),
             colour='#78ac63',
             ymax=90,
             labelActive=FALSE,
             exportPath=here('210907_PSEN2', 'pubplots', 'f61psen2_startlev2.pdf'),
             width=80,
             height=45)

ggZoomFrames(wdf=ff,
             fid='f18',
             from=(5058178-round(fps)),
             to=(5058178+3*round(fps)),
             colour='#697a87',
             ymax=90,
             labelActive=FALSE,
             exportPath=here('210907_PSEN2', 'pubplots', 'f18scr_startlev2.pdf'),
             width=80,
             height=45)



# take maximum sunsetStartle across nights --------------------------------
# assumption is that it should be closer to the real size/pigmentation of the larva in pixels
# (or at least more proportional to it)

sut1 <- read.csv('~/Dropbox/ZFAD/210907_PSEN2/pigmentBias/bhvparams/activitySunsetStartle_210907_12BEFORE.csv')
sut2 <- read.csv('~/Dropbox/ZFAD/210907_PSEN2/pigmentBias/bhvparams/activitySunsetStartle_210907_13BEFORE.csv')

sut1 <- sut1 %>%
  filter(grp!='excluded')

sut2 <- sut2 %>%
  filter(grp!='excluded')

## take the largest sunset startle of the 3 nights
sut1 <- sut1 %>%
  mutate(maxSt=apply(sut1, 1, function(row) {
    max(as.numeric(row[c('night0', 'night1', 'night2')]))
  }))

sut2 <- sut2 %>%
  mutate(maxSt=apply(sut2, 1, function(row) {
    max(as.numeric(row[c('night0', 'night1', 'night2')]))
  }))


### compare with t-tests to add on plot
t.test(maxSt ~ grp, data=sut1)
t.test(maxSt ~ grp, data=sut2)

# summary statistics mean ± sd
sut1 %>%
  group_by(grp) %>%
  summarise_at(vars(maxSt),
               list(mean=~mean(., na.rm=TRUE),
                    sd=~sd(., na.rm=TRUE),
                    n=~length(.)))

sut2 %>%
  group_by(grp) %>%
  summarise_at(vars(maxSt),
               list(mean=~mean(., na.rm=TRUE),
                    sd=~sd(., na.rm=TRUE),
                    n=~length(.)))

               


### ratio of means?
mean(sut1[which(sut1$grp=='psen2'), 'maxSt']) / mean(sut1[which(sut1$grp=='scr'), 'maxSt'])
mean(sut2[which(sut2$grp=='psen2'), 'maxSt']) / mean(sut2[which(sut2$grp=='scr'), 'maxSt'])


### 25/04/2023, created function scaleFromStartle
# scaleFromStartle also does the t-test
scaleFromStartle(ffpath='~/Dropbox/ZFAD/210907_PSEN2/210907_12_RAWs.csv',
                 genopath='~/Dropbox/ZFAD/210907_PSEN2/210907_12genotype.txt',
                 grpL='scr',
                 grpS='psen2',
                 dayduration=14)
# 0.95

scaleFromStartle(ffpath='~/Dropbox/ZFAD/210907_PSEN2/210907_13_RAWs.csv',
                 genopath='~/Dropbox/ZFAD/210907_PSEN2/210907_13genotype.txt',
                 grpL='scr',
                 grpS='psen2',
                 dayduration=14)
# 0.89



# plot sunsetStartle ------------------------------------------------------

# we can make use of ggParameter to plot maximum sunsetStartle across all nights
# just 'cheat' by calling the max sunsetStartle column nightX
colnames(sut1)[which(colnames(sut1)=='maxSt')] <- 'night9'
colnames(sut2)[which(colnames(sut2)=='maxSt')] <- 'night9'

ggParameter(pa=c('sut2', 'sut1'), # so plot clutch1, clutch2 (we swapped numbers in paper)
            grporder=c('psen2', 'scr'),
            onlyWin='night9',
            colours=c('#78ac63', '#697a87'),
            fainterExp=TRUE,
            faintMax=0.35,
            ymin=35,
            ymax=110,
            dotSize=1,
            violinWidth=0.4,
            connectMean=TRUE,
            legendOrNo=FALSE,
            xtextOrNo=FALSE,
            ynameOrNo=FALSE,
            yunitOrNo=FALSE,
            splitBy='experiment',
            titleOrNo=FALSE,
            statsOrNo=FALSE,
            width=33,
            height=55,
            exportPath=here('210907_PSEN2', 'pubplots', 'actPxMAX.pdf'))

# change back column to maxSt
colnames(sut1)[which(colnames(sut1)=='night9')] <- 'maxSt'
colnames(sut2)[which(colnames(sut2)=='night9')] <- 'maxSt'


# check approach makes sense ----------------------------------------------

# same larva should have close max deltapx between the two sunsets
# (at least closer than if you shuffle)

# will compare between sunset1 and sunset2 here
# but could also do between sunset0 and sunset1
# or between sunset0 and sunset2

# shuffle data
# to-do so, will simply shift all sunset2 one row down
# and put last value at the top

## try on box1

# first copy sunset2 column
sut1$night2shift <- sut1$night2
# then shift it all down
# take last value and put it as first & delete last value
sut1$night2shift <- sut1$night2shift[c( nrow(sut1) , 1:(nrow(sut1)-1))]

# correlation between sunsets
cor(sut1$night1, sut1$night2)
cor.test(sut1$night1, sut1$night2)

# correlation between sunsets after shuffling
cor(sut1$night1, sut1$night2shift)
cor.test(sut1$night1, sut1$night2shift)

# yes, as expected
# correlation (low but significant) before shuffling, drops after shuffling

# slope plots to illustrate

# pivot data
stl1 <- sut1 %>%
  pivot_longer(cols=-c(parameter, date, box, fish, grp),
               names_to='win',
               values_to='maxSt')

# line plot for each fish
# i.e. for each fish, deltapx sunset1 --- deltapx sunset2

# before shuffling
sloplot <- stl1 %>%
  subset(win %in% c('night1', 'night2')) %>%
  subset(grp != 'excluded') %>%
  ggplot(., aes(x=win, y=maxSt)) +
  geom_line(aes(group=fish), size=0.6) +
  coord_cartesian(ylim=c(0, 100))
sloplot
ggsave('~/Dropbox/ZFAD/210907_PSEN2/pigmentBias/slopeStartle.pdf', sloplot, width=100, height=200, units='mm')

# after shuffling
sloplot <- stl1 %>%
  subset(win %in% c('night1', 'night2shift')) %>%
  subset(grp != 'excluded') %>%
  ggplot(., aes(x=win, y=maxSt)) +
  geom_line(aes(group=fish), size=0.6) +
  coord_cartesian(ylim=c(0, 100))
sloplot
ggsave('~/Dropbox/ZFAD/210907_PSEN2/pigmentBias/slopeStartle_shuffl.pdf', sloplot, width=100, height=200, units='mm')


#### skipping a chunk from pigmentBias_v2 here
# mainly about comparing clutches: yes, clutches are different
# scale them separately with different scaling ratios


# adjust frame-by-frame data ----------------------------------------------

# function adjustPixel in FramebyFrame package
# please refer to adjustPixel.R for comments on the function
# usage adjustPixel(rawpath, grp, scale, round='down', export=TRUE)

# scaling ratios are different in means for each box
scar1 <- mean(sut1[which(sut1$grp=='psen2'), 'maxSt']) / mean(sut1[which(sut1$grp=='scr'), 'maxSt'])
scar2 <- mean(sut2[which(sut2$grp=='psen2'), 'maxSt']) / mean(sut2[which(sut2$grp=='scr'), 'maxSt'])

# box1
adjustPixel(ffpath='~/Dropbox/ZFAD/210907_PSEN2/210907_12_RAWs.csv',
            genopath='~/Dropbox/ZFAD/210907_PSEN2/210907_12genotype.txt',
            grp='scr',
            scale=scar1,
            round='down',
            export=TRUE)

# box2
adjustPixel(ffpath='~/Dropbox/ZFAD/210907_PSEN2/210907_13_RAWs.csv',
            genopath='~/Dropbox/ZFAD/210907_PSEN2/210907_13genotype.txt',
            grp='scr',
            scale=scar2,
            round='down')
# Note; saves adjusted RAWs to drive and does not return anything



# plot small timecourse to check ------------------------------------------

# can do same as the ones above as a check

# so import back adjusted RAWs
fa1 <- data.table::fread('~/Dropbox/ZFAD/210907_PSEN2/210907_12_RAWsadjusted.csv') # frame-by-frame adjusted for box1
fa2 <- data.table::fread('~/Dropbox/ZFAD/210907_PSEN2/210907_13_RAWsadjusted.csv') # frame-by-frame adjusted for box2


### plot example startle responses
ggZoomFrames(wdf=fa1,
             fid='f61', # psen2 larva
             from=(5058178-round(fps)),
             to=(5058178+3*round(fps)),
             colour='#78ac63',
             ymax=90,
             labelActive=FALSE,
             exportPath='~/Dropbox/ZFAD/210907_PSEN2/pigmentBias/f61psen2_startle_afterAdj.pdf',
             width=80,
             height=45)
# same as before, correct

ggZoomFrames(wdf=fa1,
             fid='f18',
             from=(5058178-round(fps)),
             to=(5058178+3*round(fps)),
             colour='#697a87',
             ymax=88,
             labelActive=FALSE,
             exportPath='~/Dropbox/ZFAD/210907_PSEN2/pigmentBias/f18scr_startle_afterAdj.pdf',
             width=80,
             height=45)
# yes, a tiny bit smaller


# -------------------------------------------------------------------------

# same analysis on adjusted data
# there should be no difference now

# -------------------------------------------------------------------------

## calculate activitySunsetStartle again

behaviourParameter(parameter='activitySunsetStartle',
                   ffpath=c('~/Dropbox/ZFAD/210907_PSEN2/210907_12_RAWsadjusted.csv',
                            '~/Dropbox/ZFAD/210907_PSEN2/210907_13_RAWsadjusted.csv'),
                   genopath=c('~/Dropbox/ZFAD/210907_PSEN2/210907_12genotype.txt',
                              '~/Dropbox/ZFAD/210907_PSEN2/210907_13genotype.txt'),
                   skipNight0=FALSE)
# change names to
# activitySunsetStartle_210907_12.csv > activitySunsetStartle_210907_12AFTER.csv
# activitySunsetStartle_210907_13.csv > activitySunsetStartle_210907_13AFTER.csv

LMEdaynight(pa=c('~/Dropbox/ZFAD/210907_PSEN2/pigmentBias/bhvparams/activitySunsetStartle_210907_12AFTER.csv',
                 '~/Dropbox/ZFAD/210907_PSEN2/pigmentBias/bhvparams/activitySunsetStartle_210907_13AFTER.csv'),
            grporder=c('scr', 'psen2'),
            skipNight0=TRUE,
            silent=FALSE,
            detailsOrNo=FALSE)
# now ns.

ggParameter(pa=c('~/Dropbox/ZFAD/210907_PSEN2/pigmentBias/bhvparams/activitySunsetStartle_210907_12AFTER.csv',
                 '~/Dropbox/ZFAD/210907_PSEN2/pigmentBias/bhvparams/activitySunsetStartle_210907_13AFTER.csv'),
            grporder=c('psen2', 'scr'),
            onlyDayorNight='night',
            colours=c('#78ac63', '#697a87'),
            fainterExp=TRUE,
            faintMax=0.35,
            ymin=NA,
            ymax=NA,
            dotSize=0.3,
            connectMean=TRUE,
            legendOrNo=FALSE,
            xtextOrNo=FALSE,
            ynameOrNo=FALSE,
            yunitOrNo=TRUE,
            splitBy='experiment',
            titleOrNo=FALSE,
            statsOrNo=FALSE,
            width=70,
            height=50,
            exportPath='~/Dropbox/ZFAD/210907_PSEN2/plots/actSun_adjust.pdf')
### as expected, flattened



# take maximum sunsetStartle across nights --------------------------------
# assumption is that it should be closer to the real size/pigmentation of the larva in pixels
# (or at least more proportional to it)

sut1 <- read.csv('~/Dropbox/ZFAD/210907_PSEN2/bhvparams/activitySunsetStartle_210907_12.csv')
sut2 <- read.csv('~/Dropbox/ZFAD/210907_PSEN2/bhvparams/activitySunsetStartle_210907_13.csv')

sut1 <- sut1 %>%
  filter(grp!='excluded')

sut2 <- sut2 %>%
  filter(grp!='excluded')

## take the largest sunset startle of the 3 nights
sut1 <- sut1 %>%
  mutate(maxSt=apply(sut1, 1, function(row) {
    max(as.numeric(row[c('night0', 'night1', 'night2')]))
  }))

sut2 <- sut2 %>%
  mutate(maxSt=apply(sut2, 1, function(row) {
    max(as.numeric(row[c('night0', 'night1', 'night2')]))
  }))


### compare with t-tests to add on plot
t.test(maxSt ~ grp, data=sut1)
t.test(maxSt ~ grp, data=sut2)


### ratio of means?
mean(sut1[which(sut1$grp=='psen2'), 'maxSt']) / mean(sut1[which(sut1$grp=='scr'), 'maxSt'])
mean(sut2[which(sut2$grp=='psen2'), 'maxSt']) / mean(sut2[which(sut2$grp=='scr'), 'maxSt'])



# plot sunsetStartle ------------------------------------------------------

# we can make use of ggParameter to plot maximum sunsetStartle across all nights
# just 'cheat' by calling the max sunsetStartle column nightX
colnames(sut1)[which(colnames(sut1)=='maxSt')] <- 'night9'
colnames(sut2)[which(colnames(sut2)=='maxSt')] <- 'night9'

ggParameter(pa=c('sut2', 'sut1'), # clutch1, clutch2 (nomenclature in paper)
            grporder=c('psen2', 'scr'),
            onlyWin='night9',
            colours=c('#78ac63', '#697a87'),
            fainterExp=TRUE,
            faintMax=0.35,
            ymin=35,
            ymax=110,
            dotSize=1,
            violinWidth=0.4,
            connectMean=TRUE,
            legendOrNo=FALSE,
            xtextOrNo=FALSE,
            ynameOrNo=FALSE,
            yunitOrNo=FALSE,
            splitBy='experiment',
            titleOrNo=FALSE,
            statsOrNo=FALSE,
            width=33,
            height=55,
            exportPath="~/Dropbox/ZFAD/210907_PSEN2/pubplots/actPxMAX_adj.pdf")

# change back column to maxSt
colnames(sut1)[which(colnames(sut1)=='night9')] <- 'maxSt'
colnames(sut2)[which(colnames(sut2)=='night9')] <- 'maxSt'


### Notes from pigmentBias_v2

# it seems to overshoot by 0.5 px (i.e. psen2 is now 0.5 px darker)
# I suppose it is somewhat mathematically expected from the rounding down
# I think will keep this way
# the issue with using normal rounding (i.e. 2.4 becomes 2) is that it prevents the creation of new inactive frames
# indeed, the minimum active frame is 1, so when it is multiplied by 0.9, it becomes 0.9 px
# if normal rounding: becomes 1 again
# if rounding down: becomes 0, we created a new inactive frame
# which is what we want as the same movement would have turned on 0.9 px in psen2 and hence would have been 0 in the original data


# trace plot before/after -------------------------------------------------

# I want to use ggActivityTraceByGroup
# but somewhat different here as two separated RAWs.csv files
# will try to make a fake RAWs.csv file with unprocessed/adjusted data
dir.create('~/Dropbox/ZFAD/210907_PSEN2/pigmentBias/trace/')

## import frame-by-frame data before/after adjustment
ff1 <- data.table::fread('~/Dropbox/ZFAD/210907_PSEN2/210907_12_RAWs.csv')
fa1 <- data.table::fread('~/Dropbox/ZFAD/210907_PSEN2/210907_12_RAWsadjusted.csv')

## keep only the SCR larvae
# import genotype
geno1 <- importGenotype(genopath='~/Dropbox/ZFAD/210907_PSEN2/210907_12genotype.txt')
scrfis <- geno1$scr
scrfis <- scrfis[!is.na(scrfis)]
# the fish IDs we want to keep are:
scrfis <- sprintf('f%i', scrfis)

## take only those columns from ff1/fa1
ff1 <- as.data.frame(ff1)
fa1 <- as.data.frame(fa1)

ff1scr <- ff1[, scrfis]
# change column names so f1, f2, ..., f47
colnames(ff1scr) <- sprintf('f%i', 1:ncol(ff1scr))

fa1scr <- fa1[, scrfis]
# change column names so f48, f49, ..., f94
colnames(fa1scr) <- sprintf('f%i', (ncol(fa1scr)+1):(ncol(fa1scr)+ncol(fa1scr)))

# cbind all together
ffba <- cbind(ff1[,c('fullts', 'zhrs', 'exsecs')], ff1scr, fa1scr)
# ba for before/after
# in summary: time columns; then f1 to f47 is SCR before; then f48 to f94 is SCR after
# (bit of a dirty solution but should work)
# might complain that there isn't 96 columns
ffba$f95 <- 0
ffba$f96 <- 0

data.table::fwrite(ffba, '~/Dropbox/ZFAD/210907_PSEN2/pigmentBias/trace/210907_12_RAWsSCRbeforeafter.csv', row.names=FALSE)

# now prepare manually genotype file that is
# 1 >> 47 = before
# 48 >> 95 = after

ggActivityTraceByGroup(ffpath='~/Dropbox/ZFAD/210907_PSEN2/pigmentBias/trace/210907_12_RAWsSCRbeforeafter.csv',
                       genopath='~/Dropbox/ZFAD/210907_PSEN2/pigmentBias/trace/210907_12genotypeSCRbeforeafter.txt',
                       smoothOrNo=TRUE,
                       smooth_nsecs=30*60,
                       binOrNo=TRUE,
                       bin_nsecs=10*60,
                       grporder=c('before', 'after'),
                       tracecols=c('#aab5bd', '#697986'),
                       ribboncols=c('#d3dadd', '#b2bbc3'),
                       linethick=0.4,
                       xname='',
                       yname='',
                       xtextOrNo=FALSE,
                       ytextOrNo=FALSE,
                       ymin=0,
                       ymax=35000,
                       xstart=24,
                       xstop=48,
                       trimstart=24,
                       trimstop=48,
                       nightBgOrNo=TRUE,
                       legendOrNo=FALSE,
                       exportOrNo=TRUE,
                       exportPath='~/Dropbox/ZFAD/210907_PSEN2/pubplots/traceBeforeAfter_box12.pdf',
                       width=38,
                       height=55)


#### repeat with box13

## import frame-by-frame data before/after adjustment
ff2 <- data.table::fread('~/Dropbox/ZFAD/210907_PSEN2/210907_13_RAWs.csv')
fa2 <- data.table::fread('~/Dropbox/ZFAD/210907_PSEN2/210907_13_RAWsadjusted.csv')

## keep only the SCR larvae
# import genotype
geno2 <- importGenotype(genopath='~/Dropbox/ZFAD/210907_PSEN2/210907_13genotype.txt')
scrfis <- geno2$scr
scrfis <- scrfis[!is.na(scrfis)]
# box2 so need to +96
scrfis <- scrfis+96
# the fish IDs we want to keep are:
scrfis <- sprintf('f%i', scrfis)

## take only those columns from ff1/fa1
ff2 <- as.data.frame(ff2)
fa2 <- as.data.frame(fa2)

ff2scr <- ff2[, scrfis]
# change column names so f1, f2, ..., f48
colnames(ff2scr) <- sprintf('f%i', 1:ncol(ff2scr))

fa2scr <- fa2[, scrfis]
# change column names so f49, f50, ..., f96
colnames(fa2scr) <- sprintf('f%i', (ncol(fa2scr)+1):(ncol(fa2scr)+ncol(fa2scr)))

# cbind all together
ffba <- cbind(ff2[,c('fullts', 'zhrs', 'exsecs')], ff2scr, fa2scr)
# ba for before/after
# in summary: time columns; then f1 to f47 is SCR before; then f48 to f94 is SCR after
# (bit of a dirty solution but should work)

data.table::fwrite(ffba, '~/Dropbox/ZFAD/210907_PSEN2/pigmentBias/trace/210907_13_RAWsSCRbeforeafter.csv', row.names=FALSE)

# now prepare manually genotype file that is
# 1 >> 48 = before
# 49 >> 96 = after

ggActivityTraceByGroup(ffpath='~/Dropbox/ZFAD/210907_PSEN2/pigmentBias/trace/210907_13_RAWsSCRbeforeafter.csv',
                       genopath='~/Dropbox/ZFAD/210907_PSEN2/pigmentBias/trace/210907_13genotypeSCRbeforeafter.txt',
                       smoothOrNo=TRUE,
                       smooth_nsecs=30*60,
                       binOrNo=TRUE,
                       bin_nsecs=10*60,
                       grporder=c('before', 'after'),
                       tracecols=c('#aab5bd', '#697986'),
                       ribboncols=c('#d3dadd', '#b2bbc3'),
                       linethick=0.4,
                       xname='',
                       yname='activity (sum of Δ px/10 min)',
                       xtextOrNo=FALSE,
                       ytextOrNo=TRUE,
                       ymin=0,
                       ymax=35000,
                       xstart=24,
                       xstop=48,
                       trimstart=24,
                       trimstop=48,
                       nightBgOrNo=TRUE,
                       legendOrNo=FALSE,
                       exportOrNo=TRUE,
                       exportPath='~/Dropbox/ZFAD/210907_PSEN2/pubplots/traceBeforeAfter_box13.pdf',
                       width=45,
                       height=55)


