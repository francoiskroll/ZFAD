###################################################
# ~~~ ZFAD/FramebyFrame package ~~~

# Comparison sleep parameters frame-by-frame vs one-minute activity analysis

# Francois Kroll 2022
# francois@kroll.be
###################################################

# no need to have thousands of larvae,
# just use one good experiment

# using
# 220531_14_SORL1

# v1: one this script was taking sleep parameters from the .mat file obtained from sleep_analysis2020.m written by Jason
# to compare exactly what we do in the Rihel lab
# however that introduces a new parameter which is Freeze threshold
# sleep may be *higher* in middur analysis because it misses movements below Freeze threshold

# v2: while v1 is the most accurate comparison for the Rihel lab,
# it would be clearer for the ZFAD paper not to introduce the Freeze threshold, which makes the explanations difficult
# so, will generate middur data with Freeze threshold at 0
# then measure sleep parameters, which we should be able to do using the functions written for predPharma


# packages ----------------------------------------------------------------


# library(devtools)
# install_github('francoiskroll/FramebyFrame')
library(FramebyFrame)

library(R.matlab)
library(here)
library(data.table)
library(dplyr)
library(tibble)
library(tidyr)
library(naturalsort)
library(ggplot2)
library(magrittr)

# sleep parameters
allpanm <- c('sleepHours', 'sleepNumNaps', 'sleepNapDuration', 'sleepLatency')

# directory for plots
dir.create(here('FbyFcomparison', 'pubplots'), showWarnings=FALSE)



# calculate sleep parameters by FramebyFrame ------------------------------

multiBehaviourParameter(parameters='sleep',
                        ffpath=here('FbyFcomparison', 'data', '220531_14_RAWs.csv'),
                        genopath=here('FbyFcomparison', 'data', '220531_14genotype.txt'),
                        skipNight0=FALSE,
                        dayduration=14)

# there are recorded in Environment
# make a list from there
ff1 <- list(sleepHours_220531_14, sleepNumNaps_220531_14, sleepNapDuration_220531_14, sleepLatency_220531_14)
names(ff1) <- allpanm


# calculate middur data ---------------------------------------------------
# with Freeze threshold at 0

# rawToMiddur(ffpath=here('220531_SORL1', '220531_14_RAWs.csv'),
#             zebpath=here('220531_SORL1', '220531_14_15_sorl1.xls'),
#             freezing=0, # !
#             burst=200,
#             exportOrNo=TRUE)

# ! there may be 220531_14_middur.csv in the folder already calculated for predPharma prediction
# move them in another folder to avoid re-writing
# then rename the one we just created to ..._freeze0.csv


# get sleep parameters from middur ----------------------------------------

source('~/Dropbox/predPharma/paramsFromMid_inaThr0.R') # ! only source it now as it creates conflicts with FramebyFrame
# !!!
# special version duplicated 25/08/2023 with inaThr to = 0
# it should be = 0.1 to match SCRAP.m, but here trying to make the simplest possible comparison one-minute vs FramebyFrame
# so not using 0 would add another parameter to track
# set it back to 0.1 after

# using function from predPharma project, script paramsFromMid.R
# sourced above
mid <- read.csv(here('220531_SORL1', '220531_14_middurFreeze0.csv'))
paral <- calculateParameters(mid=mid,
                             genopath=here('220531_SORL1', '220531_14genotype.txt'))


# keep only the sleep parameters
# re-create the list so we are sure about the order too
paral <- list(paral$sleep, paral$sleepBout, paral$sleepLength, paral$sleepLatency)
names(paral) <- c('sleepHours', 'sleepNumNaps', 'sleepNapDuration', 'sleepLatency')

# ! sleepHours is currently in minutes, convert it
paral$sleepHours[, c('night0', 'day1', 'night1', 'day2', 'night2')] <- paral$sleepHours[, c('night0', 'day1', 'night1', 'day2', 'night2')] / 60


# check we have the same fish IDs in mid & ff parameters ------------------

identical(paral$sleepHours$fish, ff1$sleepHours$fish)


# prepare tables for plots ------------------------------------------------

# will be a bit manual
# I think best is to have four tables

# exp1 -- day, all parameters
# exp1 -- night, all parameters
# exp2 -- day, all parameters
# exp2 -- night, all parameters

# then can do one figure for each experiment
# for each, grid of two plots
# night, with facet by parameter
# day, with facet by parameter

# EXP1: we will analyse day2 & night2

# exp1 -- day
e1day <- rbind(

  data.frame(parameter=ff1$sleepHours$parameter,
             fid=ff1$sleepHours$fish,
             grp=ff1$sleepHours$grp,
             ff_day2=ff1$sleepHours$day2,
             mid_day2=paral$sleepHours$day2),

  data.frame(parameter=ff1$sleepNumNaps$parameter,
             fid=ff1$sleepNumNaps$fish,
             grp=ff1$sleepNumNaps$grp,
             ff_day2=ff1$sleepNumNaps$day2,
             mid_day2=paral$sleepNumNaps$day2),

  data.frame(parameter=ff1$sleepNapDuration$parameter,
             fid=ff1$sleepNapDuration$fish,
             grp=ff1$sleepNapDuration$grp,
             ff_day2=ff1$sleepNapDuration$day2,
             mid_day2=paral$sleepNapDuration$day2)
)
# note, sleep latency is not defined for day in FramebyFrame

# e1day <- e1day %>%
#   mutate(param_fid=paste(parameter, fid, sep='_'), .before='parameter')

# exp1 -- night
e1night <- rbind(

  data.frame(parameter=ff1$sleepHours$parameter,
             fid=ff1$sleepHours$fish,
             grp=ff1$sleepHours$grp,
             ff_night2=ff1$sleepHours$night2,
             mid_night2=paral$sleepHours$night2),

  data.frame(parameter=ff1$sleepNumNaps$parameter,
             fid=ff1$sleepNumNaps$fish,
             grp=ff1$sleepNumNaps$grp,
             ff_night2=ff1$sleepNumNaps$night2,
             mid_night2=paral$sleepNumNaps$night2),

  data.frame(parameter=ff1$sleepNapDuration$parameter,
             fid=ff1$sleepNapDuration$fish,
             grp=ff1$sleepNapDuration$grp,
             ff_night2=ff1$sleepNapDuration$night2,
             mid_night2=paral$sleepNapDuration$night2),

  data.frame(parameter=ff1$sleepLatency$parameter,
             fid=ff1$sleepLatency$fish,
             grp=ff1$sleepLatency$grp,
             ff_night2=ff1$sleepLatency$night2,
             mid_night2=paral$sleepLatency$night2)
)



# experiment1: night plots ------------------------------------------------

######
# some common settings
dotsize <- 2
linesize <- 0.3
width <- 90
height <- 65
######

e1nightl <- e1night %>%
  filter(grp=='scr') %>% # can keep only SCR to make it less cluttered
  pivot_longer(-c(parameter, fid, grp),
               names_to='source',
               values_to='val')

e1nightl$parameter <- factor(e1nightl$parameter, levels=c('sleepNumNaps', 'sleepNapDuration', 'sleepLatency', 'sleepHours'))
e1nightl$source <- factor(e1nightl$source, levels=c('mid_night2', 'ff_night2'))

# prepare titles
titls <- c(sleepNumNaps='number of \nsleep bouts', sleepNapDuration='sleep bout \nlength', sleepLatency='sleep latency', sleepHours='total sleep')

ggzcom <- ggplot(e1nightl, aes(x=source, y=val, group=fid)) +
  facet_wrap(~parameter, scales='free_y', nrow=1,
             labeller=labeller(parameter=titls)) +
  geom_line(size=linesize, colour='#595E60') +
  geom_point(pch=21, size=dotsize, colour='white', fill='#595E60') +
  scale_x_discrete(labels=c('1-min', 'FbyF'), expand=c(0.1,0.1)) +
  theme_minimal() +
  theme(
    panel.background=element_rect(fill='#d2d2d2', colour=NA),
    strip.text.x=element_text(size=7, hjust=0.5, margin=margin(t=0, r=0, b=10, l=0)),
    axis.title.x=element_blank(),
    axis.title.y=element_blank(),
    axis.text.x=element_text(size=7),
    axis.text.y=element_text(size=7, margin=margin(t=0, r=-1, b=0, l=0)),
    panel.grid.minor.y=element_blank(),
    panel.spacing = unit(15, 'pt')
  )
ggzcom

ggsave(here('FbyFcomparison', 'pubplots', 'exp1_zCom.pdf'), ggzcom, width=width, height=height, units='mm', device=cairo_pdf)



# experiment1: night statistics -------------------------------------------

e1nitot <- e1night %>%
  filter(parameter=='sleepHours') %>%
  filter(grp=='scr')
mean(e1nitot$ff_night2 - e1nitot$mid_night2)
sd(e1nitot$ff_night2 - e1nitot$mid_night2)
# % error?
mean( (e1nitot$ff_night2 - e1nitot$mid_night2) / e1nitot$mid_night2 )
sd ((e1nitot$ff_night2 - e1nitot$mid_night2) / e1nitot$mid_night2 )


e1nico <- e1night %>%
  filter(parameter=='sleepNumNaps') %>%
  filter(grp=='scr')
mean(e1nico$ff_night2 - e1nico$mid_night2)
sd(e1nico$ff_night2 - e1nico$mid_night2)
# % error?
mean( (e1nico$ff_night2 - e1nico$mid_night2) / e1nico$mid_night2 )
sd ((e1nico$ff_night2 - e1nico$mid_night2) / e1nico$mid_night2 )
# assuming FbyF is real number of sleep bouts
# how much in % did one-min catch?
e1nico$mid_night2 / e1nico$ff_night2
# so sleep bouts missed are
mean(1 - (e1nico$mid_night2 / e1nico$ff_night2))
sd(1 - (e1nico$mid_night2 / e1nico$ff_night2))

e1nidur <- e1night %>%
  filter(parameter=='sleepNapDuration') %>%
  filter(grp=='scr')
mean(e1nidur$ff_night2 - e1nidur$mid_night2)
sd(e1nidur$ff_night2 - e1nidur$mid_night2)
# % error?
mean( (e1nidur$ff_night2 - e1nidur$mid_night2) / e1nidur$mid_night2 )
sd ((e1nidur$ff_night2 - e1nidur$mid_night2) / e1nidur$mid_night2 )

e1nilat <- e1night %>%
  filter(parameter=='sleepLatency') %>%
  filter(grp=='scr')
mean(e1nilat$ff_night2 - e1nilat$mid_night2)
sd(e1nilat$ff_night2 - e1nilat$mid_night2)
# % error?
mean( (e1nilat$ff_night2 - e1nilat$mid_night2) / e1nilat$mid_night2 )
sd ((e1nilat$ff_night2 - e1nilat$mid_night2) / e1nilat$mid_night2 )


### correlations
# using weird pipe from magrittr %$%, dplyr pipe does not work with correlation for some reason
e1night %>%
  filter(parameter=='sleepHours') %$%
  cor(.$ff_night2, .$mid_night2, use='complete.obs')

e1night %>%
  filter(parameter=='sleepNapDuration') %$%
  cor(.$ff_night2, .$mid_night2, use='complete.obs')

e1night %>%
  filter(parameter=='sleepNumNaps') %$%
  cor(.$ff_night2, .$mid_night2, use='complete.obs')

e1night %>%
  filter(parameter=='sleepLatency') %$%
  cor(.$ff_night2, .$mid_night2, use='complete.obs')




# experiment1: day plots --------------------------------------------------

e1dayl <- e1day %>%
  filter(grp=='scr') %>% # can keep only SCR to make it less cluttered
  pivot_longer(-c(parameter, fid, grp),
               names_to='source',
               values_to='val')

e1dayl$parameter <- factor(e1dayl$parameter, levels=c('sleepNumNaps', 'sleepNapDuration', 'sleepLatency', 'sleepHours'))
e1dayl$source <- factor(e1dayl$source, levels=c('mid_day2', 'ff_day2'))

# prepare titles
titls <- c(sleepNumNaps='number of \nsleep bouts', sleepNapDuration='sleep bout \nlength', sleepLatency='sleep latency', sleepHours='total sleep')

ggzcom <- ggplot(e1dayl, aes(x=source, y=val, group=fid)) +
  facet_wrap(~parameter, scales='free_y', nrow=1,
             labeller=labeller(parameter=titls)) +
  geom_line(size=linesize, colour='#595E60') +
  geom_point(pch=21, size=dotsize, colour='white', fill='#595E60') +
  scale_x_discrete(labels=c('1-min', 'FbyF'), expand=c(0.1,0.1)) +
  theme_minimal() +
  theme(
    #panel.background=element_rect(fill='#d2d2d2', colour=NA),
    strip.text.x=element_text(size=7, hjust=0.5, margin=margin(t=0, r=0, b=10, l=0)),
    axis.title.x=element_blank(),
    axis.title.y=element_blank(),
    axis.text.x=element_text(size=7),
    axis.text.y=element_text(size=7, margin=margin(t=0, r=-1, b=0, l=0)),
    panel.grid.minor.y=element_blank(),
    panel.spacing = unit(15, 'pt')
  )
ggzcom

ggsave(here('FbyFcomparison', 'pubplots', 'exp1_zComDay.pdf'), ggzcom, width=67, height=height, units='mm', device=cairo_pdf)




# experiment1: day statistics ---------------------------------------------

e1daytot <- e1day %>%
  filter(parameter=='sleepHours') %>%
  filter(grp=='scr')
mean(e1daytot$ff_day2 - e1daytot$mid_day2)
sd(e1daytot$ff_day2 - e1daytot$mid_day2)

e1dayco <- e1day %>%
  filter(parameter=='sleepNumNaps') %>%
  filter(grp=='scr')
mean(e1dayco$ff_day2 - e1dayco$mid_day2)
sd(e1dayco$ff_day2 - e1dayco$mid_day2)

e1daydur <- e1day %>%
  filter(parameter=='sleepNapDuration') %>%
  filter(grp=='scr')
mean(e1daydur$ff_day2 - e1daydur$mid_day2, na.rm=TRUE)
sd(e1daydur$ff_day2 - e1daydur$mid_day2, na.rm=TRUE)


### correlations
# using weird pipe from magrittr %$%, dplyr pipe does not work with correlation for some reason
e1day %>%
  filter(parameter=='sleepHours') %$%
  cor(.$ff_day2, .$mid_day2, use='complete.obs')

e1day %>%
  filter(parameter=='sleepNapDuration') %$%
  cor(.$ff_day2, .$mid_day2, use='complete.obs')

e1day %>%
  filter(parameter=='sleepNumNaps') %$%
  cor(.$ff_day2, .$mid_day2, use='complete.obs')


