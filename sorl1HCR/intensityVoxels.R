#####################################################
# ~ ZFAD: analysis of SORL1 HCR data ~
#
#
# Francois Kroll 2023
# francois@kroll.be
#####################################################

# script is based on activeVoxels.R which are analysing counts of active voxels within masks
# here, analysing actual voxel intensity
# calculated with script sumSignalMasks.ipynb I created by editing script from Chintan


# packages ----------------------------------------------------------------

library(here)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(ggbeeswarm)
library(lme4)

dir.create(here('sorl1HCR', 'plots'))


# function strNthSplit ----------------------------------------------------

strNthSplit <- function(stri,
                        split,
                        nth) {
  
  # confirm we are given string(s)
  stri <- as.character(unlist(stri))
  
  as.character(sapply(strsplit(stri, split=split),
                      function(s) {
                        s[nth]
                      }))
}



# import tph data ---------------------------------------------------------

# each csv file is data for one fish / one channel (sert genes or tph genes)
# for each channel, loop through files and import in a list

tph_paths <- list.files(here('sorl1HCR', 'intensityVoxels', 'tph'), full.names=TRUE)


tphL <- lapply(1:length(tph_paths), function(i) {
  
  # file name is Fx_grp_channel
  filnm <- list.files(here('sorl1HCR', 'intensityVoxels', 'tph'))[i]
  
  # get the fX component
  fx <- tolower(strNthSplit(filnm, '_', 1))
  
  # get the grp component
  grp <- tolower(strNthSplit(filnm, '_', 2))

  # import the data
  dat <- read.csv(tph_paths[i])
  
  # add fX and grp as columns
  dat <- dat %>%
    mutate(grp=grp, .before=1) %>%
    mutate(fid=fx, .before=1)
  
  return(dat)
})

# rbind all tables together
tph <- do.call(rbind, tphL)

# anatomy are e.g. Torus semicircularis6.tif, remove .tif
tph$anatomy <- gsub('.tif', '', tph$anatomy)


# import sert data --------------------------------------------------------

ser_paths <- list.files(here('sorl1HCR', 'intensityVoxels', 'sert'), full.names=TRUE)

serL <- lapply(1:length(ser_paths), function(i) {
  
  # file name is Fx_grp_channel
  filnm <- list.files(here('sorl1HCR', 'intensityVoxels', 'sert'))[i]
  
  # get the fX component
  fx <- tolower(strNthSplit(filnm, '_', 1))
  
  # get the grp component
  grp <- tolower(strNthSplit(filnm, '_', 2))

  # import the data
  dat <- read.csv(ser_paths[i])
  
  # add fX and grp as columns
  dat <- dat %>%
    mutate(grp=grp, .before=1) %>%
    mutate(fid=fx, .before=1)
  
  return(dat)
})

# rbind all tables together
ser <- do.call(rbind, serL)

# anatomy are e.g. Torus semicircularis6.tif, remove .tif
ser$anatomy <- gsub('.tif', '', ser$anatomy)


rm(serL)
rm(tphL)


# remove regions where unexpressed ----------------------------------------

# here, will use activeVoxels analysis
# make sure to run activeVoxels.R first so it writes files which list anatomical regions below threshold

# import the files with regions to exclude
ser_exclude <- read.csv(file=here('sorl1HCR', 'sert_regionstoexclude.csv'))
tph_exclude <- read.csv(file=here('sorl1HCR', 'tph_regionstoexclude.csv'))

# for ser, remove the regions to exclude
serf <- ser %>%
  filter(! anatomy %in% ser_exclude$toexclude) # i.e. we keep regions not mentioned in exclude file

# originally, 168 anatomical regions
length(unique(ser$anatomy))
# now
length(unique(serf$anatomy))

###

# for tph, remove the regions to exclude
tphf <- tph %>%
  filter(! anatomy %in% tph_exclude$toexclude) # i.e. we keep regions not mentioned in exclude file

# originally, 168 anatomical regions
length(unique(tph$anatomy))
# now
length(unique(tphf$anatomy))



# plot ser ----------------------------------------------------------------

## ser, before filtering anatomical regions
ggser <- ggplot(ser, aes(x=anatomy, y=sumVx, colour=grp)) +
  geom_quasirandom(width=0.5) +
  theme_minimal() +
  theme(
    axis.text.x=element_text(angle=45, hjust=1),
    panel.grid.minor.y=element_blank(),
    axis.title.x=element_blank()
  ) +
  geom_hline(yintercept=0.03)
ggser
ggsave(here('sorl1HCR', 'plots', 'sert_vxsum_nofilt.pdf'), ggser, width=500, height=200, unit='mm', device=cairo_pdf)

ggser <- ggplot(ser, aes(x=anatomy, y=meanVx, colour=grp)) +
  geom_quasirandom(width=0.5) +
  theme_minimal() +
  theme(
    axis.text.x=element_text(angle=45, hjust=1),
    panel.grid.minor.y=element_blank(),
    axis.title.x=element_blank()
  ) +
  geom_hline(yintercept=0.03)
ggser
ggsave(here('sorl1HCR', 'plots', 'sert_vxmean_nofilt.pdf'), ggser, width=500, height=200, unit='mm', device=cairo_pdf)

## ser, after filtering
ggserf <- ggplot(serf, aes(x=sumVx, y=anatomy, colour=grp, group=grp)) +
  geom_quasirandom(width=0.1, size=0.8, dodge.width=0.5) +
  theme_minimal() +
  theme(
    axis.text.x=element_text(size=7, margin=margin(t=-1, r=0, b=0, l=0)),
    axis.text.y=element_text(size=6, margin=margin(t=0, r=-2, b=0, l=0), colour='#020304'),
    panel.grid.minor.x=element_blank(),
    axis.title.x=element_text(size=9, margin=margin(t=0, r=0, b=0, l=0)),
    axis.title.y=element_blank(),
    legend.position='none'
  ) +
  scale_colour_manual(values=c('#697a87', '#EE7163')) +
  coord_cartesian(xlim=c(0, 600000)) +
  scale_x_continuous(breaks=c(0, 200000, 400000, 600000), labels=c('0', '200,000', '400,000', '600,000')) +
  xlab('total intensity')
ggserf
ggsave(here('sorl1HCR', 'plots', 'sert_vxsum_filt.pdf'), ggserf, width=95, height=55, unit='mm', device=cairo_pdf)

ggserf <- ggplot(serf, aes(x=meanVx, y=anatomy, colour=grp, group=grp)) +
  geom_quasirandom(width=0.1, size=0.2, dodge.width=0.5) +
  theme_minimal() +
  theme(
    axis.text.x=element_text(size=7, margin=margin(t=-1, r=0, b=0, l=0)),
    axis.text.y=element_text(size=5, margin=margin(t=0, r=-2, b=0, l=0)),
    panel.grid.minor.x=element_blank(),
    axis.title.y=element_blank(),
    legend.position='none'
  ) +
  scale_colour_manual(values=c('#697a87', '#EE7163')) +
  xlab('average intensity')
ggserf
ggsave(here('sorl1HCR', 'plots', 'sert_vxmean_filt.pdf'), ggserf, width=95, height=70, unit='mm', device=cairo_pdf)



# plot tph ----------------------------------------------------------------

## tph, before filtering anatomical regions
ggtph <- ggplot(tph, aes(x=anatomy, y=sumVx, colour=grp)) +
  geom_quasirandom(width=0.5) +
  theme_minimal() +
  theme(
    axis.text.x=element_text(angle=45, hjust=1),
    panel.grid.minor.y=element_blank(),
    axis.title.x=element_blank()
  ) +
  geom_hline(yintercept=0.03)
ggtph
ggsave(here('sorl1HCR', 'plots', 'tph_vxsum_nofilt.pdf'), ggtph, width=500, height=200, unit='mm', device=cairo_pdf)


## tph, after filtering
ggtphf <- ggplot(tphf, aes(x=sumVx, y=anatomy, colour=grp, group=grp)) +
  geom_quasirandom(width=0.1, size=0.6, dodge.width=0.5) +
  theme_minimal() +
  theme(
    axis.text.x=element_text(size=7, margin=margin(t=-1, r=0, b=0, l=0)),
    axis.text.y=element_text(size=6, margin=margin(t=0, r=-2, b=0, l=0), colour='#020304'),
    # axis.text.y=element_blank(),
    axis.title.x=element_text(size=9, margin=margin(t=0, r=0, b=0, l=0)),
    panel.grid.minor.x=element_blank(),
    axis.title.y=element_blank(),
    legend.position='none'
  ) +
  scale_colour_manual(values=c('#697a87', '#EE7163')) +
  coord_cartesian(xlim=c(0, 1300000)) +
  scale_x_continuous(breaks=c(0, 500000, 1000000), labels=c(0, '500,000', '1,000,000')) +
  xlab('total intensity')
ggtphf
ggsave(here('sorl1HCR', 'plots', 'tph_vxsum_filt.pdf'), ggtphf, width=95, height=55, unit='mm', device=cairo_pdf)

### ~ vxMean ~ ###

## tph, after filtering
ggtphf <- ggplot(tphf, aes(x=meanVx, y=anatomy, colour=grp, group=grp)) +
  geom_quasirandom(width=0.1, size=0.2, dodge.width=0.5) +
  theme_minimal() +
  theme(
    axis.text.x=element_text(size=7, margin=margin(t=-1, r=0, b=0, l=0)),
    axis.text.y=element_text(size=5, margin=margin(t=0, r=-2, b=0, l=0)),
    panel.grid.minor.x=element_blank(),
    axis.title.y=element_blank(),
    legend.position='none'
  ) +
  scale_colour_manual(values=c('#697a87', '#EE7163')) +
  xlab('total intensity')
ggtphf
ggsave(here('sorl1HCR', 'plots', 'tph_vxmean_filt.pdf'), ggtphf, width=110, height=70, unit='mm', device=cairo_pdf)

### about Ns
tphf %>%
  filter(anatomy=='Thalamus') %>%
  group_by(grp) %>%
  tally()


# note, some brains are at 0 or very low, even for regions where good expression
# not sure what to do with them
# is it always the same?

# tphf %>%
#   filter(activeVox<0.03)
# yes, only a few:
# f1 scr (10x)
# f10 sorl1 (26x, all regions)
# f2 scr (25x, all regions minus 1)
# f5 sorl1 (only once, so do not need to worry)
# f6 scr (25x, all regions minus 1)
# f7 sorl1 (25x, all regions minus 1)
# f8 scr (only once)

# so basically 4 brains where null everywhere
# 2 scr / 2 sorl1, so not related to genotypes

# f2_scr
# f6_scr
# f7_sorl1
# f10_sorl1

# e.g. tph sorl1, there is obvious expression in stack though

# will temporarily exclude them
# but need to understand what is happening

# tphf <- tphf %>%
#   mutate(fg=paste(fid, grp, sep='_'), .before=1)
# 
# tphf <- tphf %>%
#   filter(!fg %in% c('f2_scr', 'f6_scr', 'f7_sorl1', 'f10_sorl1'))


# statistics --------------------------------------------------------------

# we loop through anatomical regions,
# for each, we do a simple t-test SCR vs SORL1
sapply(unique(serf$anatomy), function(an) {
  
  # look at data for only that anatomical region
  dan <- serf %>%
    filter(anatomy==an)
  
  t.test(sumVx ~ grp, data=dan)
  
})
# all ns


sapply(unique(tphf$anatomy), function(an) {
  
  # look at data for only that anatomical region
  dan <- tphf %>%
    filter(anatomy==an)
  
  t.test(sumVx ~ grp, data=dan)
  
})
# all ns except
# Hypothalamus - Caudal Zone1 0.01210539

# effect size in Hypothalamus?
scrHyp <- tphf %>%
  filter(grp=='scr', anatomy=='Hypothalamus - Caudal Zone1') %>%
  summarise_at(vars(sumVx),
               list(
                 mean= ~ mean(., na.rm=TRUE),
                 sd= ~ sd(., na.rm=TRUE),
                 median= ~ median(., na.rm=TRUE),
                 n= ~ length(.)
               ))

sorHyp <- tphf %>%
  filter(grp=='sorl1', anatomy=='Hypothalamus - Caudal Zone1') %>%
  summarise_at(vars(sumVx),
               list(
                 mean= ~ mean(., na.rm=TRUE),
                 sd= ~ sd(., na.rm=TRUE),
                 median= ~ median(., na.rm=TRUE),
                 n= ~ length(.)
               ))
100*(1-(sorHyp$mean / scrHyp$mean))


# LME statistics ----------------------------------------------------------

# see activeVoxels.R for explanations about LME model

### tph ###
lm <- lmer(sumVx ~ grp + (1|fid) + (1|anatomy), data=tphf, REML=FALSE)
lmnull <- lmer(sumVx ~ 1 + (1|fid) + (1|anatomy), data=tphf, REML=FALSE)
summary(lm)
anova(lmnull, lm)
# sorl1 -27578
# p-val 2.154e-07 ***

## convert effect size to % decrease
# mean sumVx in SCR
scrm <- tphf %>%
  filter(grp=='scr') %>%
  group_by(anatomy) %>%
  summarise_at(vars(sumVx),
               list(
                 mean= ~ mean(., na.rm=TRUE),
                 sd= ~ sd(., na.rm=TRUE),
                 median= ~ median(., na.rm=TRUE),
                 n= ~ length(.)
               ))
# can either use the LME effect size or directly calculate percentage increase/decrease
# I think will calculate directly?
sorl1m <- tphf %>%
  filter(grp=='sorl1') %>%
  group_by(grp, anatomy) %>%
  summarise_at(vars(sumVx),
               list(
                 mean= ~ mean(., na.rm=TRUE),
                 sd= ~ sd(., na.rm=TRUE),
                 median= ~ median(., na.rm=TRUE),
                 n= ~ length(.)
               ))
mean(1 - (sorl1m$mean / scrm$mean))
# 6.77%
# which is in same range as LME as found 27,000 out of in average 250-300,000
# so 9-10%
# OR LME method:
# slope -27578
100 * mean(1- ((scrm$mean - 27578) / scrm$mean))

### sert ###
lm <- lmer(sumVx ~ grp + (1|fid) + (1|anatomy), data=serf, REML=FALSE)
lmnull <- lmer(sumVx ~ 1 + (1|fid) + (1|anatomy), data=serf, REML=FALSE)
summary(lm)
anova(lmnull, lm)
# ns



# tangent: is it true sorl1 always has lower max? -------------------------

t.test(maxVx ~ grp, data=ser)
t.test(maxVx ~ grp, data=tph)


# take max of any region for each fish
maxSer <- ser %>%
  group_by(fid, grp) %>%
  summarise_at(vars(maxVx),
               list(
                 max= ~ max(., na.rm=TRUE)
               ))
# now compare
t.test(max ~ grp, data=maxSer)
# yes, bit lower in sorl1
# scr 153.8 vs sorl1 131.1


# take max of any region for each fish
maxTph <- tph %>%
  group_by(fid, grp) %>%
  summarise_at(vars(maxVx),
               list(
                 max= ~ max(., na.rm=TRUE)
               ))
# now compare
t.test(max ~ grp, data=maxTph)
# yes, lower in sorl1
# scr 201.0 vs sorl1 171.4

# >> I think intuition was broadly correct about effect on rescaling