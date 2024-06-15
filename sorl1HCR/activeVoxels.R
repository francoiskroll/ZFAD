#####################################################
# ~ ZFAD: analysis of SORL1 HCR data ~
#
#
# Francois Kroll 2023
# francois@kroll.be
#####################################################


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

tph_paths <- list.files(here('sorl1HCR', 'activeVoxels_v3', 'tph'), full.names=TRUE)


tphL <- lapply(1:length(tph_paths), function(i) {
  
  # file name is Fx_grp_channel
  filnm <- list.files(here('sorl1HCR', 'activeVoxels_v3', 'tph'))[i]
  
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
  
  # calculate proportion active voxels
  dat <- dat %>%
    mutate(proActive=Active.Voxels/Total.Voxels) # ! R will add . to avoid the space
  
  return(dat)
})

# rbind all tables together
tph <- do.call(rbind, tphL)

# delete original 'Ratio (%)' column
tph$Ratio.... <- NULL
colnames(tph) <- c('fid', 'grp', 'anatomy', 'activeVox', 'totalVox', 'proActive')

# anatomy are e.g. Torus semicircularis6.tif, remove .tif
tph$anatomy <- gsub('.tif', '', tph$anatomy)


# import sert data --------------------------------------------------------

ser_paths <- list.files(here('sorl1HCR', 'activeVoxels_v3', 'sert'), full.names=TRUE)

serL <- lapply(1:length(ser_paths), function(i) {
  
  # file name is Fx_grp_channel
  filnm <- list.files(here('sorl1HCR', 'activeVoxels_v3', 'sert'))[i]
  
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
  
  # calculate proportion active voxels
  dat <- dat %>%
    mutate(proActive=Active.Voxels/Total.Voxels) # ! R will add . to avoid the space
  
  return(dat)
})

# rbind all tables together
ser <- do.call(rbind, serL)

# delete original 'Ratio (%)' column
ser$Ratio.... <- NULL
colnames(ser) <- c('fid', 'grp', 'anatomy', 'activeVox', 'totalVox', 'proActive')

# anatomy are e.g. Torus semicircularis6.tif, remove .tif
ser$anatomy <- gsub('.tif', '', ser$anatomy)


rm(serL)
rm(tphL)


# remove regions where unexpressed ----------------------------------------

# dilemma here is: do not want to test every region statistically as this will make significance threshold after multiple hypothesis testing correction very difficult to meet
# but at the same time do not want to miss anything interesting
# e.g. ectopic expression in sorl1? might miss if simply delete regions where no expression in scr
# I think clear from unfiltered plot that there is no such case, so will set a threshold on SCR expression
# e.g. all SCR brains must have at least 1% active voxels

# solution1: remove anatomical regions that have below 1% active voxels in at least one SCR fish AND one SORL1 fish
# edit, I do not think this is a good idea, e.g. pineal complex has a few at 0, but clear expression in many of the brains

# solution2: simply < threshold in *every* brain imaged

# for ser,
# loop through anatomical region,
# for each, gather data from all fish and ask whether all < threshold
serfilter <- sapply(unique(ser$anatomy), function(an) {
  all(ser[which(ser$grp=='scr' & ser$anatomy==an), 'proActive'] < 0.05)
})

# ! outlier f8 sorl1, abnormally high
# I think rescaling causes it somehow
# therefore, will also help to only look at SCR expression when deciding to exclude anatomical region or not

# the TRUE regions are the ones we want to exclude
# write them to drive so can be used in intensityVoxels.R
write.csv(data.frame(toexclude=sort(names(serfilter)[which(serfilter)])), file=here('sorl1HCR', 'sert_regionstoexclude.csv'), row.names=FALSE)
# 25/05/2023 new version of tiffs from Josh: was excluding 152 regions, now excluding 127 regions
# (with threshold at 0.05)

# keep the ones that are FALSE, this is the ones we want to keep
serf <- ser %>%
  filter(anatomy %in% names(serfilter)[!serfilter])
# checked on plot, filter works
# but might need to put threshold a bit higher as it keeps some regions where clearly no robust expression

# originally, 168 anatomical regions
length(unique(ser$anatomy))
# filter at 0.01 leaves 37
# filter at 0.02 leaves 26
# filter at 0.03 leaves 23
length(unique(serf$anatomy))
# looking at plot, I think 0.03 seems appropriate, 0.02 leaves some regions where no/very low expression


# for tph,
# loop through anatomical region,
# for each, gather data from all fish and ask whether all < threshold
tphfilter <- sapply(unique(tph$anatomy), function(an) {
  all(tph[which(tph$grp=='scr' & tph$anatomy==an), 'proActive'] < 0.05)
})

# the TRUE regions are the ones we want to exclude
# write them to drive so can be used in intensityVoxels.R
write.csv(data.frame(toexclude=sort(names(tphfilter)[which(tphfilter)])), file=here('sorl1HCR', 'tph_regionstoexclude.csv'), row.names=FALSE)
# 25/05/2023 new version of tiffs from Josh: was excluding 152 regions, now excluding 127 regions
# (with threshold at 0.05)

# keep the ones that are FALSE, this is the ones we want to keep
tphf <- tph %>%
  filter(anatomy %in% names(tphfilter)[!tphfilter])
# checked on plot, filter works
# but might need to put threshold a bit higher as it keeps some regions where clearly no robust expression

# originally, 168 anatomical regions
length(unique(tph$anatomy))
# filter at 0.01 leaves 34
# filter at 0.02 leaves 30
# filter at 0.03 leaves 27
length(unique(tphf$anatomy))
# looking at plot, I think 0.03 seems appropriate, 0.02 leaves some regions where no/very low expression



# plot ser ----------------------------------------------------------------

## ser, before filtering anatomical regions
ggser <- ggplot(ser, aes(x=anatomy, y=proActive, colour=grp)) +
  geom_quasirandom(width=0.5) +
  theme_minimal() +
  theme(
    axis.text.x=element_text(angle=45, hjust=1),
    panel.grid.minor.y=element_blank(),
    axis.title.x=element_blank()
  ) +
  geom_hline(yintercept=0.05)
ggser
ggsave(here('sorl1HCR', 'plots', 'sertVox_nofilt.pdf'), ggser, width=500, height=200, unit='mm', device=cairo_pdf)

## ser, after filtering
ggserf <- ggplot(serf, aes(x=proActive, y=anatomy, colour=grp, group=grp)) +
  geom_quasirandom(width=0.1, size=0.8, dodge.width=0.5) +
  theme_minimal() +
  theme(
    axis.text.x=element_text(size=7, margin=margin(t=-1, r=0, b=0, l=0)),
    axis.text.y=element_blank(),
    panel.grid.minor.x=element_blank(),
    axis.title.x=element_text(size=9, margin=margin(t=0, r=0, b=0, l=0)),
    axis.title.y=element_blank(),
    legend.position='none'
  ) +
  #geom_hline(yintercept=0.1) +
  scale_x_continuous(breaks=c(0, 0.25, 0.5, 0.75, 1), labels=c(0, 25, 50, 75, 100)) +
  scale_colour_manual(values=c('#697a87', '#EE7163')) +
  xlab('positive voxels (%)') +
  coord_cartesian(xlim=c(0.02, 0.9))
ggserf
ggsave(here('sorl1HCR', 'plots', 'serVox_h.pdf'), ggserf, width=69, height=55, unit='mm', device=cairo_pdf)




# plot tph ----------------------------------------------------------------

## tph, before filtering anatomical regions
ggtph <- ggplot(tph, aes(x=anatomy, y=proActive, colour=grp)) +
  geom_quasirandom(width=0.5) +
  theme_minimal() +
  theme(
    axis.text.x=element_text(angle=45, hjust=1),
    panel.grid.minor.y=element_blank(),
    axis.title.x=element_blank()
  ) +
  geom_hline(yintercept=0.03)
ggtph
ggsave(here('sorl1HCR', 'plots', 'tphVox_nofilt.pdf'), ggtph, width=500, height=200, unit='mm', device=cairo_pdf)


## tph, after filtering
ggtphf <- ggplot(tphf, aes(x=proActive, y=anatomy, colour=grp, group=grp)) +
  geom_quasirandom(width=0.1, size=0.6, dodge.width=0.5) +
  theme_minimal() +
  theme(
    axis.text.x=element_text(size=7, margin=margin(t=-1, r=0, b=0, l=0)),
    axis.text.y=element_blank(),
    axis.title.x=element_text(size=9, margin=margin(t=0, r=0, b=0, l=0)),
    panel.grid.minor.x=element_blank(),
    axis.title.y=element_blank(),
    legend.position='none'
  ) +
  #geom_hline(yintercept=0.1) +
  scale_x_continuous(breaks=c(0, 0.25, 0.5, 0.75, 1), labels=c(0, 25, 50, 75, 100)) +
  scale_colour_manual(values=c('#697a87', '#EE7163')) +
  xlab('positive voxels (%)') +
  coord_cartesian(xlim=c(0.02, 0.9))
ggtphf
ggsave(here('sorl1HCR', 'plots', 'tphVox_h.pdf'), ggtphf, width=69, height=55, unit='mm', device=cairo_pdf)

### about Ns
tphf %>%
  filter(anatomy=='Thalamus') %>%
  group_by(grp) %>%
  tally()


# note, some brains are at 0 or very low, even for regions where good expression
# not sure what to do with them
# is it always the same?

tphf %>%
  filter(activeVox<0.03)
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
  
  t.test(proActive ~ grp, data=dan)
  
})
# all ns
# below 0.1 is
# Rhombomere 1,2 - Ventral2


sapply(unique(tphf$anatomy), function(an) {
  
  # look at data for only that anatomical region
  dan <- tphf %>%
    filter(anatomy==an)
  
  t.test(proActive ~ grp, data=dan)
  
})
# all ns
# below 0.1 is
# Rhombomere 2 - Ventromedial
# Hypothalamus - Caudal Zone1



# LME statistics ----------------------------------------------------------

### tph ###
# alternative would be to do an LME test across regions
# logic is different anatomical regions are like repeated measures, same as different day/night for behaviour
# I think we also need to add fish ID here as one fish may be fainter across anatomical regions
# this is very close to tutorial https://arxiv.org/pdf/1308.5499.pdf
# example: pitch	~	politeness	+	sex	+	(1|subject)	+	(1|item)	+	ε
# in tutorial case, experiment was: ask people to say different things (item) in different tone (politeness) to see if their voice (pitch) change
# very similar design here: item is anatomy, as each subject (fish) is measured once per region
# but each subject may have a different pitch at baseline; same as different fish may have different brightness at baseline
lm <- lmer(proActive ~ grp + (1|fid) + (1|anatomy), data=tphf, REML=FALSE)
lmnull <- lmer(proActive ~ 1 + (1|fid) + (1|anatomy), data=tphf, REML=FALSE)
summary(lm)
anova(lmnull, lm) # ns
# sorl1 +0.0119
# p-val 0.08041

### sert ###
lm <- lmer(proActive ~ grp + (1|fid) + (1|anatomy), data=serf, REML=FALSE)
lmnull <- lmer(proActive ~ 1 + (1|fid) + (1|anatomy), data=serf, REML=FALSE)
summary(lm)
anova(lmnull, lm) # ns
# sorl1 +0.021891
# p-val 0.001333 **


# ns ----------------------------------------------------------------------

serf %>%
  filter(anatomy=='Superior raphe') %>%
  group_by(grp) %>%
  tally()
