#####################################################
# ~ ZFAD: 230214_sorl1Citalopram drug differential effect ~
#
# testing here differential effect using Euclidean distance in multidimensional space
# idea is:
#
# 1- build (imagine) multidimensional space
# where each dimension is a unique behavioural parameter
# where the origin (point at 0, 0, 0, 0, ...) is for each parameter the SCR + 0 µM mean
# 2- in this multidimensional space, every SCR larva sits at some distance from origin
# 3- taking all SCR + 1 µM larvae, we can say "adding 1 µM citalopram displaced the larvae by mean ± sd from 0 µM"
# 4- taking all SCR + 1 µM larvae, we can say "adding 1 µM citalopram displaced the larvae by mean ± sd from 0 µM"
# 5- repeat 1-, now origin is for each parameter the SORL1 + 0 µM mean
# 6- repeat 2- for SORL1 larvae
# 7- repeat 3- for SORL1 + 1 µM
# 8- repeat 4- for SORL1 + 10 µM
# 9- we can compare these distances, i.e. 3- vs 7- and 4- vs 8-
# evidence of differential effect if displacement in multidimensional space when adding drug on SORL1 is greater/lesser than displacement when adding drug on SCR

# Francois Kroll 2023
# francois@kroll.be
#####################################################


# packages ----------------------------------------------------------------

library(here)

# library(devtools)
# install_github('francoiskroll/FramebyFrame')
library(FramebyFrame)

library(dplyr)
library(tidyr)
library(tibble)
library(data.table)
library(ggplot2)
library(ggbeeswarm)


# calculate larvae fingerprints -------------------------------------------

# (step1) calculate SCR larva fingerprints
# this the fingerprint of each SCR larva in reference to the SCR + 0 µM larvae
fgscr <- calculateFingerprint(paDir=here('230214_sorl1Citalopram', 'bhvparams'),
                              controlGrp='scr_0',
                              singleFish=TRUE,
                              grporder=c('scr_0', 'scr_1', 'scr_10'),
                              skipNight0=TRUE)

# (step5) calculate SORL1 larva fingerprints
# this the fingerprint of each SORL1 larva in reference to the SORL1 + 0 µM larvae
fgso <- calculateFingerprint(paDir=here('230214_sorl1Citalopram', 'bhvparams'),
                             controlGrp='sorl1_0',
                             singleFish=TRUE,
                             grporder=c('sorl1_0', 'sorl1_1', 'sorl1_10'),
                             skipNight0=TRUE)


# calculate Euclidean distances from 0 µM origin --------------------------

# each larva is a point in the multidimensional space
# each dimension is Z-score from controls for one parameter
# e.g. larva #5 is at 2.7, 1.2, 0.1, 0.9, ...
# and we want to calculate the distance between this point and origin at 0, 0, 0, 0, ...
# origin represents the mean of controls (scr_0 or sorl1_0) by definition of the fingerprint

# copying some chunks of fingerprintSimilarity...

# pivot wider so that each column is one larva's fingerprint
# each row is a unique parameter
# or in other words, each column are all the coordinates to place this larva in the space

# add column date_box_grp
fgscr <- fgscr %>%
  mutate(date_box_grp=paste(date, box, grp, sep='_'))

fgscrw <- fgscr %>%
  group_by(date_box_grp) %>%
  pivot_wider( id_cols=uparam,
               names_from=date_box_grp_fish,
               values_from=pazm)

# ! dist() wants two rows
# top row gives the coordinates for the first point, i.e. x, y, z, ...
# bottom row gives the coordinates for the second point, i.e. x, y, z, ...
# so use rbind to achieve this
# as example, here is the distance between one larva and the origin:
dist( rbind(fgscrw$`230214_14_scr_0_f13`, rep(0, nrow(fgscrw))) , method='euclidean')
# rep(0, `number of parameters` times) to give the coordinates of the origin, i.e. 0, 0, 0, 0, ...

# (step2) now we do this for every SCR larva
discr <- apply(fgscrw[,2:ncol(fgscrw)], 2, function(fcoo) { # fcoo for set of coordinates for that fish
  dist( rbind(fcoo, rep(0, nrow(fgscrw))) , method='euclidean')
})

# repeat for SORL1 larvae
# add column date_box_grp
fgso <- fgso %>%
  mutate(date_box_grp=paste(date, box, grp, sep='_'))

fgsow <- fgso %>%
  group_by(date_box_grp) %>%
  pivot_wider( id_cols=uparam,
               names_from=date_box_grp_fish,
               values_from=pazm)

# (step6)
diso <- apply(fgsow[,2:ncol(fgsow)], 2, function(fcoo) { # fcoo for set of coordinates for that fish
  dist( rbind(fcoo, rep(0, nrow(fgsow))) , method='euclidean')
})

### re-arrange the distances into a usable dataframe
discr <- data.frame(fid=names(discr), dist=discr)
row.names(discr) <- NULL

diso <- data.frame(fid=names(diso), dist=diso)
row.names(diso) <- NULL

did <- rbind(discr, diso) # distance dataframe

# recreate some meta columns
# strNthSplit is from FramebyFrame
did <- did %>%
  mutate(exp=paste(strNthSplit(fid, split='_', nth=1),
                   strNthSplit(fid, split='_', nth=2), sep='_'), .before=dist) %>%
  mutate(grp=strNthSplit(fid, split='_', nth=3), .before=dist) %>%
  mutate(dose=strNthSplit(fid, split='_', nth=4), .before=dist) %>%
  mutate(grp_dose=paste(strNthSplit(fid, split='_', nth=3),
                        strNthSplit(fid, split='_', nth=4), sep='_'), .before=dist)

# in summary, we now have:
# for every SCR larva, its distance from the origin which is the mean of the SCR + 0 µM larvae (of its experiment)
# for every SORL1 larva, its distance from the origin which is the mean of the SORL1 + 0 µM larvae (of its experiment)
# Note, remember the fingerprint was calculated using the controls within that experiment
# so there are really 4 reference points/origins being used here: SCR + 0 µM for box14 / SORL1 + 0 µM for box14 / SCR + 0 µM for box15 / SORL1 + 0 µM for box15


# plot --------------------------------------------------------------------

### plot every group
dodgeby <- 0.7

ggdist <- ggplot(did, aes(x=grp_dose, y=dist, colour=grp_dose)) +
  geom_quasirandom(width=0.3, size=2) +
  stat_summary(aes(group=grp_dose), fun=mean, geom='point', colour='#595E60', shape=3, size=1.2, stroke=0.8, position=position_dodge(dodgeby)) +
  facet_grid(~grp, scales='free_x', space='free_x') +
  scale_colour_manual(values=c('#d9dde1', '#b3bcc3', '#697a87',
                               '#FDEFEE', '#F6B8B1', '#EE7163')) +
  theme_minimal() +
  theme(
    panel.grid.minor.y=element_blank(),
    axis.text.x=element_blank(),
    axis.text.y=element_text(size=7),
    axis.title.x=element_blank(),
    axis.title.y=element_text(size=9),
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    legend.position='none'
  ) +
  ylab('')
ggdist
ggsave(here('230214_sorl1Citalopram', 'pubplots', 'distbygrp.pdf'), ggdist, width=65, height=60, units='mm', device=cairo_pdf)

### or separate doses
ggdist2 <- ggplot(did, aes(x=grp_dose, y=dist)) +
  geom_quasirandom() +
  facet_grid(~dose, scales='free_x', space='free_x') +
  theme_minimal()
ggdist2
ggsave(here('230214_sorl1Citalopram', 'plots', 'alldist_dose.pdf'), ggdist2, width=200, height=200, units='mm')

# eventually, will probably delete the 0 µM plot


# statistics --------------------------------------------------------------

# easiest is to do pairwise t tests (with no adjustment for multiple testing) and only look at
# SCR 1 µM vs SORL1 1 µM
# SCR 10 µM vs SORL1 10 µM
pairwise.t.test(did$dist, did$grp_dose, p.adjust.method='none')
# reaches same conclusion as with t-values

# Ns
did %>%
  group_by(grp_dose) %>%
  tally()


# eLife reviews -----------------------------------------------------------
# test whether significant result would survive without the larvae with very high values in SORL1 + 10 uM
# test removing top 3 (above 20)
# test removing top 5 (above 15)

## remove top3 ##

# rows to remove:
did_rm3 <- did[-which(did$grp_dose=='sorl1_10' & did$dist>20),]

# check we removed the right ones
ggplot(did_rm3, aes(x=grp_dose, y=dist, colour=grp_dose)) +
  geom_quasirandom(width=0.3, size=2) +
  stat_summary(aes(group=grp_dose), fun=mean, geom='point', colour='#595E60', shape=3, size=1.2, stroke=0.8, position=position_dodge(dodgeby)) +
  facet_grid(~grp, scales='free_x', space='free_x') +
  scale_colour_manual(values=c('#d9dde1', '#b3bcc3', '#697a87',
                               '#FDEFEE', '#F6B8B1', '#EE7163')) +
  theme_minimal() +
  theme(
    panel.grid.minor.y=element_blank(),
    axis.text.x=element_blank(),
    axis.text.y=element_text(size=7),
    axis.title.x=element_blank(),
    axis.title.y=element_text(size=9),
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    legend.position='none'
  ) +
  ylab('')

pairwise.t.test(did_rm3$dist, did_rm3$grp_dose, p.adjust.method='none') # p = 0.18

## remove top5 ##
# rows to remove:
did_rm5 <- did[-which(did$grp_dose=='sorl1_10' & did$dist>14.5),]

# check we removed the right ones
ggplot(did_rm5, aes(x=grp_dose, y=dist, colour=grp_dose)) +
  geom_quasirandom(width=0.3, size=2) +
  stat_summary(aes(group=grp_dose), fun=mean, geom='point', colour='#595E60', shape=3, size=1.2, stroke=0.8, position=position_dodge(dodgeby)) +
  facet_grid(~grp, scales='free_x', space='free_x') +
  scale_colour_manual(values=c('#d9dde1', '#b3bcc3', '#697a87',
                               '#FDEFEE', '#F6B8B1', '#EE7163')) +
  theme_minimal() +
  theme(
    panel.grid.minor.y=element_blank(),
    axis.text.x=element_blank(),
    axis.text.y=element_text(size=7),
    axis.title.x=element_blank(),
    axis.title.y=element_text(size=9),
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    legend.position='none'
  ) +
  ylab('')

pairwise.t.test(did_rm5$dist, did_rm5$grp_dose, p.adjust.method='none') # p = 0.87
