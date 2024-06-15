#####################################################
# ~ ZFAD: small function to calculate effect size in % from parameter table ~
#
# Francois Kroll 2023
# francois@kroll.be
#####################################################


### EXAMPLE USAGE
# percEffectSize(pa=c(here('210907_PSEN2', 'bhvparams', 'sleepNumNaps_210907_13.csv'),
#                     here('210907_PSEN2', 'bhvparams', 'sleepNumNaps_210907_12.csv')),
#                lmePath=here('LMEreports/psen2_f0_LME'),
#                grporder=c('psen2', 'scr'),
#                dayornight='day')

library(dplyr)
library(tidyr)

percEffectSize <-  function(pa,
                            lmePath,
                            grporder,
                            dayornight) {
  
  ## currently support only looking up a single group
  # so grporder must be n=2
  if(length(grporder)!=2) stop('\t \t \t \t Error: grporder must be length 2 \n')
  
  par <- paramReadPivot(pa=pa,
                        grporder=grporder,
                        skipNight0=TRUE)
  
  ## calculate average of reference group for each day or night
  refgrp <- grporder[length(grporder)]
  cat('\t \t \t \t >>> Assuming', refgrp, 'is the reference group \n')
  
  parm <- par %>%
    filter(daynight==dayornight) %>%
    filter(grp==refgrp) %>%
    group_by(date_box_win, grp) %>%
    summarise_at(vars(param),
                 list(
                   mean= ~ mean(., na.rm=TRUE),
                   sd= ~ sd(., na.rm=TRUE),
                   median= ~ median(., na.rm=TRUE),
                   n= ~ length(.)
                 ))
  
  cat('\t \t \t \t >>> Averages for reference group:\n')
  print(parm)
  
  ## get the effect size from the LME report
  lme <- read.csv(lmePath)
  # which parameter are we looking at?
  param <- strNthSplit(afterLastSlash(pa[1]), '_', 1)
  
  # ! LME may have other groups than the two listed in grporder
  # so we need to also look-up beingGroup is first group of grporder
  slo <- lme[which(lme$parameter==param & lme$daynight==dayornight & lme$beingGroup==grporder[1]), 'LMEslope']
  
  cat('\t \t \t \t >>> LME slope:', slo,'\n')
  
  # now add slope to day/night averages
  # this will predict the treatment data, essentially
  trea <- parm$mean + slo
  
  cat('\t \t \t \t >>> Predicted treatment means:', trea,'\n')
  
  # calculate % increase or decrease
  # we calculate effect size for each experiment and each day or night, then average them
  # ! formula different based on direction of effect
  if(slo > 0) {
    peref <- 100 * (mean(-1+trea/parm$mean))
  } else if(slo < 0) {
    peref <- -100 * mean(1-(trea/parm$mean))
  }
  
  cat('\t \t \t \t \t >>> Effect size:', peref,'% \n')
  
  # also calculate as fold change
  # which is simply treatment means / reference means
  # then means of all fold change
  folef <- trea/parm$mean
  cat('\t \t \t \t \t >>> Fold changes:', folef,' \n')
  cat('\t \t \t \t \t >>> Average fold change:', mean(folef),' \n')
  
  invisible(peref)
  
}

