#####################################################
# ~ ZFAD, eLife reviews: histogram of inactive bout lengths for psen2 ~
#
#
# Francois Kroll 2024
# francois@kroll.be
#####################################################


# packages ----------------------------------------------------------------

library(FramebyFrame)
library(here)
library(tibble)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggbreak)

sem_narm <- function(x) {
  sd(x, na.rm=TRUE)/sqrt(length(x))
}

# import ------------------------------------------------------------------

ff1 <- importRAWs(here('210907_PSEN2/210907_12_RAWs.csv'))
ff2 <- importRAWs(here('210907_PSEN2/210907_13_RAWs.csv'))


# function transitionsToInactiveBoutLengths -------------------------------

transitionsToInactiveBoutLengths <- function(dnai,
                                             minLength=25) {
  
  ### loop fish by fish
  # for each, calculate active bout start/stop
  # above lapply to loop thru days
  # below sapply to loop thru larvae
  dnibl <- lapply(dnai, function(dni) {
    
    # ibl is a list, one slot for each fish
    # each element is a vector of inactive bout lengths
    ibl <- sapply(4:ncol(dni), function(fi) {
      sast <- transitionsToActiveBoutStartStop(dni[,..fi])
      # list of two elements: all start frames, all stop frames
      # transitionsToActiveBoutStartStop is made so that exactly same length
      # so index corresponds to active bout1, 2, 3, etc.
      # to convert to inactive start/stop
      # stop frames of active bouts give the start of each inactive bout
      # start frames of active bouts give the stop of each inactive bout
      # but we skip the first start (now stop) frame, so inactive bout is from stop1 to start2, etc.
      # at the end, the last stop frame gives the start of the last inactive bout, but we do not see its end
      # so instead we finish at the last start frame, so we skip the last stop frame
      # make a new list to be clear it is now inactive bouts
      isast <- list(saf=NA, sof=NA)
      isast$saf <- sast$sof[1 : length(sast$sof)-1] # the stop frames are now start frames, and we skip the last one
      isast$sof <- sast$saf[2 : length(sast$saf)] # the start frames are now stop frames, and we skip the first one
      
      # now to calculate the inactive bout lengths is simply:
      ina <- isast$sof - isast$saf
      
      ### ! we will only keep inactive bouts above 1 sec
      # I think not particularly interesting for the question we want to answer here
      # and they take a huge proportion of the distribution
      ina <- ina[which(ina >= minLength)]
      
      # we return this
      return(ina)
    })
  })

  return(dnibl)
  # dnibl is a list of lists
  # top level is day1, day2
  # bottom level is f1, f2, f3, etc.
  # and each element is vector of inactive bout lengths
  # I would have liked to change the names of bottom list to f1, f2, ... but it just changes the data to the names for some reason
}


### inactiveLengthsToSecBins
inactiveLengthsToFreqBins <- function(dnibl,
                                     breaks) {
  dnfr <- lapply(dnibl1, function(dni) {
    frl <- lapply(dni, function(fi) {
      his <- hist(fi, breaks=brks, plot=FALSE)
      # make a dataframe with the results
      hisd <- data.frame(breaks=his$breaks[2:length(his$breaks)],
                         counts=his$counts)
      # note column breaks is end boundary of each bin
      # i.e. breaks 25 means bin is 0 to 25 frames
      
      # add frequencies of inactive bouts of different lengths
      freq <- hisd$counts / sum(hisd$counts)
      # we will only keep frequencies
      return(freq)
    })
    # cbind all the results from above
    # this will give a dataframe with 96 columns, one per larva
    frd <- as.data.frame(do.call(cbind, frl))
    colnames(frd) <- sprintf('f%i', 1:96)
    # add the breaks as a column
    frd <- frd %>%
      mutate(breaks=brks[2:length(brks)], .before=1)
    return(frd)
  })
  
  return(dnfr)
}

# obtain active bout start/stop -------------------------------------------

# in order, need to do:
# 1- split ff by day/night to get dn
# 2- run framesToActiveInactiveTransitions, returns dnai
# 3- transitionsToActiveBoutStartStop, runs on one column at a time of dnai
# so store results in a list probably
# 4- convert active bout start/stop to inactive bout start/stop
# 5- inactive bout start/stop to inactive bout lengths
# >> will write a function for 3- and 4- and 5-

# 1- splitFramesbyDayNight
dn1 <- splitFramesbyDayNight(tc=ff1,
                            path=here('210907_PSEN2/210907_12_13_PSEN2.xls'))
dn2 <- splitFramesbyDayNight(tc=ff2,
                            path=here('210907_PSEN2/210907_12_13_PSEN2.xls'))

fps1 <- round(averageFramerate(dn1$day1$exsecs))
fps2 <- round(averageFramerate(dn2$day1$exsecs))
# they are both 25 fps...
fps <- 25

# 2- framesToActiveInactiveTransitions
dnai1 <- framesToActiveInactiveTransitions(dn1)
dnai2 <- framesToActiveInactiveTransitions(dn2)

# we are only interested in day1 & day2
dnai1$night0 <- NULL
dnai1$night1 <- NULL
dnai1$night2 <- NULL

dnai2$night0 <- NULL
dnai2$night1 <- NULL
dnai2$night2 <- NULL

# 3-, 4-, 5- transitionsToInactiveBoutLengths
dnibl1 <- transitionsToInactiveBoutLengths(dnai=dnai1, minLength=fps1)
dnibl2 <- transitionsToInactiveBoutLengths(dnai=dnai2, minLength=fps2)
# ! excluding inactive bouts below 1 sec

# for each larva & day, calculate the counts & frequency of inactive bout lengths within a bin

# bin in every second
# i.e. count inactive bout lengths which last 0-1 sec; 1-2 sec; etc.
# longest possible inactive bout is 14 hours so 14 * 60 * 60 * fps frames
fps <- round(averageFramerate(dn1$day1$exsecs))
brks <- seq(0, 14*60*60*fps, 10*fps) # gives 0, 25, 50, 75, etc.
# EDIT: will do 10 sec
# EDIT: trying 20 sec
# EDIT: trying 30 sec
# challenge is to make a readable plot, > 98% of bouts are below 20 sec
# for summary in text: do 60 sec
# was 10 sec for final figure

fibl1 <- inactiveLengthsToFreqBins(dnibl=dnibl1,
                                   breaks=brks)

fibl2 <- inactiveLengthsToFreqBins(dnibl=dnibl2,
                                   breaks=brks)
# each fibl is a list, day1 & day2
# then dataframe rows = bins, columns = fish, cells = frequency of inactive bouts of that length

# we need to calculate mean ± sd for each genotype & break

### box1
fiblsum1 <- lapply(fibl1, function(dni) {
  
  # first add genotypes to column names
  dni <- addGenotypefXColnames(df=dni,
                               genopath=here('210907_PSEN2/210907_12genotype.txt'))
  
  # I think easiest is then to pivot_long and summarise_at
  # if not slow
  # pivot_long
  dnil <- dni %>%
    pivot_longer(-breaks,
                 names_to='fid',
                 values_to='freq')
  # now split fis to get genotype column
  dnil <- dnil %>%
    mutate(fnum=strNthSplit(dnil$fid, '_', 1), .after='fid') %>%
    mutate(grp=strNthSplit(dnil$fid, '_', 2), .after='fid')
  
  # now summarise in bins
  dnisum <- dnil %>%
    group_by(breaks, grp) %>%
    summarise_at(vars(freq),
                 list(
                   mean= ~ mean(., na.rm=TRUE),
                   sd= ~ sd(., na.rm=TRUE),
                   sem= ~ sem_narm(.),
                   median= ~ median(., na.rm=TRUE),
                   n= ~ length(.)
                 ))
})


### box2
fiblsum2 <- lapply(fibl2, function(dni) {
  
  # first add genotypes to column names
  dni <- addGenotypefXColnames(df=dni,
                               genopath=here('210907_PSEN2/210907_13genotype.txt')) # ! different genotype
  
  # I think easiest is then to pivot_long and summarise_at
  # if not slow
  # pivot_long
  dnil <- dni %>%
    pivot_longer(-breaks,
                 names_to='fid',
                 values_to='freq')
  # now split fis to get genotype column
  dnil <- dnil %>%
    mutate(fnum=strNthSplit(dnil$fid, '_', 1), .after='fid') %>%
    mutate(grp=strNthSplit(dnil$fid, '_', 2), .after='fid')
  
  # now summarise in bins
  dnisum <- dnil %>%
    group_by(breaks, grp) %>%
    summarise_at(vars(freq),
                 list(
                   mean= ~ mean(., na.rm=TRUE),
                   sd= ~ sd(., na.rm=TRUE),
                   sem= ~ sem_narm(.),
                   median= ~ median(., na.rm=TRUE),
                   n= ~ length(.)
                 ))
})



# plot --------------------------------------------------------------------

# will need to do four times (2 experiments, 2 days each)
# so need a function

ggInaLenDist <- function(fiblsum_dn,
                         legendOrNo=TRUE,
                         exportPath) {
  fisum <- fiblsum_dn %>%
    filter(grp!='excluded')
  
  sumleft <- sapply(1:nrow(fisum), function(i) {
    sum(fisum$mean[i:nrow(fisum)])
  })
  # fisum <- fisum[1:which(sumleft==0)[1],]
  
  # ! only keeping below 3 min
  fisum <- fisum[which(fisum$breaks<=fps*60*3),]
  
  # add column in seconds
  fisum <- fisum %>%
    mutate(breaksec=breaks/fps, .after=breaks)
  
  fisum$breaks <- as.factor(fisum$breaks)
  fisum$breaksec <- as.factor(fisum$breaksec)
  
  fisum$grp <- factor(fisum$grp, levels=c('scr', 'psen2'))
  
  ggInaDist <- ggplot(fisum, aes(x=breaksec, y=mean, fill=grp)) + 
    geom_col(position=position_dodge(0.7), width=0.9, colour='black') +
    geom_errorbar(aes(ymin = mean-sem, ymax=mean+sem), 
                  position=position_dodge(0.7), width=0) +
    ylim(0,1) +
    #scale_y_continuous(breaks=c(0, 0.01, 0.02, 0.03, 0.90, 0.95, 1.0)) +
    scale_y_break(breaks=c(0.035, 0.90), scale=0.05) +
    scale_fill_manual(values=c('#b3bcc3', '#bbd4ae')) +
    theme_minimal() +
    {if(!legendOrNo) theme(legend.position='none')} +
    theme(
      panel.grid.minor.y=element_blank(),
      axis.title.x=element_blank(),
      axis.title.y=element_blank(),
      axis.text.y=element_blank()
    )
  
  ggsave(exportPath, ggInaDist, width=150, height=200, units='mm', device=cairo_pdf)
  return(ggInaDist)
  
}



# plot for each experiment/day --------------------------------------------

# in average, 43% of inactive bouts are < 1 sec

## box12
ggInaLenDist(fiblsum_dn=fiblsum1$day1,
             legendOrNo=FALSE,
             exportPath=here('210907_PSEN2/pubplots/210913_12_day1_inadist.pdf'))

ggInaLenDist(fiblsum_dn=fiblsum1$day2,
             legendOrNo=FALSE,
             exportPath=here('210907_PSEN2/pubplots/210913_12_day2_inadist.pdf'))

## box13
ggInaLenDist(fiblsum_dn=fiblsum2$day1,
             legendOrNo=FALSE,
             exportPath=here('210907_PSEN2/pubplots/210913_13_day1_inadist.pdf'))

ggInaLenDist(fiblsum_dn=fiblsum2$day2,
             legendOrNo=FALSE,
             exportPath=here('210907_PSEN2/pubplots/210913_13_day2_inadist.pdf'))


# ratio difference between genotypes for each bin -------------------------

### function to plot fold change

ggFoldChange <- function(fiblsum_dn,
                         exportPath) {
  ## pivot wider
  fisum <- fiblsum_dn %>%
    filter(grp!='excluded')
  
  # add column in seconds
  fisum <- fisum %>%
    mutate(breaksec=breaks/fps, .after=breaks)
  
  # ! only keeping below 3 min
  fisum <- fisum[which(fisum$breaks<=fps*60*3),]
  
  fiw <- fisum %>%
    pivot_wider(id_cols=breaksec,
                names_from='grp',
                values_from='mean')
  
  fiw <- fiw %>%
    mutate(ratio=psen2/scr) %>%
    mutate(sleep=breaksec>60, .after='breaksec')
  
  fiw$breaksec <- factor(fiw$breaksec)
  
  print(fiw)
  
  ## plot
  ggFold <- ggplot(fiw, aes(x=breaksec, y=ratio, fill=sleep)) +
    geom_point(size=1.8, pch=21, colour='black', stroke=0.4) +
    geom_hline(yintercept=1, linetype=2) +
    scale_fill_manual(values=c('white', '#a3bbdb')) +
    theme_minimal() +
    theme(
      panel.grid.minor.y=element_blank(),
      axis.title.x=element_blank(),
      axis.title.y=element_blank(),
      axis.text.x=element_text(size=7, angle=90, hjust=1, vjust=0.5, margin=margin(t=-1, r=0, b=0, l=0)),
      axis.text.y=element_text(size=7, margin=margin(t=0, r=-1, b=0, l=0)),
      legend.position='none'
    ) +
    scale_y_continuous(breaks=c(0, 1, 5, 10)) +
    coord_cartesian(ylim=c(0, 10))
  ggsave(exportPath, ggFold, width=50, height=70, units='mm', device=cairo_pdf)
  return(ggFold)
}

ggFoldChange(fiblsum_dn=fiblsum1$day1,
             exportPath=here('210907_PSEN2/pubplots/210913_12_day1_fold.pdf'))

ggFoldChange(fiblsum_dn=fiblsum1$day2,
             exportPath=here('210907_PSEN2/pubplots/210913_12_day2_fold.pdf'))

ggFoldChange(fiblsum_dn=fiblsum2$day1,
             exportPath=here('210907_PSEN2/pubplots/210913_13_day1_fold.pdf'))

ggFoldChange(fiblsum_dn=fiblsum2$day2,
             exportPath=here('210907_PSEN2/pubplots/210913_13_day2_fold.pdf'))
