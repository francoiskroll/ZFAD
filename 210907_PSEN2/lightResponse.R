#####################################################
# ~ ZFAD, eLife reviews: test how psen2 larvae react to the light transitions ~
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

library(waffle)

sem_narm <- function(x) {
  sd(x, na.rm=TRUE)/sqrt(length(x))
}

# import ------------------------------------------------------------------

ff1 <- importRAWs(here('210907_PSEN2/210907_12_RAWs.csv'))
ff2 <- importRAWs(here('210907_PSEN2/210907_13_RAWs.csv'))



# test using sunset -------------------------------------------------------

# from 210907_12_lights.csv
# sunset 6 dpf was frame #2898195
# supposed to be frame just before the light transition
sunset1 <- 2898195

ggDeltaPx(ffpath=here('210907_PSEN2/210907_12_RAWs.csv'),
          well='f96',
          frstart=sunset1-1500,
          frstop=sunset1+100,
          colour='black',
          ymax=100,
          xnsecs=1,
          width=70,
          height=60,
          exportPath=here('210907_PSEN2/pubplots/delete.pdf'))
# f1: big bout starting at 2898188
# f2: big bout starting at 2898187
# f3: big bout starting at 2898187
# f4: big bout starting at 2898188
# >> real sunset must be 2898195 - 9 frames = 2898186
# make it 10 frames to be a bit safer

# take last 1500 frames of day1
sunsetcor1 <- 2898195 - 10
lastmin1 <- ff1[(sunsetcor1-1500) : sunsetcor1, ]

sapply(lastmin1[,4:ncol(lastmin1)], function(fi) {
  sum(fi)
})
# 2 are asleep
# maybe we should look at sunrise instead?


# using SUNRISE -----------------------------------------------------------

# there should be more asleep

### EXP1 - SUNRISE1
# exsecs 65528
exp1sun1 <- 1638209 # from lights.csv

ggDeltaPx(ffpath=here('210907_PSEN2/210907_12_RAWs.csv'),
          well='f20',
          frstart=exp1sun1-10,
          frstop=exp1sun1+100,
          colour='black',
          ymax=100,
          xnsecs=1,
          width=70,
          height=60,
          exportPath=here('210907_PSEN2/pubplots/delete.pdf'))
# tried a few to find clear examples, not as robust as for sunset
# f8: big bout starting at -8
# f9: small bout starting at -9
# f10: small bout starting at -9
# f11: small bout starting at -9
# so should be same shift, 10 frames
exp1sun1 <- exp1sun1 - 10


### EXP1 - SUNRISE2
exp1sun2 <- 3798188 # from lights.csv

ggDeltaPx(ffpath=here('210907_PSEN2/210907_12_RAWs.csv'),
          well='f7',
          frstart=exp1sun2-7,
          frstop=exp1sun2+100,
          colour='black',
          ymax=100,
          xnsecs=1,
          width=70,
          height=60,
          exportPath=here('210907_PSEN2/pubplots/delete.pdf'))
# f1: bout at -6
# f2: bout at -8
# f3: bout at -5
# f4: missed transition
# f5: missed transition
# f7: bout at -7
# so again probably same shift, 10 frames
exp1sun2 <- exp1sun2 - 10


### EXP2 - SUNRISE1
exp2sun1 <- 1638207 # from lights.csv

ggDeltaPx(ffpath=here('210907_PSEN2/210907_13_RAWs.csv'),
          well='f104',
          frstart=exp2sun1-10,
          frstop=exp2sun1+100,
          colour='black',
          ymax=100,
          xnsecs=1,
          width=70,
          height=60,
          exportPath=here('210907_PSEN2/pubplots/delete.pdf'))
# tried a few to find clear examples, not as robust as for sunset
# f97: bout at -9
# f98: bout at -9
# f99: bout at -8
# f102: bout at -3
# f103: bout at -9
# f104: bout at -9
# will put -10 again
exp2sun1 <- exp2sun1 - 10


### EXP2 - SUNRISE2
exp2sun2 <- 3798189 # from lights.csv

ggDeltaPx(ffpath=here('210907_PSEN2/210907_13_RAWs.csv'),
          well='f103',
          frstart=exp2sun2,
          frstop=exp2sun2+100,
          colour='black',
          ymax=100,
          xnsecs=1,
          width=70,
          height=60,
          exportPath=here('210907_PSEN2/pubplots/delete.pdf'))
# tried a few to find clear examples, not as robust as for sunset
# f99: bout at -8
# f102: bout at -7
# f103: bout at +1
# f104: bout at -7
# will put -10 again
exp2sun2 <- exp2sun2 - 10

# >>> we have all sunrise frames we need



# function reactAsleepSunrise ---------------------------------------------

reactAsleepSunrise <- function(ff,
                               sunriseFrame) {
  # reacting in the next second seems like a good definition
  # so we know for sure larvae is reacting to the stimulus, not simply starting the day
  fps <- 25
  
  # take last minute of previous night
  lastmin <- ff[(sunriseFrame-25*60) : sunriseFrame, ]
  
  sumlastmin <- sapply(lastmin[,4:ncol(lastmin)], function(fi) {
    sum(fi)
  })
  # more are asleep, probably more interesting, good mix of both in fact
  
  asleep <- names(sumlastmin)[sumlastmin==0] # larvae that were asleep when lights switched on
  
  # now we take first 25 seconds of day
  # and sum, if >0, then there was a bout and it reacted to the lights switching on
  firstsec <- ff[(sunriseFrame+1) : (sunriseFrame+1+25), ]
  sumfirstsec <- sapply(firstsec[,4:ncol(firstsec)], function(fi) {
    sum(fi)
  })
  
  # create a dataframe with info we have
  sunrea <- data.frame(fis=names(sumlastmin), lastmin=sumlastmin, firstsec=sumfirstsec)
  # sunrea for sunrise reaction
  row.names(sunrea) <- NULL
  
  # add genotype
  # create a tiny dataframe, just to give it to addGenotypefXColnames function
  genodf <- as.data.frame(matrix(nrow=1, ncol=96))
  colnames(genodf) <- names(sumlastmin)
  
  # 
  if(colnames(ff)[4]=='f1') {
    genodf <- addGenotypefXColnames(df=genodf,
                                    genopath=here('210907_PSEN2/210907_12genotype.txt'))
  } else if(colnames(ff)[4]=='f97') {
    genodf <- addGenotypefXColnames(df=genodf,
                                    genopath=here('210907_PSEN2/210907_13genotype.txt'))
  }

  # and now we take those column names for sunrea above
  sunrea <- sunrea %>%
    mutate(fid=colnames(genodf), .after='fis')
  # now split to get the group
  sunrea <- sunrea %>%
    mutate(grp=strNthSplit(fid, '_', 2), .after='fid')
  
  # add asleep column
  # which is simply if lastmin = 0
  # and react column
  # which is simply if firstsec is > 0
  sunrea <- sunrea %>%
    mutate(asleep=(lastmin==0)) %>% 
    mutate(react=(firstsec>0))
  
  # remove excluded
  sunrea <- sunrea %>%
    filter(grp!='excluded')
  
  # return
  return(sunrea)
  
}



# plot --------------------------------------------------------------------

# easiest will be to do 4 small plots and combine them, I think

### to do each waffle plot
ggWaf <- function(sub,
                  exportPath) {
  subco <- sub %>%
    group_by(react) %>%
    tally(name='nfis')
  
  print(subco) # so we know Ns
  
  # should normalise count to 25
  # first calculate proportions
  subco <- subco %>%
    mutate(pro=nfis/sum(nfis))
  # preallocate countnorm column
  subco$countnorm <- NA
  # second, take react FALSE and calculate how many out of 25
  subco[which(!subco$react), 'countnorm'] <- round(subco[which(!subco$react), 'pro'] * 25)
  # calculate left over, it is 25 - what we just wrote
  subco[which(subco$react), 'countnorm'] <- 25 - subco[which(!subco$react), 'countnorm']
  
  # now total will always be 25
  
  # now make sure total is 35 so tiles are always the same size
  # add a blank count to fill
  # subco <- subco %>%
  #   add_row(react='blank', nfis=30-sum(subco$nfis))
  
  #subco$react <- factor(subco$react, levels=c('blank', 'TRUE', 'FALSE'))
  
  # colours:
  # if psen2, dark green / light green
  if(unique(sub$grp)=='scr') {
    colos <- c('#b4bdc4', '#5a6974')
  } else if(unique(sub$grp)=='psen2') {
    colos <- c('#d8e7d2', '#78ac63')
  }
  
  ggwaf <- ggplot(subco, aes(fill=react, values=countnorm)) +
    geom_waffle(n_rows=5, size=0.3, colour='white',
                show.legend=FALSE) +
    scale_fill_manual(name=NULL,
                      values=colos) +
    coord_equal() +
    theme_void()
  ggsave(exportPath, ggwaf, width=25, height=25, units='mm', device=cairo_pdf)
  return(ggwaf)
}



# EXP1, DAY1 --------------------------------------------------------------

# prepare data
sunrea <- reactAsleepSunrise(ff=ff1,
                             sunriseFrame=exp1sun1)

### topleft: ASLEEP / PSEN2
sub <- sunrea %>%
  filter(asleep & grp=='psen2')

ggWaf(sub=sub,
      exportPath=here('210907_PSEN2/pubplots/heat_box12_day1_1.pdf'))

### bottomleft: AWAKE / PSEN2
sub <- sunrea %>%
  filter(!asleep & grp=='psen2')

ggWaf(sub=sub,
      exportPath=here('210907_PSEN2/pubplots/heat_box12_day1_2.pdf'))

### topright: ASLEEP / SCR
sub <- sunrea %>%
  filter(asleep & grp=='scr')

ggWaf(sub=sub,
      exportPath=here('210907_PSEN2/pubplots/heat_box12_day1_3.pdf'))

### bottomright: AWAKE / SCR
sub <- sunrea %>%
  filter(!asleep & grp=='scr')

ggWaf(sub=sub,
      exportPath=here('210907_PSEN2/pubplots/heat_box12_day1_4.pdf'))


# EXP1, DAY2 --------------------------------------------------------------

# prepare data
sunrea <- reactAsleepSunrise(ff=ff1,
                             sunriseFrame=exp1sun2)

### topleft: ASLEEP / PSEN2
sub <- sunrea %>%
  filter(asleep & grp=='psen2')

ggWaf(sub=sub,
      exportPath=here('210907_PSEN2/pubplots/heat_box12_day2_1.pdf'))

### bottomleft: AWAKE / PSEN2
sub <- sunrea %>%
  filter(!asleep & grp=='psen2')

ggWaf(sub=sub,
      exportPath=here('210907_PSEN2/pubplots/heat_box12_day2_2.pdf'))

### topright: ASLEEP / SCR
sub <- sunrea %>%
  filter(asleep & grp=='scr')

ggWaf(sub=sub,
      exportPath=here('210907_PSEN2/pubplots/heat_box12_day2_3.pdf'))

### bottomright: AWAKE / SCR
sub <- sunrea %>%
  filter(!asleep & grp=='scr')

ggWaf(sub=sub,
      exportPath=here('210907_PSEN2/pubplots/heat_box12_day2_4.pdf'))


# EXP2, DAY1 --------------------------------------------------------------

# prepare data
sunrea <- reactAsleepSunrise(ff=ff2,
                             sunriseFrame=exp2sun1)

### topleft: ASLEEP / PSEN2
sub <- sunrea %>%
  filter(asleep & grp=='psen2')

ggWaf(sub=sub,
      exportPath=here('210907_PSEN2/pubplots/heat_box13_day1_1.pdf'))

### bottomleft: AWAKE / PSEN2
sub <- sunrea %>%
  filter(!asleep & grp=='psen2')

ggWaf(sub=sub,
      exportPath=here('210907_PSEN2/pubplots/heat_box13_day1_2.pdf'))

### topright: ASLEEP / SCR
sub <- sunrea %>%
  filter(asleep & grp=='scr')

ggWaf(sub=sub,
      exportPath=here('210907_PSEN2/pubplots/heat_box13_day1_3.pdf'))

### bottomright: AWAKE / SCR
sub <- sunrea %>%
  filter(!asleep & grp=='scr')

ggWaf(sub=sub,
      exportPath=here('210907_PSEN2/pubplots/heat_box13_day1_4.pdf'))



# EXP2, DAY2 --------------------------------------------------------------

# prepare data
sunrea <- reactAsleepSunrise(ff=ff2,
                             sunriseFrame=exp2sun2)

### topleft: ASLEEP / PSEN2
sub <- sunrea %>%
  filter(asleep & grp=='psen2')

ggWaf(sub=sub,
      exportPath=here('210907_PSEN2/pubplots/heat_box13_day2_1.pdf'))

### bottomleft: AWAKE / PSEN2
sub <- sunrea %>%
  filter(!asleep & grp=='psen2')

ggWaf(sub=sub,
      exportPath=here('210907_PSEN2/pubplots/heat_box13_day2_2.pdf'))

### topright: ASLEEP / SCR
sub <- sunrea %>%
  filter(asleep & grp=='scr')

ggWaf(sub=sub,
      exportPath=here('210907_PSEN2/pubplots/heat_box13_day2_3.pdf'))

### bottomright: AWAKE / SCR
sub <- sunrea %>%
  filter(!asleep & grp=='scr')

ggWaf(sub=sub,
      exportPath=here('210907_PSEN2/pubplots/heat_box13_day2_4.pdf'))




# prepare overall summaries -----------------------------------------------

# make a list of "sunreact" dataframes

sunrl <- list(reactAsleepSunrise(ff=ff1,
                                 sunriseFrame=exp1sun1),
              reactAsleepSunrise(ff=ff1,
                                 sunriseFrame=exp1sun2),
              reactAsleepSunrise(ff=ff2,
                                 sunriseFrame=exp2sun1),
              reactAsleepSunrise(ff=ff2,
                                 sunriseFrame=exp2sun2))

# for each, simply count how many reacted & total n
invisible(lapply(sunrl, function(sl) {
  
  cat('\t \t \t \t >>> total n = ', nrow(sl), '\n')
  cat('\t \t \t \t >>> reacted n = ', sum(sl$react), '\n')
  cat('\t \t \t \t >>> ratio reacted =', sum(sl$react)/nrow(sl), '\n')
  # probability to react if you are psen2
  # is proportion of psen2 which reacted
  cat('\t \t \t \t >>> psen2 n = ', length(which(sl$grp=='psen2')), '\n')
  cat('\t \t \t \t >>> ratio of psen2 which reacted =',
      length(which(sl$grp=='psen2' & sl$react)) / length(which(sl$grp=='psen2')),
      '\n')
  # probability to react if you are scr
  # is proportion of scr which reacted
  cat('\t \t \t \t >>> ratio of scr which reacted =',
      length(which(sl$grp=='scr' & sl$react)) / length(which(sl$grp=='scr')),
      '\n')
  
  # proportion asleep
  cat('\t \t \t \t >>> ratio asleep =', sum(sl$asleep)/nrow(sl), '\n')
  
  # proportion asleep which reacted
  asleepreact <- length(which(sl$asleep & sl$react)) / sum(sl$asleep)
  awakereact <- length(which(!sl$asleep & sl$react)) / sum(!sl$asleep)
  cat('\t \t \t \t >>> ratio asleep that reacted =', asleepreact, '\n')
  cat('\t \t \t \t >>> ratio awake that reacted =', awakereact, '\n')
  cat('\t \t \t \t >>> difference =', awakereact-asleepreact, '\n')
  cat('\t \t \t \t >>> foldchange =', awakereact/asleepreact, '\n')
  cat('\n')
  
  # fold change
  # proportion of asleep psen2 larvae which reacted
  # vs.
  # proportion of asleep scr larvae which reacted
  # below, each time is number of genotype larvae which were asleep & reacted / total number of genotype larvae which were asleep
  cat('\t \t \t \t >>> foldchange psen2 vs scr probability react if you are asleep =',
      ( length(which(sl$grp=='psen2' & sl$asleep & sl$react)) / length(which(sl$grp=='psen2' & sl$asleep)) ) / 
        ( length(which(sl$grp=='scr' & sl$asleep & sl$react)) / length(which(sl$grp=='scr' & sl$asleep)) )
      , '\n')
  
  cat('\t \t \t \t >>> foldchange psen2 vs scr probability react if you are awake =',
      ( length(which(sl$grp=='psen2' & !sl$asleep & sl$react)) / length(which(sl$grp=='psen2' & !sl$asleep)) ) / 
        ( length(which(sl$grp=='scr' & !sl$asleep & sl$react)) / length(which(sl$grp=='scr' & !sl$asleep)) )
      , '\n')
  
  ###
  # ns for figure
  # awake psen2
  cat('\t \t \t \t >>> N awake & psen2 = ', length(which(!sl$asleep & sl$grp=='psen2')), '\n')
  # awake scrambled
  cat('\t \t \t \t >>> N awake & scrambled = ', length(which(!sl$asleep & sl$grp=='scr')), '\n')
  # asleep psen2
  cat('\t \t \t \t >>> N asleep & psen2 = ', length(which(sl$asleep & sl$grp=='psen2')), '\n')
  # asleep scrambled
  cat('\t \t \t \t >>> N asleep & scrambled = ', length(which(sl$asleep & sl$grp=='scr')), '\n')
  
  cat('\n \t \t \t \t _________________ \n')
}))
