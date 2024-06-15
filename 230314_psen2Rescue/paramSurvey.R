#####################################################
# ~ ZFAD: paramSurvey functions ~
#
# functions for script 230314_psen2Rescue_paramSurvey
#
# Francois Kroll 2023
# francois@kroll.be
#####################################################



# function ggParamSurvey(...) ---------------------------------------------

ggParamSurvey <- function(fgcl,
                          lmePath,
                          conGrp,
                          treGrp,
                          ytextOrNo=TRUE,
                          exportPath,
                          width=60,
                          height=105) {
  
  
  # import LME statistics ---------------------------------------------------
  
  # import LMEreport
  lme <- read.csv(lmePath)
  
  # add uparam column
  lme <- lme %>%
    mutate(uparam=paste(daynight, parameter, sep='_'), .before=1) %>%
    # and keep only groups we care about
    filter(beingGroup %in% c(conGrp, treGrp))
    
  
  # pivot wider so that beingGroup are columns
  lmw <- lme %>%
    pivot_wider(id_cols=uparam,
                names_from=beingGroup,
                values_from=posthocpvalsymbol)
  
  
  # pool boxes by doing mean of fingerprints --------------------------------
  
  fgd <- fgcl %>%
    filter(grp %in% c(conGrp, treGrp))
  
  
  # to pool boxes, we do, for each uparam, mean of the two means and mean of the two SEMs
  fgpo <- fgd %>%
    group_by(grp, uparam) %>%
    summarise_at(vars(mean, sem),
                 list(
                   mean=~mean(.)
                 ))
  # fgpo for pool
  
  
  # add parameter categories ------------------------------------------------
  
  # now add one column of labels for the drug treatment
  lmw <- lmw %>%
    mutate(rescue=labelCat(.[,conGrp], .[,treGrp]))
  
  # best is now to join lmw
  fgpo <- left_join(fgpo, lmw[,c('uparam', 'rescue')], by='uparam')
  
  # get the parameters in the same order as in fingerprint (from top to bottom here)
  # alluparams is recorded in package FramebyFrame
  # but I left compressibility in alluparams by accident, keep only parameters we actually have in the data
  alluparams2 <- alluparams[which(alluparams %in% lmw$uparam)]
  fgpo$uparam <- factor(fgpo$uparam, levels=rev(alluparams2))
  
  # set order of groups
  fgpo$grp <- factor(fgpo$grp, levels=c(conGrp, treGrp))
  # order the rows so changes the plotting order, we want drug on top of untreated (DMSO grey)
  fgpo <- fgpo[order(fgpo$grp),]
  
  
  # colour only the drug treatment
  # we can switch DMSO to something else
  fgpo[which(fgpo$grp==conGrp), 'rescue'] <- 'dmso'
  
  fgpo$rescue <- factor(fgpo$rescue,
                        levels=c('notSideEffect', 'sideEffect', 'missedRescue', 'rescue', 'dmso'))
  
  # ! there may not be every category represented so set the colours, then take those that are needed
  colod <- data.frame(cat=c('notSideEffect', 'sideEffect', 'missedRescue', 'rescue', 'dmso'),
                      colo=c('#cacacc', '#fcb505', '#EE7163', '#78ac63', '#697a87'))
  colod <- colod %>%
    filter(cat %in% unique(fgpo$rescue))
  labcolours <- colod$colo
  
  
  ### add colour dots for each label
  
  # will prepare a column to use as geom_point
  # where X is the position of the dot
  xstart <- -1.5
  xstop <- 1.5
  
  fgpo <- fgpo %>%
    mutate(labeldot=xstart)
  # but ! should only have one per uparam, so turn to NA the SCR + DMSO ones
  fgpo[which(fgpo$grp==conGrp), 'labeldot'] <- NA
  
  # # for the linerange, we want black for drug / grey for DMSO
  # # add column
  # fgpo <- fgpo %>%
  #   mutate(isDMSO=(grp=='PSEN2.DMSO'))
  
  ####
  
  #library(ggnewscale)
  
  # clean-up Y axis so give simple parameter names
  # could use function param2Title from FramebyFrame but still gives relatively long names
  # probably need to do manually
  # then will add night background in Illustrator
  # rev() because wrote it while looking at plot top to bottom
  yparam <- rev(c('% time active', 'total px', 'slope', 'transition delta', 'fractal dimension',
                  'active bout length', 'active bout mean', 'active bout sum', 'active bout std', 'active bout min', 'active bout max', 'active bout num',
                  'sleep total', 'sleep bout num', 'sleep bout length',
                  
                  '% time active', 'total px', 'startle', 'slope', 'transition delta', 'fractal dimension',
                  'active bout length', 'active bout mean', 'active bout sum', 'active bout std', 'active bout min', 'active bout max', 'active bout num',
                  'sleep total', 'sleep bout num', 'sleep bout length', 'sleep latency'))
  
  # be careful, if order of parameter changes, Y axis will not as putting labels manually
  # checked carefully 07/04/2023
  # check again by commenting out scale_y_discrete
  
  dodg <- 0.6
  
  ggpasurv <- ggplot(data=fgpo, aes(x=mean_mean, y=uparam, fill=rescue)) +
    
    geom_col(position=position_dodge(dodg), width=0.8) +
    scale_fill_manual(values=labcolours) +
    
    geom_linerange(aes(xmin=mean_mean - sem_mean, xmax=mean_mean + sem_mean, colour=rescue),
                   position=position_dodge(dodg)) +
    scale_colour_manual(values=labcolours) +
    
    #new_scale_colour() +
    
    geom_point(aes(x=labeldot, colour=rescue), size=1.7) +
    #scale_colour_manual(values=labcolours) +
    scale_y_discrete(labels=yparam) +
    theme_minimal() +
    theme(
      axis.title.x=element_text(size=9),
      axis.text.x=element_text(size=7),
      axis.text.y=element_text(size=5),
      axis.title.y=element_blank(),
      legend.position='none',
      legend.title=element_blank(),
      panel.grid.minor.x=element_blank()
    ) +
    
    {if(!ytextOrNo) theme(axis.text.y=element_blank())} +
    
    coord_cartesian(xlim=c(xstart, xstop)) +
    xlab('z-score\nfrom scrambled + DMSO')
  
  ggsave(exportPath, width=width, height=height, units='mm', device=cairo_pdf)
  
  return(ggpasurv)
  
}


# function labelCat(...) --------------------------------------------------

# small function to label parameters with categories based on significance in LME
# need to be given column of untreated pvalAsterisk & column of treated pvalAsterisk
# expect both KO, i.e. untreated is e.g. KO + DMSO and treated is KO + drug
# but statistics should be in reference to WT or SCR
# so can assess rescue

# there are four categories of parameters
# 1/ "not a side effect" -- parameter was not significant in PSEN2 alone and drug does not affect it
# 2/ "side effect" -- parameter was not significant in PSEN2 alone but drug moved it

# 3/ "missed rescue" -- parameter was significant in PSEN2 alone and drug did not improve it
# 4/ "rescue" -- parameter was significant in PSEN2 alone and drug brought it to ns.

labelCat <- function(colControl, colDrug) {
  
  # make sure we are receiving simple vectors
  colControl <- as.vector(unlist(colControl))
  colDrug <- as.vector(unlist(colDrug))
  
  cats <- sapply(1:length(colControl), function(i) {
    
    # category1 -- not a side effect
    # this is a parameter which was not modified in KO and drug did not move it either
    if(colControl[i]=='ns' & colDrug[i]=='ns'){
      return('notSideEffect') 
      
      # category2 -- side effect
      # this is a parameter which was not modified in KO and drug did move it
    } else if (colControl[i]=='ns' & colDrug[i]!='ns') {
      return('sideEffect')
      
      # category3 -- missed rescue
      # this is a parameter which was modified in KO and drug did not change that (also includes if drug made it worse)
    } else if (colControl[i]!='ns' & colDrug[i]!='ns') {
      return('missedRescue')
      
      # category4 -- rescue
      # this is a parameter which was modified in KO and drug rescued it, i.e. made it non-significant
    } else if (colControl[i]!='ns' & colDrug[i]=='ns') {
      return('rescue')
    }
    
  })
  
  return(cats)
  
}