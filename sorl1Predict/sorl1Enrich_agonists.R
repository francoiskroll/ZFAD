#####################################################
# ~ ZFAD: agonist/antagonist analysis of SORL1 x serotonin ~
#
# Francois Kroll 2023
# francois@kroll.be
#####################################################

# predPharma predicted serotonin for SORL1 f0 mean
# experiment with Citalopram seems to confirm the prediction
# what was the predicted direction?
# i.e. do agonists give a more similar phenotype and antagonist a more different phenotype; or no bias either way?

# I think we can test this by analysing the ranked list of drugs
# 1/ rank drugs in comparison with SORL1 f0 mean
# 2/ split list in two: all positive cos and all negative cos
# 3/ in the positive cos; how many are antagonists/how many are agonists?
# 4/ in the negative cos; how many are antagonists/how many are agonists?
# 5/ Chi-sq test on each side


# packages & functions ----------------------------------------------------

library(here)
library(tictoc)
library(ggplot2)

# excuse the absolute paths!
source('~/Dropbox/predPharma/drawEnrich_v5.R')
source('~/Dropbox/predPharma/legacyFingerprint.R')
source('~/Dropbox/predPharma/gglegacyFingerprint.R')
source('~/Dropbox/predPharma/ggEnrich.R')
source('~/Dropbox/predPharma/paramsFromMid.R')



# calculate fingerprint ---------------------------------------------------

#### f0 mean ###
lfp <- legacyFingerprintMEAN(matPaths=c(here('220531_SORL1', 'legacyMiddur', '220531_14', '220531_14.mat'),
                                        here('220531_SORL1', 'legacyMiddur', '220531_15', '220531_15.mat')),
                             conGrp='scr',
                             treGrp='sorl1',
                             days=c(2,3),
                             nights=c(2,3))

### correct way to proceed
mid1 <- read.csv('~/Dropbox/predPharma_dataDev/220531_14_middur.csv')
lfp1 <- legacyFingerprintMid(mid=mid1,
                     genopath='~/Dropbox/predPharma_dataDev/220531_14genotype.txt',
                     treGrp='sorl1',
                     conGrp='scr',
                     nights=c('night1', 'night2'),
                     days=c('day1', 'day2'))

mid2 <- read.csv('~/Dropbox/predPharma_dataDev/220531_15_middur.csv')
lfp2 <- legacyFingerprintMid(mid=mid2,
                             genopath='~/Dropbox/predPharma_dataDev/220531_15genotype.txt',
                             treGrp='sorl1',
                             conGrp='scr',
                             nights=c('night1', 'night2'),
                             days=c('day1', 'day2'))

### prepare mean
tmp <- data.frame(lfp1$`220531_14_sorl1`, lfp2$`220531_15_sorl1`)
lfp <- lfp1
lfp$`220531_14_sorl1` <- apply(tmp, 1, mean)
colnames(lfp)[4] <- 'zsco'


# rank drugs --------------------------------------------------------------

vdbr <- rankDrugDb(legacyFgp=lfp,
                   dbPath='~/Dropbox/predPharma/drugDb.csv', 
                   metric='cosine')


# plot drugs that interact with serotonin transporters --------------------

# for illustration, ggBarcode of drugs that interact with serotonin transporter

# from Therapeutic Target Database
ggBarcode(vdbr=vdbr,
          namesPath='~/Dropbox/predPharma/compounds.csv',
          annotationPath=here('annotateDrugDb', 'TTDtargets.csv'),
          annotation='TTDtargets',
          testAnnotation='T27812', # human SLC6A4
          minScore=NA,
          exportPath=here('sorl1Predict/plots/sorl1f0mean_TTDslc6a4.pdf'),
          width=500,
          height=80)

# from zebrafish STITCH
ggBarcode(vdbr=vdbr,
          namesPath='~/Dropbox/predPharma/compounds.csv',
          annotationPath=here('annotateDrugDb', 'zebrafishSTITCH.csv'),
          annotation='zebrafishSTITCH',
          testAnnotation='ENSDARP00000081336', # zebrafish slc6a4a
          minScore=900,
          exportPath=here('sorl1Predict/plots/sorl1f0mean_zSTITCHslc6a4a.pdf'),
          width=500,
          height=80)

ggBarcode(vdbr=vdbr,
          namesPath='~/Dropbox/predPharma/compounds.csv',
          annotationPath=here('annotateDrugDb', 'zebrafishSTITCH.csv'),
          annotation='zebrafishSTITCH',
          testAnnotation='ENSDARP00000081315', # zebrafish slc6a4b
          minScore=900,
          exportPath=here('sorl1Predict/plots/sorl1f0mean_zSTITCHslc6a4b.pdf'),
          width=500,
          height=80)

# from human STITCH
ggBarcode(vdbr=vdbr,
          namesPath='~/Dropbox/predPharma/compounds.csv',
          annotationPath=here('annotateDrugDb', 'humanSTITCH.csv'),
          annotation='humanSTITCH',
          testAnnotation='ENSP00000261707', # human SLC6A4
          minScore=900,
          exportPath=here('sorl1Predict/plots/sorl1f0mean_hSTITCHslc6a4a.pdf'),
          width=500,
          height=80)

# there may be too few examples of negative cosine to do the Chi-sq on both sides
# would being less stringent on minScore help?
ggBarcode(vdbr=vdbr,
          namesPath='~/Dropbox/predPharma/compounds.csv',
          annotationPath=here('annotateDrugDb', 'humanSTITCH.csv'),
          annotation='humanSTITCH',
          testAnnotation='ENSP00000261707', # human SLC6A4
          minScore=800,
          exportPath=here('sorl1Predict/plots/sorl1f0mean_hSTITCHslc6a4aminScore800.pdf'),
          width=100,
          height=25)
# yes, many more

# TTD & zebrafish STITCH & human STITCH look fairly similar, which is reassuring
# i.e. drugs interacting with serotonin transporter are concentrated in the first half (so positive correlation > ~ 0.4)
# in zebrafish STITCH, I think slc6a4a and slc6a4b results are exactly the same, i.e. probably exactly the same drugs are annotated

# below takes a while!
# den <- drugEnrichment(vdbr=vdbr,
#                       namesPath='~/Dropbox/predPharma/compounds.csv',
#                       annotationPath=here('annotateDrugDb', 'humanSTITCH.csv'),
#                       annotation='humanSTITCH',
#                       whichRank='rank0',
#                       minScore=800,
#                       minNex=3,
#                       ndraws=10000,
#                       alphaThr=0.05,
#                       maxPval=0.1,
#                       statsExport=NA)


# human STITCH: annotate SLC6A4 inhibition/activation ---------------------

# import human STITCH
hsti <- read.csv(here('annotateDrugDb', 'humanSTITCH.csv'))

# filter on scores
hsti <- hsti %>%
  filter(score >= 800)

# check values of action column
unique(hsti$action) # activation or inhibition or ''
# replace empty strings by NA
hsti[which(hsti$action==''), 'action'] <- NA

# in hsti, there are some annotated interactions where compound_is_acting is FALSE
# should we keep them or no?
# read notes here http://stitch.embl.de/download/README about a_is_acting column (in hsti: compound_is_acting)
# I do not understand e.g. mode = activation / action = activation / compound_is_acting = FALSE
# how can one know it is activation if not sure about whether the compound is active partner?
# asked by email

# will keep for now

# go through compounds in ranked list one by one
# for each, look whether they interact with SLC6A4 in human STITCH
# use CID as query
slcCol <- sapply(vdbr$cid, function(ci) {
  
  ## keep only interactions this compound -- SLC6A4
  sli <- hsti %>%
    filter(cid==ci) %>%
    filter(gene_symbol=='SLC6A4')
  # sli for interactions with SLC6A4
  
  if(nrow(sli)==0) {
    # if no interaction this compound -- SLC6A4, return NA
    return(NA)
    
  } else if(nrow(sli)==1) {
    # if one interaction this compound -- SLC6A4, return the action of this interaction
    # note, this may be NA
    return(as.character(sli$action))
    
  } else if(nrow(sli)>1) {
    # if more than one interaction this compound -- SLC6A4...
    
    # if all interactions have the same action, report that action
    # i.e. possible cases: all NA or all activation or all inhibition
    if(length(unique(sli$action))==1) {
      return(unique(sli$action))
      
      # if only one action is annotated and the other ones are NA, report that action
    } else if(length(unique(sli$action)[!is.na(sli$action)])==1) { # remove NA from unique actions, is there only one left?
      # only case where this does not work is if there is both activation & inhibition annotated
      cat('\t \t \t \t >>> Some NA and one action, returning that action. \n')
      return( unique(sli$action)[!is.na(sli$action)] )
    } else if( length(unique(sli$action)[!is.na(sli$action)])==2 ) {
      cat('\t \t \t \t >>> activation & inhibition both reported, returning both. \n')
      return('both') # we will duplicate rows after in data, if needed
    }
  }
  
})


# add column to ranked list of drugs
nrow(vdbr) == length(slcCol)

vdbr <- vdbr %>%
  mutate(interactSLC=slcCol, .after='cor')

# minScore=900: I think there is 0 activation...
which(slcCol=='activation')
# in a way, correct to say all are inhibition but bit of a strawman if activation was not really "tested"
# this is with minScore=900
# repeating with minScore=700, 'high' confidence on STITCH website

# better

# deal with 'both'
# i.e. both means that for that compound, there is one (or more) 'activation' interaction(s) and one (or more) 'inhibition' interaction(s)
# we will duplicate the rows in vdbr and have one 'activation' and one 'inhibition'
# copy those rows aside
inhAdd <- vdbr[which(vdbr$interactSLC=='both'),]
# change both to inhibition
inhAdd[,'interactSLC'] <- 'inhibition'
# change the original rows to activation
vdbr[which(vdbr$interactSLC=='both'), 'interactSLC'] <- 'activation'
# add the rows we copied
vdbr <- rbind(vdbr, inhAdd)
# order again according to cos
vdbr <- vdbr[order(rev(vdbr$cos)),]
# looks correct


# split ranked list in positive & negative cos ----------------------------

pos <- vdbr[which(vdbr$cos>0),]
neg <- vdbr[which(vdbr$cos<0),]

# in each, how many activation/inhibition?
pot <- pos %>%
  group_by(interactSLC) %>%
  tally()

net <- neg %>%
  group_by(interactSLC) %>%
  tally()

# place in a small dataframe
tal <- as.data.frame(full_join(pot, net, by='interactSLC')) # tal for tally
colnames(tal) <- c('interactSLC', 'posCos', 'negCos')
# replace NA by unannotated
tal[which(is.na(tal$interactSLC)),'interactSLC'] <- 'unannotated'


# small plot

tall <- tal %>%
  pivot_longer(-interactSLC,
               names_to='side',
               values_to='count')

tall$side <- factor(tall$side, levels=c('posCos', 'negCos'))

# remove unannotated
tall <- tall %>%
  filter(interactSLC!='unannotated')

ggplot(tall, aes(x=side, y=count, fill=interactSLC)) +
  geom_col()

# normalise
vdbr <- vdbr %>%
  mutate(posCos=(cos>0), .after='cor')

# make contingency table
cont <- table(vdbr$posCos, vdbr$interactSLC)

chi <- chisq.test(cont, correct=FALSE)

chi$expected
chi$observed
# in positive (TRUE): we observed more inhibition than expected
# in negative (FALSE): we observed less inhibition than expected

# only N=6 activation in negative cos (with minScore=700)
# would be good to have more examples, but I would not want to play around with minScore until having a significant p-value

# still N=6 with minScore=400

# N=7 with minScore=150

# with minScore=150, observed counts is basically = expected

# likely conclusion is no strong bias, but hard to be confident without more examples of activation in negative cosines

# need to think about the assumptions

# repeat analysis with zebrafish STITCH
# but likely we'll have even fewer activation

# could we add serotonin receptors? probably not, as they were not significant in f0mean_zSTITCH?