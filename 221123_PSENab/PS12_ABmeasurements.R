###################################################
# ~~~ ZFAD: AB measurements on PSEN1 / PSEN2 F0 knockouts ~~~

# Francois Kroll 2023
# francois@kroll.be
###################################################

# samples were processed on 12/12/2022

# packages ----------------------------------------------------------------

library(here)
library(ggplot2)
library(dplyr)
library(tibble)
library(stringr)
library(ggbeeswarm)



# functions ---------------------------------------------------------------

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


# import ------------------------------------------------------------------

ab <- read.csv(here('221123_PSENab', 'ABmeasurements', '221212_PS1PS2_AB.csv'))

# remove comments
ab <- ab[-(which(ab[,1]=='Comments'):nrow(ab)),]
# it creates empty rows
ab <- ab[-which(ab$PlateName==''),]


# some operations ---------------------------------------------------------

# create composite column sample_assay to have unique ID
# assay is which AB species was tested: AB40 or AB42 or AB38
ab <- ab %>%
  mutate(id=paste(Sample, Assay, Replicate, sep='_'), .before=1)

# we can use the diluted samples as additional technical replicates by multiplying their concentrations by two
ab <- ab %>%
  mutate(CalcConcentration2=CalcConcentration, .after='CalcConcentration')

ab[str_detect(ab$Sample, 'dil'), 'CalcConcentration2'] <- ab[str_detect(ab$Sample, 'dil'), 'CalcConcentration'] * 2

# CalcConcentration is in pg/mL
# buffer volumes are constant (source: Guliz) so one way to normalise would be to divide by number of animals
# concentrations in pg/mL are calculated from calibration curve
# based on protocol https://www.mesoscale.com/~/media/files/product%20inserts/abeta%20peptide%20panel%20-1.pdf
# 25 µL of each sample is added to plate
# I think original buffer added was 100 µL (source: Guliz)
# so assuming no added water in samples (which is inaccurate, but hoping it diluted each sample by roughly the same factor)
# and assuming 100% recovery from the sample (i.e. that 100% of the AB escapes from the mashed tissue into the buffer)
# then total mass of each AB species should be CalcConcentration divided by 10
# because we go from mass in one mL (1000 µL) to mass in 100 µL, which we assume was total mass of AB in all the animals pooled
ab <- ab %>%
  mutate(totalMass=CalcConcentration2/10)

# now that we have total mass, we can divide by number of animals
# to get average mass of each AB species in one animal
ab <- ab %>%
  mutate(massPerFish=totalMass/nfish)


# see emails with Mesoscale 17/01/2023
# my understanding is limit of detection* is simply concentration of standard8 (S008)
# which is calculated by them and included in the data, see Guliz_221212_all_data.csv, column concentration

# * note, this is not how they call it,
# but I do not understand the logic of defining limit of detection as some threshold under which you can still calculate a concentration
# so my definition of limit of detection here is simply: what is the minimum concentration under which Mesoscale returns NA
# I will not worry about the other thresholds they set

ab$LOD <- NA
ab[which(ab$Assay=='AB38'), 'LOD'] <- 28.80859375 # *
ab[which(ab$Assay=='AB40'), 'LOD'] <- 39.79492188 # *
ab[which(ab$Assay=='AB42'), 'LOD'] <- 4.565429688 # *
# * copied from Guliz_221212_all_data, column concentration, row S008
# in pg/mL

# how do we convert those into massPerFish though?
# calibrators are prepared in the same way so I think we can also convert those values into totalMass
ab$LOD_totalMass <- ab$LOD/10

# and to convert into massPerFish?
# we could divide by nfish to have a rough estimate
ab$LOD_massPerFish <- ab$LOD_totalMass/ab$nfish

# to draw a line on plot, might do the highest LOD for each assay
lods <- ab %>%
  group_by(Assay) %>%
  summarise_at(vars(LOD_massPerFish),
               list(
                 max= ~ max(., na.rm=TRUE)
               ))
lod <- lods$max
names(lod) <- lods$Assay


# in plot, we will want to pool technical replicates together
# to do so, create a new column sample2
# it comes down to simply removing the dilution information, if any
ab <- ab %>%
  mutate(sample2=strNthSplit(Sample, split='_', nth=1), .after='Sample')


# plot --------------------------------------------------------------------

# order conditions
ab$sample2 <- factor(ab$sample2, levels=c('uninj', 'SCR4x', 'PSEN1clutch1', 'PSEN1clutch2', 'PSEN2', 'uninjlarvae', 'PSEN12doublelarvae'))

###

ab38 <- ab[which(ab$Assay=='AB38'),]

gg38 <- ggplot(ab38, aes(x=sample2, y=massPerFish)) +
  geom_point() + # geom_quasirandom throws an error because all points are NA
  geom_rect(xmin=0.95, xmax=1.05, ymin=0, ymax=lod['AB38'], fill='#595E60') +
  geom_rect(xmin=1.95, xmax=2.05, ymin=0, ymax=lod['AB38'], fill='#595E60') +
  geom_rect(xmin=2.95, xmax=3.05, ymin=0, ymax=lod['AB38'], fill='#595E60') +
  geom_rect(xmin=3.95, xmax=4.05, ymin=0, ymax=lod['AB38'], fill='#595E60') +
  geom_rect(xmin=4.95, xmax=5.05, ymin=0, ymax=lod['AB38'], fill='#595E60') +
  geom_rect(xmin=5.95, xmax=6.05, ymin=0, ymax=lod['AB38'], fill='#595E60') +
  geom_rect(xmin=6.95, xmax=7.05, ymin=0, ymax=lod['AB38'], fill='#595E60') +
  geom_hline(yintercept=lod['AB38']) +
  theme_minimal() +
  theme(
    axis.text.x=element_text(angle=45, hjust=1),
    axis.title.x=element_blank(),
  ) +
  ylab('AB38 per juvenile (pg)') + 
  coord_cartesian(ylim=c(0, 6)) +
  scale_x_discrete(labels=c('uninjected', 'non-targeting', 'psen1 F0', 'psen1 F0', 'psen2', 'uninjected larvae', 'psen1/2 F0 larvae'))
gg38

ggsave(here('figures', 'ab38.pdf'), gg38, width=80, height=70, unit='mm')

###

ab40 <- ab[which(ab$Assay=='AB40'),]

gg40 <- ggplot(ab40, aes(x=sample2, y=massPerFish)) +
  geom_quasirandom(width=0.1) +
  geom_rect(xmin=2.95, xmax=3.05, ymin=0, ymax=lod['AB40'], fill='#595E60') +
  geom_rect(xmin=3.95, xmax=4.05, ymin=0, ymax=lod['AB40'], fill='#595E60') +
  geom_rect(xmin=5.95, xmax=6.05, ymin=0, ymax=lod['AB40'], fill='#595E60') +
  geom_rect(xmin=6.95, xmax=7.05, ymin=0, ymax=lod['AB40'], fill='#595E60') +
  geom_hline(yintercept=lod['AB40']) +
  theme_minimal() +
  theme(
    axis.text.x=element_text(angle=45, hjust=1),
    axis.title.x=element_blank(),
  ) +
  ylab('AB40 per juvenile (pg)') + 
  coord_cartesian(ylim=c(0, 6)) +
  scale_x_discrete(labels=c('uninjected', 'non-targeting', 'psen1 F0', 'psen1 F0', 'psen2', 'uninjected larvae', 'psen1/2 F0 larvae'))
gg40

ggsave(here('figures', 'ab40.pdf'), gg40, width=80, height=70, unit='mm')

##
# version without larvae

ab40 <- ab %>%
  filter(Assay=='AB40') %>%
  filter(! sample2 %in% c('uninjlarvae', 'PSEN12doublelarvae'))

gg40juv <- ggplot(ab40, aes(x=sample2, y=massPerFish)) +
  geom_quasirandom(width=0.1) +
  geom_rect(xmin=2.95, xmax=3.05, ymin=0, ymax=lod['AB40'], fill='#595E60') +
  geom_rect(xmin=3.95, xmax=4.05, ymin=0, ymax=lod['AB40'], fill='#595E60') +
  geom_rect(xmin=5.95, xmax=6.05, ymin=0, ymax=lod['AB40'], fill='#595E60') +
  geom_rect(xmin=6.95, xmax=7.05, ymin=0, ymax=lod['AB40'], fill='#595E60') +
  geom_hline(yintercept=lod['AB40']) +
  theme_minimal() +
  theme(
    axis.text.x=element_text(angle=45, hjust=1),
    axis.title.x=element_blank(),
  ) +
  ylab('AB40 per juvenile (pg)') + 
  coord_cartesian(ylim=c(0, 6)) +
  scale_x_discrete(labels=c('uninjected', 'non-targeting', 'psen1 F0', 'psen1 F0', 'psen2 F0'))
gg40juv
ggsave(here('figures', 'ab40juv.pdf'), gg40juv, width=80, height=80, unit='mm')

###

ab42 <- ab[which(ab$Assay=='AB42'),]

gg42 <- ggplot(ab42, aes(x=sample2, y=massPerFish)) +
  geom_quasirandom(width=0.1) +
  geom_rect(xmin=2.95, xmax=3.05, ymin=0, ymax=lod['AB42'], fill='#595E60') +
  geom_rect(xmin=3.95, xmax=4.05, ymin=0, ymax=lod['AB42'], fill='#595E60') +
  geom_rect(xmin=5.95, xmax=6.05, ymin=0, ymax=lod['AB42'], fill='#595E60') +
  geom_rect(xmin=6.95, xmax=7.05, ymin=0, ymax=lod['AB42'], fill='#595E60') +
  geom_hline(yintercept=lod['AB42']) +
  theme_minimal() +
  theme(
    axis.text.x=element_text(angle=45, hjust=1),
    axis.title.x=element_blank(),
  ) +
  ylab('AB42 per juvenile (pg)') + 
  coord_cartesian(ylim=c(0, 0.5)) +
  scale_x_discrete(labels=c('uninjected', 'non-targeting', 'psen1 F0', 'psen1 F0', 'psen2', 'uninjected larvae', 'psen1/2 F0 larvae'))
gg42

ggsave(here('figures', 'ab42.pdf'), gg42, width=80, height=80, unit='mm')


# version without larvae
ab42 <- ab %>%
  filter(Assay=='AB42') %>%
  filter(! sample2 %in% c('uninjlarvae', 'PSEN12doublelarvae'))

gg42 <- ggplot(ab42, aes(x=sample2, y=massPerFish)) +
  geom_quasirandom(width=0.1) +
  geom_rect(xmin=2.95, xmax=3.05, ymin=0, ymax=lod['AB42'], fill='#595E60') +
  geom_rect(xmin=3.95, xmax=4.05, ymin=0, ymax=lod['AB42'], fill='#595E60') +
  geom_rect(xmin=5.95, xmax=6.05, ymin=0, ymax=lod['AB42'], fill='#595E60') +
  geom_rect(xmin=6.95, xmax=7.05, ymin=0, ymax=lod['AB42'], fill='#595E60') +
  geom_hline(yintercept=lod['AB42']) +
  theme_minimal() +
  theme(
    axis.text.x=element_text(angle=45, hjust=1),
    axis.title.x=element_blank(),
  ) +
  ylab('AB42 per juvenile (pg)') + 
  coord_cartesian(ylim=c(0, 0.5)) +
  scale_x_discrete(labels=c('uninjected', 'non-targeting', 'psen1 F0', 'psen1 F0', 'psen2 F0'))
gg42
ggsave(here('figures', 'ab42juv.pdf'), gg42, width=80, height=80, unit='mm')



# note 19/04/2023 ---------------------------------------------------------

# from Bradford: seems like very low recovery
# so I do not think analysis above is OK as pg / animal assumes 100% recovery
# while we probably have < 20% recovery from sample
# better to do normalisation with Bradford


# normalise by protein concentration from Bradford ------------------------

brad <- read.csv(here('221123_PSENab', 'ABmeasurements', 'bradford_results.csv'))

# change column names so matches correctly when joining
colnames(brad) <- c('sample2', 'bradfordConc')

# add concentrations to AB measurements
ab <- left_join(ab, brad, by='sample2')

# Bradford concentration (column bradfordConc) are in µg/mL
# Mesoscale concentrations (column CalcConcentration) are in pg/mL

# so we can simply divide Mesoscale / Bradford, will give pg AB / µg of total protein
# low recovery (see bradford.R) should not be a big issue as same buffer for Mesoscale/Bradford
# so say recovery was 100%, we would simply multiply numerator (AB concentration) and denominator (total protein) by same ratio
ab <- ab %>%
  mutate(pgPerugProt=CalcConcentration2/bradfordConc)

# paper e.g. https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2663406/ uses ng/g so will convert to that
# slightly more intuitive unit, perhaps
# it is * 1,000,000, / 1,000, so * 1,000
ab <- ab %>%
  mutate(ngPergProt=pgPerugProt*1000)


### sanity check for Bradford, check if it correlates roughly with number of fish
bracheck <- ab[!duplicated(ab$sample2),]
# remove larvae
bracheck <- bracheck %>%
  filter(! sample2 %in% c('PSEN12doublelarvae', 'uninjlarvae'))

ggplot(bracheck, aes(x=nfish, y=bradfordConc)) +
  geom_point()
# too many points to really say, but seems reasonable
cor(bracheck$nfish, bracheck$bradfordConc) # r=0.85


# deal with LODs ----------------------------------------------------------

# e.g. AB40, lowest concentration we could have detected is 39.79 pg/mL

# what is the average sample concentration?
# do not take the diluted samples
brad2 <- brad %>%
  filter(!str_detect(sample2, '_dil15'))
avconc <- mean(brad2$bradfordConc)

# so if sample was AB40 = 39.79 pg/mL and total protein concentration was 3788 µg/mL
# then would be 0.01050422 pg/µg = 10.50422 ng/g

# so, for each sample, we can tell what the limit of detection is
# calculate each sample's limit of detection
# formula is (LOD / total protein concentration) * 1000
# LODs are already in data
# add column sampleLOD
ab <- ab %>%
  mutate(sampleLOD_ngperg=1000*LOD/bradfordConc)

write.csv(ab, here('221123_PSENab', '221212_PS1PS2_AB_processed.csv'))


# ZFAD publication plots --------------------------------------------------

dir.create(here('221123_PSENab', 'pubplots'), showWarnings=FALSE)

### AB38 is all below detection, so forget about it
# can say concentrations below...
ab38 <- ab %>%
  filter(Assay=='AB38')
unique(ab38$sampleLOD_ngperg)
mean(unique(ab38$sampleLOD_ngperg), na.rm=TRUE)


### AB40
ab40 <- ab %>%
  filter(Assay=='AB40')
# PSEN1 clutch1: all NA, 4 replicates (2 of which diluted 1:2)
# PSEN1 clutch2: all NA, 4 replicates (2 of which diluted 1:2)
# so simple, can mark a cross below LOD
# PSEN12doublelarvae & uninjlarvae all NA, each 2 replicates (no diluted), will skip
# 3 'sample2' left:
# PSEN2; SCR4x; uninj, all have 4 replicates (2 of which diluted 1:2)

# average those in the plot:
ab40av <- ab40 %>%
  group_by(sample2) %>%
  summarise_at(vars(ngPergProt),
               list(
                 ABngperg= ~ mean(.)
               )) %>%
  # add back sample LOD etc.
  left_join(.,
            ab40[!duplicated(ab40$sample2), c('sample2', 'bradfordConc', 'sampleLOD_ngperg')],
            by='sample2') %>%
  # delete uninjlarvae & PSEN12doublelarvae
  filter(! sample2 %in% c('PSEN12doublelarvae', 'uninjlarvae'))

# order the samples
ab40av$sample2 <- factor(ab40av$sample2, levels=c('uninj', 'SCR4x', 'PSEN1clutch1', 'PSEN1clutch2', 'PSEN2'))


# prepare data to place the small lines showing the LODs
xsta <- (1:nrow(ab40av)) - 0.4
xsto <- (1:nrow(ab40av)) + 0.4

# add a column which is simply to add a small cross below LOD for samples which could not be measured
# e.g. LOD / 3, so it places the cross below
ab40av$cross <- ab40av$sampleLOD_ngperg/3
# ! only for samples which could not be measured, so switch cross to NA for the others
ab40av[which(!is.na(ab40av$ABngperg)), 'cross'] <- NA 

# ready to plot
gg40 <- ggplot(ab40av, aes(x=sample2, y=ABngperg, colour=sample2)) +
  geom_point(size=1) +
  geom_segment(aes(x=xsta, xend=xsto, y=sampleLOD_ngperg, yend=sampleLOD_ngperg),
               linewidth=0.35, colour='black') +
  geom_point(aes(x=sample2, y=cross), shape=4, size=0.8) + # divided by 3 is arbitrary, just placing a cross below LOD to show sample could not be measured
  scale_colour_manual(values=c('#697a87', '#697a87', '#db5072', '#db5072', '#78ac63')) +
  theme_minimal() +
  theme(
    panel.grid.minor.y=element_blank(),
    #axis.text.x=element_text(angle=45, hjust=1),
    axis.text.x=element_blank(),
    axis.title.x=element_blank(),
    axis.title.y=element_text(size=9, margin=margin(t=0, r=0, b=0, l=0)),
    axis.text.y=element_text(size=7, margin=margin(t=0, r=-2, b=0, l=0)),
    legend.position='none'
  ) +
  coord_cartesian(ylim=c(0, 245)) +
  ylab('') # annoying with beta symbol, will do manually
gg40

ggsave(here('221123_PSENab', 'pubplots', 'ab40.pdf'), gg40, width=30, height=35, units='mm', device=cairo_pdf)



### AB42
ab42 <- ab %>%
  filter(Assay=='AB42')
# PSEN1 clutch1: all NA, 4 replicates (2 of which diluted 1:2)
# PSEN1 clutch2: all NA, 4 replicates (2 of which diluted 1:2)
# so simple, can mark a cross below LOD
# PSEN12doublelarvae & uninjlarvae all NA, each 2 replicates (no diluted), will skip
# 3 'sample2' left:

# PSEN2:
# both diluted give NA
# undiluted, one give a measure, other give NA; what should we do?

# SCR4x:
# diluted, one NA, one OK
# undiluted: both OK

# uninj:
# diluted, both OK
# undiluted, one NA, one OK (surprisingly)

# solution could be to replace NA by LOD and do average
# if all samples are NA, do nothing

ab42[which(ab42$sample2 %in% c('PSEN2', 'SCR4x', 'uninj') & is.na(ab42$CalcConcentration)), 'ngPergProt'] <- 
  ab42[which(ab42$sample2 %in% c('PSEN2', 'SCR4x', 'uninj') & is.na(ab42$CalcConcentration)), 'sampleLOD_ngperg']

# now average
ab42av <- ab42 %>%
  group_by(sample2) %>%
  summarise_at(vars(ngPergProt),
               list(
                 ABngperg= ~ mean(.)
               )) %>%
  # add back sample LOD etc.
  left_join(.,
            ab42[!duplicated(ab40$sample2), c('sample2', 'bradfordConc', 'sampleLOD_ngperg')],
            by='sample2') %>%
  # delete uninjlarvae & PSEN12doublelarvae
  filter(! sample2 %in% c('PSEN12doublelarvae', 'uninjlarvae'))

# order the samples
ab42av$sample2 <- factor(ab42av$sample2, levels=c('uninj', 'SCR4x', 'PSEN1clutch1', 'PSEN1clutch2', 'PSEN2'))


# prepare data to place the small lines showing the LODs
xsta <- (1:nrow(ab42av)) - 0.4
xsto <- (1:nrow(ab42av)) + 0.4

# add a column which is simply to add a small cross below LOD for samples which could not be measured
# e.g. LOD / 3, so it places the cross below
ab42av$cross <- ab42av$sampleLOD_ngperg/3
# ! only for samples which could not be measured, so switch cross to NA for the others
ab42av[which(!is.na(ab42av$ABngperg)), 'cross'] <- NA 

# ready to plot
gg42 <- ggplot(ab42av, aes(x=sample2, y=ABngperg, colour=sample2)) +
  geom_point(size=1) +
  geom_segment(aes(x=xsta, xend=xsto, y=sampleLOD_ngperg, yend=sampleLOD_ngperg),
               linewidth=0.35, colour='black') +
  geom_point(aes(x=sample2, y=cross), shape=4, size=0.8) + # divided by 3 is arbitrary, just placing a cross below LOD to show sample could not be measured
  scale_colour_manual(values=c('#697a87', '#697a87', '#db5072', '#db5072', '#78ac63')) +
  theme_minimal() +
  theme(
    panel.grid.minor.y=element_blank(),
    #axis.text.x=element_text(angle=45, hjust=1),
    axis.text.x=element_blank(),
    axis.title.x=element_blank(),
    axis.title.y=element_text(size=9, margin=margin(t=0, r=0, b=0, l=0)),
    axis.text.y=element_text(size=7, margin=margin(t=0, r=-2, b=0, l=0)),
    legend.position='none'
  ) +
  coord_cartesian(ylim=c(0, 14.8)) +
  ylab('')
gg42

ggsave(here('221123_PSENab', 'pubplots', 'ab42.pdf'), gg42, width=30, height=35, units='mm', device=cairo_pdf)
