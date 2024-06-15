#####################################################
# ~ ZFAD: main script for MiSeq run 25/05/2022 ~
#
# note, also plots data from 220915_miseq
#
# Francois Kroll 2023
# francois@kroll.be
#####################################################

# expecting fastq in many folders in /data/220524_miseq/


# packages & functions ----------------------------------------------------

library(FramebyFrame) # to use function gatherFiles from pathUtilities
library(here)
library(tibble)

source(here('utilities', 'ggFrameshift_v2.R'))


# gather all files --------------------------------------------------------
# want all fastq files in one folder

# use gatherFiles function, from FramebyFrame (pathUtilities)
# see there for comments

# gatherFiles(parent='~/Dropbox/phd/220524_miseq/data/220524_miseq/',
#             output='~/Dropbox/phd/220524_miseq/data/reads/')

# at this point, can delete folder 220524_miseq


# -------------------------------------------------------------------------

# up until ampliCan: see README.md file -----------------------------------

# -------------------------------------------------------------------------


# AMPLICAN ----------------------------------------------------------------

# packages & functions for ampliCan ---------------------------------------

library(BiocManager)
# BiocManager::install('amplican', force=TRUE)
# some issues with dependencies, I said NO to "Do you want to install from sources...?" and that seemed to work
library(amplican)

# if re-run ampliCan to solve an issue; it cannot overwrite files
# so if output folders exist already; delete them before running again ampliCan
# below, small function to do so
resetampliCan <- function() {
  if (dir.exists(paste(output_folder, 'alignments', sep='/'))) {
    unlink((paste(output_folder, 'alignments', sep='/')), recursive=TRUE)
  }
  
  if (dir.exists(paste(output_folder, 'reports', sep='/'))) {
    unlink((paste(output_folder, 'reports', sep='/')), recursive=TRUE)
  }
  
  if (file.exists(paste(output_folder, 'barcode_reads_filters.csv', sep='/'))) {
    unlink((paste(output_folder, 'barcode_reads_filters.csv', sep='/')), recursive=TRUE)
  }
  
  if (file.exists(paste(output_folder, 'config_summary.csv', sep='/'))) {
    unlink((paste(output_folder, 'config_summary.csv', sep='/')), recursive=TRUE)
  }
  
  if (file.exists(paste(output_folder, 'RunParameters.txt', sep=''))) {
    unlink((paste(output_folder, 'RunParameters.txt', sep='/')), recursive=TRUE)
  }
}




# run ampliCan on may 2022 MiSeq ------------------------------------------

config <- here('220524_miseq', '220524_config.csv')
fastq_folder <- here('220524_miseq', 'data', 'reads', 'filterfastq')

dir.create(here('220524_miseq', 'amplican_output2'))
output_folder <- here('220524_miseq', 'amplican_output2')

resetampliCan()
amplicanPipeline(config=config,
                 fastq_folder=fastq_folder,
                 results_folder=output_folder,
                 min_freq=0.005, cut_buffer=12,
                 event_filter=FALSE)
# see below for cut_buffer explanation; here just to be consistent
# event_filter=FALSE, I am doing filtering already so I think OK to turn off

# normalisation to controls is OFF, see README.md

# frameshift plot -- all --------------------------------------------------

amp <- read.csv(here('220524_miseq', 'amplican_output', 'config_summary.csv'))

### add 220915 MiSeq ###
amp2 <- read.csv(here('220915_miseq', 'amplican_output', 'config_summary.csv'))


write.csv(rbind(amp, amp2), file=here('220524_miseq', 'amplican_output', 'configBoth.csv'), row.names=FALSE)

# path to pooled config

# create folders for plots
dir.create(here('220524_miseq', 'plots'), showWarnings=FALSE)
dir.create(here('220524_miseq', 'pubplots'), showWarnings=FALSE)

amp <- ggFrameshift(amplican=here('220524_miseq', 'amplican_output', 'config_summary.csv'),
                    onlygene='all',
                    onlysource='all',
                    covFilter=TRUE,
                    mincov=20,
                    mutated_col='#5a6974',
                    frameshift_col='#f1876b',
                    exportOrNo=TRUE,
                    width=159.1,
                    height=82.03,
                    exportfull=here('220524_miseq', 'plots', 'all_framestack.pdf'))

# pooled
amp <- ggFrameshift(amplican=here('220524_miseq', 'amplican_output', 'configBoth.csv'),
                    onlygene='all',
                    onlysource='all',
                    covFilter=TRUE,
                    mincov=20,
                    mutated_col='#5a6974',
                    frameshift_col='#f1876b',
                    exportOrNo=TRUE,
                    width=159.1,
                    height=82.03,
                    exportfull=here('220524_miseq', 'plots', 'all_framestackAll.pdf'))

amp %>%
  filter(startsWith(grp, 'ko')) %>%
  summarise_at(vars(edit),
               list(
                 min= ~ min(., na.rm=TRUE),
                 max= ~ max(., na.rm=TRUE),
                 mean= ~ mean(., na.rm=TRUE),
                 sd= ~ sd(., na.rm=TRUE),
                 n= ~ length(.)
               ))

amp %>%
  filter(startsWith(grp, 'ko')) %>%
  summarise_at(vars(frameshift),
               list(
                 min= ~ min(., na.rm=TRUE),
                 max= ~ max(., na.rm=TRUE),
                 mean= ~ mean(., na.rm=TRUE),
                 sd= ~ sd(., na.rm=TRUE),
                 n= ~ length(.)
               ))


# frameshift plots by gene ------------------------------------------------

### sorl1
ggFrameshift(amplican=amp,
             onlygene='sorl1',
             covFilter=TRUE,
             mincov=20,
             mutated_col='#5a6974',
             frameshift_col='#f1876b',
             xtext=FALSE,
             annotateOrNo=TRUE,
             annotateSize=2,
             exportOrNo=TRUE,
             width=65,
             height=44,
             exportfull=here('220524_miseq', 'pubplots', 'sorl1_framestack.pdf'))


### clu
# there are two samples clu.3, but clu.3 was re-ran as clu.3v2, should use those

amp %>%
  filter(locus!='clu.3') %>%
  ggFrameshift(amplican=.,
               onlygene='clu',
               covFilter=TRUE,
               mincov=20,
               mutated_col='#5a6974',
               frameshift_col='#f1876b',
               xtext=FALSE,
               annotateOrNo=TRUE,
               annotateSize=2,
               exportOrNo=TRUE,
               width=68,
               height=42,
               exportfull=here('220524_miseq', 'pubplots', 'clu_framestack.pdf'))


### cd2ap
# there are three samples cd2ap.1, but locus was re-ran as cd2ap.1v2, should use the v2 data
amp %>%
  filter(locus!='cd2ap.1') %>%
  ggFrameshift(amplican=.,
               onlygene='cd2ap',
               covFilter=TRUE,
               mincov=20,
               mutated_col='#5a6974',
               frameshift_col='#f1876b',
               xtext=FALSE,
               annotateOrNo=TRUE,
               annotateSize=2,
               exportOrNo=TRUE,
               width=68,
               height=42,
               exportfull=here('220524_miseq', 'pubplots', 'cd2ap_framestack.pdf'))


### apoea
amp %>%
  ggFrameshift(amplican=.,
               onlygene='apoea',
               covFilter=TRUE,
               mincov=20,
               mutated_col='#5a6974',
               frameshift_col='#f1876b',
               xtext=FALSE,
               annotateOrNo=TRUE,
               annotateSize=2,
               exportOrNo=TRUE,
               width=46,
               height=42,
               exportfull=here('220524_miseq', 'pubplots', 'apoea_framestack.pdf'))


### apoeb
amp %>%
  ggFrameshift(amplican=.,
               onlygene='apoeb',
               covFilter=TRUE,
               mincov=20,
               mutated_col='#5a6974',
               frameshift_col='#f1876b',
               xtext=FALSE,
               ytextOrNo=FALSE,
               ynameOrNo=FALSE,
               annotateOrNo=TRUE,
               annotateSize=2,
               exportOrNo=TRUE,
               width=38,
               height=42,
               exportfull=here('220524_miseq', 'pubplots', 'apoeb_framestack.pdf'))

# summary stats for apoea/apoeb
# can get them by generating frameshift plot for both
tmp <- amp %>%
  ggFrameshift(amplican=.,
               onlygene=c('apoea', 'apoeb'),
               covFilter=TRUE,
               mincov=20,
               mutated_col='#5a6974',
               frameshift_col='#f1876b',
               xtext=FALSE,
               ytextOrNo=FALSE,
               ynameOrNo=FALSE,
               annotateOrNo=TRUE,
               annotateSize=2,
               exportOrNo=FALSE,
               width=NA,
               height=NA,
               exportfull=NA)



### appa
amp %>%
  ggFrameshift(amplican=.,
               onlygene='appa',
               covFilter=TRUE,
               mincov=20,
               mutated_col='#5a6974',
               frameshift_col='#f1876b',
               xtext=FALSE,
               annotateOrNo=TRUE,
               annotateSize=2,
               exportOrNo=TRUE,
               width=49,
               height=42,
               exportfull=here('220524_miseq', 'pubplots', 'appa_framestack.pdf'))

### appb
amp %>%
  ggFrameshift(amplican=.,
               onlygene='appb',
               covFilter=TRUE,
               mincov=20,
               mutated_col='#5a6974',
               frameshift_col='#f1876b',
               xtext=FALSE,
               ytextOrNo=FALSE,
               ynameOrNo=FALSE,
               annotateOrNo=TRUE,
               annotateSize=2,
               exportOrNo=TRUE,
               width=45,
               height=42,
               exportfull=here('220524_miseq', 'pubplots', 'appb_framestack.pdf'))

# for appa: 96.59 ± 8.39% mutated reads, 70.45 ± 23.23% of all reads had a frameshift mutation
# for appb: 95.01 ± 5.65% mutated reads, 82.56 ± 19.05% of all reads had a frameshift mutation

# summary stats for appa/appb
# can get them by generating frameshift plot for both
tmp <- amp %>%
  ggFrameshift(amplican=.,
               onlygene=c('appa', 'appb'),
               covFilter=TRUE,
               mincov=20,
               mutated_col='#5a6974',
               frameshift_col='#f1876b',
               xtext=FALSE,
               ytextOrNo=FALSE,
               ynameOrNo=FALSE,
               annotateOrNo=TRUE,
               annotateSize=2,
               exportOrNo=FALSE,
               width=NA,
               height=NA,
               exportfull=NA)


# for ZFAD text: average of all larvae sequenced --------------------------

amp <- read.csv(here('220524_miseq', 'amplican_output', 'configBoth.csv'))

# need to add octoberMiSeq
# we only need plate1
# plate2 is appa samples, which we do not use here (new gRNAs/samples)
ampoct <- read.csv(here('october2021_MiSeq', 'plate1_amplican', 'config_summary.csv'))

# pool both
all <- rbind(amp, ampoct)

# need to remove some samples

# bit of a loophole, by running ggFrameshift already, we can split ID column in gene, locus, etc.
all2 <- ggFrameshift(amplican=all,
                     covFilter=TRUE,
                     mincov=20,
                     mutated_col='#5a6974',
                     frameshift_col='#f1876b',
                     xtext=FALSE,
                     ytextOrNo=FALSE,
                     ynameOrNo=FALSE,
                     annotateOrNo=TRUE,
                     annotateSize=2,
                     exportOrNo=FALSE,
                     width=NA,
                     height=NA,
                     exportfull=NA)
# remember it will also exclude coverage < mincov here

# now remove some samples:
# csnk1db: not part of ZFAD project
# there are two samples clu.3, but clu.3 was re-ran as clu.3v2, should use those
# there are three samples cd2ap.1, but locus was re-ran as cd2ap.1v2, should use the v2 data
# I think that is all

all2 %>%
  filter(gene != 'csnk1db') %>%
  filter(! locus %in% c('clu.3', 'cd2ap.1')) %>%
  ggFrameshift(amplican=.,
               covFilter=TRUE,
               mincov=20,
               mutated_col='#5a6974',
               frameshift_col='#f1876b',
               xtext=FALSE,
               annotateOrNo=TRUE,
               annotateSize=2,
               exportOrNo=FALSE,
               width=NA,
               height=NA,
               exportfull=NA)

# it says n = 126
# counting manually on figures:
# psen1: 14
# psen2: 15
# appa: 12
# appb: 12
# sorl1: 16
# cd2ap: 18
# apoea: 12
# apoeb: 9
# clu: 18
# total 126


### also mutation rate of double F0-knockout experiment:
all2 %>%
  filter(gene %in% c('psen1', 'psen2', 'apoea', 'apoeb', 'appa', 'appb')) %>%
  # amyloid beta experiment double psen1/2 KOs used psen1 #1 & #3 // psen2 #1 & #2
  filter(! locus %in% c('psen1.2', 'psen2.3')) %>%
  ggFrameshift(amplican=.,
               covFilter=TRUE,
               mincov=20,
               mutated_col='#5a6974',
               frameshift_col='#f1876b',
               xtext=TRUE,
               annotateOrNo=TRUE,
               annotateSize=2,
               exportOrNo=FALSE,
               width=NA,
               height=NA,
               exportfull=NA)
