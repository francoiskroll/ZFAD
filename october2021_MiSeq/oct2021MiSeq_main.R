# main script to follow for September 2021 (csnk1db samples) + October 2021 MiSeq

# unzip oct2021MiSeq_alldata.zip
# rename folder oct2021MiSeq_alldata

# in sept2021_csnk1db/, unzip sept2021_allMiSeqdata.zip
# rename folder sept2021_allMiSeqdata

# (can delete folder oct2021MiSeq_alldata &  and keep only zip after running gatherFiles below)


# packages & functions ----------------------------------------------------

library(here)
library(tibble)
library(FramebyFrame) # for function gatherFiles (in pathUtilities)

source(here('utilities', 'ggFrameshift_v2.R'))


# gather all files --------------------------------------------------------
# want all fastq files in one folder

# use gatherFiles function, sourced from pathUtilities.R
# see there for comments

# for sept2021 run
gatherFiles(parent='~/Dropbox/phd/october2021_MiSeq/sept2021_csnk1db/sept2021MiSeq_alldata/PayneMiseqSept21-292781502/',
            output='~/Dropbox/phd/october2021_MiSeq/sept2021_csnk1db/allfastq/')

# for oct2021 run
gatherFiles(parent='~/Dropbox/phd/october2021_MiSeq/oct2021MiSeq_alldata/PayneMiseqOct21-295623328/',
            output='~/Dropbox/phd/october2021_MiSeq/allfastq/')


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



# run ampliCan on september 2021 MiSeq ------------------------------------

config <- here('sept2021_csnk1db/sept2021_config.csv')
fastq_folder <- here('sept2021_csnk1db/reads/filterfastq')

dir.create(here('sept2021_csnk1db/plate1_amplican'))
output_folder <- here('sept2021_csnk1db/sept2021_amplican')

resetampliCan()
amplicanPipeline(config=config, fastq_folder=fastq_folder, results_folder=output_folder, min_freq=0.005, average_quality=25, cut_buffer=12)
# see below for cut_buffer explanation; here just to be consistent


# run ampliCan on plate1 --------------------------------------------------

config <- here('plate1_config.csv')
fastq_folder <- here('plate1_data', 'filterfastq')

dir.create(here('plate1_amplican'))
output_folder <- here('plate1_amplican')

resetampliCan()
amplicanPipeline(config=config, fastq_folder=fastq_folder, results_folder=output_folder, min_freq=0.005, average_quality=25, cut_buffer=12)
# cut_buffer = how far it looks on either side of the PAM, default is 5
# I extended it after noticing a sample with a 1-bp deletion at 8 bp of the PAM (clearly made by Cas9 as not present in 3 SCR larvae)


# run ampliCan on plate2 --------------------------------------------------

config <- here('plate2_config.csv')
fastq_folder <- here('plate2_data', 'filterfastq')

dir.create(here('plate2_amplican'))
output_folder <- here('plate2_amplican')

resetampliCan()
amplicanPipeline(config=config, fastq_folder=fastq_folder, results_folder=output_folder, min_freq=0.005, average_quality=25, cut_buffer=12)


# pool together config files ----------------------------------------------
# i.e. pool together september 2021 + october 2021 plate 1 + october 2021 plate 2

# import sept2021 MiSeq and add source column
sept2021 <- read.csv(here('sept2021_csnk1db', 'sept2021_amplican', 'config_summary.csv'))
sept2021 <- sept2021 %>%
  add_column(source='sept2021_plate1', .before='Barcode')

# import oct2021 plate1 MiSeq and add source column
plate1 <- read.csv(here('plate1_amplican', 'config_summary.csv'))
plate1 <- plate1 %>%
  add_column(source='oct2021_plate1', .before='Barcode')

# import oct2021 plate2 MiSeq and add source column
plate2 <- read.csv(here('plate2_amplican', 'config_summary.csv'))
plate2 <- plate2 %>%
  add_column(source='oct2021_plate2', .before='Barcode')

# rbind them all
amplican_all <- rbind(sept2021, plate1, plate2)

# ! some csnk1db samples sent twice, so will have same ID
# add full_ID column, which is composite source_sample
amplican_all <- amplican_all %>%
  mutate(full_ID=paste(amplican_all$source, amplican_all$ID, sep='_'), .after=ID)
  
write.csv(amplican_all, file=here('amplican_all.csv'), row.names=FALSE)

rm(sept2021)
rm(plate1)
rm(plate2)


# deal with csnk1db duplicates --------------------------------------------

### ### ### CHECKPOINT ### ### ###

# if want to start from here;
amplican_all <- read.csv(file=here('october2021_MiSeq', 'amplican_all.csv'))

# some csnk1db were sent twice, might get messy if we do not deal now with duplicates sample ID

# below rows of samples sent twice
# so can choose which duplicate to keep
amplican_all[which(amplican_all$ID %in% amplican_all$ID[duplicated(amplican_all$ID)]) ,
             c('ID', 'full_ID', 'Reads_Filtered', 'Reads_Del', 'Reads_In', 'Reads_Edited', 'Reads_Frameshifted')]

# csnk1db.B_ko_3
# sept2021: coverage = 29x
# oct2021: coverage = 0x DELETE

# csnk1db.B_ko_4
# sept2021: coverage = 383x DELETE
# oct2021: coverage = 3047x
# best is to keep oct2021

# csnk1db.B_ko_4: Are results replicated?
# sept2021: 353/383 = 92% del; 55/383 = 14% in; 353/383 = 92% edited; 353/383 = 92% frameshift
# oct2021: 3044/3047 = 99.9% del; 2/3047 ~ 0% in; 3044/3047 ~ 100% edited; 2287/3047 = 75% frameshift
# could be differences in PCR or sequencing (probably former)

# csnk1db.D_ko_4
# sept2021: coverage = 24x
# oct2021: coverage = 0x DELETE

# csnk1db.D_ko_5
# sept2021: coverage = 20x DELETE
# oct2021: coverage = 3549x 

# ! delete by full_ID so unique
todelete <- c('oct2021_plate1_csnk1db.B_ko_3',
              'sept2021_plate1_csnk1db.B_ko_4',
              'oct2021_plate1_csnk1db.D_ko_4',
              'sept2021_plate1_csnk1db.D_ko_5')

amp <- amplican_all %>%
  subset(! full_ID %in% todelete)

# now all ID should be unique
if(sum(duplicated(amp$ID) != 0)) stop('\t \t \t \t >>> Stop! Still some duplicate samples \n')


# frameshift plot -- all --------------------------------------------------

# create folder for plots
dir.create(here('october2021_MiSeq', 'plots'), showWarnings=FALSE)
dir.create(here('october2021_MiSeq', 'pubplots'), showWarnings=FALSE)

ggFrameshift(amplican=amp,
             onlygene='all',
             onlysource='all',
             covFilter=TRUE,
             mincov=20,
             mutated_col='#5a6974',
             frameshift_col='#f1876b',
             exportOrNo=TRUE,
             width=159.1,
             height=82.03,
             exportfull=here('plots', 'all_framestack.pdf'))

# frameshift plots by gene ------------------------------------------------

# csnk1db
ggFrameshift(amplican=amp,
             onlygene='csnk1db',
             covFilter=TRUE,
             mincov=20,
             mutated_col='#5a6974',
             frameshift_col='#f1876b',
             xtext=FALSE,
             exportOrNo=FALSE,
             width=57,
             height=45,
             exportfull=here('plots', 'csnk1db_framestack.pdf'))

# psen1
ggFrameshift(amplican=amp,
             onlygene='psen1',
             covFilter=TRUE,
             mincov=20,
             mutated_col='#5a6974',
             frameshift_col='#f1876b',
             xtext=FALSE,
             annotateOrNo=TRUE,
             annotateSize=2,
             exportOrNo=TRUE,
             width=62,
             height=42,
             exportfull=here('october2021_MiSeq', 'pubplots', 'psen1_framestack.pdf'))

# psen2
amp_ps2 <- ggFrameshift(amplican=amp,
                        onlygene='psen2',
                        covFilter=TRUE,
                        mincov=20,
                        mutated_col='#5a6974',
                        frameshift_col='#f1876b',
                        xtext=FALSE,
                        annotateOrNo=TRUE,
                        annotateSize=2,
                        exportOrNo=TRUE,
                        width=65,
                        height=42,
                        exportfull=here('october2021_MiSeq', 'pubplots', 'psen2_framestack.pdf'))

# appa
amp2 <- ggFrameshift(amplican=amp,
                     onlygene='appa',
                     covFilter=TRUE,
                     mincov=20,
                     mutated_col='#5a6974',
                     frameshift_col='#f1876b',
                     xtext=FALSE,
                     exportOrNo=FALSE,
                     width=55,
                     height=51.5,
                     exportfull=here('plots', 'appa_framestack.pdf'))

# summary statistics of appa for stats
# need mean ± std of loci other than #2
appasum <- amp2 %>%
  subset(locus=='appa.1' | locus=='appa.3') %>%
  subset(grp=='ko') %>%
  select(ID, gene, locus, locusnum, grp, samplenum, edit, frameshift, nonframeshift) %>%
  pivot_longer(cols=-c(ID, gene, locus, locusnum, grp, samplenum),
               names_to='type',
               values_to='pro') %>%
  group_by(grp, type) %>%
  summarise_at(vars(pro),
               list(
                 mean= ~ mean(.),
                 sd= ~ sd(.),
                 sem= ~ sem(.),
                 nspl= ~ length(.)
               ))
print(appasum)

appasum$mean
appasum$sd