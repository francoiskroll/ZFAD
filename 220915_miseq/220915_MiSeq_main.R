# main script for MiSeq 15/09/2022

# packages & functions ----------------------------------------------------

library(here)
library(tibble)

source('~/Dropbox/phd/framebyframe/pathUtilities.R')
source('~/Dropbox/phd/utilities/ggFrameshift.R')


# gather all files --------------------------------------------------------
# want all fastq files in one folder

# use gatherFiles function, sourced from pathUtilities.R
# see there for comments

gatherFiles(parent='~/Dropbox/phd/220915_miseq/data/',
            output='~/Dropbox/phd/220915_miseq/data/reads/')

# at this point, can delete folder 220915_miseq


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





# run ampliCan ------------------------------------------------------------

config <- here('220915_config.csv')
fastq_folder <- here('data/reads/filterfastq/')

dir.create(here('amplican_output'))
output_folder <- here('amplican_output')

resetampliCan()
amplicanPipeline(config=config, fastq_folder=fastq_folder, results_folder=output_folder,
                 min_freq=0.005, cut_buffer=12,
                 event_filter=FALSE)
# see below for cut_buffer explanation; here just to be consistent
# event_filter=FALSE, I am doing filtering already so I think OK to turn off

# normalisation to controls is OFF, see README.md

# frameshift plot -- all --------------------------------------------------

amp <- here('amplican_output', 'config_summary.csv')

# create folder for plots
dir.create(here('plots'))

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